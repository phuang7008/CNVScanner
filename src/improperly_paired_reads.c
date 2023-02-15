/*
 * =====================================================================================
 *
 *       Filename:  improperly_paired_reads.c
 *
 *    Description:  The detailed implementation of improperly_paired_reads.h methods
 *
 *        Version:  1.0
 *        Created:  01/23/2023 08:09:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#include "improperly_paired_reads.h"

void NotProperlyPairedReadsInit(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking){
    uint32_t i;

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        improperly_PR_array[i] = calloc(1, sizeof(Not_Properly_Paired_Reads_Array));
        improperly_PR_array[i]->chrom_id = strdup(chrom_tracking->chromosome_ids[i]);

        improperly_PR_array[i]->num_of_groups = -1;
        improperly_PR_array[i]->capacity = INIT_SIZE;
        improperly_PR_array[i]->grouped_improperly_PRs = 
            calloc(improperly_PR_array[i]->capacity, sizeof(Grouped_Not_Properly_Paired_Reads));
    }
}

void getAllUnmappedReads(khash_t(khStrInt) *unmapped_read_hash, hts_idx_t *sfh_idx, bam_hdr_t* header, samFile* sfh) {
    uint64_t total_unmapped_reads=0;
    hts_itr_t *iter_umpd = sam_itr_querys(sfh_idx, header, "*");
    int res, absent;

    if (iter_umpd) {
        bam1_t *b = bam_init1();
        while ((res = sam_itr_next(sfh, iter_umpd, b)) >= 0) {
            khiter_t iter = kh_put(khStrInt, unmapped_read_hash, bam_get_qname(b), &absent);
            if (absent) 
                kh_key(unmapped_read_hash, iter) = strdup(bam_get_qname(b));
            kh_value(unmapped_read_hash, iter) = 1;
            total_unmapped_reads++;
        }
        fprintf(stderr, "total unmapped reads with *\t%"PRIu64"\n", total_unmapped_reads);
        bam_destroy1(b);
        hts_itr_destroy(iter_umpd);
    }
}

uint32_t fetchImproperPRArrayChrIndex(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking* chrom_tracking, uint32_t chrom_index) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (strcmp(chrom_tracking->chromosome_ids[chrom_index], improperly_PR_array[i]->chrom_id) == 0) {
            return i;
        }
    }

    fprintf(stderr, "Cann't locate the chrom index for the input chromosome %s in function fetchImproperPRArrayChrIndex\n", chrom_tracking->chromosome_ids[chrom_index]);
    exit(EXIT_FAILURE);
}

int32_t getMateMatchLengthFromMCTag(char *mate_cigar) {
    char *c = mate_cigar;
    int32_t match_length = 0;

    while (*c && *c != '*') {
        long num = 0;

        // get the value
        //
        if (isdigit((int)*c)) {
            num = strtol(c, &c, 10);
        } else {
            num = 1;
        }

        // get the operator
        //
        if (*c == 'M') {
            match_length += num;
        }

        c++;
    }

    return match_length;
}

void processImproperlyPairedReads(Not_Properly_Paired_Reads_Array* improperly_PR_array, khash_t(khStrInt)* unmapped_read_hash, bam1_t *rec) {
    // Here MAPQ == 0 will be eliminated, no matter the situation
    //
    if (rec->core.qual == 0)
        return;

    // return if mate is unmapped
    //
    if (rec->core.flag & BAM_FMUNMAP)
        return;

    // we are not going to handle INS on different chromosomes
    // so return if paired reads mapped to different chromosomes
    //
    if (rec->core.tid != rec->core.mtid)
        return;

    // check if the current reads is NOT perfect match or 97% match, RETURN
    //
    uint32_t * cigar = bam_get_cigar(rec);      // get cigar info
    uint32_t i;
    int length=0;
    for (i=0; i<rec->core.n_cigar; ++i) {
        int cop = cigar[i] & BAM_CIGAR_MASK;    // operation
        int cln = cigar[i] >> BAM_CIGAR_SHIFT;  // length

        /*
         * The “M” CIGAR op (BAM_CMATCH) is forbidden in PacBio BAM files. 
         * PacBio BAM files use the more explicit ops “X” (BAM_CDIFF) and “=” (BAM_CEQUAL). 
         * PacBio software will abort if BAM_CMATCH is found in a CIGAR field.
         */
        if (cop == BAM_CMATCH || cop == BAM_CEQUAL || cop == BAM_CDIFF)
            length += cln;
    }

    if (length < rec->core.l_qseq - 5)
        return;

    // check its mate cigar
    //
    uint8_t *m_cigar = bam_aux_get(rec,"MC");
    char *mc_tag = bam_aux2Z(m_cigar);

    if (getMateMatchLengthFromMCTag(mc_tag) < rec->core.l_qseq - 5)
        return;

    bool new_group=false;
    int32_t g_idx = improperly_PR_array->num_of_groups;   // need to be signed

    if (g_idx >= 0 && (rec->core.pos - improperly_PR_array->grouped_improperly_PRs[g_idx].group_start_end) > 150) {
        // Start a new group! However,
        // before moving forward, need to handle the mate_ends_hash from the previous group
        // the purpose is to cluster TLEN based on the length info as some group have very
        // different set of TLEN from 1234 to 1234567890 in length. So need to elimiate singletons
        // Here we choose the end positions is because all the starts are closed by (~150bp)
        //
        /*uint16_t size = improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000;
        if (size >= 2) {
            uint32_t * ends = calloc(size, sizeof(uint32_t));
            khiter_t iter;
            uint32_t j=0;
            for (iter=kh_begin(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash); \
                    iter!=kh_end(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash); iter++) {
                if (kh_exist(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash, iter)) {
                    // Note, all ends at the mate_end_hash have the TLEN >= 1000bp when compare to the corresponding start
                    //
                    uint32_t i;
                    for (i=0; i<kh_value(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash, iter); i++) {
                        ends[j] = kh_key(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash, iter);
                        j++;
                    }
                }
            }

            // sort the ends array
            //
            qsort(ends, j, sizeof(uint32_t), compare);

            // group ends into bins so that data within a group have sizes differ < 1000
            // eliminate all singletons and update the num_of_pairs_TLEN_ge_1000 to the 
            // one with the group size of the highest number of members
            //
            uint32_t i, group_mate_end;
            uint16_t cur_group_size=1, max=0;
            for (i=1; i<j; i++) {
                if (ends[i] - ends[i-1] < 300) {
                    // same group
                    //
                    cur_group_size++;
                } else {
                    if (cur_group_size > max && cur_group_size >= 2) {
                        max = cur_group_size;
                        cur_group_size = 1;
                        group_mate_end = ends[i-1];
                    }
                }
            }

            // update the num_of_pairs_TLEN_ge_1000 and reset the group_mate_end
            //
            improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000 = max;
            improperly_PR_array->grouped_improperly_PRs[g_idx].group_mate_end = group_mate_end;

            if (ends != NULL) {
                free(ends);
                ends = NULL;
            }
        }*/

        improperly_PR_array->num_of_groups++;
        g_idx++;
        new_group = true;

    }

    if (g_idx == -1 || new_group) {
        if (g_idx == -1) {
            g_idx = 0;
            improperly_PR_array->num_of_groups = 0;
        }

        improperly_PR_array->grouped_improperly_PRs[g_idx].group_start = rec->core.pos;
        improperly_PR_array->grouped_improperly_PRs[g_idx].group_mate_end = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].total_paired_reads = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000 = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash = kh_init(m32);
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_unmapped_reads = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_diff_chrom = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_same_chrom = 0;
    }

    improperly_PR_array->grouped_improperly_PRs[g_idx].group_start_end = rec->core.pos;
    improperly_PR_array->grouped_improperly_PRs[g_idx].total_paired_reads++;

    khiter_t iter_g = kh_get(khStrInt, unmapped_read_hash, bam_get_qname(rec));
    if (iter_g == kh_end(unmapped_read_hash)) {
        // both are mapped reads
        //
        if (rec->core.tid == rec->core.mtid) {
            improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_same_chrom++;
            
            // we are only take care of right-hand side as the other side will be left-hand size with negative value
            //
            if (rec->core.isize >= 1000) {
                improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000++;
                improperly_PR_array->grouped_improperly_PRs[g_idx].group_mate_end = improperly_PR_array->grouped_improperly_PRs[g_idx].group_start + rec->core.isize;

                addValueToKhashBucket32(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash, rec->core.mpos, 1);
                // we also need to track the TLEN array here
                // TODO
            }
        } else {
            improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_diff_chrom++;
        }
    } else {
        // one of the paired reads is unmapped reads
        //
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_unmapped_reads++;
    }

    // resize the group size of improperly_PR_array dynamically
    //
    if (improperly_PR_array->num_of_groups + 3 > improperly_PR_array->capacity) {
        improperly_PR_array->capacity += 50;
        improperly_PR_array->grouped_improperly_PRs = realloc(improperly_PR_array->grouped_improperly_PRs, improperly_PR_array->capacity * sizeof(Grouped_Not_Properly_Paired_Reads));
        failureExit(improperly_PR_array->grouped_improperly_PRs, "improperly_PR_array->grouped_improperly_PRs");
    }
}

void outputGroupedImproperlyPairedReads(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    FILE *fh = fopen ("Grouped_Improperly_Paired_Reads.txt", "w");
    fileOpenError(fh, "Grouped_Improperly_Paired_Reads.txt");

    fprintf(fh, "chr\tgroup_start\tgroup_start_end\tgroup_mate_end\tlength\tcumulative_TLEN\ttotal_paired_reads\tnum_of_pairs_TLEN_ge_1000\tnum_of_unmapped_reads\tnum_of_mapped_reads_on_diff_chrom\tnum_of_mapped_reads_on_same_chrom\n");

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        int32_t g;
        for (g=0; g<improperly_PR_array[i]->num_of_groups; g++) {
            uint32_t cumulative_TLEN = 0;
            if (improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end > 0)
                cumulative_TLEN = improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end - \
                                  improperly_PR_array[i]->grouped_improperly_PRs[g].group_start;

            fprintf(fh, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", improperly_PR_array[i]->chrom_id,
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start,  \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start_end,    \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end,     \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start_end -   \
                            improperly_PR_array[i]->grouped_improperly_PRs[g].group_start,   \
                    cumulative_TLEN,
                    improperly_PR_array[i]->grouped_improperly_PRs[g].total_paired_reads,    \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_pairs_TLEN_ge_1000,  \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_unmapped_reads,      \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_mapped_reads_on_diff_chrom, \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_mapped_reads_on_same_chrom
            );
        }
    }

    fclose(fh);
}
