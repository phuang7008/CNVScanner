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

void NotProperlyPairedReadsDestroy(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking){
    uint32_t i;

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (improperly_PR_array[i]->chrom_id) {
            free(improperly_PR_array[i]->chrom_id);
            improperly_PR_array[i]->chrom_id = NULL;
        }

        if (improperly_PR_array[i]->grouped_improperly_PRs) {
            if (improperly_PR_array[i]->num_of_groups >= 0) {
                // should be signed value
                //
                int32_t g_idx = improperly_PR_array[i]->num_of_groups;
                int32_t j;
                for (j=0; j<g_idx; j++) {
                    if (improperly_PR_array[i]->grouped_improperly_PRs[j].mate_ends_hash)
                        kh_destroy(m32, improperly_PR_array[i]->grouped_improperly_PRs[j].mate_ends_hash);
                }
            }

            if (improperly_PR_array[i]->grouped_improperly_PRs != NULL) {
                free(improperly_PR_array[i]->grouped_improperly_PRs);
                improperly_PR_array[i]->grouped_improperly_PRs = NULL;
            }
        }

        if (improperly_PR_array[i]) {
            free(improperly_PR_array[i]);
            improperly_PR_array[i] = NULL;
        }
    }

    if (improperly_PR_array) {
        free(improperly_PR_array);
        improperly_PR_array = NULL;
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

void processPairedReadsWithinTheSameGroup(Not_Properly_Paired_Reads_Array* improperly_PR_array) {
    int32_t g_idx = improperly_PR_array->num_of_groups;   // need to be signed
    if (g_idx < 0) return;

    uint16_t size = improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000;
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
        uint32_t i, group_mate_end=0;
        uint16_t cur_group_size=1, max=0;
        for (i=1; i<j; i++) {
            if (ends[i] - ends[i-1] < 300) {
                // same group
                //
                cur_group_size++;
                group_mate_end = ends[i];
            } else {
                if (cur_group_size > max && cur_group_size >= 2) {
                    max = cur_group_size;
                    cur_group_size = 1;
                }
            }
        }

        if (cur_group_size > max)
            max = cur_group_size;

        // update the num_of_pairs_TLEN_ge_1000 and reset the group_mate_end
        //
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000 = max;
        improperly_PR_array->grouped_improperly_PRs[g_idx].group_mate_end = group_mate_end;

        if (ends != NULL) {
            free(ends);
            ends = NULL;
        }
    }

    improperly_PR_array->num_of_groups++;
}

void processImproperlyPairedReads(Not_Properly_Paired_Reads_Array* improperly_PR_array, bam1_t *rec) {
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

    // Kim suggested to use 140 as cut off
    //
    //if (length < rec->core.l_qseq - 5)
    if (length < rec->core.l_qseq - 10)
        return;

    // check its mate cigar
    //
    uint8_t *m_cigar = bam_aux_get(rec,"MC");
    // some bam file doesn't have the 'MC' tag produced, 
    // in this case, we simply don't use this info
    //
    if (m_cigar == NULL)
        return;

    char *mc_tag = bam_aux2Z(m_cigar);

    // Kim suggested to use 140 as cutoff
    //
    //if (getMateMatchLengthFromMCTag(mc_tag) < rec->core.l_qseq - 5)
    if (getMateMatchLengthFromMCTag(mc_tag) < rec->core.l_qseq - 10)
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
        processPairedReadsWithinTheSameGroup(improperly_PR_array);
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
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_diff_chrom = 0;
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_same_chrom = 0;
    }

    improperly_PR_array->grouped_improperly_PRs[g_idx].group_start_end = rec->core.pos;
    improperly_PR_array->grouped_improperly_PRs[g_idx].total_paired_reads++;

    // both are mapped reads
    //
    improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_mapped_reads_on_same_chrom++;
            
    // we are only take care of right-hand side as the other side will be left-hand size with negative value
    //
    if (rec->core.isize >= 1000) {
        improperly_PR_array->grouped_improperly_PRs[g_idx].num_of_pairs_TLEN_ge_1000++;
        improperly_PR_array->grouped_improperly_PRs[g_idx].group_mate_end = improperly_PR_array->grouped_improperly_PRs[g_idx].group_start + rec->core.isize;

        setValueToKhashBucket32(improperly_PR_array->grouped_improperly_PRs[g_idx].mate_ends_hash, rec->core.mpos, 1);
    }

    // resize the group size of improperly_PR_array dynamically
    //
    if (improperly_PR_array->num_of_groups + 3 > improperly_PR_array->capacity) {
        improperly_PR_array->capacity += 50;
        improperly_PR_array->grouped_improperly_PRs = realloc(improperly_PR_array->grouped_improperly_PRs, improperly_PR_array->capacity * sizeof(Grouped_Not_Properly_Paired_Reads));
        failureExit(improperly_PR_array->grouped_improperly_PRs, "improperly_PR_array->grouped_improperly_PRs");
    }
}

// There are many special cases that might interfere with the downstream analysis
// 1). Two intervals have the same ends, but starts are different more than 500bps
//          eithr merge them or use different ends by 5bps (choose the latter one)
// 2). Two intervals overlaps for more than 50% of each other (not handled)
//
void organizeImproperlyPairedReadArray(Not_Properly_Paired_Reads_Array* improperly_PR_array) {
    
    khash_t(m32) *seen_ends_hash  = kh_init(m32);

    int32_t i;
    for (i=0; i<improperly_PR_array->num_of_groups; i++) {
        uint32_t end = improperly_PR_array->grouped_improperly_PRs[i].group_mate_end;
        khiter_t iter = kh_get(m32, seen_ends_hash, end);
        if (iter == kh_end(seen_ends_hash)) {
            setValueToKhashBucket32(seen_ends_hash, end, 1);
        } else {
            // added to resolve the issues that might be 3 or 4 or more endings with the same value
            //
            improperly_PR_array->grouped_improperly_PRs[i].group_mate_end += 5 * kh_value(seen_ends_hash, iter);
            addValueToKhashBucket32(seen_ends_hash, end, 1);    // increment 1
        }
    }

    kh_destroy(m32, seen_ends_hash);
}

void outputGroupedImproperlyPairedReads(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    FILE *fh = fopen ("Grouped_Improperly_Paired_Reads.txt", "w");
    fileOpenError(fh, "Grouped_Improperly_Paired_Reads.txt");

    fprintf(fh, "chr\tgroup_start\tgroup_start_end\tgroup_mate_end\tlength\tcumulative_TLEN\ttotal_paired_reads\tnum_of_pairs_TLEN_ge_1000\tnum_of_mapped_reads_on_diff_chrom\tnum_of_mapped_reads_on_same_chrom\n");

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        int32_t g;
        for (g=0; g<improperly_PR_array[i]->num_of_groups; g++) {
            uint32_t cumulative_TLEN = 0;
            if (improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end > 0)
                cumulative_TLEN = improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end - \
                                  improperly_PR_array[i]->grouped_improperly_PRs[g].group_start;

            fprintf(fh, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", improperly_PR_array[i]->chrom_id,
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start,  \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start_end,    \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_mate_end,     \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].group_start_end -   \
                            improperly_PR_array[i]->grouped_improperly_PRs[g].group_start,   \
                    cumulative_TLEN,
                    improperly_PR_array[i]->grouped_improperly_PRs[g].total_paired_reads,    \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_pairs_TLEN_ge_1000,  \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_mapped_reads_on_diff_chrom, \
                    improperly_PR_array[i]->grouped_improperly_PRs[g].num_of_mapped_reads_on_same_chrom
            );
        }
    }

    fclose(fh);
}
