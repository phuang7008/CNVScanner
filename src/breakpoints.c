/*
 * =====================================================================================
 *
 *       Filename:  breakpoints.c
 *
 *    Description:  the detail implementation of breakpoints header file
 *
 *        Version:  1.0
 *        Created:  10/05/2021 02:04:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX, USA
 *
 * =====================================================================================
 */

#include "breakpoints.h"

void BreakpointArrayInit(Breakpoint_Array *breakpoint_array, Chromosome_Tracking *chrom_tracking) {
    breakpoint_array->size  = chrom_tracking->number_of_chromosomes;
    breakpoint_array->chrom_ids = calloc(breakpoint_array->size, sizeof(char*));
    breakpoint_array->bpts_per_chr = calloc(breakpoint_array->size, sizeof(Breakpoints_Per_Chromosome));

    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        breakpoint_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);

        breakpoint_array->bpts_per_chr[i].capacity = INIT_SIZE;
        breakpoint_array->bpts_per_chr[i].size = 0;
        breakpoint_array->bpts_per_chr[i].breakpoints = 
                calloc(breakpoint_array->bpts_per_chr[i].capacity, sizeof(Breakpoint));
    }
}

void BreakpointArrayDestroy(Breakpoint_Array *breakpoint_array) {
    uint32_t i, j;
    for (i=0; i<breakpoint_array->size; i++) {
        // free chromosome ids
        //
        if (breakpoint_array->chrom_ids[i] != NULL) {
            free(breakpoint_array->chrom_ids[i]);
            breakpoint_array->chrom_ids[i]=NULL;
        }

        // free Breakpoint_Array->bpts_per_chr
        //
        for (j=0; j<breakpoint_array->bpts_per_chr[i].size; j++) {
            if (breakpoint_array->bpts_per_chr[i].breakpoints[j].read_name != NULL) {
                free(breakpoint_array->bpts_per_chr[i].breakpoints[j].read_name);
                breakpoint_array->bpts_per_chr[i].breakpoints[j].read_name = NULL;
            }
        }

        if (breakpoint_array->bpts_per_chr[i].breakpoints != NULL) {
            free(breakpoint_array->bpts_per_chr[i].breakpoints);
            breakpoint_array->bpts_per_chr[i].breakpoints = NULL;
        }
    }

    if (breakpoint_array->chrom_ids != NULL) {
        free(breakpoint_array->chrom_ids);
        breakpoint_array->chrom_ids = NULL;
    }

    if (breakpoint_array->bpts_per_chr != NULL) {
        free(breakpoint_array->bpts_per_chr);
        breakpoint_array->bpts_per_chr = NULL;
    }

    if (breakpoint_array != NULL) {
        free(breakpoint_array);
        breakpoint_array = NULL;
    }
}

uint32_t fetchBreakpointArrayChrIndex(Breakpoint_Array *breakpoint_array, char * chrom_id) {
    uint32_t i;
    for (i=0; i<breakpoint_array->size; i++) {
        if (strcmp(chrom_id, breakpoint_array->chrom_ids[i]) == 0) {
            return i;
        }
    }
    
    return breakpoint_array->size + 10;
}

// for type: 1 -> soft clipping while 2 -> hard clipping
//
void storeCurrentReadBreakpointInfo(uint32_t current_ref_pos, bam1_t *rec, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, khash_t(khStrInt) *breakpoint_pairs_hash, int type) {

    // if paired reads are not on the same chromosome, we won't process it as it is translocation
    //
    if (rec->core.tid != rec->core.mtid) return;

    if (bpt_arr == NULL) return;

    uint32_t bpt_arr_ind = bpt_arr->bpts_per_chr[bpt_chr_idx].size;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].type = type;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_index = bpt_arr_ind;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].breakpoint_position = current_ref_pos;

    // according to sam.h,  @field  pos: 0-based leftmost coordinate of the current read
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_read_start = rec->core.pos;

    // according to sam.h, @field  mpos: 0-based leftmost coordinate of next read in template
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_read_start = rec->core.mpos;

    // according to sam.h, @field  isize: observed template length ("insert size")
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].gap_distance_TLEN = rec->core.isize;

    // according to sam.h, @field  tid: chromosome ID, defined by sam_hdr_t
    //
    //bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_chr_id = rec->core.tid;

    // according to sam.h, @field  mtid: chromosome ID of next read in template, defined by sam_hdr_t
    //
    //bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_chr_id = rec->mtid;

    // @field  data: all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    // the read name should be accessed using bam_get_qname
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].read_name = strdup(bam_get_qname(rec));

    // now we need to record this read name in the hash table for quick lookup
    //
    khiter_t iter = kh_get(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec));
    if (iter == kh_end(breakpoint_pairs_hash)) {
        // first time
        //
        int absent;
        iter = kh_put(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec), &absent);
        if (absent) {
            kh_key(breakpoint_pairs_hash, iter) = strdup(bam_get_qname(rec));
            kh_value(breakpoint_pairs_hash, iter) = bpt_arr_ind;
        }
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_index = -1;
    } else {
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_index = kh_value(breakpoint_pairs_hash, iter);

        // now need to go the mate data and reset its mate index
        //
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[kh_value(breakpoint_pairs_hash, iter)].mate_index = bpt_arr_ind;
    }

    bpt_arr->bpts_per_chr[bpt_chr_idx].size++;

    if (bpt_arr->bpts_per_chr[bpt_chr_idx].size + 10 >
            bpt_arr->bpts_per_chr[bpt_chr_idx].capacity) {

        // need to dynamically increase the array size
        //
        dynamicBreakpointPerChrArraySizeIncrease(&bpt_arr->bpts_per_chr[bpt_chr_idx]);
    }
}

void dynamicBreakpointPerChrArraySizeIncrease(Breakpoints_Per_Chromosome *bpts_per_chr) {
    if (bpts_per_chr == NULL) return;

    bpts_per_chr->capacity += INIT_SIZE;

    if (bpts_per_chr->breakpoints) {
        bpts_per_chr->breakpoints = 
            realloc(bpts_per_chr->breakpoints, bpts_per_chr->capacity * sizeof(Breakpoint));

        failureExit(bpts_per_chr->breakpoints, "Breakpoints_Per_Chromosome->breakpoints");
    } else {
        fprintf(stderr, "Error: The bpts_per_chr.breakpoints is NULL\n");
        exit(EXIT_FAILURE);
    }
}

void outputBreakpointArray(Breakpoint_Array *bpt_arr) {
    FILE *fp = fopen("Breakpoint_Details_for_Debugging.txt", "w");
    uint32_t i, j;
    for (i=0; i<bpt_arr->size; i++) {
        fprintf(fp, "For chromosome: %s\n", bpt_arr->chrom_ids[i]);

        for (j=0; j<bpt_arr->bpts_per_chr[i].size ; j++) {
            fprintf(fp, "  Type: %d\n", bpt_arr->bpts_per_chr[i].breakpoints[j].type);
            fprintf(fp, "  breakpoint_position: %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].breakpoint_position);
            fprintf(fp, "  read_name:     %s\n", bpt_arr->bpts_per_chr[i].breakpoints[j].read_name);
            fprintf(fp, "  current_index: %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].current_index);
            fprintf(fp, "  mate_index:    %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].mate_index);
            fprintf(fp, "  current_read_start: %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].current_read_start);
            fprintf(fp, "  mate_read_start:    %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].mate_read_start);
            fprintf(fp, "  gap_distance_TLEN:  %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].gap_distance_TLEN);
        }
    }
    fclose(fp);
}

void PairedReadsAcrossBreakpointsArrayInit(Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array, Chromosome_Tracking *chrom_tracking) {
    pread_x_bpts_array->size = chrom_tracking->number_of_chromosomes;
    pread_x_bpts_array->chrom_ids = calloc(pread_x_bpts_array->size, sizeof(char*));
    pread_x_bpts_array->preads_x_bpts_per_chr_arr = calloc(pread_x_bpts_array->size, sizeof(khash_t(khIntPrArray)*));

    uint32_t i;
    for (i=0; i<pread_x_bpts_array->size; i++) {
        pread_x_bpts_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);
                                            
        pread_x_bpts_array->preads_x_bpts_per_chr_arr[i] = kh_init(khIntPrArray);

    }
}

void PairedReadsAcrossBreakpointsArrayDestroy(Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array) {
    uint32_t i;
    for (i=0; i< pread_x_bpts_array->size; i++) {
        if (pread_x_bpts_array->chrom_ids[i] != NULL) {
            free(pread_x_bpts_array->chrom_ids[i]);
            pread_x_bpts_array->chrom_ids[i] = NULL;
        }

        cleanKhashIntPrArray(pread_x_bpts_array->preads_x_bpts_per_chr_arr[i]);
    }

    if (pread_x_bpts_array->chrom_ids != NULL) {
        free(pread_x_bpts_array->chrom_ids);
        pread_x_bpts_array->chrom_ids = NULL;
    }

    if (pread_x_bpts_array->preads_x_bpts_per_chr_arr != NULL) {
        free(pread_x_bpts_array->preads_x_bpts_per_chr_arr);
        pread_x_bpts_array->preads_x_bpts_per_chr_arr = NULL;
    }
}

uint32_t fetchPReadsXBreakpointArrayChrIndex(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, char * chrom_id) {
    uint32_t i;
    for (i=0; i<preads_x_bpt_arr->size; i++) {
        if (strcmp(preads_x_bpt_arr->chrom_ids[i], chrom_id) == 0) {
            return i;
        }
    }

    return preads_x_bpt_arr->size + 10;
}

void storePairedReadsAcrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array, uint32_t pr_chr_ind, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh) {
    // need to store processed breakpoints
    //
    khash_t(m32) *seen_breakpoints_hash = kh_init(m32);

    // walk through breakpoints of the chromosome with the bpt_chr_idx
    //
    uint32_t k;
    for (k=0; k<bpt_arr->bpts_per_chr[bpt_chr_idx].size; k++) {
        // breakpoint position is 0-based
        // now form a region for search
        //
        char region[200];
        uint32_t bpt_pos = bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].breakpoint_position;

        // check to see if we have seen this breakpoint position
        //
        khiter_t iter_h = kh_get(m32, seen_breakpoints_hash, bpt_pos);
        int absent;
        if (iter_h == kh_end(seen_breakpoints_hash)) {
            iter_h = kh_put(m32, seen_breakpoints_hash, bpt_pos, &absent);
            if (absent) {
                kh_key(seen_breakpoints_hash, iter_h) = bpt_pos;
                kh_value(seen_breakpoints_hash, iter_h) = bpt_pos;
            }
        } else {
            // need to update the pread_x_a_bpt.current_breakpoint_count
            //
            iter_h = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], bpt_pos);
            if (iter_h == kh_end(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind])) {
                fprintf(stderr, "Error: Failed to find the breakpoint key at %"PRIu32"\n", bpt_pos);
                exit(EXIT_FAILURE);
            }
            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->current_breakpoint_count++;

            if (bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].type == 1) {
                kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->num_of_soft_clipping++;
            } else {
                kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->num_of_hard_clipping++;
            }
            continue;
        }

        sprintf(region, "%s:%"PRIu32"-%"PRIu32"", bpt_arr->chrom_ids[bpt_chr_idx], bpt_pos, bpt_pos+1);
        hts_itr_t *hts_itr = sam_itr_querys(sfh_idx, header, region);

        if (hts_itr == NULL) {
            fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", region);
            exit(EXIT_FAILURE);
        }

        // need to record the query read name for quick lookup as we are only need to record paired reads once
        //
        khash_t(khStrInt) *paired_read_name_hash = kh_init(khStrInt);

        bam1_t *b = bam_init1();
        while (sam_itr_next(sfh, hts_itr, b) >= 0) {
            // if paired reads are not on the same chromosome, don't process it, just skip it
            //
            if (b->core.tid != b->core.mtid) continue;

            // check if we have seen this read name before, skip it
            //
            khiter_t iter_rn = kh_get(khStrInt, paired_read_name_hash, bam_get_qname(b));
            if (iter_rn == kh_end(paired_read_name_hash)) {
                // first time
                //
                iter_rn = kh_put(khStrInt, paired_read_name_hash, bam_get_qname(b), &absent);
                if (absent) {
                    kh_key(paired_read_name_hash, iter_rn) = strdup(bam_get_qname(b));
                    kh_value(paired_read_name_hash, iter_rn) = bpt_pos;
                }
            } else {
                continue;
            }

            // now add the current read to the preads_x_bpts_per_chr_arr[pr_chr_idx]
            //
            iter_h = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], bpt_pos);
            if (iter_h == kh_end(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind])) {
                iter_h = kh_put(khIntPrArray, pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], bpt_pos, &absent);
                if (absent) {
                    kh_key(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h) = bpt_pos;
                    kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h) = \
                                                        calloc(1, sizeof(Paired_Reads_Across_A_Breakpoint_Array));
                    kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->capacity = PR_INIT_SIZE;
                    kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->size = 0;
                    kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->pread_x_a_bpt = \
                                                        calloc(PR_INIT_SIZE, sizeof(Paired_Reads_Across_A_Breakpoint));

                    // only record once for each breakpoint info
                    //
                    kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->current_breakpoint_count++;
                    if (bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].type == 1) {
                        kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->num_of_soft_clipping++;
                    } else {
                        kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->num_of_hard_clipping++;
                    }
                }
            }

            uint32_t ind = kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->size;
            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->pread_x_a_bpt[ind].current_start_position = b->core.pos;
            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->pread_x_a_bpt[ind].mate_start_position = b->core.mpos;
            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->pread_x_a_bpt[ind].read_name = strdup(bam_get_qname(b));
            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->pread_x_a_bpt[ind].tlen = abs(b->core.isize);
            if (abs(b->core.isize) >= 1000)
                kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->num_TLEN_ge_1000++;

            kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->size++;


            if (kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->size + 10 >
                        kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h)->capacity)
                dynamicPairedReadsAcrossABreakpointArraySizeIncrease(kh_value(pread_x_bpts_array->preads_x_bpts_per_chr_arr[pr_chr_ind], iter_h));
        }
        bam_destroy1(b);
        cleanKhashStrInt(paired_read_name_hash);
        hts_itr_destroy(hts_itr);
    }

    cleanKhashInt(seen_breakpoints_hash);
}

void dynamicPairedReadsAcrossABreakpointArraySizeIncrease(Paired_Reads_Across_A_Breakpoint_Array *preads_x_bpts_arr) {
    if (preads_x_bpts_arr == NULL) return;

    preads_x_bpts_arr->capacity += PR_INIT_SIZE;

    if (preads_x_bpts_arr->pread_x_a_bpt) {
        preads_x_bpts_arr->pread_x_a_bpt = realloc(preads_x_bpts_arr->pread_x_a_bpt, preads_x_bpts_arr->capacity * sizeof(Paired_Reads_Across_A_Breakpoint));
        failureExit(preads_x_bpts_arr, "Paired_Reads_Across_A_Breakpoint_Array *preads_x_bpts_arr->pread_x_a_bpt");
    } else {
        fprintf(stderr, "Error: The preads_x_bpts_arr is NULL\n");
        exit(EXIT_FAILURE);
    }
}

void outputPairedReadsAcrossBreakpointsArray(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr) {
    FILE *fp = fopen("Paired_Reads_Across_Breakpoints_Array.txt" ,"w");

    uint32_t i, j;
    for (i=0; i<preads_x_bpt_arr->size; i++) {      // at the top array level
        fprintf(fp, "For chromosome: %s\n", preads_x_bpt_arr->chrom_ids[i]);

        khint_t k;
        for (k=kh_begin(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i]); \
                k!=kh_end(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i]); ++k) {   // at the hash array per chr
            if (kh_exist(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)) {
                // print Breakpoint key first
                //
                fprintf(fp, "At_Breakpoint\t%"PRIu32"\n", kh_key(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k));
                fprintf(fp, "  Associated_Number_of_Paired_Reads\t%"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->size);
                fprintf(fp, "  Num_of_Current_Breakpoint_only\t%d\n", \
                        kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->current_breakpoint_count);
                fprintf(fp, "  Num_of_soft_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->num_of_soft_clipping);
                fprintf(fp, "  Num_of_hard_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->num_of_hard_clipping);
                fprintf(fp, "  Number_of_Paired_Reads_with_Insertion_size >= 1000:\t%d\n", \
                        kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->num_TLEN_ge_1000);

                for (j=0; j<kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->size; j++) {    // at each breakpoint array
                    fprintf(fp, "    Read Name: %s\n", kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->pread_x_a_bpt[j].read_name);
                    fprintf(fp, "    Start Position 1: %"PRIu32"\n", \
                            kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->pread_x_a_bpt[j].current_start_position);
                    fprintf(fp, "    Start Position 2: %"PRIu32"\n", \
                            kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->pread_x_a_bpt[j].mate_start_position);
                    fprintf(fp, "    TLEN: %"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_bpts_per_chr_arr[i], k)->pread_x_a_bpt[j].tlen);
                }
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void BreakpointStatsArrayInit(Breakpoint_Stats_Array *bpt_stats_array, Chromosome_Tracking *chrom_tracking) {
    bpt_stats_array->size = chrom_tracking->number_of_chromosomes;
    bpt_stats_array->chrom_ids = calloc(bpt_stats_array->size, sizeof(char*));
    bpt_stats_array->bp_stats_per_chr = calloc(bpt_stats_array->size, sizeof(Breakpoint_Stats_Per_Chromosome));

    uint32_t i;
    for (i=0; i<bpt_stats_array->size; i++) {
        bpt_stats_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);

        bpt_stats_array->bp_stats_per_chr[i].size = 0;
        bpt_stats_array->bp_stats_per_chr[i].capacity = INIT_SIZE;
        bpt_stats_array->bp_stats_per_chr[i].bpts_stats =
            calloc(bpt_stats_array->bp_stats_per_chr[i].capacity, sizeof(Breakpoint_Stats));
    }
}

void BreakpointStatsArrayDestroy(Breakpoint_Stats_Array *bpt_stats_array) {
    uint32_t i;
    for (i=0; i< bpt_stats_array->size; i++) {
        if (bpt_stats_array->chrom_ids[i] != NULL) {
            free(bpt_stats_array->chrom_ids[i]);
            bpt_stats_array->chrom_ids[i] = NULL;
        }

        if (bpt_stats_array->bp_stats_per_chr[i].bpts_stats != NULL) {
            free(bpt_stats_array->bp_stats_per_chr[i].bpts_stats);
            bpt_stats_array->bp_stats_per_chr[i].bpts_stats = NULL;
        }
    }

    if (bpt_stats_array->chrom_ids != NULL) {
        free(bpt_stats_array->chrom_ids);
        bpt_stats_array->chrom_ids = NULL;
    }

    if (bpt_stats_array->bp_stats_per_chr != NULL) {
        free(bpt_stats_array->bp_stats_per_chr);
        bpt_stats_array->bp_stats_per_chr = NULL;
    }
    
    //if (bpt_stats_array->bp_stats_per_chr[i].breakpoints !=NULL) { }  // should be handled by BreakpointArrayDestroy()
}

void dynamicBreakpointStatsPerChrSizeIncrease(Breakpoint_Stats_Per_Chromosome *bp_stats_per_chr) {
    if (bp_stats_per_chr == NULL) return;

    bp_stats_per_chr->capacity += INIT_SIZE;

    if (bp_stats_per_chr->bpts_stats) {
        bp_stats_per_chr->bpts_stats = 
            realloc(bp_stats_per_chr->bpts_stats, bp_stats_per_chr->capacity * sizeof(Breakpoint_Stats));
        failureExit(bp_stats_per_chr->bpts_stats, "Breakpoint_Stats_Per_Chromosome->bpts_stats");
    } else {
        fprintf(stderr, "Error: The Breakpoint_Stats_Per_Chromosome->bpts_stats is NULL\n");
        exit(EXIT_FAILURE);
    }
}
