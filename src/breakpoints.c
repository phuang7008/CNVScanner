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
    pread_x_bpts_array->preads_x_per_anchor_bpt_arr = calloc(pread_x_bpts_array->size, sizeof(khash_t(khIntPrArray)*));

    uint32_t i;
    for (i=0; i<pread_x_bpts_array->size; i++) {
        pread_x_bpts_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);
                                            
        pread_x_bpts_array->preads_x_per_anchor_bpt_arr[i] = kh_init(khIntPrArray);

    }
}

void PairedReadsAcrossBreakpointsArrayDestroy(Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array) {
    uint32_t i;
    for (i=0; i< pread_x_bpts_array->size; i++) {
        if (pread_x_bpts_array->chrom_ids[i] != NULL) {
            free(pread_x_bpts_array->chrom_ids[i]);
            pread_x_bpts_array->chrom_ids[i] = NULL;
        }

        PairedReadsAcrossABreakpointPerAnchorArrayDestroy(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[i]);
        cleanKhashIntPrArray(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[i]);
        kh_destroy(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_arr[i]);
    }

    if (pread_x_bpts_array->chrom_ids != NULL) {
        free(pread_x_bpts_array->chrom_ids);
        pread_x_bpts_array->chrom_ids = NULL;
    }

    if (pread_x_bpts_array->preads_x_per_anchor_bpt_arr != NULL) {
        free(pread_x_bpts_array->preads_x_per_anchor_bpt_arr);
        pread_x_bpts_array->preads_x_per_anchor_bpt_arr = NULL;
    }
}

void PairedReadsAcrossABreakpointPerAnchorArrayDestroy(khash_t(khIntPrArray) * preads_x_per_anchor_bpt_arr) {
    khint_t k;
    for (k=kh_begin(preads_x_per_anchor_bpt_arr); k!=kh_end(preads_x_per_anchor_bpt_arr); ++k) {
        if (kh_exist(preads_x_per_anchor_bpt_arr, k)) {
            cleanKhashStrInt(kh_value(preads_x_per_anchor_bpt_arr, k)->seen_paired_read_hash);
            // the followings will be handled at cleanKhashIntPrArray()
            //free(kh_value(preads_x_per_anchor_bpt_arr, k)->pread_x_a_bpt->read_name);
            //free(kh_value(preads_x_per_anchor_bpt_arr, k)->pread_x_a_bpt);
        }
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

// This method will process all breakpoints of a chromosome 
// breakpoints will be grouped if they are 5 bp away
//
void storePairedReadsAcrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array, uint32_t pr_chr_ind, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh) {

    // obtain breakpoint coordinate array and sort it
    //
    fprintf(stderr, "Total number of breakpoints: %"PRIu32"\n", bpt_arr->bpts_per_chr[bpt_chr_idx].size);
    uint32_t *sorted_breakpoints = calloc(bpt_arr->bpts_per_chr[bpt_chr_idx].size, sizeof(uint32_t));
    getSortedBreakpointArray(sorted_breakpoints, bpt_arr, bpt_chr_idx);

    // need to find out the anchor breakpoints for grouping neighboring breakpoints <= 5 bps away
    //
    khash_t(m32) *anchor_breakpoints_hash = kh_init(m32);
    uint32_t num_of_anchors = recordAnchorBreakpoints(sorted_breakpoints, anchor_breakpoints_hash, bpt_arr, bpt_chr_idx);

    // need to store the processed breakpoints, so that the breakpoint will be handled only once
    //
    khash_t(m32) *seen_breakpoints_hash = kh_init(m32);

    // walk through breakpoints of the chromosome with the index == bpt_chr_idx
    //
    int absent;
    uint32_t k;
    for (k=0; k<bpt_arr->bpts_per_chr[bpt_chr_idx].size; k++) {
        // breakpoint position is 0-based
        //
        uint32_t bpt_pos = bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].breakpoint_position;
        //if (bpt_pos == 245975) {
        //    printf("stop\n");
        // }

        // get the anchor breakpoint for this breakpoint
        //
        uint32_t current_anchor_pos=0;
        khiter_t iter_anchor = kh_get(m32, anchor_breakpoints_hash, bpt_pos);
        if (iter_anchor == kh_end(anchor_breakpoints_hash)) {
            continue;
        } else {
            current_anchor_pos = kh_value(anchor_breakpoints_hash, iter_anchor);

            // check if current anchor position exist
            //
            iter_anchor = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], current_anchor_pos);
            if (iter_anchor == kh_end(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind])) {
                // Starts a new group with the current anchor position
                // and initialize current anchor breakpoint bucket
                //
                iter_anchor = kh_put(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], current_anchor_pos, &absent);
                if (absent) {
                    kh_key(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor) = current_anchor_pos;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor) = \
                                    calloc(1, sizeof(Paired_Reads_Across_Per_Anchor_Breakpoint_Array));
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->capacity = PR_INIT_SIZE;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->size = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->my_group_size = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash = kh_init(khStrInt);
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->pread_x_a_bpt = \
                                    calloc(PR_INIT_SIZE, sizeof(Paired_Reads_Across_A_Breakpoint));
                }
            }
        }
            
        // Now get the iterator for the current anchor breakpoint group
        //
        iter_anchor = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], current_anchor_pos);
        if (iter_anchor == kh_end(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind])) {
            fprintf(stderr, "Error: Failed to find the current anchor breakpoint at %"PRIu32"\t2a\n", current_anchor_pos);
            exit(EXIT_FAILURE);
        }

        // check to see if we have seen this breakpoint position
        //
        khiter_t iter_h = kh_get(m32, seen_breakpoints_hash, bpt_pos);
        if (iter_h == kh_end(seen_breakpoints_hash)) {      // first time processing this breakpoint
            iter_h = kh_put(m32, seen_breakpoints_hash, bpt_pos, &absent);
            if (absent) {
                kh_key(seen_breakpoints_hash, iter_h) = bpt_pos;
                kh_value(seen_breakpoints_hash, iter_h) = bpt_pos;
            }

            // only record once for each breakpoint info
            //
            int counter = kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->my_group_size;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->current_breakpoint_count++;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->my_breakpoint_group[counter] = bpt_pos;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->my_group_size++;
            if (bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].type == 1) {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->num_of_soft_clipping++;
            } else {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->num_of_hard_clipping++;
            }
        } else {
            // the breakpoint has been processed before, therefore, we only need to update the pread_x_a_bpt.current_breakpoint_count
            //
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->current_breakpoint_count++;

            if (bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].type == 1) {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->num_of_soft_clipping++;
            } else {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->num_of_hard_clipping++;
            }

            continue;
        }

        // now declare a string 'region' for search, something like chr3:2-3 (breakpoint only 1 bp long)
        // For search region, we need to use bpt_pos not the current anchor position
        //
        char region[200];
        sprintf(region, "%s:%"PRIu32"-%"PRIu32"", bpt_arr->chrom_ids[bpt_chr_idx], bpt_pos, bpt_pos+1);
        hts_itr_t *hts_itr = sam_itr_querys(sfh_idx, header, region);

        if (hts_itr == NULL) {
            fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", region);
            exit(EXIT_FAILURE);
        }

        bam1_t *b = bam_init1();
        while (sam_itr_next(sfh, hts_itr, b) >= 0) {
            // if paired reads are not on the same chromosome, don't process it, just skip it
            //
            if (b->core.tid != b->core.mtid) continue;

            // check if we have seen this read name before. If so, skip it
            //
            khiter_t iter_rn = kh_get(khStrInt, \
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash, bam_get_qname(b));
            if (iter_rn == kh_end(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash)) {
                // first time
                //
                iter_rn = kh_put(khStrInt, \
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash, bam_get_qname(b), &absent);
                if (absent) {
                    kh_key(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash, iter_rn) = strdup(bam_get_qname(b));
                    kh_value(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->seen_paired_read_hash, iter_rn) = bpt_pos;
                }
            } else {
                continue;
            }

            uint32_t ind = kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->size;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->pread_x_a_bpt[ind].current_start_position = b->core.pos;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->pread_x_a_bpt[ind].mate_start_position = b->core.mpos;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->pread_x_a_bpt[ind].read_name = strdup(bam_get_qname(b));
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->pread_x_a_bpt[ind].tlen = abs(b->core.isize);
            if (abs(b->core.isize) >= 1000)
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->num_TLEN_ge_1000++;

            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->size++;


            if (kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->size + 10 >
                        kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor)->capacity)
                dynamicPairedReadsAcrossABreakpointArraySizeIncrease(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_arr[pr_chr_ind], iter_anchor));
        }
        bam_destroy1(b);
        hts_itr_destroy(hts_itr);
    }
    free(sorted_breakpoints);
    kh_destroy(m32, anchor_breakpoints_hash);

    eliminateUnwantedBreakpoints(bpt_arr->chrom_ids[bpt_chr_idx], pread_x_bpts_array, pr_chr_ind, num_of_anchors);
}

// this method will remove those breakpoints where:
// 1). only associated with a single breakpoint AND
// 2). the paired reads across breakpoint with tlen < 1000
//
void eliminateUnwantedBreakpoints(char *chr_id,  Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t pr_chr_ind, uint32_t num_of_anchors) {
    char filename[50];
    sprintf(filename, "%s_breakpoints_sorted.bed", chr_id);

    FILE *fp = fopen(filename, "w");

    uint32_t *anchor_breakpoints = calloc(num_of_anchors, sizeof(uint32_t));
    uint32_t counter=0;

    khint_t k;
    for (k=kh_begin(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind]); \
            k!=kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind]); ++k) {   // at the hash array per chr
        if (kh_exist(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)) {
            anchor_breakpoints[counter] = kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k);
            counter++;
        }
    }

    if (counter != num_of_anchors) {
        fprintf(stderr, "ERROR: counter %"PRIu32"doesn't match to the number of anchors %"PRIu32"\n", counter, num_of_anchors);
        exit(EXIT_FAILURE);
    }

    qsort(anchor_breakpoints, counter, sizeof(uint32_t), compare);

    uint32_t i, j;
    for (j=0; j<counter; j++) {
        k = kh_get(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], anchor_breakpoints[j]);
        if (k == kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind])) {
            fprintf(stderr, "ERROR: can't find the anchor key at position: %"PRIu32"\n", anchor_breakpoints[j]);
            exit(EXIT_FAILURE);
        } else {
            if (!(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->current_breakpoint_count > 1 &&
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->num_TLEN_ge_1000 > 0)) {
                // free memories allocated for the pointer variables in the current bucket
                //
                for (i=0; i<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->size; i++) {
                    if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt[i].read_name != NULL) {
                        free(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt[i].read_name);
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt[i].read_name = NULL;
                    }
                }

                if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt != NULL) {
                    free(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt);
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->pread_x_a_bpt = NULL;
                }

                cleanKhashStrInt(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->seen_paired_read_hash);

                // remove the key-value pair
                //
                kh_del(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k);
            } else {
                // print breakpoint in bedfile format
                // chr_id   start   end   #breakpoints   #tlen>1000
                //
                //fprintf(fp, "%s", preads_x_bpt_arr->chrom_ids[pr_chr_ind]);
                //fprintf(fp, "\t%"PRIu32, kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k));
                //fprintf(fp, "\t%"PRIu32, kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)+1);
                //fprintf(fp, "\t%"PRIu8, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->current_breakpoint_count);
                //fprintf(fp, "\t%"PRIu8"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->num_TLEN_ge_1000);
                fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu8"\t%"PRIu8"\n", preads_x_bpt_arr->chrom_ids[pr_chr_ind], kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k), kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)+1, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->current_breakpoint_count, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[pr_chr_ind], k)->num_TLEN_ge_1000);
            }
        }
    }
    fclose(fp);
    free(anchor_breakpoints);
}

void dynamicPairedReadsAcrossABreakpointArraySizeIncrease(Paired_Reads_Across_Per_Anchor_Breakpoint_Array *preads_x_bpts_arr) {
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
        for (k=kh_begin(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i]); \
                k!=kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i]); ++k) {   // at the hash array per chr
            if (kh_exist(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)) {
                // print Breakpoint key first
                //
                fprintf(fp, "At_Breakpoint\t%"PRIu32"\n", kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k));
                fprintf(fp, "  Breakpoint Group Info");
                int p;
                for (p=0; p<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->my_group_size; p++) {
                    fprintf(fp, "\t%"PRIu32, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->my_breakpoint_group[p]);
                }
                fprintf(fp, "\n");
                fprintf(fp, "  Associated_Number_of_Paired_Reads\t%"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->size);
                fprintf(fp, "  Num_of_Current_Breakpoint_only\t%d\n", \
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->current_breakpoint_count);
                fprintf(fp, "  Num_of_soft_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->num_of_soft_clipping);
                fprintf(fp, "  Num_of_hard_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->num_of_hard_clipping);
                fprintf(fp, "  Number_of_Paired_Reads_with_Insertion_size >= 1000:\t%d\n", \
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->num_TLEN_ge_1000);

                for (j=0; j<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->size; j++) {    // at each breakpoint array
                    fprintf(fp, "    Read Name: %s\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->pread_x_a_bpt[j].read_name);
                    fprintf(fp, "    Start Position 1: %"PRIu32"\n", \
                            kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->pread_x_a_bpt[j].current_start_position);
                    fprintf(fp, "    Start Position 2: %"PRIu32"\n", \
                            kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->pread_x_a_bpt[j].mate_start_position);
                    fprintf(fp, "    TLEN: %"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_arr[i], k)->pread_x_a_bpt[j].tlen);
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

void getSortedBreakpointArray(uint32_t *sorted_breakpoints, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx) {
    uint32_t i;
    for (i=0; i<bpt_arr->bpts_per_chr[bpt_chr_idx].size; i++) {
        sorted_breakpoints[i] = bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[i].breakpoint_position;       
    }

    qsort(sorted_breakpoints, bpt_arr->bpts_per_chr[bpt_chr_idx].size, sizeof(uint32_t), compare);

    // output for debugging
    //
    FILE *fp = fopen("sorted_breakpoints.txt", "w");
    for (i=0; i<bpt_arr->bpts_per_chr[bpt_chr_idx].size; i++) {
        fprintf(fp, "%"PRIu32"\n", sorted_breakpoints[i]);
    }
}

// this method is to set the anchor point for each breakpoint in anchor_breakpoints_hash
// key: current breakpoint; value: current_anchor
//
uint32_t recordAnchorBreakpoints(uint32_t *sorted_breakpoints, khash_t(m32) *anchor_breakpoints_hash, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx) {
    uint32_t i, current_anchor=0, current_breakpoint=0, num_of_anchors=0;
    khash_t(m32) *tmp_anchor_count_hash = kh_init(m32);
    khiter_t iter_h;
    int absent;
    
    for (i=0; i<bpt_arr->bpts_per_chr[bpt_chr_idx].size; i++) {
        //current_breakpoint = bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[i].breakpoint_position;
        current_breakpoint = sorted_breakpoints[i];
        if (current_breakpoint == 59118983) {
            printf("stop\n");
        }

        if (current_anchor == 0 || (current_breakpoint - current_anchor > 5)) {
            // starts a new group, need to initialize current breakpoint bucket
            //
            iter_h = kh_get(m32, anchor_breakpoints_hash, current_breakpoint);
            if (iter_h == kh_end(anchor_breakpoints_hash)) {
                iter_h = kh_put(m32, anchor_breakpoints_hash, current_breakpoint, &absent);
                if (absent) {
                    kh_key(anchor_breakpoints_hash, iter_h) = current_breakpoint;
                    kh_value(anchor_breakpoints_hash, iter_h) = current_breakpoint;
                }
            }
            current_anchor = current_breakpoint;

            // record count information for current anchor breakpoint
            //
            iter_h = kh_get(m32, tmp_anchor_count_hash, current_anchor);
            if (iter_h == kh_end(tmp_anchor_count_hash)) {
                iter_h = kh_put(m32, tmp_anchor_count_hash, current_anchor, &absent);
                if (absent) {
                    kh_key(tmp_anchor_count_hash, iter_h) = current_anchor;
                    kh_value(tmp_anchor_count_hash, iter_h) = 1;
                    if (current_anchor == 59118983) printf("59118983 is added!\n");
                }
            }
        } else {
            // check to see if it is closer to the previous breakpoint (i.e. current_anchor)
            // here the key is current breakpoint, while the value is the current_anchor
            //
            if ( (current_anchor > 0) && (current_breakpoint - current_anchor <= 5) ) {
                // the current_anchor should be the anchor (as value, not key)
                //
                iter_h = kh_get(m32, anchor_breakpoints_hash, current_breakpoint);
                if (iter_h == kh_end(anchor_breakpoints_hash)) {
                    iter_h = kh_put(m32, anchor_breakpoints_hash, current_breakpoint, &absent);
                    if (absent) {
                        kh_key(anchor_breakpoints_hash, iter_h) = current_breakpoint;
                        kh_value(anchor_breakpoints_hash, iter_h) = current_anchor;
                    }
                }

                iter_h = kh_get(m32, tmp_anchor_count_hash, current_anchor);
                if (iter_h == kh_end(tmp_anchor_count_hash)) {
                    fprintf(stderr, "ERROR: current anchor %"PRIu32" doesn't exist at 1!\n", current_anchor);
                    exit(EXIT_FAILURE);
                } else {
                    kh_value(tmp_anchor_count_hash, iter_h)++;
                }
            } else {
                // something is wrong, output the warning message
                //
                fprintf(stderr, "Warning: Need to check at breakpoint: %"PRIu32" with current anchor %"PRIu32"\n", current_breakpoint, current_anchor);
            }
        }
    }

    // output anchor breakpoint grouping info for debugging purpose
    // and get rid of anchors with only single breakpoint
    //
    FILE *fp = fopen("Anchor_Breakpoint_Details_Unique.txt", "w");
    khint_t k;
    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            fprintf(fp, "%d\t%d\n", kh_key(anchor_breakpoints_hash, k), kh_value(anchor_breakpoints_hash, k));

            iter_h = kh_get(m32, tmp_anchor_count_hash, kh_value(anchor_breakpoints_hash, k));
            if (iter_h == kh_end(tmp_anchor_count_hash)) {
                fprintf(stderr, "ERROR: current anchor %"PRIu32" doesn't exist at 2!\n", current_anchor);
                exit(EXIT_FAILURE);
            } else {
                if (kh_value(tmp_anchor_count_hash, iter_h) == 1)
                    kh_del(m32, anchor_breakpoints_hash, k);
            }
        }
    }
    fclose(fp);
    kh_destroy(m32, tmp_anchor_count_hash);

    // now output shrinked anchor hash
    //
    fp = fopen("Anchor_Breakpoint_Used_For_CNV_Detection.txt", "w");
    khash_t(m32) *seen_anchors_hash = kh_init(m32);

    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            fprintf(fp, "%d\t%d\n", kh_key(anchor_breakpoints_hash, k), kh_value(anchor_breakpoints_hash, k));
            
            iter_h = kh_get(m32, seen_anchors_hash, kh_value(anchor_breakpoints_hash, k));
            if (iter_h == kh_end(seen_anchors_hash)) {
                iter_h = kh_put(m32, seen_anchors_hash, kh_value(anchor_breakpoints_hash, k), &absent);
                if (absent) {
                    kh_key(seen_anchors_hash, iter_h) = kh_value(anchor_breakpoints_hash, k);
                    kh_value(seen_anchors_hash, iter_h) = 1;
                    num_of_anchors++;
                }
            }
        }
    }
    fclose(fp);
    kh_destroy(m32, seen_anchors_hash);

    return num_of_anchors;
}
