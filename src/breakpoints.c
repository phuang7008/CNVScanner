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

void BreakpointArrayInit(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        breakpoint_array[i] = calloc(1, sizeof(Breakpoint_Array));
        breakpoint_array[i]->chrom_id = strdup(chrom_tracking->chromosome_ids[i]);

        breakpoint_array[i]->size = 0;
        breakpoint_array[i]->capacity = INIT_SIZE;
        breakpoint_array[i]->breakpoints = calloc(breakpoint_array[i]->capacity, sizeof(Breakpoint));
    }
}

void BreakpointArrayDestroy(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i, j;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (breakpoint_array[i]->chrom_id) {
            free(breakpoint_array[i]->chrom_id);
            breakpoint_array[i]->chrom_id = NULL;
        }

        for (j=0; j<breakpoint_array[i]->size; j++) {
            // free Breakpoint_Array->breakpoints array
            //
            //if (breakpoint_array[i]->breakpoints[j].read_name != NULL) {
            //    free(breakpoint_array[i]->breakpoints[j].read_name);
            //    breakpoint_array[i]->breakpoints[j].read_name = NULL;
            //}
        }

        if (breakpoint_array[i]->breakpoints != NULL) {
            free(breakpoint_array[i]->breakpoints);
            breakpoint_array[i]->breakpoints = NULL;
        }

        if (breakpoint_array[i] != NULL) {
            free(breakpoint_array[i]);
            breakpoint_array[i] = NULL;
        }
    }

    if (breakpoint_array != NULL) {
        free(breakpoint_array);
        breakpoint_array = NULL;
    }
}

uint32_t fetchBreakpointArrayChrIndex(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking, uint32_t chrom_idx) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (strcmp(chrom_tracking->chromosome_ids[chrom_idx], breakpoint_array[i]->chrom_id) == 0) {
            return i;
        }
    }
    
    fprintf(stderr, "Cann't locate the chrom index for the input chromosome %s in function fetchBreakpointArrayChrIndex\n", chrom_tracking->chromosome_ids[chrom_idx]);
    exit(EXIT_FAILURE);
}

// for type: 1 -> soft clipping while 2 -> hard clipping
//
void storeCurrentReadBreakpointInfo(uint32_t current_ref_pos, bam1_t *rec, Breakpoint_Array *bpt_arr, int type) {

    // if paired reads are not on the same chromosome, we won't process it as it is translocation
    //
    if (rec->core.tid != rec->core.mtid) return;

    if (bpt_arr == NULL) return;

    bpt_arr->breakpoints[bpt_arr->size].type = type;
    //bpt_arr->breakpoints[bpt_arr->size].current_index = bpt_arr->size;
    bpt_arr->breakpoints[bpt_arr->size].breakpoint_position = current_ref_pos;

    // according to sam.h,  @field  pos: 0-based leftmost coordinate of the current read
    //
    bpt_arr->breakpoints[bpt_arr->size].current_read_start = rec->core.pos;

    // according to sam.h, @field  mpos: 0-based leftmost coordinate of next read in template
    //
    bpt_arr->breakpoints[bpt_arr->size].mate_read_start = rec->core.mpos;

    // according to sam.h, @field  isize: observed template length ("insert size")
    //
    bpt_arr->breakpoints[bpt_arr->size].gap_distance_TLEN = rec->core.isize;

    // according to sam.h, @field  tid: chromosome ID, defined by sam_hdr_t
    //
    //bpt_arr->breakpoints[bpt_arr->size].current_chr_id = rec->core.tid;

    // according to sam.h, @field  mtid: chromosome ID of next read in template, defined by sam_hdr_t
    //
    //bpt_arr->breakpoints[bpt_arr->size].mate_chr_id = rec->mtid;

    // @field  data: all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    // the read name should be accessed using bam_get_qname
    //
    //bpt_arr->breakpoints[bpt_arr->size].read_name = strdup(bam_get_qname(rec));
    //bpt_arr->breakpoints[bpt_arr->size].read_name = NULL;

    // now we need to record this read name in the hash table for quick lookup
    //
    /*khiter_t iter = kh_get(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec));
    if (iter == kh_end(breakpoint_pairs_hash)) {
        // first time
        //
        int absent;
        iter = kh_put(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec), &absent);
        if (absent) {
            kh_key(breakpoint_pairs_hash, iter) = strdup(bam_get_qname(rec));
            kh_value(breakpoint_pairs_hash, iter) = bpt_arr->size;  // recorded at current index which is the current size
        }
        bpt_arr->breakpoints[bpt_arr->size].mate_index = -1;
    } else {
        bpt_arr->breakpoints[bpt_arr->size].mate_index = kh_value(breakpoint_pairs_hash, iter);

        // now need to go the mate data and reset its mate index
        //
        bpt_arr->breakpoints[kh_value(breakpoint_pairs_hash, iter)].mate_index = bpt_arr->size;
    }*/

    bpt_arr->size++;

    if (bpt_arr->size + 10 > bpt_arr->capacity) {

        // need to dynamically increase the array size
        //
        dynamicBreakpointPerChrArraySizeIncrease(bpt_arr);
    }
}

void dynamicBreakpointPerChrArraySizeIncrease(Breakpoint_Array *cur_bpts_array) {
    if (cur_bpts_array == NULL) return;

    cur_bpts_array->capacity += INIT_SIZE;

    if (cur_bpts_array->breakpoints) {
        cur_bpts_array->breakpoints = 
            realloc(cur_bpts_array->breakpoints, cur_bpts_array->capacity * sizeof(Breakpoint));

        failureExit(cur_bpts_array->breakpoints, "Breakpoints_Array->breakpoints");
    } else {
        fprintf(stderr, "Error: The cur_bpts_array.breakpoints is NULL\n");
        exit(EXIT_FAILURE);
    }
}

void outputBreakpointArray(Breakpoint_Array *bpt_arr) {
    char filename[100];
    sprintf(filename, "Breakpoint_Details_for_Debugging_%s.txt", bpt_arr->chrom_id);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "For chromosome: %s\n", bpt_arr->chrom_id);

    uint32_t j;
    for (j=0; j<bpt_arr->size; j++) {
        fprintf(fp, "  Type: %d\n", bpt_arr->breakpoints[j].type);
        fprintf(fp, "  breakpoint_position: %"PRIu32"\n", bpt_arr->breakpoints[j].breakpoint_position);
        //fprintf(fp, "  read_name:     %s\n", bpt_arr->breakpoints[j].read_name);
        //fprintf(fp, "  current_index: %"PRId32"\n", bpt_arr->breakpoints[j].current_index);
        //fprintf(fp, "  mate_index:    %"PRId32"\n", bpt_arr->breakpoints[j].mate_index);
        fprintf(fp, "  current_read_start: %"PRIu32"\n", bpt_arr->breakpoints[j].current_read_start);
        fprintf(fp, "  mate_read_start:    %"PRIu32"\n", bpt_arr->breakpoints[j].mate_read_start);
        fprintf(fp, "  gap_distance_TLEN:  %"PRId32"\n", bpt_arr->breakpoints[j].gap_distance_TLEN);
    }
    fclose(fp);
}

void PairedReadsAcrossBreakpointsArrayInit(Paired_Reads_Across_Breakpoints_Array **pread_x_bpts_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        pread_x_bpts_array[i] = calloc(1, sizeof(Paired_Reads_Across_Breakpoints_Array));
        pread_x_bpts_array[i]->chrom_id = strdup(chrom_tracking->chromosome_ids[i]);
                                            
        pread_x_bpts_array[i]->preads_x_per_anchor_bpt_hash = kh_init(khIntPrArray);
    }
}

void PairedReadsAcrossBreakpointsArrayDestroy(Paired_Reads_Across_Breakpoints_Array **pread_x_bpts_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // handle each chromosome
        //
        if (pread_x_bpts_array[i]->chrom_id != NULL) {
            free(pread_x_bpts_array[i]->chrom_id);
            pread_x_bpts_array[i]->chrom_id = NULL;
        }

        cleanKhashIntPrArray(pread_x_bpts_array[i]->preads_x_per_anchor_bpt_hash);
        kh_destroy(khIntPrArray, pread_x_bpts_array[i]->preads_x_per_anchor_bpt_hash);

        if (pread_x_bpts_array[i])
            free(pread_x_bpts_array[i]);
    }
}

uint32_t fetchPReadsXBreakpointArrayChrIndex(Paired_Reads_Across_Breakpoints_Array **preads_x_bpt_arr, Chromosome_Tracking *chrom_tracking, uint32_t chrom_idx) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (strcmp(preads_x_bpt_arr[i]->chrom_id, chrom_tracking->chromosome_ids[chrom_idx]) == 0) {
            return i;
        }
    }

    fprintf(stderr, "Cann't locate the chrom index for the input chromosome %s in function fetchPReadsXBreakpointArrayChrIndex\n", chrom_tracking->chromosome_ids[chrom_idx]);
    exit(EXIT_FAILURE);
}

// This method will process all breakpoints of a chromosome 
// breakpoints will be grouped if they are 5 bp away
//
void storePairedReadsAcrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh, User_Input *user_inputs) {

    // obtain breakpoint coordinate array and sort it
    //
    //fprintf(stderr, "Total number of breakpoints: %"PRIu32"\n", bpt_arr->size);
    uint32_t *sorted_breakpoints = calloc(bpt_arr->size, sizeof(uint32_t));
    getSortedBreakpointArray(sorted_breakpoints, bpt_arr, user_inputs);

    // need to find out the anchor breakpoints to group neighboring breakpoints <= 5 bps away
    //
    khash_t(m32) *anchor_breakpoints_hash = kh_init(m32);
    uint32_t num_of_anchors = recordAnchorBreakpoints(sorted_breakpoints, anchor_breakpoints_hash, bpt_arr, user_inputs);

    // need to store the processed breakpoints, so that the breakpoint will be handled only once
    //
    khash_t(m32) *seen_breakpoints_hash = kh_init(m32);

    // walk through breakpoints of the chromosome with the index == bpt_chr_idx
    //
    int absent;
    uint32_t k;
    uint32_t prev_anchor_pos=0;
    for (k=0; k<bpt_arr->size; k++) {
        // breakpoint position is 0-based
        //
        uint32_t bpt_pos = sorted_breakpoints[k];
        //uint32_t bpt_pos = bpt_arr->breakpoints[k].breakpoint_position;
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
            //printf("%"PRIu32"\n", current_anchor_pos);

            if (prev_anchor_pos > 0 && current_anchor_pos > prev_anchor_pos && current_anchor_pos - prev_anchor_pos > 5) {
                // clean-up the seen_paired_read_hash of the previous anchor and remove it to safe some space
                //
                iter_anchor = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_hash, prev_anchor_pos);
                if (kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash) {
                    cleanKhashStrInt(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash);
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash = NULL;
                }
            }

            prev_anchor_pos = current_anchor_pos;

            // check if current anchor position exist
            //
            iter_anchor = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_hash, current_anchor_pos);
            if (iter_anchor == kh_end(pread_x_bpts_array->preads_x_per_anchor_bpt_hash)) {
                // Starts a new group with the current anchor position
                // and initialize current anchor breakpoint bucket
                //
                iter_anchor = kh_put(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_hash, current_anchor_pos, &absent);
                if (absent) {
                    // create anchor key/value pair
                    //
                    kh_key(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor) = current_anchor_pos;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor) = \
                                    calloc(1, sizeof(Paired_Reads_Across_Per_Anchor_Breakpoint_Array));

                    // setup breakpoint info (grouping and each individual breakpoint)
                    //
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->my_group_size = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->total_paired_reads = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_soft_clipping = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_hard_clipping = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->current_breakpoint_count = 0;

                    // setup paired reads across this anchor breakpoint group
                    //
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->size = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->capacity = PR_INIT_SIZE;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->pread_x_a_bpt = \
                                    calloc(PR_INIT_SIZE, sizeof(Paired_Reads_Across_A_Breakpoint));
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_TLEN_ge_1000 = 0;
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash = kh_init(khStrInt);
                }
            }
        }
            
        // Now get the iterator for the current anchor breakpoint group
        //
        iter_anchor = kh_get(khIntPrArray, pread_x_bpts_array->preads_x_per_anchor_bpt_hash, current_anchor_pos);
        if (iter_anchor == kh_end(pread_x_bpts_array->preads_x_per_anchor_bpt_hash)) {
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
            int counter = kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->my_group_size;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->current_breakpoint_count++;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->my_breakpoint_group[counter] = bpt_pos;
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->my_group_size++;
            if (bpt_arr->breakpoints[k].type == 1) {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_soft_clipping++;
            } else {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_hard_clipping++;
            }
        } else {
            // the breakpoint has been processed before, therefore, we only need to update the pread_x_a_bpt.current_breakpoint_count
            //
            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->current_breakpoint_count++;

            if (bpt_arr->breakpoints[k].type == 1) {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_soft_clipping++;
            } else {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_of_hard_clipping++;
            }

            continue;
        }

        // now declare a string 'region' for search, something like chr3:2-3 (breakpoint only 1 bp long)
        // For search region, we need to use bpt_pos not the current anchor position
        //
        char region[200];
        sprintf(region, "%s:%"PRIu32"-%"PRIu32"", bpt_arr->chrom_id, bpt_pos-300, bpt_pos+300);
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
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash, bam_get_qname(b));
            if (iter_rn == kh_end(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash)) {
                // first time
                //
                iter_rn = kh_put(khStrInt, \
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash, bam_get_qname(b), &absent);
                if (absent) {
                    kh_key(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash, iter_rn) = strdup(bam_get_qname(b));
                    kh_value(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->seen_paired_read_hash, iter_rn) = bpt_pos;
                }
            } else {
                continue;
            }

            if (abs(b->core.isize) >= 1000) {
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_TLEN_ge_1000++;
                uint32_t ind = kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->size;
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->pread_x_a_bpt[ind].current_start_position = b->core.pos;
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->pread_x_a_bpt[ind].mate_start_position = b->core.mpos;
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->pread_x_a_bpt[ind].read_name = strdup(bam_get_qname(b));
                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->pread_x_a_bpt[ind].tlen = abs(b->core.isize);

                kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->size++;

                if (kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->size + 5 >
                            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->capacity)
                    dynamicPairedReadsAcrossABreakpointArraySizeIncrease(kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor));
            }

            kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->total_paired_reads++;

            // check point: to prevent cases where the paired reads across a breakpoint > 10000 (in chr16)
            // this doesn't make sense. So I will stop looping if the paired reads here > 100
            //
            if ((strcmp(bpt_arr->chrom_id, "16") == 0 || strcmp(bpt_arr->chrom_id, "chr16") == 0) && 
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->total_paired_reads > 100 &&
                    kh_value(pread_x_bpts_array->preads_x_per_anchor_bpt_hash, iter_anchor)->num_TLEN_ge_1000 >= 3)
                break;
        }
        bam_destroy1(b);
        hts_itr_destroy(hts_itr);
    }
    free(sorted_breakpoints);
    kh_destroy(m32, seen_breakpoints_hash);
    kh_destroy(m32, anchor_breakpoints_hash);

    eliminateUnwantedBreakpoints(bpt_arr->chrom_id, pread_x_bpts_array, num_of_anchors);
}

// this method will remove those breakpoints where:
// 1). only associated with a single breakpoint AND
// 2). the paired reads across breakpoint with tlen < 1000
//
void eliminateUnwantedBreakpoints(char *chr_id, Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t num_of_anchors) {
    char filename[50];
    sprintf(filename, "breakpoints_sorted_%s.txt", chr_id);

    FILE *fp = fopen(filename, "w");
    fileOpenError(fp, filename);

    uint32_t *anchor_breakpoints = calloc(num_of_anchors, sizeof(uint32_t));
    uint32_t counter=0;

    khint_t k;
    for (k=kh_begin(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); \
            k!=kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); ++k) {   // at the hash array per chr
        if (kh_exist(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)) {
            anchor_breakpoints[counter] = kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k);
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
        k = kh_get(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, anchor_breakpoints[j]);
        if (k == kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash)) {
            fprintf(stderr, "ERROR: can't find the anchor key at position: %"PRIu32"\n", anchor_breakpoints[j]);
            exit(EXIT_FAILURE);
        } else {
            if (!(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count > 1 &&
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000 > 0)) {
                // free memories allocated for the pointer variables in the current bucket
                //
                for (i=0; i<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->size; i++) {
                    if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[i].read_name != NULL) {
                        free(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[i].read_name);
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[i].read_name = NULL;
                    }
                }

                if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt != NULL) {
                    free(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt);
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt = NULL;
                }

                if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->seen_paired_read_hash) {
                    cleanKhashStrInt(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->seen_paired_read_hash);
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->seen_paired_read_hash = NULL;
                }

                if (kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k) != NULL) {
                    free(kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k));
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k) = NULL;
                }

                // remove the key-value pair
                //
                kh_del(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k);
            } else {
                // print breakpoint in bedfile format
                // chr_id   start   end   #breakpoints   #tlen>1000
                //
                //fprintf(fp, "%s", preads_x_bpt_arr->chrom_id);
                //fprintf(fp, "\t%"PRIu32, kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k));
                //fprintf(fp, "\t%"PRIu32, kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)+1);
                //fprintf(fp, "\t%"PRIu8, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count);
                //fprintf(fp, "\t%"PRIu8"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000);
                fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu8"\t%"PRIu8"\n", preads_x_bpt_arr->chrom_id, kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k), kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)+1, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000);
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
    char filename[100];
    sprintf(filename, "Paired_Reads_Across_Breakpoints_Array_%s.txt", preads_x_bpt_arr->chrom_id);
    FILE *fp = fopen(filename ,"w");
    fileOpenError(fp, filename);

    //fprintf(fp, "For chromosome: %s\n", preads_x_bpt_arr->chrom_id);

    khint_t k;
    for (k=kh_begin(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); \
            k!=kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); ++k) {   // at the hash array per chr
        if (kh_exist(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)) {
            // print Breakpoint key first
            //
            fprintf(fp, "At_Breakpoint\t%"PRIu32"\n", kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k));
            fprintf(fp, "  Breakpoint Group Info");
            int p;
            for (p=0; p<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->my_group_size; p++) {
                fprintf(fp, "\t%"PRIu32, kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->my_breakpoint_group[p]);
            }
            fprintf(fp, "\n");
            fprintf(fp, "  Associated_Number_of_Paired_Reads\t%"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->total_paired_reads);
            fprintf(fp, "  Num_of_Current_Breakpoint_only\t%d\n", \
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count);
            fprintf(fp, "  Num_of_soft_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_of_soft_clipping);
            fprintf(fp, "  Num_of_hard_clipping:\t%d\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_of_hard_clipping);
            fprintf(fp, "  Number_of_Paired_Reads_with_Insertion_size >= 1000:\t%d\n", \
                    kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000);

            uint32_t j;
            for (j=0; j<kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->size; j++) {    // at each breakpoint array
                fprintf(fp, "    Read Name: %s\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[j].read_name);
                fprintf(fp, "    Start Position 1: %"PRIu32"\n", \
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[j].current_start_position);
                fprintf(fp, "    Start Position 2: %"PRIu32"\n", \
                        kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[j].mate_start_position);
                fprintf(fp, "    TLEN: %"PRIu32"\n", kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->pread_x_a_bpt[j].tlen);
            }
        }
    }
    fclose(fp);
}

void getSortedBreakpointArray(uint32_t *sorted_breakpoints, Breakpoint_Array *bpt_arr, User_Input *user_inputs) {
    uint32_t i;
    for (i=0; i<bpt_arr->size; i++) {
        sorted_breakpoints[i] = bpt_arr->breakpoints[i].breakpoint_position;       
    }

    qsort(sorted_breakpoints, bpt_arr->size, sizeof(uint32_t), compare);

    // output for debugging
    //
    if (user_inputs->debug_ON) {
        /*char filename[100];
        sprintf(filename, "sorted_breakpoints_%s.txt", bpt_arr->chrom_id);
        FILE *fp = fopen(filename, "w");
        fileOpenError(fp, filename);

        for (i=0; i<bpt_arr->size; i++) {
            fprintf(fp, "%"PRIu32"\n", sorted_breakpoints[i]);
        }
        fclose(fp);*/
    }
}

// this method is to set the anchor point for each breakpoint in anchor_breakpoints_hash
// key: current breakpoint; value: current_anchor
//
uint32_t recordAnchorBreakpoints(uint32_t *sorted_breakpoints, khash_t(m32) *anchor_breakpoints_hash, Breakpoint_Array *bpt_arr, User_Input *user_inputs) {
    uint32_t i, current_anchor=0, current_breakpoint=0, num_of_anchors=0;
    khash_t(m32) *tmp_anchor_count_hash = kh_init(m32);
    khiter_t iter_h;
    int absent;
    
    for (i=0; i<bpt_arr->size; i++) {
        //current_breakpoint = bpt_arr->breakpoints[i].breakpoint_position;
        current_breakpoint = sorted_breakpoints[i];

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
    char *filename = calloc(strlen(bpt_arr->chrom_id) + 100, sizeof(char));
    FILE *fp = NULL;

    if (user_inputs->debug_ON) {
        sprintf(filename, "Anchor_Breakpoint_Details_Unique_%s.txt", bpt_arr->chrom_id);
        fp = fopen(filename, "w");
        fileOpenError(fp, filename);
    }

    khint_t k;
    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            if (user_inputs->debug_ON)
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
    if (user_inputs->debug_ON) fclose(fp);
    kh_destroy(m32, tmp_anchor_count_hash);

    // now output shrinked anchor hash
    //
    if (user_inputs->debug_ON) {
        sprintf(filename, "Anchor_Breakpoint_Used_For_CNV_Detection_%s.txt", bpt_arr->chrom_id);
        fp = fopen(filename, "w");
        fileOpenError(fp, filename);
    }
        
    khash_t(m32) *seen_anchors_hash = kh_init(m32);

    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            if (user_inputs->debug_ON)
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
    if (user_inputs->debug_ON) fclose(fp);
    kh_destroy(m32, seen_anchors_hash);
    if (filename) free(filename);

    return num_of_anchors;
}
