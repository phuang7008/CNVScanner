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

void AnchorBreakpointsHashArrayInit(khash_t(m32) **anchor_breakpoints_hash_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        anchor_breakpoints_hash_array[i] = kh_init(m32);
    }
}

void AnchorBreakpointsHashArrayDestroy(khash_t(m32) **anchor_breakpoints_hash_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        kh_destroy(m32, anchor_breakpoints_hash_array[i]);
    }

    if (anchor_breakpoints_hash_array)
        free(anchor_breakpoints_hash_array);
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

void processBreakpoints(Breakpoint_Array *bpt_arr, khash_t(m32) *anchor_breakpoints_hash, User_Input *user_inputs) {
    // obtain breakpoint coordinate array and sort it
    //
    // fprintf(stderr, "Total number of breakpoints: %"PRIu32"\n", bpt_arr->size);
    //
    uint32_t *sorted_breakpoints = calloc(bpt_arr->size, sizeof(uint32_t));
    getSortedBreakpointArray(sorted_breakpoints, bpt_arr, user_inputs);

    // need to find out the anchor breakpoints to group neighboring breakpoints <= 5 bps away
    //
    recordAnchorBreakpoints(sorted_breakpoints, anchor_breakpoints_hash, bpt_arr, user_inputs);

    free(sorted_breakpoints);
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
    khiter_t iter_h;
    int absent;
    
    for (i=0; i<bpt_arr->size; i++) {
        //current_breakpoint = bpt_arr->breakpoints[i].breakpoint_position;
        current_breakpoint = sorted_breakpoints[i];

        if (current_anchor == 0 || (current_breakpoint - current_anchor > BREAKPOINT_DISTANCE_TO_GROUP)) {
            // starts a new group, need to initialize current breakpoint bucket
            //
            iter_h = kh_get(m32, anchor_breakpoints_hash, current_breakpoint);
            if (iter_h == kh_end(anchor_breakpoints_hash)) {
                iter_h = kh_put(m32, anchor_breakpoints_hash, current_breakpoint, &absent);
                if (absent) {
                    kh_key(anchor_breakpoints_hash, iter_h) = current_breakpoint;
                    kh_value(anchor_breakpoints_hash, iter_h) = 1;
                }
            }
            current_anchor = current_breakpoint;

        } else {
            // check to see if it is closer to the previous breakpoint (i.e. current_anchor)
            // here the key is current_anchor, while the value is the current_anchor count
            //
            if ( (current_anchor > 0) && (current_breakpoint - current_anchor <= BREAKPOINT_DISTANCE_TO_GROUP) ) {
                // the current_anchor should be the anchor (as key)
                //
                iter_h = kh_get(m32, anchor_breakpoints_hash, current_anchor);
                if (iter_h == kh_end(anchor_breakpoints_hash)) {
                    fprintf(stderr, "ERROR: can't find current anchor %"PRIu32"\n", current_anchor);
                    exit(EXIT_FAILURE);
                } else {
                    kh_value(anchor_breakpoints_hash, iter_h)++;
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

            if (kh_value(anchor_breakpoints_hash, k) == 1)
                kh_del(m32, anchor_breakpoints_hash, k);
        }
    }
    if (user_inputs->debug_ON) fclose(fp);

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
