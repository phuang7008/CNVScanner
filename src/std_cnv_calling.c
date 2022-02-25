/*
 * =====================================================================================
 *
 *       Filename:  std_cnv_calling.c
 *
 *    Description:  Implementing details using Standard Deviation approach for CNV calling
 *
 *        Version:  1.0
 *        Created:  01/20/2022 09:23:38 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#include "std_cnv_calling.h"

#include <omp.h>

void generateCNVs(CNV_Array **equal_bin_cnv_array, Binned_Data_Wrapper **equal_size_window_wrappers, Binned_Data_Wrapper **raw_bin_data_wrappers, Paired_Reads_Across_Breakpoints_Array **preads_x_bpts_array, uint32_t number_of_chromosomes, Simple_Stats *the_stats, User_Input *user_inputs) {
#pragma omp parallel shared(the_stats) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        uint32_t cnv_array_index;
        for (cnv_array_index=0; cnv_array_index<number_of_chromosomes; ++cnv_array_index) {
#pragma omp task
          {
            int thread_id = omp_get_thread_num();
            printf("Current thread id in generating CNV calls: %d\n", thread_id);

            // find the corresponding index in equal_size_window_wrappers, raw_bin_data_wrappers and preads_x_bpts_array
            //
            uint32_t equal_bin_index, raw_bin_index, pr_x_bpts_arr_index;
            for (equal_bin_index=0; equal_bin_index<number_of_chromosomes; equal_bin_index++) {
                if (strcmp(equal_size_window_wrappers[equal_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (raw_bin_index=0; raw_bin_index<number_of_chromosomes; raw_bin_index++) {
                if (strcmp(raw_bin_data_wrappers[raw_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (pr_x_bpts_arr_index=0; pr_x_bpts_arr_index<number_of_chromosomes; pr_x_bpts_arr_index++) {
                if (strcmp(preads_x_bpts_array[pr_x_bpts_arr_index]->chrom_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            mergeNeighboringBinsBasedOnZscore(equal_bin_cnv_array[cnv_array_index], equal_size_window_wrappers[equal_bin_index], the_stats, equal_bin_cnv_array[cnv_array_index]->chromosome_id, user_inputs, 2);

            expandMergedCNVWithRawBins(raw_bin_data_wrappers[raw_bin_index], equal_bin_cnv_array[cnv_array_index], the_stats);

            checkBreakpointForEachCNV(equal_bin_cnv_array[cnv_array_index], preads_x_bpts_array[pr_x_bpts_arr_index]);

            outputCNVArray(equal_bin_cnv_array[cnv_array_index], equal_bin_cnv_array[cnv_array_index]->chromosome_id, 2);

          } // end task 
        } // end for loop
#pragma omp taskwait

      } // end single
    } // end parallel
}

// Type 1: raw varying sized bins;  Type 2: Equal window bins
// it will do the merge for both equal bins and raw bins
// but currently, it is used for the equal bin
// I could use the raw bin merging for debugging
//
void mergeNeighboringBinsBasedOnZscore(CNV_Array *cnv_array, Binned_Data_Wrapper *equal_size_window_wrapper, Simple_Stats *equal_window_stats, char *chrom_id, User_Input *user_inputs, int type) {
    FILE *fp = fopen("bin_merging_info_for_debugging.txt", "w");

    uint32_t j;
    double hap_cutoff = equal_window_stats->average_coverage - equal_window_stats->zScore;
    double dup_cutoff = equal_window_stats->average_coverage + equal_window_stats->zScore;

    // flag used to indicate which category current bin belongs to: 
    // 0 -> excluded; 1 -> hap; 2 -> dip; 3 -> dup;
    //
    uint8_t prev_flag=0, cur_flag=0;     
    uint32_t cnv_counter=0, bin_counter=0, bin_capacity=0;
    uint32_t prev_start=0, prev_end=0, prev_len;
    double total_coverage=0;
    Equal_Window_Bin *merged_equal_bin_array = NULL;

    for (j=0; j<equal_size_window_wrapper->size; j++) {
        if ( equal_size_window_wrapper->data[j].start == 20621000) {
            printf("stop\n");
        }
        if (equal_size_window_wrapper->data[j].length == 0) {
            if (prev_flag == 0) {
                continue;
            } else {
                // need to store the previous merged CNV
                // it will go to the following else condition
                //
                cur_flag = 0;
            }
        }

        if (prev_flag == 0) {
            // starts a new CNV where the length > 0
            //
            if (equal_size_window_wrapper->data[j].ave_coverage <= hap_cutoff) {
                prev_flag = 1;
            } else if (equal_size_window_wrapper->data[j].ave_coverage >= dup_cutoff) {
                prev_flag = 3;
            } else {
                //prev_flag = 2;
                continue;
            }

            prev_start = equal_size_window_wrapper->data[j].start;
            prev_end   = equal_size_window_wrapper->data[j].end;
            prev_len   = equal_size_window_wrapper->data[j].length;
            total_coverage = equal_size_window_wrapper->data[j].ave_coverage * prev_len;

            // store current bin info (clean it before using it)
            //
            if (merged_equal_bin_array) {
                free(merged_equal_bin_array);
                merged_equal_bin_array=NULL;
                bin_counter = 0;
                bin_capacity = 0;
            }

            bin_capacity = EQUAL_BIN_SIZE;
            merged_equal_bin_array = calloc(bin_capacity, sizeof(Equal_Window_Bin));
            merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper->data[j].start;
            merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper->data[j].end;
            merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper->data[j].length;
            merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper->data[j].ave_coverage;
            bin_counter++;
        } else {
            if (equal_size_window_wrapper->data[j].length == 0) {
                cur_flag = 0;
            } else if (equal_size_window_wrapper->data[j].ave_coverage <= hap_cutoff) {
                cur_flag = 1;
            } else if (equal_size_window_wrapper->data[j].ave_coverage >= dup_cutoff) {
                cur_flag = 3;
            } else {
                cur_flag = 2;
            }

            // combine neighboring bins together
            // For raw bins, if the distance (or gap) between neighboring bins is <= 300, merge them
            // This value of 300 is based on 2 x seq-length = 2 x 150 = 300 (Qiaoyan's suggestion)
            // This value is defined in the terms.h
            //
            if (cur_flag == prev_flag) {
                if ((type == 2 && equal_size_window_wrapper->data[j].start == prev_end) || 
                    (type == 1 && equal_size_window_wrapper->data[j].start - prev_end <= DISTANCE_CUTOFF)) {
                    // elongate bin size by resetting the end position and length
                    //
                    prev_end  = equal_size_window_wrapper->data[j].end;
                    prev_len += equal_size_window_wrapper->data[j].length;
                    total_coverage += equal_size_window_wrapper->data[j].ave_coverage * equal_size_window_wrapper->data[j].length;

                    // add current window bin info
                    //
                    merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper->data[j].start;
                    merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper->data[j].end;
                    merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper->data[j].length;
                    merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper->data[j].ave_coverage;
                    bin_counter++;

                    if (bin_counter == bin_capacity - 3) {
                        // using 3 here, so that I don't have to check the capacity when doing the extension using
                        // excluded regions
                        //
                        bin_capacity = bin_capacity * 2;
                        dynamicIncreaseBinArraySize(&merged_equal_bin_array, bin_capacity);
                    }
                } else {

                }
            } else {
                // the previous bin flag differs from current bin flag. 
                // Let's store the previous merged bin CNV first if its total length > user specified length
                // NOTE: we won't count this as CNV if the length used for CNV call is < 500bp
                //   OR: if the ratio between prev_len / (prev_end - prev_start) < 0.25, skip it
                //
                if (prev_end - prev_start >= 1000) {
                    if (prev_len >= 500 || ((double)prev_len/((double)prev_end - (double)prev_start)) >= 0.25) {
                        double coverage = total_coverage / prev_len;
                        // further extend the CNV at both ends by one bin (500bp) each 
                        // if both ends connected with excluded regions (Ns regions, repeatmaskers etc.)
                        //
                        storeCurrentCNVtoArray(cnv_array, prev_start, prev_end, prev_len, \
                                                    coverage, merged_equal_bin_array, bin_counter, cnv_counter, prev_flag);
                        extendBothEndsByOneBin(cnv_array, equal_size_window_wrapper, cnv_counter, bin_counter, j);
                        cnv_counter += combineNeighboringCNVs(cnv_array, cnv_counter);
                        cnv_counter++;
                    } else {
                        if (user_inputs->debug_ON)
                            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.f\texclude_0.25\n", chrom_id, prev_start, \
                                prev_end, prev_len, (double)prev_len/((double)prev_end - (double)prev_start));
                    }
                } else {
                    if (user_inputs->debug_ON)
                        fprintf(fp,"%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\ttoo_short\n", chrom_id, prev_start, prev_end, prev_len);
                }

                // clean-up
                //
                //fprintf(stderr, "bin counter is %"PRIu32"\n", bin_counter);
                //fprintf(stderr, "cnv counter is %"PRIu32"\n", cnv_counter);
                if (merged_equal_bin_array != NULL) {
                    free(merged_equal_bin_array);
                    merged_equal_bin_array=NULL;
                }

                // reset the following variable
                //
                prev_start = 0;
                prev_end   = 0;
                prev_len   = 0;
                bin_counter  = 0;

                if (cur_flag == 2 || cur_flag == 0) {
                    // if current bin is diploid (normal) or it is excluded regions, we just skip it
                    // otherwise, we need to record it for the next CNV
                    //
                    prev_flag = 0;
                    continue;
                }

                // current window is not diploid, but a different CNV type. So record values;
                //
                prev_flag = cur_flag;
                prev_start = equal_size_window_wrapper->data[j].start;
                prev_end   = equal_size_window_wrapper->data[j].end;
                prev_len   = equal_size_window_wrapper->data[j].length;
                total_coverage = equal_size_window_wrapper->data[j].ave_coverage * prev_len;

                bin_capacity = EQUAL_BIN_SIZE;
                merged_equal_bin_array = calloc(bin_capacity, sizeof(Equal_Window_Bin));
                merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper->data[j].start;
                merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper->data[j].end;
                merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper->data[j].length;
                merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper->data[j].ave_coverage;
                bin_counter++;
            }
        }
    }

    cnv_array->size = cnv_counter;

    if (merged_equal_bin_array) {
        free(merged_equal_bin_array);
        merged_equal_bin_array=NULL;
    }

    fclose(fp);
}

// extend current CNV on both ends if they are excluded regions
// this is to help merge two separate neighboring CNVs into one CNV
//
void extendBothEndsByOneBin(CNV_Array *cnv_array, Binned_Data_Wrapper *equal_size_window_wrapper, uint32_t cnv_index, uint32_t num_of_bin_used, uint32_t end_index) {
    // **NOTE**: the start_index need to be set here as the num_of_bin_used might be increased 
    //
    int32_t start_index = end_index - num_of_bin_used - 1;      // using signed int as it might go to negative
    if (cnv_array->cnvs[cnv_index].equal_bin_start == 20713500)
        printf("debug stopped\n");

    // check the bin just after the CNV
    //
    if (equal_size_window_wrapper->data[end_index].length == 0 && 
            cnv_array->cnvs[cnv_index].equal_bin_end == equal_size_window_wrapper->data[end_index+1].start) {
        cnv_array->cnvs[cnv_index].equal_bin_end = equal_size_window_wrapper->data[end_index+1].end;
        cnv_array->cnvs[cnv_index].equal_bin_array[num_of_bin_used].start = equal_size_window_wrapper->data[end_index].start;
        cnv_array->cnvs[cnv_index].equal_bin_array[num_of_bin_used].end   = equal_size_window_wrapper->data[end_index].end;
        cnv_array->cnvs[cnv_index].equal_bin_array[num_of_bin_used].length = equal_size_window_wrapper->data[end_index].length;
        cnv_array->cnvs[cnv_index].equal_bin_array[num_of_bin_used].type   = 'X';
        cnv_array->cnvs[cnv_index].equal_bin_array[num_of_bin_used].ave_coverage = \
                                        equal_size_window_wrapper->data[end_index+1].ave_coverage;
        num_of_bin_used++;
    }

    // check the bin just before the CNV
    //
    if (start_index >= 0 && equal_size_window_wrapper->data[start_index].length == 0 &&
            cnv_array->cnvs[cnv_index].equal_bin_start == equal_size_window_wrapper->data[start_index].end) {
        cnv_array->cnvs[cnv_index].equal_bin_start = equal_size_window_wrapper->data[start_index].start;

        // re-order th equal bin array
        //
        int32_t i;      // it might go to negative
        for (i=num_of_bin_used+1; i>0; i--) {
            cnv_array->cnvs[cnv_index].equal_bin_array[i].start  = cnv_array->cnvs[cnv_index].equal_bin_array[i-1].start;
            cnv_array->cnvs[cnv_index].equal_bin_array[i].end    = cnv_array->cnvs[cnv_index].equal_bin_array[i-1].end;
            cnv_array->cnvs[cnv_index].equal_bin_array[i].length = cnv_array->cnvs[cnv_index].equal_bin_array[i-1].length;
            cnv_array->cnvs[cnv_index].equal_bin_array[i].type   = cnv_array->cnvs[cnv_index].equal_bin_array[i-1].type;
            cnv_array->cnvs[cnv_index].equal_bin_array[i].ave_coverage = \
                            cnv_array->cnvs[cnv_index].equal_bin_array[i-1].ave_coverage;
        }

        // now add the before bin
        //
        cnv_array->cnvs[cnv_index].equal_bin_array[0].start  = equal_size_window_wrapper->data[start_index].start;
        cnv_array->cnvs[cnv_index].equal_bin_array[0].end    = equal_size_window_wrapper->data[start_index].end;
        cnv_array->cnvs[cnv_index].equal_bin_array[0].length = equal_size_window_wrapper->data[start_index].length;
        cnv_array->cnvs[cnv_index].equal_bin_array[0].type   = 'X';
        cnv_array->cnvs[cnv_index].equal_bin_array[0].ave_coverage  = \
                            equal_size_window_wrapper->data[start_index].ave_coverage;
    }
}

// The method will first check to see if the previous CNV should be merged with the current CNV 
// if their distance is <= 1000 bp away; if not, it will add the new one
//
void storeCurrentCNVtoArray(CNV_Array *cnv_array, uint32_t start, uint32_t end, uint32_t length, double coverage, Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_size, uint32_t cnv_index, uint8_t cnv_flag) {
    uint32_t k;

    cnv_array->cnvs[cnv_index].equal_bin_start  = start;
    cnv_array->cnvs[cnv_index].equal_bin_end    = end; 
    cnv_array->cnvs[cnv_index].length = length;
    cnv_array->cnvs[cnv_index].ave_coverage = coverage;

    if (cnv_flag == 1) {
        cnv_array->cnvs[cnv_index].cnv_type = 'L';
    } else if (cnv_flag == 3) {
        cnv_array->cnvs[cnv_index].cnv_type = 'P';
    } else {
        fprintf(stderr, "ERROR: cnv type %"PRIu8" is niether dup nor del\n", cnv_flag);
    }

    // store all equal window bin info to the current CNV
    //
    cnv_array->cnvs[cnv_index].size = bin_size;
    cnv_array->cnvs[cnv_index].capacity = bin_size + 3;
    cnv_array->cnvs[cnv_index].equal_bin_array = calloc(cnv_array->cnvs[cnv_index].capacity, sizeof(Equal_Window_Bin));

    for (k=0; k<cnv_array->cnvs[cnv_index].size; k++) {
        cnv_array->cnvs[cnv_index].equal_bin_array[k].start  = merged_equal_bin_array[k].start;
        cnv_array->cnvs[cnv_index].equal_bin_array[k].end    = merged_equal_bin_array[k].end;
        cnv_array->cnvs[cnv_index].equal_bin_array[k].length = merged_equal_bin_array[k].length;
        cnv_array->cnvs[cnv_index].equal_bin_array[k].type   = 'E';           // adding type is delayed to this step for equal bin
        cnv_array->cnvs[cnv_index].equal_bin_array[k].ave_coverage = merged_equal_bin_array[k].ave_coverage;
    }

    // initialize other variables. It is only need to be done once!
    //
    cnv_array->cnvs[cnv_index].raw_bin_start = 0; 
    cnv_array->cnvs[cnv_index].raw_bin_end   = 0; 
    cnv_array->cnvs[cnv_index].left_breakpoint = 0;
    cnv_array->cnvs[cnv_index].left_num_of_TLEN_ge_1000 = 0;
    cnv_array->cnvs[cnv_index].left_num_of_breakpoints = 0;
    //cnv_array->cnvs[cnv_index].final_start = 0;
    cnv_array->cnvs[cnv_index].right_breakpoint = 0;
    cnv_array->cnvs[cnv_index].right_num_of_TLEN_ge_1000 = 0;
    cnv_array->cnvs[cnv_index].right_num_of_breakpoints = 0;
    //cnv_array->cnvs[cnv_index].final_end = 0;

}

int combineNeighboringCNVs(CNV_Array *cnv_array, uint32_t cnv_index) {
    if (cnv_index > 0) {
        if (cnv_array->cnvs[cnv_index].equal_bin_start <= cnv_array->cnvs[cnv_index-1].equal_bin_end + 1000 && 
                cnv_array->cnvs[cnv_index-1].equal_bin_end + 1000 <= cnv_array->cnvs[cnv_index].equal_bin_end) {
            // merge with the previous one and re-calculate the coverage
            //
            cnv_array->cnvs[cnv_index-1].equal_bin_end = cnv_array->cnvs[cnv_index].equal_bin_end;
            cnv_array->cnvs[cnv_index-1].ave_coverage = cnv_array->cnvs[cnv_index].length * cnv_array->cnvs[cnv_index].ave_coverage + \
                        cnv_array->cnvs[cnv_index-1].length * cnv_array->cnvs[cnv_index-1].ave_coverage;
            cnv_array->cnvs[cnv_index-1].length += cnv_array->cnvs[cnv_index].length;
            cnv_array->cnvs[cnv_index-1].ave_coverage /= cnv_array->cnvs[cnv_index-1].length;

            // add the current CNV's merged_equal_bin_array
            //
            uint32_t tmp_size = cnv_array->cnvs[cnv_index-1].size + cnv_array->cnvs[cnv_index].size;
            cnv_array->cnvs[cnv_index-1].capacity = tmp_size + 3;   // give extra 3 for excluded bins on the CNV's both sides
            cnv_array->cnvs[cnv_index-1].equal_bin_array = realloc(cnv_array->cnvs[cnv_index-1].equal_bin_array, \
                                                    cnv_array->cnvs[cnv_index-1].capacity*sizeof(Equal_Window_Bin));
            failureExit(cnv_array->cnvs[cnv_index-1].equal_bin_array, "cnv_array[cnv_index-1].equal_bin_array memory realloc failed\n");

            uint32_t i=0, k=0;
            for (k=cnv_array->cnvs[cnv_index-1].size; k<tmp_size; k++) {
                i = k - cnv_array->cnvs[cnv_index-1].size;
                fprintf(stderr, "i value is %"PRIu32" while k value is %"PRIu32"\n", i, k);
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].start  = cnv_array->cnvs[cnv_index].equal_bin_array[i].start;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].end    = cnv_array->cnvs[cnv_index].equal_bin_array[i].end;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].length = cnv_array->cnvs[cnv_index].equal_bin_array[i].length;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].type   = cnv_array->cnvs[cnv_index].equal_bin_array[i].type;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].ave_coverage = \
                                                cnv_array->cnvs[cnv_index].equal_bin_array[i].ave_coverage;
            }

            // **NOTE**: the size needs to be reset after the for loop!!!
            //
            cnv_array->cnvs[cnv_index-1].size = tmp_size;

            return -1;
        }
    }

    return 0;
}

//                     s1  e1  s2   e2     s3     e3       s4   e4       s5         e5
// CNV array             ===     =====       ======          =====         ===========
// Raw bins     -- --- ---  --- --      ----------       -----          ----  ----- ---- ------ ----  ---
//                                 |
//                              current-CNV
// looping?:      to continue                             to break
//
void expandMergedCNVWithRawBins(Binned_Data_Wrapper *binned_data_wrapper, CNV_Array *cnv_array, Simple_Stats *equal_window_stats) {
    double hap_cutoff = equal_window_stats->average_coverage - equal_window_stats->zScore;
    double dup_cutoff = equal_window_stats->average_coverage + equal_window_stats->zScore;

    // raw_bin_start is used to prevent looping from the start point
    //
    uint32_t i, j, raw_bin_start=0;

    for (i=0; i<cnv_array->size;i++) {
        if (cnv_array->cnvs[i].equal_bin_start == 11231000) {
            printf("stop3\n");
        }
        for (j=raw_bin_start; j<binned_data_wrapper->size; j++) {
            if (binned_data_wrapper->data[j].start == 11230812) {
                printf("stop4\n");
            }
            // check if the distance is within 300bp away
            // the value of 300 bp is based on seq-length (150bp) x 2 (Qiaoyan's suggestion)
            // This is defined in terms.h
            //
            if (cnv_array->cnvs[i].equal_bin_start > binned_data_wrapper->data[j].end + DISTANCE_CUTOFF) { 
                // no intersect, record restart position and then skip
                //
                continue;
            } else if (cnv_array->cnvs[i].equal_bin_end + DISTANCE_CUTOFF < binned_data_wrapper->data[j].start) {
                // no intersect. The current raw bin has pass the current CNV bin
                // no need to continue, just break out
                //
                break;
            } else {
                // they intersect
                //
                raw_bin_start = j;
                uint32_t p;

                // because we added DISTANCE_CUTOFF to the ends of the CNVs
                // check which raw bin index is the closest to the current CNV's start position
                //
                if (binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_start) { 
                    // Handle left hand side to find the left most raw bin for the current CNV
                    //
                    for (p=j+1; p<binned_data_wrapper->size; p++) {
                        // raw bin passes the start position of current CNV, so stop looping
                        //
                        if (cnv_array->cnvs[i].equal_bin_start < binned_data_wrapper->data[p].start)
                            break;
                    }
                    j = p - 1;      // we find the left most raw bin index
                    raw_bin_start = j;

                    // need to use a tmp variable tmp_raw_bin_start for the testing and looping
                    //
                    uint32_t tmp_raw_bin_start = cnv_array->cnvs[i].equal_bin_start;

                    if (binned_data_wrapper->data[j].start <= tmp_raw_bin_start &&
                            tmp_raw_bin_start <= binned_data_wrapper->data[j].end) {
                        // they intersect, record it
                        //
                        cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                        tmp_raw_bin_start = cnv_array->cnvs[i].raw_bin_start;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                        j--;
                    }
                    
                    // the while loop will handle all the left raw bins that don't intersect
                    // with the current CNV
                    //
                    while ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                           (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff 
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        if (tmp_raw_bin_start - binned_data_wrapper->data[j].end <= DISTANCE_CUTOFF) {
                            cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                            tmp_raw_bin_start = cnv_array->cnvs[i].raw_bin_start;
                            addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                            j--;
                        } else {
                            break;
                        }
                    }
                    j = raw_bin_start;
                    continue;
                }

                // Now handle right hand side, 
                // the current j will be the raw bin index that is the closest to the current CNV
                // first, check if they intersect
                //
                uint32_t tmp_raw_bin_end = cnv_array->cnvs[i].equal_bin_end;

                if (binned_data_wrapper->data[j].start <= tmp_raw_bin_end &&
                        tmp_raw_bin_end <= binned_data_wrapper->data[j].end) {
                    cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].end;
                    tmp_raw_bin_end = cnv_array->cnvs[i].raw_bin_end;
                    addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                    j++;
                }

                // now let's extend the CNV further at the right end
                //
                while ((binned_data_wrapper->data[j].ave_coverage <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                       (binned_data_wrapper->data[j].ave_coverage >= dup_cutoff
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                    // check if the distance is within 300bp away
                    // This value of 300 is based on 2 x seq-length = 2 x 150 = 300 (Qiaoyan's suggestion)
                    //
                    if (binned_data_wrapper->data[j].start - tmp_raw_bin_end <= DISTANCE_CUTOFF) {
                        cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].end;
                        tmp_raw_bin_end = cnv_array->cnvs[i].raw_bin_end;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                        j++;
                        raw_bin_start = j;
                    } else {
                        break;
                    }
                }
            }
        }
    }
}

void addRawBinToCNV(Binned_Data_Wrapper *binned_data_wrapper, uint32_t raw_bin_index, CNV *cnv, uint32_t cnv_index) {
    // need to dynamically increase the size of bins used for the current CNV calling
    //
    if (cnv[cnv_index].size + 3 >= cnv[cnv_index].capacity) {
        cnv[cnv_index].capacity += 10;
        cnv[cnv_index].equal_bin_array = realloc(cnv[cnv_index].equal_bin_array, \
                                         cnv[cnv_index].capacity * sizeof(Equal_Window_Bin));
        failureExit(cnv[cnv_index].equal_bin_array, "cnv[cnv_index].equal_bin_array memory realloc failed\n");
    }

    cnv[cnv_index].equal_bin_array[cnv[cnv_index].size].start  = binned_data_wrapper->data[raw_bin_index].start;
    cnv[cnv_index].equal_bin_array[cnv[cnv_index].size].end    = binned_data_wrapper->data[raw_bin_index].end;
    cnv[cnv_index].equal_bin_array[cnv[cnv_index].size].length = binned_data_wrapper->data[raw_bin_index].length;
    cnv[cnv_index].equal_bin_array[cnv[cnv_index].size].ave_coverage = binned_data_wrapper->data[raw_bin_index].ave_coverage;
    cnv[cnv_index].equal_bin_array[cnv[cnv_index].size].type    = 'R';
    cnv[cnv_index].size++;
}

void checkBreakpointForEachCNV(CNV_Array *cnv_array, Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr) {
    // first, we need to get the sorted anchor breakpoints
    //
    uint32_t capacity = PR_INIT_SIZE*10;    // the initial size is set to 2000
    uint32_t *anchor_breakpoints = calloc(capacity, sizeof(uint32_t));
    uint32_t counter=0;

    khint_t k;
    for (k=kh_begin(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); \
            k!=kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash); ++k) {   // at the hash array per chr
        if (kh_exist(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)) {
            anchor_breakpoints[counter] = kh_key(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k);
            counter++;

            // dynamically increase the capacity of anchor breakpoint array
            //
            if (counter + 5 >= capacity) {
                capacity += PR_INIT_SIZE*10;
                anchor_breakpoints = realloc(anchor_breakpoints, capacity * sizeof(uint32_t));
                failureExit(anchor_breakpoints, "anchor_breakpoints memory realloc failed\n");
            }
        }
    }

    qsort(anchor_breakpoints, counter, sizeof(uint32_t), compare);

    // then walk through the CNV array
    //
    uint32_t i, j, breakpoint_start=0;

    for (i=0; i<cnv_array->size;i++) {
        for (j=breakpoint_start; j<counter; j++) {
            if (anchor_breakpoints[j] == 12698930) {
                printf("stop\n");
            }
            uint32_t tmp_start = (cnv_array->cnvs[i].raw_bin_start > 0) ? \
                                    cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
            uint32_t tmp_end = (cnv_array->cnvs[i].raw_bin_end > 0) ? \
                                    cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;

            if (anchor_breakpoints[j] + DISTANCE_CUTOFF < tmp_start) {
                continue;
            } else if (tmp_end + DISTANCE_CUTOFF < anchor_breakpoints[j]) {
                break;
            } else {
                // we consider they intersect.
                // Now we need to find out the breakpoint on the current CNV
                //
                if (anchor_breakpoints[j] == tmp_start) {
                    addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 1);
                } else if (anchor_breakpoints[j] == tmp_end) {
                    addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 2);
                } else if (tmp_start < anchor_breakpoints[j] && anchor_breakpoints[j] < tmp_end) {
                    if (anchor_breakpoints[j] - tmp_start < tmp_end - anchor_breakpoints[j]) {
                        // closer to the start position
                        //
                        addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 1);
                    } else {
                        addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 2);
                    }
                } else if ((anchor_breakpoints[j] >= tmp_start - DISTANCE_CUTOFF) && 
                            (anchor_breakpoints[j] <= tmp_start + DISTANCE_CUTOFF)) {
                    addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 1);
                } else if ((tmp_end + DISTANCE_CUTOFF >= anchor_breakpoints[j]) &&
                            (tmp_end - DISTANCE_CUTOFF <= anchor_breakpoints[j])) {
                    addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j], 2);
                } else {
                    fprintf(stderr, "Warning: anchor breakpoint %"PRIu32" shouldn't be here\n", anchor_breakpoints[j]);
                }
            }
        }
    }
}

// pos_type: 1 left side (start), while 2 right side (end)
//
void addBreakpointInfo(CNV_Array *cnv_array, uint32_t cnv_index, Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t anchor_breakpoint, int pos_type) {
    khint_t k = kh_get(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, anchor_breakpoint);
    if (k != kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash)) {
        // add breakpoint info
        //
        if (pos_type == 1) {
            cnv_array->cnvs[cnv_index].left_breakpoint  = anchor_breakpoint;
            //cnv_array->cnvs[cnv_index].final_start  = anchor_breakpoint;
            cnv_array->cnvs[cnv_index].left_num_of_TLEN_ge_1000 = kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000;
            cnv_array->cnvs[cnv_index].left_num_of_breakpoints  = kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count;
        } else {
            cnv_array->cnvs[cnv_index].right_breakpoint = anchor_breakpoint;
            //cnv_array->cnvs[cnv_index].final_end = anchor_breakpoint;
            cnv_array->cnvs[cnv_index].right_num_of_TLEN_ge_1000 = kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000;
            cnv_array->cnvs[cnv_index].right_num_of_breakpoints  = kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count;
        }
    } else {
        fprintf(stderr, "Warning: something is wrong, the anchor breakpoint %"PRIu32" doesn't exist\n", anchor_breakpoint);
    }

}

void outputCNVArray(CNV_Array *cnv_array, char *chrom_id, int type) {
    char filename[100];
    FILE *fp;
    if (type == 1) {        // For raw varying size bin output
        sprintf(filename, "%s_CNV_Array_raw_bin.txt", chrom_id);
        fp = fopen(filename, "w");
    } else {                // For equal size window output
        sprintf(filename, "%s_CNV_Array_equal_bin.txt", chrom_id);
        fp = fopen(filename, "w");
    }
    fprintf(fp, "chr\tstart\tend\tlength\tCNV-call\tAvg_Cov\tbkpt_L\t#bkpt_L\t#Insert_size_L>=1000\t");
    fprintf(fp, "bkpt_R\t#bkpt_R\t#Insert_Size_R>=1000\tStart_Eq\tEnd_Eq\tLen_Eq\tLen_Eq_Based_used\tavg_cov_Eq\tBin_Type\t");
    fprintf(fp, "Start_Raw\tEnd_Raw\tLen_Raw\tavg_cov_Raw\tBin_Type\n");

    uint32_t j;
    for (j=0; j<cnv_array->size; j++) {
        if (cnv_array->cnvs[j].left_breakpoint > 0 || cnv_array->cnvs[j].right_breakpoint > 0) {
            if (cnv_array->cnvs[j].left_num_of_breakpoints <=2 && cnv_array->cnvs[j].left_num_of_TLEN_ge_1000 <=1 &&
                    cnv_array->cnvs[j].right_num_of_breakpoints == 0 && cnv_array->cnvs[j].right_num_of_TLEN_ge_1000 == 0)
                continue;

            if (cnv_array->cnvs[j].left_num_of_breakpoints == 0 && cnv_array->cnvs[j].left_num_of_TLEN_ge_1000 == 0 &&
                cnv_array->cnvs[j].right_num_of_breakpoints <=2 && cnv_array->cnvs[j].right_num_of_TLEN_ge_1000 <= 1)
                continue;

            uint32_t cnv_start = (cnv_array->cnvs[j].left_breakpoint > 0) ? cnv_array->cnvs[j].left_breakpoint : \
                        (cnv_array->cnvs[j].raw_bin_start > 0) ? cnv_array->cnvs[j].raw_bin_start : cnv_array->cnvs[j].equal_bin_start;
            uint32_t cnv_end = (cnv_array->cnvs[j].right_breakpoint > 0) ? cnv_array->cnvs[j].right_breakpoint : \
                        (cnv_array->cnvs[j].raw_bin_end > 0) ? cnv_array->cnvs[j].raw_bin_end : cnv_array->cnvs[j].equal_bin_end;

            char CNV[10];
            (cnv_array->cnvs[j].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "DUP");

            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%s\t%.2f\t", chrom_id, \
                    cnv_start, cnv_end, cnv_end-cnv_start, CNV, cnv_array->cnvs[j].ave_coverage);

            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", cnv_array->cnvs[j].left_breakpoint, \
                    cnv_array->cnvs[j].left_num_of_breakpoints, cnv_array->cnvs[j].left_num_of_TLEN_ge_1000);

            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", cnv_array->cnvs[j].right_breakpoint, \
                    cnv_array->cnvs[j].right_num_of_breakpoints, cnv_array->cnvs[j].right_num_of_TLEN_ge_1000);

            fprintf(fp, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t", cnv_array->cnvs[j].equal_bin_start, \
                    cnv_array->cnvs[j].equal_bin_end, cnv_array->cnvs[j].equal_bin_end - cnv_array->cnvs[j].equal_bin_start,\
                    cnv_array->cnvs[j].length, cnv_array->cnvs[j].ave_coverage);

            fprintf(fp, "\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array->cnvs[j].raw_bin_start, \
                    cnv_array->cnvs[j].raw_bin_end, cnv_array->cnvs[j].raw_bin_end - cnv_array->cnvs[j].raw_bin_start,\
                    cnv_array->cnvs[j].length, cnv_array->cnvs[j].ave_coverage);

            uint32_t k;
            for (k=0; k<cnv_array->cnvs[j].size; k++) {
                if (cnv_array->cnvs[j].equal_bin_array[k].type == 'E') {
                    fprintf(fp, "\t\t\t\t\t\t\t\t\t\t\t\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%c\n", 
                        cnv_array->cnvs[j].equal_bin_array[k].start, cnv_array->cnvs[j].equal_bin_array[k].end, \
                        cnv_array->cnvs[j].equal_bin_array[k].end - cnv_array->cnvs[j].equal_bin_array[k].start,
                        cnv_array->cnvs[j].equal_bin_array[k].length, cnv_array->cnvs[j].equal_bin_array[k].ave_coverage, \
                        cnv_array->cnvs[j].equal_bin_array[k].type);
                } else {
                    fprintf(fp, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%c\n",
                            cnv_array->cnvs[j].equal_bin_array[k].start, cnv_array->cnvs[j].equal_bin_array[k].end, \
                            cnv_array->cnvs[j].equal_bin_array[k].length, cnv_array->cnvs[j].equal_bin_array[k].ave_coverage, \
                            cnv_array->cnvs[j].equal_bin_array[k].type);
                }
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

void dynamicIncreaseBinArraySize(Equal_Window_Bin **merged_equal_bin_array, uint32_t bin_capacity) {
    if (*merged_equal_bin_array) {
        *merged_equal_bin_array = realloc(*merged_equal_bin_array, bin_capacity * sizeof(Equal_Window_Bin));
        failureExit(merged_equal_bin_array, "Equal_Window_Bin *merged_equal_bin_array memory realloc failed\n");
    } else {
        return;
    }
}

void cnvArrayInit(CNV_Array **cnv_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        cnv_array[i] = calloc(1, sizeof(CNV_Array));
        cnv_array[i]->chromosome_id = strdup(chrom_tracking->chromosome_ids[i]);

        // initialize CNVs on each chromosome
        //
        cnv_array[i]->capacity = 5*PR_INIT_SIZE;
        cnv_array[i]->cnvs = calloc(cnv_array[i]->capacity, sizeof(CNV));
        cnv_array[i]->size = 0;
    }
}

void cnvArrayDestroy(CNV_Array **cnv_array, uint32_t number_of_chromosomes) {
    uint32_t i, j;
    for (i=0; i<number_of_chromosomes; i++) {
        // handle chrom id
        //
        if (cnv_array[i]->chromosome_id) {
            free(cnv_array[i]->chromosome_id);
            cnv_array[i]->chromosome_id = NULL;
        }

        for (j=0; j<cnv_array[i]->size;j++) {
            if (cnv_array[i]->cnvs[j].equal_bin_array) {
                free(cnv_array[i]->cnvs[j].equal_bin_array);
                cnv_array[i]->cnvs[j].equal_bin_array = NULL;
            }
        }

        if (cnv_array[i]->cnvs) {
            free(cnv_array[i]->cnvs);
            cnv_array[i]->cnvs = NULL;
        }
    }
    if (cnv_array)
        free(cnv_array);
}
