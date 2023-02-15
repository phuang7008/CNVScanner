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

void generateCNVs(CNV_Array **equal_bin_cnv_array, Binned_Data_Wrapper **equal_size_window_wrappers, Binned_Data_Wrapper **raw_bin_data_wrappers, Paired_Reads_Across_Breakpoints_Array **preads_x_bpts_array, Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking, Simple_Stats *the_stats, User_Input *user_inputs) {
#pragma omp parallel shared(the_stats) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        uint32_t cnv_array_index;
        for (cnv_array_index=0; cnv_array_index<chrom_tracking->number_of_chromosomes; ++cnv_array_index) {
#pragma omp task
          {
            int thread_id = omp_get_thread_num();
            printf("Current thread id in generating CNV calls: %d\n", thread_id);

            // find the corresponding index in equal_size_window_wrappers, raw_bin_data_wrappers and preads_x_bpts_array
            //
            uint32_t equal_bin_index, raw_bin_index, pr_x_bpts_arr_index, improper_array_index;
            for (equal_bin_index=0; equal_bin_index<chrom_tracking->number_of_chromosomes; equal_bin_index++) {
                if (strcmp(equal_size_window_wrappers[equal_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (raw_bin_index=0; raw_bin_index<chrom_tracking->number_of_chromosomes; raw_bin_index++) {
                if (strcmp(raw_bin_data_wrappers[raw_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (pr_x_bpts_arr_index=0; pr_x_bpts_arr_index<chrom_tracking->number_of_chromosomes; pr_x_bpts_arr_index++) {
                if (strcmp(preads_x_bpts_array[pr_x_bpts_arr_index]->chrom_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (improper_array_index=0; improper_array_index<chrom_tracking->number_of_chromosomes; improper_array_index++) {
                if (strcmp(improperly_PR_array[improper_array_index]->chrom_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            mergeNeighboringBinsBasedOnZscore(equal_bin_cnv_array[cnv_array_index], equal_size_window_wrappers[equal_bin_index], the_stats, equal_bin_cnv_array[cnv_array_index]->chromosome_id, user_inputs, 2);

            expandMergedCNVWithRawBins(raw_bin_data_wrappers[raw_bin_index], equal_bin_cnv_array[cnv_array_index], the_stats);

            checkBreakpointForEachCNV(equal_bin_cnv_array[cnv_array_index], preads_x_bpts_array[pr_x_bpts_arr_index]);

            setLeftRightCNVBreakpoints(equal_bin_cnv_array[cnv_array_index]);

            checkImproperlyPairedReadsForEachCNV(equal_bin_cnv_array[cnv_array_index], improperly_PR_array[improper_array_index]);

            outputCNVArray(equal_bin_cnv_array[cnv_array_index], equal_bin_cnv_array[cnv_array_index]->chromosome_id, 2);
            outputCNVArray(equal_bin_cnv_array[cnv_array_index], equal_bin_cnv_array[cnv_array_index]->chromosome_id, 3);

          } // end task 
        } // end for loop
#pragma omp taskwait

      } // end single
    } // end parallel

    // produce vcf file
    //
    FILE *fh = fopen(user_inputs->vcf_output_file, "w");
    generateVCF_MetaData(user_inputs, chrom_tracking, fh);
    generateVCFresults(equal_bin_cnv_array, chrom_tracking, fh);
    fclose(fh);
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
        if ( equal_size_window_wrapper->data[j].start == 44338500) {
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
                    fprintf(stderr, "WARNING: failed to elongate CNV using equal bins");
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
                        if (cnv_counter + 3 >= cnv_array->capacity) {
                            cnv_array->capacity += PR_INIT_SIZE * 2;
                            cnv_array->cnvs = realloc(cnv_array->cnvs, cnv_array->capacity * sizeof(CNV));
                            failureExit(cnv_array->cnvs, "mergeNeighboringBinsBasedOnZscore cnv_array->cnvs memory realloc failed\n");
                        }
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
    //if (cnv_array->cnvs[cnv_index].equal_bin_start == 20713500)
    //    printf("debug stopped\n");

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

void storeCurrentCNVtoArray(CNV_Array *cnv_array, uint32_t start, uint32_t end, uint32_t length, double coverage, Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_size, uint32_t cnv_index, uint8_t cnv_flag) {
    uint32_t k;

    cnv_array->cnvs[cnv_index].equal_bin_start  = start;
    cnv_array->cnvs[cnv_index].equal_bin_end    = end; 
    cnv_array->cnvs[cnv_index].length = length;
    cnv_array->cnvs[cnv_index].ave_coverage = coverage;

    if (cnv_flag == 1) {
        cnv_array->cnvs[cnv_index].cnv_type = 'L';      // for DEL
    } else if (cnv_flag == 3) {
        cnv_array->cnvs[cnv_index].cnv_type = 'P';      // for DUP
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

    // set the cnv_breakpoints to NULL, and initialization other varilabes here
    // The cnv_breakpoints_size, left_start_index, right_end_index MUST be set here even though the cnv_breakpoints == NULL
    //
    cnv_array->cnvs[cnv_index].cnv_breakpoints = NULL;
    cnv_array->cnvs[cnv_index].cnv_breakpoints_size = 0;
    cnv_array->cnvs[cnv_index].left_start_index = -1;   // 0-index is valid, so need to use -1 for it
    cnv_array->cnvs[cnv_index].right_end_index  = -1;   // 0-index is valid, so need to use -1 for it
}

// The method will first check to see if the previous CNV should be merged with the current CNV 
// if their distance is <= 1000 bp away; if not, keep everything as is
//
int combineNeighboringCNVs(CNV_Array *cnv_array, uint32_t cnv_index) {
    if (cnv_index > 0) {
        if (cnv_array->cnvs[cnv_index].equal_bin_start <= cnv_array->cnvs[cnv_index-1].equal_bin_end + 1000 && 
                cnv_array->cnvs[cnv_index-1].equal_bin_end + 1000 <= cnv_array->cnvs[cnv_index].equal_bin_end &&
                cnv_array->cnvs[cnv_index-1].cnv_type == cnv_array->cnvs[cnv_index].cnv_type) {
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
                //fprintf(stderr, "i value is %"PRIu32" while k value is %"PRIu32"\n", i, k);
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].start  = cnv_array->cnvs[cnv_index].equal_bin_array[i].start;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].end    = cnv_array->cnvs[cnv_index].equal_bin_array[i].end;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].length = cnv_array->cnvs[cnv_index].equal_bin_array[i].length;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].type   = cnv_array->cnvs[cnv_index].equal_bin_array[i].type;
                cnv_array->cnvs[cnv_index-1].equal_bin_array[k].ave_coverage = \
                                                cnv_array->cnvs[cnv_index].equal_bin_array[i].ave_coverage;
            }

            // need to clean-up current CNV
            //
            if (cnv_array->cnvs[cnv_index].equal_bin_array) {
                free(cnv_array->cnvs[cnv_index].equal_bin_array);
                cnv_array->cnvs[cnv_index].equal_bin_array = NULL;
            }
                
            cnv_array->cnvs[cnv_index].equal_bin_start = 0;
            cnv_array->cnvs[cnv_index].equal_bin_end = 0;
            cnv_array->cnvs[cnv_index].raw_bin_start = 0;
            cnv_array->cnvs[cnv_index].raw_bin_end = 0;
            cnv_array->cnvs[cnv_index].length = 0;
            cnv_array->cnvs[cnv_index].ave_coverage = 0;
            cnv_array->cnvs[cnv_index].cnv_type = 'D';
            cnv_array->cnvs[cnv_index].size = 0;

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
        if (cnv_array->cnvs[i].equal_bin_start == 92160000) {
            printf("stop3\n");
        }
        for (j=raw_bin_start; j<binned_data_wrapper->size; j++) {
            // Note the code will handle the extension of extra raw-bins later on
            //
            if (binned_data_wrapper->data[j].end < cnv_array->cnvs[i].equal_bin_start) { 
                // no intersect, just skip
                // ------=----------------------------------------=---------
                //  Raw-bin-end               <               CNV-start
                continue;
            } else if (cnv_array->cnvs[i].equal_bin_end  < binned_data_wrapper->data[j].start) {
                // no intersect. The current raw bin has pass the current CNV bin,  just break out
                // ------=----------------------------------------=---------
                //    CNV-end                 <              Raw_bin_start
                //
                break;
            } else {
                // they intersect
                //
                raw_bin_start = j;

                // Handle left hand side raw bin for the current CNV
                // 
                if (binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_start &&
                        cnv_array->cnvs[i].equal_bin_start <= binned_data_wrapper->data[j].end) {
                    // ------=-----------------=------------------------=---------------
                    //   Raw-bin-start     CNV-start               Raw-bin-end
                    //
                    // next, need to check the average coverage between overlapped raw-bin and equal-bin's CNV
                    //
                    if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff && 
                                cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                        (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff && 
                                cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) { 
                        cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                    } else {
                        // intersected Raw-bin
                        //      ---------------- 25X
                        //              ------------------------------------ equal-bin CNV 15X
                        //                      ============================ adjusted CNV  15X
                        // the same schema is true for the INS
                        //
                        if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized > hap_cutoff && 
                                cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                            (binned_data_wrapper->data[j].ave_cov_map_gc_normalized < dup_cutoff &&
                                cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                            cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].end + 1;
                            addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                            j++;
                        }
                    }
            
                    // the while loop will handle all the left raw bins that don't intersect
                    // to extend the current CNV on the left-hand side
                    //
                    j--;
                    while ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                           (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff 
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        // raw-bin at j    raw bin at j+1 (intersected)
                        // -------------- -----------------------------
                        //          CNV-start ------------------------------------- CNV-end
                        // ======================================================== extended CNV
                        //
                        cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                        j--;
                    }
                    j = raw_bin_start+1;
                    continue;
                }

                // Now handle right hand side, 
                // the current j will be the raw bin index that is the closest to the current CNV
                // first, check if they intersect
                //
                if (binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_end &&
                       cnv_array->cnvs[i].equal_bin_end  <= binned_data_wrapper->data[j].end) {
                    // ------=----------------=-----------------=---------
                    // Raw-bin-start    CNV-bin-end       Raw-bin-end
                    //
                    if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff &&
                                cnv_array->cnvs[i].ave_coverage <= hap_cutoff) ||
                        (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff &&
                                cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].end;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                    } else {
                        if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized > hap_cutoff &&
                                    cnv_array->cnvs[i].ave_coverage <= hap_cutoff) ||
                            (binned_data_wrapper->data[j].ave_cov_map_gc_normalized < dup_cutoff &&
                                    cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                            // ------=----------------=-----------------=------------------
                            // --------------------------------- CNV-bin-end (15X)
                            //      Raw-bin-start ------------------------- Raw-bin-end (33X)
                            // ================== new CNV-bin-end
                            // The same is true for INS
                            //
                            cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].start - 1;
                            addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                            j--;
                        }
                    }

                    // now let's extend the CNV further at the right end
                    //
                    j++;
                    while ((binned_data_wrapper->data[j].ave_coverage <= hap_cutoff &&
                                cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                           (binned_data_wrapper->data[j].ave_coverage >= dup_cutoff &&
                                cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        //      raw-bin at j-1 (intersected)                  raw bin at j
                        //      --------------------------------------------  -----------------
                        //  CNV-start 
                        //  ...------------------------------------- CNV-end
                        //  ...================================================================ extended CNV
                        //
                        cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].end;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                        j++;
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
    uint32_t capacity = PR_INIT_SIZE*50;    // the initial size is set to 10000
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
                capacity += PR_INIT_SIZE*20;
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
        if (i == 1100) {
            printf("stop\n");
        }
        uint32_t tmp_start = (cnv_array->cnvs[i].raw_bin_start>0 && \
                                        cnv_array->cnvs[i].raw_bin_start < cnv_array->cnvs[i].equal_bin_start) ? \
                                        cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
        uint32_t tmp_end = (cnv_array->cnvs[i].raw_bin_end > cnv_array->cnvs[i].equal_bin_end) ? \
                                        cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;

        for (j=breakpoint_start; j<counter; j++) {
            if (anchor_breakpoints[j] + DISTANCE_CUTOFF < tmp_start) {
                continue;
            } else if (tmp_end + DISTANCE_CUTOFF < anchor_breakpoints[j]) {
                break;
            } else {
                // we consider they intersect.
                // Now we need to find out the breakpoint on the current CNV
                //
                addBreakpointInfo(cnv_array, i, preads_x_bpt_arr, anchor_breakpoints[j]);
            }
        }
    }

    if (anchor_breakpoints) { free(anchor_breakpoints); }
}

// The breakpoints will be stored in an array
//
void addBreakpointInfo(CNV_Array *cnv_array, uint32_t cnv_index, Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t anchor_breakpoint) {
    khint_t k = kh_get(khIntPrArray, preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, anchor_breakpoint);
    if (k != kh_end(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash)) {
        // add breakpoint info
        //
        if (cnv_array->cnvs[cnv_index].cnv_breakpoints == NULL) {
            cnv_array->cnvs[cnv_index].cnv_breakpoints_capacity = 25;
            cnv_array->cnvs[cnv_index].cnv_breakpoints =
                            calloc(cnv_array->cnvs[cnv_index].cnv_breakpoints_capacity, sizeof(CNV_Breakpints));
            cnv_array->cnvs[cnv_index].cnv_breakpoints_size = 0;
            cnv_array->cnvs[cnv_index].left_start_index = -1;   // 0-index is valid, so need to use -1 for it
            cnv_array->cnvs[cnv_index].right_end_index  = -1;   // 0-index is valid, so need to use -1 for it
        }

        uint16_t index = cnv_array->cnvs[cnv_index].cnv_breakpoints_size;
        //printf("start position: %"PRIu32" with size: %"PRIu16" at anchor: %"PRIu32"\n", cnv_array->cnvs[cnv_index].equal_bin_start, index, anchor_breakpoint);
        //if (cnv_array->cnvs[cnv_index].equal_bin_start == 44338500) {
        //    printf("stop\n");
        //}
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].breakpoint = anchor_breakpoint;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].num_of_TLEN_ge_1000 = 
                            kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->num_TLEN_ge_1000;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].num_of_breakpoints =
                            kh_value(preads_x_bpt_arr->preads_x_per_anchor_bpt_hash, k)->current_breakpoint_count;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].orientation = 0;
        cnv_array->cnvs[cnv_index].cnv_breakpoints_size++;

        // resize the cnv_breakpoints array if needed
        //
        if ((cnv_array->cnvs[cnv_index].cnv_breakpoints_size + 3) > cnv_array->cnvs[cnv_index].cnv_breakpoints_capacity) {
            cnv_array->cnvs[cnv_index].cnv_breakpoints_capacity += 25;
            cnv_array->cnvs[cnv_index].cnv_breakpoints = realloc(cnv_array->cnvs[cnv_index].cnv_breakpoints, \
                            (cnv_array->cnvs[cnv_index].cnv_breakpoints_capacity)*sizeof(CNV_Breakpints));
            failureExit(cnv_array->cnvs[cnv_index].cnv_breakpoints, "cnv_array->cnvs[cnv_index].cnv_breakpoints");
        }
    }
}

void setLeftRightCNVBreakpoints(CNV_Array *cnv_array) {
    // When using breakpoints to set starts (left-end) and ends (right-end) of CNVs, we need the following
    // 1). number of breakpoints associated with the current anchor breakpoint >= 2
    // 2). num_of_TLEN_ge_1000 >= 1
    // 3). left and right breakpoints should be separated by >= 1000
    //
    //uint16_t num_left_breakpoints = 0, num_right_breakpoints = 0;

    uint32_t i;
    for (i=0; i<cnv_array->size; i++) {
        /*if (cnv_array->cnvs[i].equal_bin_start == 18613500) {
            printf("stop 6\n");
            printf("%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array->cnvs[i].equal_bin_start, \
                    cnv_array->cnvs[i].equal_bin_end, cnv_array->cnvs[i].equal_bin_end - cnv_array->cnvs[i].equal_bin_start,\
                    cnv_array->cnvs[i].length, cnv_array->cnvs[i].ave_coverage);

            printf("\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array->cnvs[i].raw_bin_start, \
                    cnv_array->cnvs[i].raw_bin_end, cnv_array->cnvs[i].raw_bin_end - cnv_array->cnvs[i].raw_bin_start,\
                    cnv_array->cnvs[i].length, cnv_array->cnvs[i].ave_coverage);
            printf("cnv_breakpoints_size is %"PRIu8"\n", cnv_array->cnvs[i].cnv_breakpoints_size);
        }*/

        // some of the CNVs don't have breakpoints associated them, so skip
        //
        if (cnv_array->cnvs[i].cnv_breakpoints == NULL)
            continue;

        bool left_most_set = false, right_most_set = false;
        uint32_t j;
        for (j=0; j<cnv_array->cnvs[i].cnv_breakpoints_size; j++) {
            /*printf("breakpoint: %"PRIu32"\tnum_of_breakpoints: %"PRIu8"\tnum_of_TLEN_ge_1000: %"PRIu8"\n", \
                    cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint,
                    cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints, 
                    cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000);
            */
            // All anchor breakpoints associated with the current CNV should be ordered from the left most to the right most.
            //
            if (cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints <= 3)          // singleton
                continue;
            if (cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000 <= 2)         // small indel
                continue;

            // pass the first 2 criteria
            // 
            uint32_t cur_start = ((cnv_array->cnvs[i].raw_bin_start > 0 && 
                                cnv_array->cnvs[i].raw_bin_start < cnv_array->cnvs[i].equal_bin_start)) ? \
                                        cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
            uint32_t cur_end = (cnv_array->cnvs[i].raw_bin_end > cnv_array->cnvs[i].equal_bin_end) ? \
                                        cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;

            //printf("%"PRIu32"\t%"PRIu32"\n", cur_start, cur_end);
            // because values are uint32_t, so the subtraction will cause overflow as it won't be negative
            //
            if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_start)) < \
                    abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_end))) {
                // the current anchor breakpoint closer to the left-hand side
                //
                if (!left_most_set) {
                    cnv_array->cnvs[i].cnv_breakpoints[j].orientation = 1;
                    cnv_array->cnvs[i].left_start_index = j;
                    left_most_set = true;
                } else {
                    // belongs to the right-hand breakpoint
                    //
                    if (!right_most_set) {
                        if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].left_start_index].breakpoint - \
                                    cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint)) >= 1000) {
                            cnv_array->cnvs[i].cnv_breakpoints[j].orientation = 2;
                            cnv_array->cnvs[i].right_end_index = j;
                            right_most_set = true;
                        }
                    } else {
                        // remove the previous right-hand breakpoint
                        // TODO: Here we might want to check which one is better right breakpoint
                        //
                        cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].right_end_index].orientation = 0;

                        // Now set the new right hand breakpoint
                        //
                        cnv_array->cnvs[i].cnv_breakpoints[j].orientation = 2;
                        cnv_array->cnvs[i].right_end_index = j;
                    }
                }
            } else {
                // the current anchor breakpoint closer to the right-hand side
                //
                if (!right_most_set) {
                    cnv_array->cnvs[i].cnv_breakpoints[j].orientation = 2;
                    cnv_array->cnvs[i].right_end_index = j;
                    right_most_set = true;
                } else {
                    if (!left_most_set) {
                        // if both breakpoints are closer to the right-hand side
                        //
                        if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].right_end_index].breakpoint - \
                                    cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint)) >= 1000) {
                            // set previous right hand breakpoint as the left hand breakpoint
                            //
                            cnv_array->cnvs[i].left_start_index = cnv_array->cnvs[i].right_end_index;
                            cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].left_start_index].orientation = 1;
                            left_most_set = true;
                        }

                        cnv_array->cnvs[i].right_end_index = j;
                        cnv_array->cnvs[i].cnv_breakpoints[j].orientation = 2;
                    }
                }
            }
        }
    }
}

void checkImproperlyPairedReadsForEachCNV(CNV_Array *cnv_array, Not_Properly_Paired_Reads_Array *improperly_PR_array) {
    // setup an array with all starts and ends from improperly paired reads and all CNVs
    //
    uint32_t capacity = (improperly_PR_array->num_of_groups + cnv_array->size) * 2;
    uint32_t * all_starts_ends = calloc(capacity, sizeof(uint32_t));
    int32_t count = 0;      // need to use signed for negative value checking

    khash_t(m32) *imp_PR_start_hash  = kh_init(m32);
    khash_t(m32) *imp_PR_end_hash    = kh_init(m32);
    khash_t(m32) *imp_PR_end_start_lookup = kh_init(m32);
    khash_t(m32) *imp_PR_start_end_lookup = kh_init(m32);

    // loop through the improperly paired reads and find all starts and ends
    //
    int32_t g;
    for (g=0; g<improperly_PR_array->num_of_groups; g++) {
        uint32_t num_TLEN_1000 = improperly_PR_array->grouped_improperly_PRs[g].num_of_pairs_TLEN_ge_1000;

        if (num_TLEN_1000 >= 2 && (improperly_PR_array->grouped_improperly_PRs[g].group_mate_end - 
                improperly_PR_array->grouped_improperly_PRs[g].group_start <= 100000)) {
            // TODO Once I add the condition using perfect matched paired reads, I will remove the 2nd condition checking
            // Store TLEN count info as value so that it is easy for lookup
            //
            all_starts_ends[count] = improperly_PR_array->grouped_improperly_PRs[g].group_start;
            addValueToKhashBucket32(imp_PR_start_hash, all_starts_ends[count], num_TLEN_1000);
            count++;

            all_starts_ends[count] = improperly_PR_array->grouped_improperly_PRs[g].group_mate_end;
            addValueToKhashBucket32(imp_PR_end_hash, all_starts_ends[count], num_TLEN_1000);
            count++;
            fprintf(stderr, "IMP\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", all_starts_ends[count-2], all_starts_ends[count-1], num_TLEN_1000);

            // add values to the lookup table
            //
            addValueToKhashBucket32(imp_PR_end_start_lookup, improperly_PR_array->grouped_improperly_PRs[g].group_mate_end, improperly_PR_array->grouped_improperly_PRs[g].group_start);
            addValueToKhashBucket32(imp_PR_start_end_lookup, improperly_PR_array->grouped_improperly_PRs[g].group_start, improperly_PR_array->grouped_improperly_PRs[g].group_mate_end);
        }
    }

    // walk through CNV array and get all its starts and ends
    //
    khash_t(m32) *cnv_start_hash  = kh_init(m32);
    khash_t(m32) *cnv_end_hash    = kh_init(m32);

    uint32_t i;
    for (i=0; i<cnv_array->size;i++) {
        // Store the CNV_Array index into the hash tables so that it is easy to lookup and set values
        //
        uint32_t tmp_start = cnv_array->cnvs[i].raw_bin_start>0 ?
                             cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
        addValueToKhashBucket32(cnv_start_hash, tmp_start, i);
        all_starts_ends[count] = tmp_start;
        count++;

        uint32_t tmp_end = cnv_array->cnvs[i].raw_bin_end > 0 ?
                           cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;
        addValueToKhashBucket32(cnv_end_hash, tmp_end, i);
        all_starts_ends[count] = tmp_end;
        count++;
        fprintf(stderr, "CNV\t%"PRIu32"\t%"PRIu32"\n", all_starts_ends[count-2], all_starts_ends[count-1]);

        // need to set all placeholders to 0 for the improperly paired reads
        //
        cnv_array->cnvs[i].imp_PR_start = 0;
        cnv_array->cnvs[i].imp_PR_end   = 0;
        cnv_array->cnvs[i].num_of_imp_RP_TLEN_1000 = 0;
    }

    // Here the capacity won't be equal to the count value as the latter one is filtered, 
    // not all improperly paired reads will be added into the array
    //
    fprintf(stderr, "The capacity is %"PRId32" while the count is %"PRId32" in the checkImproperlyPairedReadsForEachCNV\n\n", capacity, count);

    // now we need to resize the all_starts_ends array to the size of the 'count'
    //
    capacity = count;
    all_starts_ends = realloc(all_starts_ends, capacity * sizeof(uint32_t));

    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    // for debugging purpose
    //
    //for (g=0; g<count;g++)
    //    fprintf(stderr, "pos %"PRIu32"\n", all_starts_ends[g]);

    // now do the intersect. 
    // If they intersect, set the improperly paired reads support to true
    //
    int32_t cnv_index = -1;         // need to use signed value as sometimes, no value found
    int32_t imp_PR_start = -1;      // singed value for testing
    count = 0;
    for (i=0; i<capacity; i++) {
        if (all_starts_ends[i] == 92160217) {
            printf("I am here\n");
        }
        if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) ||
                checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i])) {

            // always decrease count if it is the end position
            //
            count--;

            if (count < 0) {
                fprintf(stderr, "Error: the count %"PRId16" should NOT be negative", count);
                fprintf(stderr, "current pos: %"PRIu32" with prev pos %"PRIu32"\n", all_starts_ends[i], all_starts_ends[i-1]);
                exit(EXIT_FAILURE);
            }

            // when current position is an end and count == 1, there is an intersect 
            // there might be multiple groups of improperly paired reads
            // need to find the one with most number of TLEN >= 1000
            //
            if (count == 1) {
                uint32_t cur_TLEN = getValueFromKhash32(imp_PR_start_hash, imp_PR_start);
                if (cnv_array->cnvs[cnv_index].num_of_imp_RP_TLEN_1000 < cur_TLEN) {
                    cnv_array->cnvs[cnv_index].imp_PR_start = imp_PR_start;
                    cnv_array->cnvs[cnv_index].imp_PR_end =
                        getValueFromKhash32(imp_PR_start_end_lookup, imp_PR_start);;
                    cnv_array->cnvs[cnv_index].num_of_imp_RP_TLEN_1000 = cur_TLEN;
                }
            }

            // it is possible that some starts and ends are the same value
            // we need to delete the end position as it goes first in both this loop and in the data
            //
            khiter_t iter;
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, cnv_end_hash, all_starts_ends[i]);
                if (iter != kh_end(cnv_end_hash)) {
                    kh_del(m32, cnv_end_hash, iter);
                }
            } else if (checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, imp_PR_end_hash, all_starts_ends[i]);
                if (iter != kh_end(imp_PR_end_hash)) {
                    kh_del(m32, imp_PR_end_hash, iter);        // this deletes the key second
                }
            }
        } else {
            // it should be in the 'start' position, always increment count
            //
            count++;

            // get current CNV start positioin
            //
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                cnv_index = getSignedValueFromKhash32(cnv_start_hash, all_starts_ends[i]);
                if (cnv_index == -1) {
                    fprintf(stderr, "Something went wrong with CNV start at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }
            }

            // get improperly paired-read start position
            //
            if (checkm32KhashKey(imp_PR_start_hash, all_starts_ends[i]))
                imp_PR_start = all_starts_ends[i];
        }
    }

    // clean-ups
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, imp_PR_end_start_lookup);
    kh_destroy(m32, imp_PR_start_end_lookup);
    kh_destroy(m32, imp_PR_start_hash);
    kh_destroy(m32, imp_PR_end_hash);
    kh_destroy(m32, cnv_start_hash);
    kh_destroy(m32, cnv_end_hash);
}

void outputCNVArray(CNV_Array *cnv_array, char *chrom_id, int type) {
    char filename[200];
    FILE *fp;
    if (type == 1) {                // For raw varying size bin output
        sprintf(filename, "%s_CNV_Array_raw_bin.txt", chrom_id);
        fp = fopen(filename, "w");
    } else if (type == 2) {         // For filtered merged CNVs from equal size window output
        sprintf(filename, "%s_CNV_Array_final_filtered.txt", chrom_id);
        //sprintf(filename, "%s_CNV_Array_equal_bin_filtered.txt", chrom_id);
        fp = fopen(filename, "w");
    } else {                        // For unfiltered merged CNVs from equal size window output
        sprintf(filename, "%s_CNV_Array_final_all.txt", chrom_id);
        //sprintf(filename, "%s_CNV_Array_equal_bin_all.txt", chrom_id);
        fp = fopen(filename, "w");
    }
    fileOpenError(fp, filename);

    fprintf(fp, "chr\tstart\tend\tlength\tCNV-call\tAvg_Cov\tbkpt_L\t#bkpt_L\t#Insert_size_L>=1000\t");
    fprintf(fp, "bkpt_R\t#bkpt_R\t#Insert_Size_R>=1000\tImp_PR_start\tImp_PR_end\tImp_PR_TLEN>=1000\t");
    fprintf(fp, "Start_Eq\tEnd_Eq\tLen_Eq\tLen_Eq_Based_used\tavg_cov_Eq\tBin_Type\t");
    fprintf(fp, "Start_Raw\tEnd_Raw\tLen_Raw\tavg_cov_Raw\tBin_Type\n");

    uint32_t j;
    for (j=0; j<cnv_array->size; j++) {
        int16_t left_idx  = cnv_array->cnvs[j].left_start_index;
        int16_t right_idx = cnv_array->cnvs[j].right_end_index;

        uint32_t left_breakpoint=0, left_num_bpoint=0, left_num_geTLEN=0;
        uint32_t right_breakpoint=0, right_num_bpoint=0, right_num_geTLEN=0;

        if (left_idx >= 0 || right_idx >= 0 || type == 3 || cnv_array->cnvs[j].num_of_imp_RP_TLEN_1000 >= 2) {
            if (cnv_array->cnvs[j].cnv_breakpoints != NULL) {
                left_breakpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].breakpoint : 0;
                left_num_bpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_breakpoints : 0;
                left_num_geTLEN  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_TLEN_ge_1000 : 0;

                right_breakpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].breakpoint : 0;
                right_num_bpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_breakpoints : 0;
                right_num_geTLEN = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_TLEN_ge_1000 : 0;
            }

            if (type != 3) {
                // filter CNVs based on the breakpoint info and improperly paired reads info
                // based need at least 2 support for all categories
                //
                if (cnv_array->cnvs[j].num_of_imp_RP_TLEN_1000 <= 1) {
                    if (left_idx >= 0 && right_idx < 0 && left_num_bpoint <=1 && left_num_geTLEN <=1)
                        continue;

                    if (left_idx < 0 && right_idx >= 0 && right_num_bpoint <=1 && right_num_geTLEN <= 1)
                        continue;
                }
            }

            /*if (cnv_array->cnvs[j].equal_bin_start == 18613500) {
                printf("stop 7\n");
            }*/
            uint32_t cnv_start = (left_idx >= 0) ? left_breakpoint : (cnv_array->cnvs[j].raw_bin_start > 0) ? \
                                    cnv_array->cnvs[j].raw_bin_start : cnv_array->cnvs[j].equal_bin_start;
            uint32_t cnv_end = (right_idx >= 0) ? right_breakpoint : (cnv_array->cnvs[j].raw_bin_end > 0) ? \
                                    cnv_array->cnvs[j].raw_bin_end : cnv_array->cnvs[j].equal_bin_end;

            char CNV[10];
            (cnv_array->cnvs[j].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "INS");

            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%s\t%.2f\t", chrom_id, \
                    cnv_start, cnv_end, cnv_end-cnv_start, CNV, cnv_array->cnvs[j].ave_coverage);


            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", left_breakpoint, left_num_bpoint, left_num_geTLEN);

            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", right_breakpoint, right_num_bpoint, right_num_geTLEN);

            fprintf(fp, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t", cnv_array->cnvs[j].imp_PR_start, cnv_array->cnvs[j].imp_PR_end, cnv_array->cnvs[j].num_of_imp_RP_TLEN_1000);

            fprintf(fp, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t", cnv_array->cnvs[j].equal_bin_start, \
                    cnv_array->cnvs[j].equal_bin_end, cnv_array->cnvs[j].equal_bin_end - cnv_array->cnvs[j].equal_bin_start,\
                    cnv_array->cnvs[j].length, cnv_array->cnvs[j].ave_coverage);

            fprintf(fp, "\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array->cnvs[j].raw_bin_start, \
                    cnv_array->cnvs[j].raw_bin_end, abs((signed)(cnv_array->cnvs[j].raw_bin_end - cnv_array->cnvs[j].raw_bin_start)),\
                    cnv_array->cnvs[j].length, cnv_array->cnvs[j].ave_coverage);

            if (type != 3) {
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
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

void generateVCF_MetaData(User_Input *user_inputs, Chromosome_Tracking *chrom_tracking, FILE *fh) {
    fprintf(fh, "##fileformat=VCFv4.2\n");
    fprintf(fh, "##source=%s\n", SOURCE_);
    fprintf(fh, "##reference=%s\n", user_inputs->reference_file);

    // for chromosome ids
    //
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        fprintf(fh, "##contig=<ID=%s,length=%"PRIu32">\n", chrom_tracking->chromosome_ids[i], chrom_tracking->chromosome_lengths[i]);
    }

    fprintf(fh, "##ALT=<ID=CNV,Description=\"Copy Number Variant\">\n");
    fprintf(fh, "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n");
    fprintf(fh, "##ALT=<ID=DUP,Description=\"Insertion/Duplication relative to the reference\">\n");
    fprintf(fh, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
    fprintf(fh, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV:DEL=Deletion, DUP=Duplication\">\n");
    fprintf(fh, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
    fprintf(fh, "##INFO=<ID=AVGCOV,Number=.,Type=Float,Description=\"Average Coverage of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTL,Number=1,Type=Integer,Description=\"Breakpoint position at the left side of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTR,Number=1,Type=Integer,Description=\"Breakpoint position at the right side of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTLCOUNT,Number=1,Type=Integer,Description=\"Number of breakpoints at the left side of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTRCOUNT,Number=1,Type=Integer,Description=\"Number of breakpoints at the right side of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTLTLEN,Number=1,Type=Integer,Description=\"Number of Reads with Insertion size >= 1000bp across the breakpoints at the end of left side of the CNV\">\n");
    fprintf(fh, "##INFO=<ID=BPTRTLEN,Number=1,Type=Integer,Description=\"Number of Reads with Insertion size >= 1000bp across the breakpoints at the end of right side of the CNV\">\n");
    fprintf(fh, "##FILTER=<ID=noBreakpoint,Description=\"CNV without breakpoint support\">\n");
    fprintf(fh, "##FILTER=<ID=littleBreakpointSupport,Description=\"CNV has breakpoint support, but the support is too little to be useful\">\n");
    fprintf(fh, "##FILTER=<ID=PASS,Description=\"CNV pass the filter\">\n");
    fprintf(fh, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", user_inputs->sample_name);
    //fprintf(fh, "");
}

void generateVCFresults(CNV_Array **equal_bin_cnv_array, Chromosome_Tracking *chrom_tracking, FILE *fp) {
    // walk through chromosome list
    //
    uint32_t i, j, cnv_start, cnv_end;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        CNV_Array *cnv_array = equal_bin_cnv_array[i];      // pointer assignment, don't need to be free-ed
        for (j=0; j<cnv_array->size; j++) {
            int16_t left_idx  = cnv_array->cnvs[j].left_start_index;
            int16_t right_idx = cnv_array->cnvs[j].right_end_index;

            uint32_t left_breakpoint=0, left_num_bpoint=0, left_num_geTLEN=0;
            uint32_t right_breakpoint=0, right_num_bpoint=0, right_num_geTLEN=0;

            if (cnv_array->cnvs[j].cnv_breakpoints != NULL) {
                left_breakpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].breakpoint : 0;
                left_num_bpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_breakpoints : 0;
                left_num_geTLEN  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_TLEN_ge_1000 : 0;

                right_breakpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].breakpoint : 0;
                right_num_bpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_breakpoints : 0;
                right_num_geTLEN = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_TLEN_ge_1000 : 0;
            }

            cnv_start = (left_breakpoint > 0) ? left_breakpoint : (cnv_array->cnvs[j].raw_bin_start > 0) ? \
                            cnv_array->cnvs[j].raw_bin_start : cnv_array->cnvs[j].equal_bin_start;

            cnv_end = (right_breakpoint > 0) ? right_breakpoint : (cnv_array->cnvs[j].raw_bin_end > 0) ? \
                            cnv_array->cnvs[j].raw_bin_end : cnv_array->cnvs[j].equal_bin_end;

            char CNV[10];
            (cnv_array->cnvs[j].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "DUP");

            char FILTER[30];
            strcpy(FILTER, "PASS");

            char GT[10];
            strcpy(GT, "./1");

            if (left_breakpoint > 0 || right_breakpoint > 0) {
                if (left_num_bpoint <=2 && left_num_geTLEN <=1 && right_num_bpoint == 0 && right_num_geTLEN == 0) {
                    strcpy(FILTER, "littleBreakpointSupport");
                    strcpy(GT, "./.");
                }

                if (left_num_bpoint == 0 && right_num_geTLEN == 0 && right_num_bpoint <=2 && right_num_geTLEN <= 1) {
                    strcpy(FILTER, "littleBreakpointSupport");
                    strcpy(GT, "./.");
                }
            } else {
                strcpy(FILTER, "noBreakpoint");;
                strcpy(GT, "./.");
            }

            fprintf(fp, "%s\t%"PRIu32"\t.\tN\t%s\t.\t%s\tEND=%"PRIu32";SVLEN=%"PRIu32";SVTYPE=%s;AVGCOV=%.2f;BPTL=%"PRIu32";BPTLCOUNT=%"PRIu8";BPTLTLEN=%"PRIu8";BPTR=%"PRIu32";BPTRCOUNT=%"PRIu8";BPTRTLEN=%"PRIu8"\tGT\t%s\n", \
                    chrom_tracking->chromosome_ids[i], cnv_start, CNV, FILTER, cnv_end, cnv_end-cnv_start, CNV, cnv_array->cnvs[j].ave_coverage, left_breakpoint, left_num_bpoint, left_num_geTLEN, right_breakpoint, right_num_bpoint, right_num_geTLEN, GT);
        } // end equal_bin_cnv_array
    } // end chromosome list
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

            if (cnv_array[i]->cnvs[j].cnv_breakpoints) {
                free(cnv_array[i]->cnvs[j].cnv_breakpoints);
                cnv_array[i]->cnvs[j].cnv_breakpoints = NULL;
            }
        }

        if (cnv_array[i]->cnvs) {
            free(cnv_array[i]->cnvs);
            cnv_array[i]->cnvs = NULL;
        }

        if (cnv_array[i] != NULL) {
            free(cnv_array[i]);
            cnv_array[i]=NULL;
        }
    }
    if (cnv_array)
        free(cnv_array);
}
