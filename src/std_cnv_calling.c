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

void generateCNVs(CNV_Array **equal_bin_cnv_array, Binned_Data_Wrapper **equal_size_window_wrappers, Binned_Data_Wrapper **raw_bin_data_wrappers, khash_t(m32) **anchor_breakpoints_hash_array, Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking, Simple_Stats *the_stats, Stats_Info *stats_info, User_Input *user_inputs, bam_hdr_t **headers, hts_idx_t **sfh_idx, samFile **sfh) {

#pragma omp parallel shared(the_stats) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        uint32_t cnv_array_index;
        for (cnv_array_index=0; cnv_array_index<chrom_tracking->number_of_chromosomes; ++cnv_array_index) {
#pragma omp task
          {
            int thread_id = omp_get_thread_num();
            printf("Current thread id in generating CNV calls: %d with chr %s\n", thread_id, chrom_tracking->chromosome_ids[cnv_array_index]);

            // find the corresponding index in equal_size_window_wrappers, raw_bin_data_wrappers and preads_x_bpts_array
            //
            uint32_t equal_bin_index, raw_bin_index, improper_array_index;
            for (equal_bin_index=0; equal_bin_index<chrom_tracking->number_of_chromosomes; equal_bin_index++) {
                if (strcmp(equal_size_window_wrappers[equal_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (raw_bin_index=0; raw_bin_index<chrom_tracking->number_of_chromosomes; raw_bin_index++) {
                if (strcmp(raw_bin_data_wrappers[raw_bin_index]->chromosome_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0)
                    break;
            }

            for (improper_array_index=0; improper_array_index<chrom_tracking->number_of_chromosomes; improper_array_index++) {
                if (strcmp(improperly_PR_array[improper_array_index]->chrom_id, equal_bin_cnv_array[cnv_array_index]->chromosome_id) == 0) {
                    // need to process the last group of improperly paired reads
                    //
                    processPairedReadsWithinTheSameGroup(improperly_PR_array[improper_array_index]);
                    break;
                }
            }

            mergeNeighboringBinsBasedOnZscore(equal_bin_cnv_array[cnv_array_index], equal_size_window_wrappers[equal_bin_index], the_stats, equal_bin_cnv_array[cnv_array_index]->chromosome_id, user_inputs, 2);

            expandMergedCNVWithRawBins(raw_bin_data_wrappers[raw_bin_index], equal_bin_cnv_array[cnv_array_index], the_stats, user_inputs);

            checkBreakpointForEachCNV(equal_bin_cnv_array[cnv_array_index], anchor_breakpoints_hash_array[cnv_array_index], headers[thread_id], sfh_idx[thread_id], sfh[thread_id], user_inputs);

            checkImproperlyPairedReadsForEachCNV(equal_bin_cnv_array[cnv_array_index], improperly_PR_array[improper_array_index]);

            processPairedReadsAcrossABreakpointTlenInfo(equal_bin_cnv_array[cnv_array_index]);

            setLeftRightCNVBreakpoints(equal_bin_cnv_array[cnv_array_index], user_inputs);

            fprintf(stderr, "%s\n", equal_bin_cnv_array[cnv_array_index]->chromosome_id);
            cleanupOverlappingCNVs(equal_bin_cnv_array[cnv_array_index], the_stats);

            outputCNVArray(equal_bin_cnv_array[cnv_array_index], equal_bin_cnv_array[cnv_array_index]->chromosome_id, user_inputs, 2);
            outputCNVArray(equal_bin_cnv_array[cnv_array_index], equal_bin_cnv_array[cnv_array_index]->chromosome_id, user_inputs, 3);

          } // end task 
        } // end for loop
#pragma omp taskwait

      } // end single
    } // end parallel

    // produce vcf file
    //
    FILE *fh = fopen(user_inputs->vcf_output_file, "w");
    FILE *smfh = fopen(user_inputs->simple_vcf_output_file, "w");
    generateVCF_MetaData(user_inputs, chrom_tracking, fh);
    generateVCFresults(equal_bin_cnv_array, chrom_tracking, the_stats, stats_info, user_inputs, fh, smfh);
    //outputLog2Ratio(equal_size_window_wrappers, chrom_tracking, user_inputs);
    fclose(fh);
    fclose(smfh);
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
        if ( equal_size_window_wrapper->data[j].start == 83071600) {
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
                if (prev_end - prev_start >= user_inputs->min_cnv_length) {
                    int used_length=500;
                    if (user_inputs->min_cnv_length < 1000)
                        used_length = user_inputs->min_cnv_length;

                    if (prev_len >= used_length || ((double)prev_len/((double)prev_end - (double)prev_start)) >= 0.25) {
                        double coverage = total_coverage / prev_len;
                        // further extend the CNV at both ends by one bin (500bp) each 
                        // if both ends connected with excluded regions (Ns regions, repeatmaskers etc.)
                        //
                        storeCurrentCNVtoArray(cnv_array, prev_start, prev_end, prev_len, \
                                                    coverage, merged_equal_bin_array, bin_counter, cnv_counter, prev_flag);
                        //extendBothEndsByOneBin(cnv_array, equal_size_window_wrapper, cnv_counter, bin_counter, j);
                        cnv_counter += combineNeighboringCNVs(cnv_array, cnv_counter, user_inputs);
                        cnv_counter++;
                        if (cnv_counter + 3 >= cnv_array->capacity) {
                            cnv_array->capacity += PR_INIT_SIZE * 2;
                            if (cnv_array->cnvs == NULL) {
                                cnv_array->cnvs = calloc(cnv_array->capacity, sizeof(CNV));
                            } else {
                                cnv_array->cnvs = realloc(cnv_array->cnvs, cnv_array->capacity * sizeof(CNV));
                                failureExit(cnv_array->cnvs, "mergeNeighboringBinsBasedOnZscore cnv_array->cnvs memory realloc failed\n");
                            }
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
        cnv_array->cnvs[cnv_index].cnv_type = 'P';      // for DUP/INS
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
    cnv_array->cnvs[cnv_index].combined = false;        // initialize to false

    // set the cnv_breakpoints to NULL, and initialization other varilabes here
    // The cnv_breakpoints_size, left_start_index, right_end_index MUST be set here even though the cnv_breakpoints == NULL
    //
    cnv_array->cnvs[cnv_index].cnv_breakpoints = NULL;
    cnv_array->cnvs[cnv_index].cnv_breakpoints_size = 0;
    cnv_array->cnvs[cnv_index].left_start_index = -1;   // 0-index is valid, so need to use -1 for it
    cnv_array->cnvs[cnv_index].right_end_index  = -1;   // 0-index is valid, so need to use -1 for it

    // for improperly paired-reads info
    //
    cnv_array->cnvs[cnv_index].imp_PR_start = 0;
    cnv_array->cnvs[cnv_index].imp_PR_end   = 0;
    cnv_array->cnvs[cnv_index].num_of_imp_PR_TLEN_1000 = 0;

    // for inner_CNV (used in CNV segmentation merging step)
    //
    cnv_array->cnvs[cnv_index].inner_cnv.start = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.end = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.passed = false;
    cnv_array->cnvs[cnv_index].inner_cnv.qual = 0.0;
    cnv_array->cnvs[cnv_index].inner_cnv.ave_coverage = 0.0;
    cnv_array->cnvs[cnv_index].inner_cnv.cnv_type = ' ';
    cnv_array->cnvs[cnv_index].inner_cnv.valid_cnv = false;

    cnv_array->cnvs[cnv_index].inner_cnv.imp_PR_start = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.imp_PR_end = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.num_larger_imp_PR_TLEN = 0;

    cnv_array->cnvs[cnv_index].inner_cnv.left_breakpoint  = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.left_breakpoint_count  = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.num_larger_TLEN_left  = 0;

    cnv_array->cnvs[cnv_index].inner_cnv.right_breakpoint = 0;  // the real left breakpoint should start at 1
    cnv_array->cnvs[cnv_index].inner_cnv.right_breakpoint_count = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.num_larger_TLEN_right = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.last_right_breakpoint_count = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.last_num_larger_TLEN_right = 0;

    cnv_array->cnvs[cnv_index].inner_cnv.breakpoint_size = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.breakpoint_capacity = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.evidence_count = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.num_merged_CNVs = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.low_mapp_length = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.gc_lt25pct_length = 0;
    cnv_array->cnvs[cnv_index].inner_cnv.gc_gt85pct_length = 0;
}

// The method will first check to see if the previous CNV should be merged with the current CNV 
// if their distance is <= 1000 bp away; if not, keep everything as is
//
int combineNeighboringCNVs(CNV_Array *cnv_array, uint32_t cnv_index, User_Input *user_inputs) {
    if (cnv_index > 0) {
        if (cnv_array->cnvs[cnv_index].equal_bin_start <= cnv_array->cnvs[cnv_index-1].equal_bin_end + user_inputs->min_cnv_length && 
                cnv_array->cnvs[cnv_index-1].equal_bin_end + user_inputs->min_cnv_length <= cnv_array->cnvs[cnv_index].equal_bin_end &&
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

            // set the 'combined' flag
            //
            cnv_array->cnvs[cnv_index-1].combined = true;

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
void expandMergedCNVWithRawBins(Binned_Data_Wrapper *binned_data_wrapper, CNV_Array *cnv_array, Simple_Stats *equal_window_stats, User_Input *user_inputs) {
    double hap_cutoff = equal_window_stats->average_coverage - equal_window_stats->zScore;
    double dup_cutoff = equal_window_stats->average_coverage + equal_window_stats->zScore;

    // raw_bin_start is used to prevent looping from the start point
    //
    uint32_t i, j, raw_bin_start=0;

    for (i=0; i<cnv_array->size;i++) {
        if (cnv_array->cnvs[i].equal_bin_end == 16293000) {
            printf("stop3\n");
        }

        for (j=raw_bin_start; j<binned_data_wrapper->size; j++) {
            if (j==39913)
                printf("stop\n");
            // Note the code will handle the extension of extra raw-bins later on
            //
            if (binned_data_wrapper->data[j].end < cnv_array->cnvs[i].equal_bin_start) { 
                // no intersect, just skip
                // ------xxx----------------------------------------==-----------
                //  Raw-bin-end               <               CNV-start
                continue;
            } else if (cnv_array->cnvs[i].equal_bin_end  < binned_data_wrapper->data[j].start) {
                // no intersect. The current raw bin has pass the current CNV bin,  just break out
                // ------====------------------------------------------x---------
                //    CNV-end                 <              Raw_bin_start
                //
                break;
            } else {
                // they intersect
                //
                raw_bin_start = j;

                // Handle left hand side raw bin for the current CNV
                // 
                if ((binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_start &&
                        cnv_array->cnvs[i].equal_bin_start - binned_data_wrapper->data[j].start <= user_inputs->min_cnv_length &&
                        cnv_array->cnvs[i].equal_bin_start <= binned_data_wrapper->data[j].end) ||
                        (cnv_array->cnvs[i].equal_bin_start <= binned_data_wrapper->data[j].start &&
                         binned_data_wrapper->data[j].start - cnv_array->cnvs[i].equal_bin_start <= user_inputs->min_cnv_length &&
                         binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_end)) {
                    // ------x-----------------==------------------------xxx---------------
                    //   Raw-bin-start     CNV-start               Raw-bin-end
                    //
                    // OR
                    // -----==----------------x---------------------xxx--------------------
                    //   CNV-start      Raw-bin-start           Raw-bin-end
                    //
                    // next, need to check the average coverage between overlapped raw-bin and equal-bin's CNV
                    // same CNV type
                    //
                    if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff && 
                                cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                        (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff && 
                                cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) { 
                        cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                        addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                    } else {
                        // different copy number categories
                        //
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
                            // need to check to avoid cnv start > cnv end
                            //
                            if (binned_data_wrapper->data[j].end < cnv_array->cnvs[i].equal_bin_end) {
                                cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].end;
                                addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                                j++;
                            }
                        }
                    }
            
                    // the while loop will handle all the left raw bins that don't intersect
                    // to extend the current CNV on the left-hand side, 
                    // Note, there are many raw bins for a CNV to check as CNV size >= 1000, while raw bin sizes vary.
                    // j-- will handle previous raw bin
                    //
                    j--;
                    while (j>=raw_bin_start && j<binned_data_wrapper->size) {
                        if (((cnv_array->cnvs[i].equal_bin_start > binned_data_wrapper->data[j].end) &&
                            ((cnv_array->cnvs[i].equal_bin_start - binned_data_wrapper->data[j].end) <= user_inputs->min_cnv_length)) ||
                            ((cnv_array->cnvs[i].raw_bin_start > binned_data_wrapper->data[j].end) &&   // raw_bin has assigned to the CNV
                             ((cnv_array->cnvs[i].raw_bin_start - binned_data_wrapper->data[j].end) <= user_inputs->min_cnv_length))) {
                            // raw-bin at j-1       
                            // ------------------
                            //                  |         | the distance is <= user_inputs->min_cnv_length (used to be 1000)
                            //                  CNV-start ------------------------------------- CNV-end
                            //                       -------------------- raw bin at j (intersected)
                            // ================================================================ extended CNV
                            //
                            if ((binned_data_wrapper->data[j].ave_cov_map_gc_normalized <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                                (binned_data_wrapper->data[j].ave_cov_map_gc_normalized >= dup_cutoff 
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                                cnv_array->cnvs[i].raw_bin_start = binned_data_wrapper->data[j].start;
                                addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                                j--;
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                    j = raw_bin_start+1;
                    continue;
                }

                // Now handle right hand side, 
                // the current j will be the raw bin index that is the closest to the current CNV
                // first, check if they intersect
                //
                if ( (binned_data_wrapper->data[j].start <= cnv_array->cnvs[i].equal_bin_end &&
                        cnv_array->cnvs[i].equal_bin_end  <= binned_data_wrapper->data[j].end &&
                        binned_data_wrapper->data[j].end - cnv_array->cnvs[i].equal_bin_end <= user_inputs->min_cnv_length) ||
                        (cnv_array->cnvs[i].equal_bin_start <= binned_data_wrapper->data[j].end &&
                         binned_data_wrapper->data[j].end <= cnv_array->cnvs[i].equal_bin_end &&
                         cnv_array->cnvs[i].equal_bin_end - binned_data_wrapper->data[j].end <= user_inputs->min_cnv_length) ) {
                    // ------x----------------==-----------------xxx---------
                    // Raw-bin-start    CNV-bin-end       Raw-bin-end
                    //
                    // OR
                    // -------=----------------xxx--------------==-----------
                    // CNV-bin-start        Raw-bin-end     CNV-bin-end
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
                            // also need to check so that we don't have cnv_start > cnv_end
                            //
                            if (binned_data_wrapper->data[j].start -1 > cnv_array->cnvs[i].equal_bin_start) {
                                cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].start - 1;
                                addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                                j--;
                            }
                        }
                    }

                    // now let's extend the CNV further at the right end
                    // but when we do it on the right hand-side, we need to ensure the next ones
                    // are within 1000bp (2 equal-bin blocks) to each other
                    //
                    j++;
                    while (j>=raw_bin_start && j<binned_data_wrapper->size) {
                        if (j==39913)
                            printf("stop\n");
                        if (((binned_data_wrapper->data[j].start > cnv_array->cnvs[i].equal_bin_end) &&
                                (binned_data_wrapper->data[j].start - cnv_array->cnvs[i].equal_bin_end <= user_inputs->min_cnv_length)) ||
                             ((binned_data_wrapper->data[j].start > cnv_array->cnvs[i].raw_bin_end) && 
                              (binned_data_wrapper->data[j].start - cnv_array->cnvs[i].raw_bin_end <= user_inputs->min_cnv_length))) {
                            //      raw-bin at j-1 (intersected)                  raw bin at j
                            //      --------------------------------------------  -----------------
                            //  CNV-start 
                            //  ...------------------------------------- CNV-end
                            //  ...================================================================ extended CNV
                            //
                            if ((binned_data_wrapper->data[j].ave_coverage <= hap_cutoff &&
                                    cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                                (binned_data_wrapper->data[j].ave_coverage >= dup_cutoff &&
                                    cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                                cnv_array->cnvs[i].raw_bin_end = binned_data_wrapper->data[j].end;
                                addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                                j++;
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
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

void checkBreakpointForEachCNV(CNV_Array *cnv_array, khash_t(m32) *anchor_breakpoints_hash, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh, User_Input *user_inputs) {
    // first, we need to get the sorted anchor breakpoints
    //
    int32_t capacity = PR_INIT_SIZE*500;       // the initial size is set to 100,000
    uint32_t *all_starts_ends = calloc(capacity, sizeof(uint32_t));

    khash_t(m32) *breakpoint_start_hash = kh_init(m32);
    khash_t(m32) *breakpoint_end_hash   = kh_init(m32);
    khash_t(m32) *breakpoint_start_end_lookup = kh_init(m32);
    khash_t(m32) *breakpoint_end_start_lookup = kh_init(m32);

    int32_t counter=0;

    khint_t k;
    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            // Only process breakpoint with count >= 2 and number of TLEN >=2
            //
            if (kh_value(anchor_breakpoints_hash, k) >= 2) {
                uint32_t breakpoint_pos = kh_key(anchor_breakpoints_hash, k);

                if ((int32_t)breakpoint_pos - user_inputs->breakpoint_distance < 0) {
                    // prevent negative value or value overflow
                    //
                    all_starts_ends[counter] = 1;
                } else {
                    all_starts_ends[counter] = breakpoint_pos - user_inputs->breakpoint_distance;
                }
                setValueToKhashBucket32(breakpoint_start_hash, all_starts_ends[counter], breakpoint_pos);
                counter++;

                if (breakpoint_pos + user_inputs->breakpoint_distance > cnv_array->chrom_length) {
                    all_starts_ends[counter] = cnv_array->chrom_length - 1;
                } else {
                    all_starts_ends[counter] = breakpoint_pos + user_inputs->breakpoint_distance;
                }
                setValueToKhashBucket32(breakpoint_end_hash, all_starts_ends[counter], breakpoint_pos);
                counter++;

                //fprintf(stderr, "Breakpoint\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", all_starts_ends[counter-2], all_starts_ends[counter-1], breakpoint_pos);

                setValueToKhashBucket32(breakpoint_start_end_lookup, all_starts_ends[counter-2], all_starts_ends[counter-1]);
                setValueToKhashBucket32(breakpoint_end_start_lookup, all_starts_ends[counter-1], all_starts_ends[counter-2]);

                // dynamically increase the capacity of anchor breakpoint array
                //
                if (counter + 5 >= capacity) {
                    capacity += PR_INIT_SIZE*25;
                    all_starts_ends = realloc(all_starts_ends, capacity * sizeof(uint32_t));
                    failureExit(all_starts_ends, "all_starts_ends in intersecting CNV w/ breakpoint memory realloc failed\n");
                }
            }
        }
    }

    //fprintf(stderr, "Pass 1 The capacity is %"PRId32" and counter is %"PRId32"\n", capacity, counter);

    // reset the size of anchor_breakpoints array to the size of counter
    // plus the cnv_array->size
    //
    capacity = counter + cnv_array->size*2;
    all_starts_ends = realloc(all_starts_ends, capacity * sizeof(uint32_t));
    failureExit(all_starts_ends, "all_starts_ends in intersecting CNV/breakpoint memory realloc failed\n");

    // walk through the CNV array and save all CNV coordinates
    //
    khash_t(m32) *cnv_start_hash = kh_init(m32);
    khash_t(m32) *cnv_end_hash   = kh_init(m32);
    khash_t(m32) *cnv_start_end_lookup = kh_init(m32);
    khash_t(m32) *cnv_end_start_lookup = kh_init(m32);

    uint32_t i;
    for (i=0; i<cnv_array->size;i++) {
        uint32_t tmp_start = cnv_array->cnvs[i].raw_bin_start>0 ? \
                                        cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
        setValueToKhashBucket32(cnv_start_hash, tmp_start, i);      // store index here
        all_starts_ends[counter] = tmp_start;
        counter++;
        
        uint32_t tmp_end = cnv_array->cnvs[i].raw_bin_end > 0 ? \
                                        cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;
        setValueToKhashBucket32(cnv_end_hash, tmp_end, i);
        all_starts_ends[counter] = tmp_end;
        counter++;

        //fprintf(stderr, "CNV\t%"PRId32"\t%"PRId32"\n", all_starts_ends[counter-2], all_starts_ends[counter-1]);

        setValueToKhashBucket32(cnv_start_end_lookup, tmp_start, tmp_end);
        setValueToKhashBucket32(cnv_end_start_lookup, tmp_end, tmp_start);
    }

    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    //fprintf(stderr, "Pass 2 The capacity is %"PRId32" and counter is %"PRId32"\n", capacity, counter);

    // do intersect
    //
    khash_t(m32) *live_cnv_start_hash = kh_init(m32);
    khash_t(m32) *live_bpt_start_hash = kh_init(m32);

    khash_t(m32) *seen_cnv_start_hash = kh_init(m32);
    khash_t(m32) *seen_bpt_start_hash = kh_init(m32);

    khiter_t iter;
    counter = 0;
    int32_t cnv_index = -1;             // need to use signed value as sometimes, no value found
    for (i=0; i<(uint32_t)capacity; i++) {
        if (all_starts_ends[i] == 31732652 || all_starts_ends[i] == 31731653)
            printf("here it is\n");

        // if a key is present in both start and end position, we have to process start first
        //
        bool start_first = false;
        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_start_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i]) && checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_bpt_start_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]) 
                && !checkm32KhashKey(seen_cnv_start_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_bpt_start_hash, all_starts_ends[i]) )
            start_first = true;

        if (!start_first && (checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) ||
                checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]))) {

            // always decrease count if it is the end position
            //
            counter--;

            if (counter < 0) {
                if (i == 0) {
                    counter = 0;    // this happens when breakpoint_pos > user_inputs->breakpoint_distance
                } else {
                    fprintf(stderr, "Error: chr id is %s and index i is %"PRIu32"\n", cnv_array->chromosome_id, i);
                    fprintf(stderr, "Error: the counter %"PRId16" should NOT be negative in breakpoint/cnv intersection\n", counter);
                    fprintf(stderr, "current pos: %"PRIu32" and the prev pos %"PRIu32"\n", all_starts_ends[i], all_starts_ends[i-1]);
                    exit(EXIT_FAILURE);
                }
            }

            /*
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "%s\tCNV_End\t%"PRIu32"\t%"PRId32"\n", cnv_array->chromosome_id, all_starts_ends[i], counter);
            } else {
                fprintf(stderr, "%s\tBreakpoint_End\t%"PRIu32"\t%"PRId32"\n", cnv_array->chromosome_id, all_starts_ends[i], counter);
            }*/

            // when current position is an end and counter >= 1, there is an intersect
            // There might be multiple breakpoints associated with this CNV 
            // and the breakpoint intervals are not Merged. Will store all of them
            //
            if (counter >= 1) {
                uint32_t cnv_start = 0;

                // because the max length of breakpoint interval used for intersection is only 600
                // so we don't have to worry that breakpoint interval completely engulf a CNV
                //
                uint32_t cur_anchor_breakpoint = 0;
                if (checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i])) {
                    // find the anchor breakpoint using the breakpoint end hash with the current end
                    //              ================================= CNV
                    //   ----------x---------- breakpoint interval 1
                    //        -----------x---------- breakpoint interval 2
                    //   bs1  bs2   cs3      be1   be2              ce3
                    //    1    2     3        2     1                0
                    //
                    // Using the current end position to get the anchor breakpoint
                    //
                    cur_anchor_breakpoint = getValueFromKhash32(breakpoint_end_hash, all_starts_ends[i]);

                    // loop through current live_cnv_start_hash
                    //
                    for (iter=kh_begin(live_cnv_start_hash); iter!=kh_end(live_cnv_start_hash); ++iter) {
                        if (kh_exist(live_cnv_start_hash, iter)) {
                            cnv_index = kh_value(live_cnv_start_hash, iter);
                            addBreakpointInfo(cnv_array, cnv_index, anchor_breakpoints_hash, cur_anchor_breakpoint, header, sfh_idx, sfh, user_inputs);
                        }
                    }

                    // delete the breakpoint start in the live_bpt_start_hash
                    //
                    iter = kh_get(m32, live_bpt_start_hash, cur_anchor_breakpoint - user_inputs->breakpoint_distance);
                    kh_del(m32, live_bpt_start_hash, iter);
                } else {
                    // find the anchor breakpoint using current breakpoint_start
                    // ================================== CNV
                    //      ---------x--------- breakpoint interval 1
                    //                   -----------x--------- breakpoint interval 2
                    // cs1  bs2          bs3  be1       ce1  be2
                    //  1    2            3    2         1    0
                    //                        ==> use breakpoint end potition to get the anchor breakpoint
                    //                                  ==> use the stored Breakpoint_start to get the anchor breakpoint
                    //
                    cnv_start = getValueFromKhash32(cnv_end_start_lookup, all_starts_ends[i]);
                    cnv_index = getValueFromKhash32(live_cnv_start_hash, cnv_start);

                    // loop through the live_bpt_start_hash
                    //
                    for (iter=kh_begin(live_bpt_start_hash); iter!=kh_end(live_bpt_start_hash); ++iter) {
                        if (kh_exist(live_bpt_start_hash, iter)) {
                            uint32_t breakpoint_start = kh_key(live_bpt_start_hash, iter);
                            cur_anchor_breakpoint = getValueFromKhash32(breakpoint_start_hash, breakpoint_start);
                            addBreakpointInfo(cnv_array, cnv_index, anchor_breakpoints_hash, cur_anchor_breakpoint, header, sfh_idx, sfh, user_inputs);
                        }
                    }

                    // now delete current cnv start at the live_cnv_start_hash
                    //
                    iter = kh_get(m32, live_cnv_start_hash, cnv_start);
                    kh_del(m32, live_cnv_start_hash, iter);
                }
            } else {
                // clean-up both live_cnv_start_hash and live_bpt_start_hash
                //
                for (iter=kh_begin(live_cnv_start_hash); iter!=kh_end(live_cnv_start_hash); ++iter) { 
                    if (kh_exist(live_cnv_start_hash, iter)) 
                        kh_del(m32, live_cnv_start_hash, iter);
                }

                for (iter=kh_begin(live_bpt_start_hash); iter!=kh_end(live_bpt_start_hash); ++iter) {
                    if (kh_exist(live_bpt_start_hash, iter))
                        kh_del(m32, live_bpt_start_hash, iter);
                }
            }

            // as some start and end positions might be the same, it is quite confusion
            // also some times both ends have the same value
            // so we need to delete the end ones once we have processed it
            //
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, cnv_end_hash, all_starts_ends[i]);
                if (iter != kh_end(cnv_end_hash)) {
                    kh_del(m32, cnv_end_hash, iter);
                }
            } else if (checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, breakpoint_end_hash, all_starts_ends[i]);
                if (iter != kh_end(breakpoint_end_hash)) {
                    kh_del(m32, breakpoint_end_hash, iter);        // this deletes the key
                }
            }
        } else {
            // start position
            //
            counter++;

            // get current CNV start position in cnv_index
            //
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_cnv_start_hash, all_starts_ends[i], i);
                cnv_index = getSignedValueFromKhash32(cnv_start_hash, all_starts_ends[i]);
                if (cnv_index == -1) {
                    fprintf(stderr, "Something went wrong with CNV/Breakpoint start at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }
                setValueToKhashBucket32(live_cnv_start_hash, all_starts_ends[i], cnv_index);
                //fprintf(stderr, "%s\tCNV_Start\t%"PRIu32"\t%"PRId32"\n", cnv_array->chromosome_id, all_starts_ends[i], counter);
            }

            // get current breakpoint start position
            //
            if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_bpt_start_hash, all_starts_ends[i], i);
                setValueToKhashBucket32(live_bpt_start_hash, all_starts_ends[i], i);
                //fprintf(stderr, "%s\tBreakpoint_Start\t%"PRIu32"\t%"PRId32"\n", cnv_array->chromosome_id, all_starts_ends[i], counter);
            }

            // Note, we don't want to delete the hash keys here as we might need to use them later one
            // instead, we will use the seen* hashtable for the 'skip' purpose
            /*
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                iter = kh_get(m32, cnv_start_hash, all_starts_ends[i]);
                if (iter != kh_end(cnv_start_hash)) {
                    kh_del(m32, cnv_start_hash, iter);
                }
            } else if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i])) {
                iter = kh_get(m32, breakpoint_start_hash, all_starts_ends[i]);
                if (iter != kh_end(breakpoint_start_hash)) {
                    kh_del(m32, breakpoint_start_hash, iter);        // this deletes the key
                }
            }*/
        }
    }

    // clean-up
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, breakpoint_end_start_lookup);
    kh_destroy(m32, breakpoint_start_end_lookup);
    kh_destroy(m32, breakpoint_start_hash);
    kh_destroy(m32, breakpoint_end_hash);
    kh_destroy(m32, cnv_start_end_lookup);
    kh_destroy(m32, cnv_end_start_lookup);
    kh_destroy(m32, cnv_start_hash);
    kh_destroy(m32, cnv_end_hash);
    kh_destroy(m32, live_cnv_start_hash);
    kh_destroy(m32, live_bpt_start_hash);
    kh_destroy(m32, seen_bpt_start_hash);
    kh_destroy(m32, seen_cnv_start_hash);
}

// The breakpoints will be stored in an array
//
void addBreakpointInfo(CNV_Array *cnv_array, uint32_t cnv_index, khash_t(m32) *anchor_breakpoints_hash, uint32_t anchor_breakpoint, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh, User_Input *user_inputs) {
    khint_t k = kh_get(m32, anchor_breakpoints_hash, anchor_breakpoint);
    if (k != kh_end(anchor_breakpoints_hash)) {
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
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].num_of_breakpoints = kh_value(anchor_breakpoints_hash, k);
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].num_of_TLEN_ge_1000=0;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].num_of_paired_reads = 0;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].capacity = 10;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].paired_read_starts = \
            calloc(cnv_array->cnvs[cnv_index].cnv_breakpoints[index].capacity, sizeof(uint32_t));
        cnv_array->cnvs[cnv_index].cnv_breakpoints[index].paired_read_ends = \
            calloc(cnv_array->cnvs[cnv_index].cnv_breakpoints[index].capacity, sizeof(uint32_t));

        storePairedReadsAcrossABreakpoint(cnv_array, cnv_index, anchor_breakpoint, header, sfh_idx, sfh, user_inputs);

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

void storePairedReadsAcrossABreakpoint(CNV_Array *cnv_array, uint32_t cnv_index, uint32_t anchor_breakpoint, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh, User_Input *user_inputs) {
    // skip if it is a singleton
    // Note, the cnv_breakpoints_size is the latest index to the anchor breakpoint array
    //
    uint32_t cnv_breakpoints_size = cnv_array->cnvs[cnv_index].cnv_breakpoints_size;
    if (cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].num_of_breakpoints <= 1)
        return;

    //if (cnv_index == 2524 && anchor_breakpoint==208664075)
    //    printf("anchor_breakpoint=208664075\n");

    uint32_t cur_start = cnv_array->cnvs[cnv_index].raw_bin_start > 0 ? \
                         cnv_array->cnvs[cnv_index].raw_bin_start : cnv_array->cnvs[cnv_index].equal_bin_start;
    uint32_t cur_end   = cnv_array->cnvs[cnv_index].raw_bin_end > 0 ? \
                         cnv_array->cnvs[cnv_index].raw_bin_end : cnv_array->cnvs[cnv_index].equal_bin_end;
    uint32_t cnv_length = cur_end - cur_start;

    // now declare a string 'region' for search, something like chr3:2-3 (breakpoint only 1 bp long)
    // For search region, we need to use bpt_pos not the current anchor position
    // Because we group nearby breakpoint within 5 bps away together, so the search will add addition 10bp both side
    //
    char region[1500];
    sprintf(region, "%s:%"PRIu32"-%"PRIu32"", cnv_array->chromosome_id, anchor_breakpoint-user_inputs->breakpoint_distance, anchor_breakpoint+user_inputs->breakpoint_distance);
    //sprintf(region, "%s:%"PRIu32"-%"PRIu32"", cnv_array->chromosome_id, anchor_breakpoint-user_inputs->breakpoint_distance-10, anchor_breakpoint+user_inputs->breakpoint_distance+10);
    hts_itr_t *hts_itr = sam_itr_querys(sfh_idx, header, region);

    if (hts_itr == NULL) {
        fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", region);
        exit(EXIT_FAILURE);
    }

    khash_t(m32) *seen_starts_hash = kh_init(m32);
    khash_t(m32) *seen_ends_hash = kh_init(m32);
    khiter_t iter;

    bam1_t *b = bam_init1();
    while (sam_itr_next(sfh, hts_itr, b) >= 0) {
        if(b->core.tid != b->core.mtid) continue;       // paired reads are not on the same chromosome

        if(b->core.flag & BAM_FDUP) continue;           // duplicated reads

        if(b->core.flag & BAM_FUNMAP) continue;         // unmapped reads

        if(b->core.flag & BAM_FSECONDARY) continue;     // not the primary reads

        if(b->core.flag & BAM_FQCFAIL) continue;        // Fails Vendor Quality Check

        if(b->core.qual < user_inputs->min_map_quality) continue;       // doesn't pass MAPQ

        //(abs(b->core.isize) >= 1000 && (abs(b->core.isize)*3 >= cnv_length) && (abs(b->core.isize) <= 3*cnv_length))  => wrong
        if (abs(b->core.isize) >= 1000 && (abs(b->core.isize) <= 3*cnv_length)) {

            uint32_t start = (b->core.pos <= b->core.mpos) ? b->core.pos : b->core.mpos;
            uint32_t end = start + abs(b->core.isize);
            uint32_t start_cutoff = start - 500;
            uint32_t end_cutoff = end + 500;

            // first, need to check if the start pos and end pos exist, if so, increase end or shink start by 1
            //
            while (1) {
                iter = kh_get(m32, seen_starts_hash, start);
                if (iter != kh_end(seen_starts_hash)) {
                    start--;
                    if (start <= start_cutoff) break;
                } else {
                    addValueToKhashBucket32(seen_starts_hash, start, 1);
                    break;
                }
            }

            while (1) {
                iter = kh_get(m32, seen_ends_hash, end);
                if (iter != kh_end(seen_ends_hash)) {
                    end++;
                    if (end >= end_cutoff) break;
                } else {
                    addValueToKhashBucket32(seen_ends_hash, end, 1);
                    break;
                }
            }

            // add to the current CNV,
            //
            uint16_t num_of_paired_reads = cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].num_of_paired_reads;
                
            //if (cnv_index == 2524 && anchor_breakpoint==208664075)
            //fprintf(stderr, "cnv_index %"PRIu32"; cnv_breakpoints_size %"PRIu32"; num_of_paired_reads %"PRIu16"; anchor_breakpoint %"PRIu32"; start %"PRIu32"; end %"PRIu32"\n", cnv_index, cnv_breakpoints_size, num_of_paired_reads, anchor_breakpoint, cnv_array->cnvs[cnv_index].equal_bin_start, cnv_array->cnvs[cnv_index].equal_bin_end );

            cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_starts[num_of_paired_reads] = start;
            cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_ends[num_of_paired_reads] = end;
            cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].num_of_paired_reads++;

            // dynamic increase the array size
            //
            if (cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].num_of_paired_reads + 3 >= cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].capacity) {
                cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].capacity += 10;
                cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_starts = 
                    realloc(cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_starts, \
                        cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].capacity * sizeof(uint32_t));
                failureExit(cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_starts, "cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_starts");

                cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_ends =
                    realloc(cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_ends, \
                        cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].capacity * sizeof(uint32_t)); 
                failureExit(cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_ends, "cnv_array->cnvs[cnv_index].cnv_breakpoints[cnv_breakpoints_size].paired_read_ends");
            }
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(hts_itr);

    kh_destroy(m32, seen_starts_hash);
    kh_destroy(m32, seen_ends_hash);
}

void processPairedReadsAcrossABreakpointTlenInfo(CNV_Array *cnv_array) {
    uint16_t i;
    for (i=0; i<cnv_array->size; i++) {
        uint16_t j;
        for (j=0; j<cnv_array->cnvs[i].cnv_breakpoints_size; j++) {
            if (cnv_array->cnvs[i].cnv_breakpoints[j].num_of_paired_reads == 0)
                continue;

            if (cnv_array->cnvs[i].cnv_breakpoints[j].num_of_paired_reads == 1) {
                cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000 = 1;
                continue;
            }

            uint16_t tmp_size = cnv_array->cnvs[i].cnv_breakpoints[j].num_of_paired_reads;
            uint32_t *all_starts_ends_array = calloc(tmp_size*2, sizeof(uint32_t));
            uint16_t k, count=0;
            khash_t(m32) *seen_starts_hash = kh_init(m32);
            khash_t(m32) *seen_ends_hash = kh_init(m32);

            for (k=0; k<tmp_size; k++) {
                all_starts_ends_array[count] = cnv_array->cnvs[i].cnv_breakpoints[j].paired_read_starts[k];
                addValueToKhashBucket32(seen_starts_hash, all_starts_ends_array[count], 1);
                count++;
                all_starts_ends_array[count] = cnv_array->cnvs[i].cnv_breakpoints[j].paired_read_ends[k];
                addValueToKhashBucket32(seen_ends_hash, all_starts_ends_array[count], 1);
                count++;
            }

            qsort(all_starts_ends_array, tmp_size, sizeof(uint32_t), compare);

            // do intersect
            //
            count = 0;
            for (k=0; k<tmp_size; k++) {
                khiter_t iter = kh_get(m32, seen_ends_hash, all_starts_ends_array[k]);
                if (k>0 && iter != kh_end(seen_ends_hash)) {
                    if (count > 0) count--;

                    // delete the end key so that to avoid collide with start
                    //
                    kh_del(m32, seen_ends_hash, iter);
                } else {
                    count++;
                    if (count > cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000) 
                        cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000 = count;
                }
            }

            //fprintf(stderr, "breakpoint: %"PRIu32"\n", cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint);

            // clean-up
            //
            kh_destroy(m32, seen_starts_hash);
            kh_destroy(m32, seen_ends_hash);
            if (all_starts_ends_array) free(all_starts_ends_array);
            seen_starts_hash = NULL;
            seen_ends_hash = NULL;
            all_starts_ends_array=NULL;
        }
    }
}

void setLeftRightCNVBreakpoints(CNV_Array *cnv_array, User_Input *user_inputs) {
    // When using breakpoints to set starts (left-end) and ends (right-end) of CNVs, we need the following
    // 1). number of breakpoints associated with the current anchor breakpoint >= 2 (handled during the intersection)
    // 2). num_of_TLEN_ge_1000 >= 2 (handled during the intersection)
    // 3). left and right breakpoints should be separated by >= 1000
    //

    uint32_t i;
    for (i=0; i<cnv_array->size; i++) {
        // Note: when checking breakpoints, we need to ensure it doesn't cross into neighboring CNVs 
        //
        int32_t prev_cnv_end = -1;      // need to be signed for checking
        int32_t next_cnv_start = -1;    // need to be signed for checking
        if (i>0) {
            if (cnv_array->cnvs[i-1].right_end_index > 0) {
                prev_cnv_end = cnv_array->cnvs[i-1].cnv_breakpoints[cnv_array->cnvs[i-1].right_end_index].breakpoint; 
            } else {
                prev_cnv_end = cnv_array->cnvs[i-1].raw_bin_end > 0 ? \
                    cnv_array->cnvs[i-1].raw_bin_end : cnv_array->cnvs[i-1].equal_bin_end;
            }
        }

        if (i < (cnv_array->size - 1)) {
            if (cnv_array->cnvs[i+1].left_start_index > 0) {
                next_cnv_start = cnv_array->cnvs[i+1].cnv_breakpoints[cnv_array->cnvs[i+1].left_start_index].breakpoint;
            } else {
                next_cnv_start = cnv_array->cnvs[i+1].raw_bin_start > 0 ? \
                    cnv_array->cnvs[i+1].raw_bin_start : cnv_array->cnvs[i+1].equal_bin_start;
            }
        }

        //if (cnv_array->cnvs[i].equal_bin_start == 15492447 || cnv_array->cnvs[i].equal_bin_end == 16293000)
        //    printf("stop 6\n");

        // Get the CNV start and end here as we need them for the checking
        //
        uint32_t cur_start = cnv_array->cnvs[i].raw_bin_start > 0 ? \
                             cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
        uint32_t cur_end   = cnv_array->cnvs[i].raw_bin_end > 0 ? \
                             cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;

        uint32_t left_distance_cutoff = 0, right_distance_cutoff = 0;

        if (prev_cnv_end >= cur_start) {
            left_distance_cutoff = 150;
        } else {
            if (cur_start - prev_cnv_end < 50) {
                left_distance_cutoff = 150;
            } else if (cur_start - prev_cnv_end < 900) {
                left_distance_cutoff = cur_start - prev_cnv_end + 100;
            } else {
                left_distance_cutoff = user_inputs->breakpoint_distance;
            }
        }

        if (next_cnv_start <= cur_end) {
            right_distance_cutoff = 150;
        } else {
            if (cur_end - next_cnv_start <= 50) {
                right_distance_cutoff = 150;
            } else if (cur_end - next_cnv_start < 900) {
                right_distance_cutoff = cur_end - next_cnv_start + 100;
            } else {
                right_distance_cutoff = user_inputs->breakpoint_distance;
            }
        }

        // Need to reset the improperly paired reads TLEN if the following condition is true
        //
        uint32_t cnv_length = cur_end - cur_start;
        uint32_t imp_PR_length = cnv_array->cnvs[i].imp_PR_end - cnv_array->cnvs[i].imp_PR_start;
        //if (imp_PR_length > 5*cnv_length)
        if (imp_PR_length > 3*cnv_length || imp_PR_length*3 < cnv_length)
            cnv_array->cnvs[i].num_of_imp_PR_TLEN_1000 = 0;

        // some of the CNVs don't have breakpoints associated them, so skip
        // the check has to be here as it is used to set the breakpoint
        // don't worry about improperly paired reads checking as it will be handled at the output
        //
        if (cnv_array->cnvs[i].cnv_breakpoints == NULL)
            continue;

        bool left_most_set = false, right_most_set = false;
        uint32_t j;
        for (j=0; j<cnv_array->cnvs[i].cnv_breakpoints_size; j++) {
            /*printf("breakpoint: %"PRIu32"\tnum_of_breakpoints: %"PRIu8"\tnum_of_TLEN_ge_1000: %"PRIu8"\n", \
                    cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint,
                    cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints, 
                    cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000);*/
            
            // All anchor breakpoints associated with the current CNV should be ordered from the left most to the right most.
            // no need to check it here as this is checking during the intersect
            //
            /*if (cnv_array->cnvs[i].cnv_type == 'L' && !(cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints >= 2
                    && cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000 >= 1))
                continue;*/

            bool left=false;

            if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_start)) < \
                    abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_end))) {
                if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_start)) > left_distance_cutoff) {
                    continue;
                } else {
                    left = true;
                }
            } else {
                if ((cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint > cur_end) && 
                        (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_end) > right_distance_cutoff) {
                    // no need to go any further
                    //
                    break;
                }

                if (abs((signed) (cur_end - cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint)) > right_distance_cutoff)
                    continue;
            }

            uint32_t prev_bpts=0, prev_tlen=0, cur_bpts=0, cur_tlen=0;

            if (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint < cur_start) {
                left = true;
            } else if (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint > cur_end) {
                left = false;
            } else {
                /*if (cnv_array->cnvs[i].num_of_imp_PR_TLEN_1000 >= 2) {
                    // need to locate the intersected region with raw-bin or equal-bin
                    // (  -------------------------------- raw-bin
                    //       =============================== improperly paired reads
                    //         xx                         breakpoint
                    //
                    if (cur_start < cnv_array->cnvs[i].imp_PR_start && cnv_array->cnvs[i].imp_PR_start < cur_end) {
                        if (cnv_array->cnvs[i].imp_PR_start - cur_start < cur_end - cnv_array->cnvs[i].imp_PR_start \
                                && cnv_array->cnvs[i].imp_PR_start - cur_start < user_inputs->breakpoint_distance)
                            cur_start = cnv_array->cnvs[i].imp_PR_start;
                    }

                    if (cur_end > cnv_array->cnvs[i].imp_PR_end && cnv_array->cnvs[i].imp_PR_end > cur_start) {
                        if (cur_end - cnv_array->cnvs[i].imp_PR_end > cnv_array->cnvs[i].imp_PR_end - cur_start \
                                && cur_end - cnv_array->cnvs[i].imp_PR_end < user_inputs->breakpoint_distance)
                            cur_end = cnv_array->cnvs[i].imp_PR_end;
                    }
                }*/
                
                if (abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_start)) < \
                    abs((signed) (cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint - cur_end)))
                    left = true;

                // NOTE: here I won't check IF (cur_start > cnv_array->cnvs[i].imp_PR_start) to reset
                // the start position as the length of the improperly paired reads might be quite long
                // i.e., 3 times longer. It will make CNV too long. The same is true for the end position
                //            -----------------------          raw-bin
                //  ==...==============================...==== improperly paired reads
                //  Do not use this case to reset start and end positions
            }

            // however, cur_start and cur_end might change based on the above code changes
            // if the left-hand position > right-hand position after assign left-hand = breakpoint, 
            // we will turn left off as false
            //
            if (left && cnv_array->cnvs[i].cnv_breakpoints[j].breakpoint > cur_end)
                left = false;

            if (left) {
                // the current anchor breakpoint closer to the left-hand side
                //
                if (!left_most_set) {
                    cnv_array->cnvs[i].left_start_index = j;
                    left_most_set = true;
                } else {
                    // left most already set, let's do the comparison
                    //
                    prev_bpts = cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].left_start_index].num_of_breakpoints;
                    prev_tlen = cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].left_start_index].num_of_TLEN_ge_1000;
                    cur_bpts = cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints;
                    cur_tlen = cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000;

                    //if ((cur_bpts >prev_bpts) && (cur_tlen > prev_tlen))
                    if ((cur_bpts + cur_tlen) > (prev_bpts + prev_tlen))
                        cnv_array->cnvs[i].left_start_index = j;
                }
            } else {
                // the current anchor breakpoint closer to the right-hand side
                //
                if (!right_most_set) {
                    cnv_array->cnvs[i].right_end_index = j;
                    right_most_set = true;
                } else {
                    // right hand side already set, let's do the comparison
                    //
                    prev_bpts = cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].right_end_index].num_of_breakpoints;
                    prev_tlen = cnv_array->cnvs[i].cnv_breakpoints[cnv_array->cnvs[i].right_end_index].num_of_TLEN_ge_1000;
                    cur_bpts = cnv_array->cnvs[i].cnv_breakpoints[j].num_of_breakpoints;
                    cur_tlen = cnv_array->cnvs[i].cnv_breakpoints[j].num_of_TLEN_ge_1000;

                    //if ((cur_bpts >prev_bpts) && (cur_tlen > prev_tlen))
                    if ((cur_bpts + cur_tlen) > (prev_bpts + prev_tlen))
                        cnv_array->cnvs[i].right_end_index = j;
                }
            }
        }
    }
}

void checkImproperlyPairedReadsForEachCNV(CNV_Array *cnv_array, Not_Properly_Paired_Reads_Array *improperly_PR_array) {
    // some bam/cram files don't have 'MC' tag, so skip the function
    // num_of_groups was originally set to -1
    //
    if (improperly_PR_array->num_of_groups <= 0)
        return;

    organizeImproperlyPairedReadArray(improperly_PR_array);

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

        /*if (improperly_PR_array->grouped_improperly_PRs[g].group_start == 6028590) {
            printf("stopped\n");
        }*/

        if (num_TLEN_1000 >= 2 ) {
            // TODO Once I add the condition using perfect matched paired reads, I will remove the 2nd condition checking
            // Store TLEN count info as value so that it is easy for lookup
            //
            all_starts_ends[count] = improperly_PR_array->grouped_improperly_PRs[g].group_start;
            setValueToKhashBucket32(imp_PR_start_hash, all_starts_ends[count], num_TLEN_1000);
            count++;

            all_starts_ends[count] = improperly_PR_array->grouped_improperly_PRs[g].group_mate_end;
            setValueToKhashBucket32(imp_PR_end_hash, all_starts_ends[count], num_TLEN_1000);
            count++;
            //fprintf(stderr, "IMP\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", all_starts_ends[count-2], all_starts_ends[count-1], num_TLEN_1000);

            // add values to the lookup table
            //
            setValueToKhashBucket32(imp_PR_end_start_lookup, improperly_PR_array->grouped_improperly_PRs[g].group_mate_end, \
                    improperly_PR_array->grouped_improperly_PRs[g].group_start);
            setValueToKhashBucket32(imp_PR_start_end_lookup, improperly_PR_array->grouped_improperly_PRs[g].group_start, \
                    improperly_PR_array->grouped_improperly_PRs[g].group_mate_end);
        }
    }

    // walk through CNV array and get all its starts and ends
    //
    khash_t(m32) *cnv_start_hash  = kh_init(m32);
    khash_t(m32) *cnv_end_hash    = kh_init(m32);
    khash_t(m32) *cnv_start_end_lookup = kh_init(m32);
    khash_t(m32) *cnv_end_start_lookup = kh_init(m32);

    uint32_t i;
    for (i=0; i<cnv_array->size;i++) {
        // Store the CNV_Array index into the hash tables so that it is easy to lookup and set values
        //
        uint32_t tmp_start = cnv_array->cnvs[i].raw_bin_start>0 ?
                             cnv_array->cnvs[i].raw_bin_start : cnv_array->cnvs[i].equal_bin_start;
        setValueToKhashBucket32(cnv_start_hash, tmp_start, i);
        all_starts_ends[count] = tmp_start;
        count++;

        uint32_t tmp_end = cnv_array->cnvs[i].raw_bin_end > 0 ?
                           cnv_array->cnvs[i].raw_bin_end : cnv_array->cnvs[i].equal_bin_end;
        setValueToKhashBucket32(cnv_end_hash, tmp_end, i);
        all_starts_ends[count] = tmp_end;
        count++;
        //fprintf(stderr, "CNV\t%"PRIu32"\t%"PRIu32"\n", all_starts_ends[count-2], all_starts_ends[count-1]);
        setValueToKhashBucket32(cnv_start_end_lookup, tmp_start, tmp_end);
        setValueToKhashBucket32(cnv_end_start_lookup, tmp_end, tmp_start);

        // need to set all placeholders to 0 for the improperly paired reads
        //
        cnv_array->cnvs[i].imp_PR_start = 0;
        cnv_array->cnvs[i].imp_PR_end   = 0;
        cnv_array->cnvs[i].num_of_imp_PR_TLEN_1000 = 0;
    }

    // Here the capacity won't be equal to the count value as the latter one is filtered, 
    // not all improperly paired reads will be added into the array
    //
    //fprintf(stderr, "The capacity is %"PRId32" while the count is %"PRId32" in the checkImproperlyPairedReadsForEachCNV\n\n", capacity, count);

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
    khash_t(m32) *live_cnv_start_hash  = kh_init(m32);
    //khash_t(m32) *live_cnv_end_hash  = kh_init(m32);
    khash_t(m32) *live_imp_start_hash  = kh_init(m32);
    //khash_t(m32) *live_imp_end_hash  = kh_init(m32);

    khash_t(m32) *seen_cnv_starts_hash = kh_init(m32);
    khash_t(m32) *seen_cnv_ends_hash   = kh_init(m32);
    khash_t(m32) *seen_imp_starts_hash = kh_init(m32);
    khash_t(m32) *seen_imp_ends_hash   = kh_init(m32);

    khiter_t iter;
    count = 0;
    for (i=0; i<capacity; i++) {
        /*if (all_starts_ends[i] == 6028000)
            printf("I am here\n");

        uint32_t num = getValueFromKhash32(imp_PR_start_hash, 6028590);
        fprintf(stderr, "IMPr: %"PRIu32"\n", num);*/

        bool start_first = false;
        // checking cnv-start vs cnv-end
        //
        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        // checking imp-start vs imp-end
        //
        if (checkm32KhashKey(imp_PR_start_hash, all_starts_ends[i]) && checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_imp_starts_hash, all_starts_ends[i]) )
            start_first = true;

        // checking cnv-start vs imp-end
        //
        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        // checking imp-start vs cnv-end
        //
        if (checkm32KhashKey(imp_PR_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_imp_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (!start_first && (checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) ||
                checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i]))) {

            // always decrease count first if it is an end position
            //
            count--;
            /*
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "CNV_End\t%"PRIu32"\t%"PRId32"\n", all_starts_ends[i], count);
            } else {
                fprintf(stderr, "IMP_End\t%"PRIu32"\t%"PRId32"\n", all_starts_ends[i], count);
            }*/

            if (count < 0) {
                if (i == 0) {
                    count = 0;
                } else {
                    fprintf(stderr, "Error: the count %"PRId16" should NOT be negative in imp/cnv intersection", count);
                    fprintf(stderr, "current pos: %"PRIu32" with prev pos %"PRIu32"\n", all_starts_ends[i], all_starts_ends[i-1]);
                    exit(EXIT_FAILURE);
                }
            }

            // when current position is an end and count >= 1, there is an intersect 
            // there might be multiple groups of improperly paired reads
            // need to combine all of them together
            //
            uint32_t imp_start=0, imp_end=0, cnv_start=0, cur_TLEN = 0;
            int32_t cnv_index = -1;         // signed value needed

            if (count >= 1) {

                if (checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i])) {
                    // an improperly paired read end here
                    // loop through all the current live_cnv_start_hash and combine them
                    //
                    cur_TLEN = getValueFromKhash32(imp_PR_end_hash, all_starts_ends[i]);
                    imp_start = getValueFromKhash32(imp_PR_end_start_lookup, all_starts_ends[i]);

                    for (iter=kh_begin(live_cnv_start_hash); iter!=kh_end(live_cnv_start_hash); iter++) {
                        if (kh_exist(live_cnv_start_hash, iter)) {
                            cnv_index = getValueFromKhash32(live_cnv_start_hash, kh_key(live_cnv_start_hash, iter));
                            if ( (all_starts_ends[i] - imp_start) > (3*cnv_array->cnvs[cnv_index].length) )
                                continue;

                            if ( (all_starts_ends[i] - imp_start)*3 < cnv_array->cnvs[cnv_index].length )
                                continue;

                            cnv_array->cnvs[cnv_index].num_of_imp_PR_TLEN_1000 += cur_TLEN;

                            if (cnv_array->cnvs[cnv_index].imp_PR_start == 0 ||
                                    cnv_array->cnvs[cnv_index].imp_PR_start > imp_start)
                                cnv_array->cnvs[cnv_index].imp_PR_start = imp_start;

                            if (cnv_array->cnvs[cnv_index].imp_PR_end == 0 ||
                                    cnv_array->cnvs[cnv_index].imp_PR_end > all_starts_ends[i])
                                cnv_array->cnvs[cnv_index].imp_PR_end = all_starts_ends[i];
                        }
                    }

                    // now delete the corresponding imp_start from live_imp_start_hash
                    //
                    iter = kh_get(m32, live_imp_start_hash, imp_start);
                    kh_del(m32, live_imp_start_hash, iter);
                } else {
                    // a CNV end here
                    // first, get the cnv_start from the cnv_end_start_lookup using all_starts_ends[i] as cnv_end
                    // next, get cnv_index from live_cnv_start_hash using the cnv_start as key
                    //
                    cnv_start = getValueFromKhash32(cnv_end_start_lookup, all_starts_ends[i]);
                    cnv_index = getValueFromKhash32(live_cnv_start_hash, cnv_start);

                    // loop through all the current live_imp_start_hash and add all counts together
                    //
                    for (iter=kh_begin(live_imp_start_hash); iter!=kh_end(live_imp_start_hash); iter++) {
                        if (kh_exist(live_imp_start_hash, iter)) {
                            // get imp_start and imp_end
                            //
                            imp_start = kh_key(live_imp_start_hash, iter);
                            imp_end = getValueFromKhash32(imp_PR_start_end_lookup, imp_start);
                            if ( (imp_end - imp_start) > (3*cnv_array->cnvs[cnv_index].length) )
                                continue;

                            if ( (imp_end - imp_start)*3 < cnv_array->cnvs[cnv_index].length )
                                continue;

                            cur_TLEN = getValueFromKhash32(imp_PR_start_hash, imp_start);
                            cnv_array->cnvs[cnv_index].num_of_imp_PR_TLEN_1000 += cur_TLEN;

                            // elongate the imp_paired_reads on both sides
                            //
                            if (cnv_array->cnvs[cnv_index].imp_PR_start == 0 || 
                                    cnv_array->cnvs[cnv_index].imp_PR_start > imp_start)
                                cnv_array->cnvs[cnv_index].imp_PR_start = imp_start;
                            if (cnv_array->cnvs[cnv_index].imp_PR_end < imp_end )
                                cnv_array->cnvs[cnv_index].imp_PR_end = imp_end;
                        }
                    }

                    // now delete the current CNV start in live_cnv_start_hash
                    //
                    iter = kh_get(m32, live_cnv_start_hash, cnv_start);
                    kh_del(m32, live_cnv_start_hash, iter);
                }
            } else {
                // clean all the values in live_cnv_start_hash and the live_imp_start_hash
                //
                for (iter=kh_begin(live_cnv_start_hash); iter!=kh_end(live_cnv_start_hash); ++iter) {
                    if (kh_exist(live_cnv_start_hash, iter)) {
                        kh_del(m32, live_cnv_start_hash, iter);
                    }
                }

                for (iter=kh_begin(live_imp_start_hash); iter!=kh_end(live_imp_start_hash); ++iter) {
                    if (kh_exist(live_imp_start_hash, iter)) {
                        kh_del(m32, live_imp_start_hash, iter);
                    }   
                }
                        
                /*if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                    cnv_start = getValueFromKhash32(cnv_end_start_lookup, all_starts_ends[i]);
                    if (getSignedValueFromKhash32(live_cnv_start_hash, cnv_start) >= 0)
                        kh_del(m32, live_cnv_start_hash, iter);
                } else if (checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i])) {
                    imp_start = getValueFromKhash32(imp_PR_end_start_lookup, all_starts_ends[i]);
                    if (getSignedValueFromKhash32(live_imp_start_hash, imp_start) >= 0)
                        kh_del(m32, live_imp_start_hash, iter);
                }*/
            }

            // it is possible that some starts and ends are the same value
            // delete the processed positions, so it won't interfere with next position
            //
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, cnv_end_hash, all_starts_ends[i]);
                if (iter != kh_end(cnv_end_hash)) {
                    setValueToKhashBucket32(seen_cnv_ends_hash, all_starts_ends[i], 1);
                }
            } else if (checkm32KhashKey(imp_PR_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, imp_PR_end_hash, all_starts_ends[i]);
                if (iter != kh_end(imp_PR_end_hash)) {
                    setValueToKhashBucket32(seen_imp_ends_hash, all_starts_ends[i], 1);
                }
            }
        } else {
            // it should be in the 'start' position, always increment count
            //
            count++;

            // get current CNV start positioin
            //
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                int32_t cnv_index = getSignedValueFromKhash32(cnv_start_hash, all_starts_ends[i]);
                if (cnv_index == -1) {
                    fprintf(stderr, "Something went wrong with CNV start at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }
                //fprintf(stderr, "%s\tCNV_Start\t%"PRIu32"\t%"PRId32"\n", cnv_array->chromosome_id, all_starts_ends[i], count);
                setValueToKhashBucket32(live_cnv_start_hash, all_starts_ends[i], cnv_index);
            }

            // get improperly paired-read start position
            //
            if (checkm32KhashKey(imp_PR_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(live_imp_start_hash, all_starts_ends[i], 1);
                //fprintf(stderr, "IMP_Start\t%"PRIu32"\t%"PRId32"\n", all_starts_ends[i], count);
            }

            // now delete all the processed starts
            //
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_cnv_starts_hash, all_starts_ends[i], 1);
            } else if (checkm32KhashKey(imp_PR_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_imp_starts_hash, all_starts_ends[i], 1);
            }
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
    kh_destroy(m32, cnv_start_end_lookup);
    kh_destroy(m32, cnv_end_start_lookup);
    kh_destroy(m32, cnv_start_hash);
    kh_destroy(m32, cnv_end_hash);
    kh_destroy(m32, live_cnv_start_hash);
    kh_destroy(m32, live_imp_start_hash);
    kh_destroy(m32, seen_cnv_starts_hash);
    kh_destroy(m32, seen_imp_starts_hash);
    kh_destroy(m32, seen_cnv_ends_hash);
    kh_destroy(m32, seen_imp_ends_hash);
}

void cleanupOverlappingCNVs(CNV_Array *cnv_array, Simple_Stats *equal_window_stats) {
    char prev_type=' ', cur_type=' ';
    double prev_qual=0.0;
    int prev_supporting_evidences = 0;
    int16_t prev_left_idx=0, prev_right_idx=0;
    uint32_t prev_start=0, prev_end=0, cur_start=0, cur_end=0;
    uint32_t prev_left_breakpoint=0, prev_right_breakpoint=0;
    uint32_t prev_left_num_bpoint=0, prev_right_num_bpoint=0;
    uint32_t prev_left_num_geTLEN=0, prev_right_num_geTLEN=0;
    uint32_t prev_num_of_imp_PR_TLEN_1000=0;
    uint32_t j;
    for (j=0; j<cnv_array->size; j++) {
        int16_t left_idx  = cnv_array->cnvs[j].left_start_index;
        int16_t right_idx = cnv_array->cnvs[j].right_end_index;
        //if (right_idx == 65535)
        //    printf("inside\n");

        uint32_t left_breakpoint=0, left_num_bpoint=0, left_num_geTLEN=0;
        uint32_t right_breakpoint=0, right_num_bpoint=0, right_num_geTLEN=0;
        uint32_t num_of_imp_PR_TLEN_1000 = cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000;

        if (left_idx >= 0 || right_idx >= 0 || num_of_imp_PR_TLEN_1000 >= 2) {
            if (cnv_array->cnvs[j].cnv_breakpoints != NULL) {
                left_breakpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].breakpoint : 0;
                left_num_bpoint  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_breakpoints : 0;
                left_num_geTLEN  = (left_idx >= 0)  ? cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_TLEN_ge_1000 : 0;

                right_breakpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].breakpoint : 0;
                right_num_bpoint = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_breakpoints : 0;
                right_num_geTLEN = (right_idx >= 0) ? cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_TLEN_ge_1000 : 0;
            }
        }

        cur_start = (left_idx >= 0) ? left_breakpoint : (cnv_array->cnvs[j].raw_bin_start > 0) ? \
                            cnv_array->cnvs[j].raw_bin_start : cnv_array->cnvs[j].equal_bin_start;
        cur_end   = (right_idx >= 0) ? right_breakpoint : (cnv_array->cnvs[j].raw_bin_end > 0) ? \
                            cnv_array->cnvs[j].raw_bin_end : cnv_array->cnvs[j].equal_bin_end;
        cur_type  = cnv_array->cnvs[j].cnv_type;

        double qual = (cnv_array->cnvs[j].ave_coverage - equal_window_stats->average_coverage) / equal_window_stats->stdev;

        int supporting_evidences = 0;
        if (cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 >= 2) {
            supporting_evidences += 2;

            if (left_idx >= 0 || right_idx >= 0) { 
                if (left_num_bpoint >= 3) supporting_evidences++;
                if (left_num_geTLEN >= 2) supporting_evidences++;
                if (right_num_bpoint >=3) supporting_evidences++;
                if (right_num_geTLEN >=2) supporting_evidences++;
            }

            /*if (cnv_array->cnvs[j].cnv_type == 'P') {
                if (left_num_bpoint >= 25 || right_num_bpoint >= 25)
                    supporting_evidences++;
            }*/
        }

        if (j > 0 && prev_right_breakpoint > 0 && prev_left_breakpoint > 0) {
            if (prev_start < cur_start && cur_start < prev_end) {
                // overlapping
                //
                if (prev_type == cur_type) {
                    // same type, so merge it. Here we only check the previous values > current value
                    //
                    if ( prev_left_breakpoint >= 0 && left_idx >= 0 && 
                            (prev_left_num_bpoint + prev_left_num_geTLEN) > (left_num_bpoint + left_num_geTLEN)) {
                        cnv_array->cnvs[j].cnv_breakpoints[left_idx].breakpoint = prev_left_breakpoint;
                        cnv_array->cnvs[j].cnv_breakpoints[left_idx].num_of_breakpoints = prev_left_num_bpoint;
                    }

                    if ( prev_right_breakpoint >= 0 && right_idx >= 0 && 
                            (prev_right_num_bpoint+prev_right_num_geTLEN) > (right_num_bpoint+right_num_geTLEN)) {
                        cnv_array->cnvs[j].cnv_breakpoints[right_idx].breakpoint = prev_right_breakpoint;
                        cnv_array->cnvs[j].cnv_breakpoints[right_idx].num_of_breakpoints = prev_right_num_bpoint;
                    }
                
                    cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 = prev_num_of_imp_PR_TLEN_1000;
                    //fprintf(stderr, "first: %"PRIu32"\n", j);
                    voidCNVFromList(cnv_array, j-1, prev_left_idx, prev_right_idx);
                } else { // not the same type
                    if (prev_supporting_evidences > supporting_evidences) {
                        //fprintf(stderr, "second: %"PRIu32"\n", j);
                        voidCNVFromList(cnv_array, j, left_idx, right_idx);           
                        right_breakpoint = 0;
                        left_breakpoint  = 0;
                    } else if (supporting_evidences > prev_supporting_evidences) {
                        //fprintf(stderr, "third: %"PRIu32"\n", j);
                        voidCNVFromList(cnv_array, j-1, prev_left_idx, prev_right_idx);
                    }

                    if (abs(prev_qual) >= abs(qual)) {
                        //fprintf(stderr, "fourth %"PRIu32"\t%d\n", j, right_idx);
                        voidCNVFromList(cnv_array, j, left_idx, right_idx);
                        right_breakpoint = 0;
                        left_breakpoint  = 0;
                    } else {
                        //fprintf(stderr, "fifth %"PRIu32"\n", j);
                        voidCNVFromList(cnv_array, j-1, prev_left_idx, prev_right_idx);
                    }
                }
            }
        }
        
        prev_start = cur_start;
        prev_end   = cur_end;
        prev_type  = cur_type;

        prev_left_idx = left_idx;
        prev_left_breakpoint = left_breakpoint;
        prev_left_num_bpoint = left_num_bpoint;
        prev_left_num_geTLEN = left_num_geTLEN;

        prev_right_idx = right_idx;
        prev_right_breakpoint = right_breakpoint;
        prev_right_num_bpoint = right_num_bpoint;
        prev_right_num_geTLEN = right_num_geTLEN;

        prev_num_of_imp_PR_TLEN_1000 = num_of_imp_PR_TLEN_1000;
        prev_supporting_evidences = supporting_evidences;
        prev_qual = qual;
    }
}

void voidCNVFromList(CNV_Array *cnv_array, uint32_t cnv_index, int16_t left_breakpoint_index, int16_t right_breakpoint_index) {
    if (left_breakpoint_index >= 0) {
        cnv_array->cnvs[cnv_index].cnv_breakpoints[left_breakpoint_index].breakpoint = 0;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[left_breakpoint_index].num_of_breakpoints = 0;
        cnv_array->cnvs[cnv_index].left_start_index = -1;
    }

    if (right_breakpoint_index >= 0) {
        cnv_array->cnvs[cnv_index].cnv_breakpoints[right_breakpoint_index].breakpoint = 0;
        cnv_array->cnvs[cnv_index].cnv_breakpoints[right_breakpoint_index].num_of_breakpoints = 0;
        cnv_array->cnvs[cnv_index].right_end_index = -1;
    }

    cnv_array->cnvs[cnv_index].num_of_imp_PR_TLEN_1000 = 0;
    cnv_array->cnvs[cnv_index].raw_bin_start = 0;
    cnv_array->cnvs[cnv_index].raw_bin_end = 0;
    cnv_array->cnvs[cnv_index].equal_bin_start = 0;
    cnv_array->cnvs[cnv_index].equal_bin_end = 0;
}

void outputCNVArray(CNV_Array *cnv_array, char *chrom_id, User_Input *user_inputs, int type) {
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
        //if (cnv_array->cnvs[j].raw_bin_start == 190154072)
        //    printf("Shopping\n");
        int16_t left_idx  = cnv_array->cnvs[j].left_start_index;
        int16_t right_idx = cnv_array->cnvs[j].right_end_index;

        uint32_t left_breakpoint=0, left_num_bpoint=0, left_num_geTLEN=0;
        uint32_t right_breakpoint=0, right_num_bpoint=0, right_num_geTLEN=0;

        if (left_idx >= 0 || right_idx >= 0 || type == 3 || cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 >= 2) {
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
                int supporting_evidences = 0;
                if (cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 >= 2) {
                    supporting_evidences += 2;

                    if (left_idx >= 0 || right_idx >= 0) { 
                        if (left_num_bpoint >= 3) supporting_evidences++;
                        if (left_num_geTLEN >= 2) supporting_evidences++;
                        if (right_num_bpoint >=3) supporting_evidences++;
                        if (right_num_geTLEN >=2) supporting_evidences++;
                    }

                    // For DUP, there will have no left_num_geTLEN and right_num_geTLEN
                    // I will relax the checking for left_num_bpoint and right_num_bpoint
                    // if left_num_bpoint or right_num_bpoint >= 15, I will PASS it
                    // Note: the VCF output need 25, this is just for debugging
                    //
                    //if (cnv_array->cnvs[j].cnv_type == 'P') {
                    //    if (left_num_bpoint >= 25 || right_num_bpoint >= 25)
                    //        supporting_evidences++;
                    //}
                }
            }

            /*if (cnv_array->cnvs[j].equal_bin_start == 18613500) {
                printf("stop 7\n");
            }*/
            uint32_t cnv_start = (left_idx >= 0) ? left_breakpoint : (cnv_array->cnvs[j].raw_bin_start > 0) ? \
                                    cnv_array->cnvs[j].raw_bin_start : cnv_array->cnvs[j].equal_bin_start;
            uint32_t cnv_end = (right_idx >= 0) ? right_breakpoint : (cnv_array->cnvs[j].raw_bin_end > 0) ? \
                                    cnv_array->cnvs[j].raw_bin_end : cnv_array->cnvs[j].equal_bin_end;

            if (cnv_end - cnv_start < user_inputs->min_cnv_length)
                continue;

            char CNV[10];
            (cnv_array->cnvs[j].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "INS");

            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%s\t%.2f\t", chrom_id, \
                    cnv_start, cnv_end, cnv_end-cnv_start, CNV, cnv_array->cnvs[j].ave_coverage);


            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", left_breakpoint, left_num_bpoint, left_num_geTLEN);

            fprintf(fp,"%"PRIu32"\t%"PRIu8"\t%"PRIu8"\t", right_breakpoint, right_num_bpoint, right_num_geTLEN);

            fprintf(fp, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t", cnv_array->cnvs[j].imp_PR_start, cnv_array->cnvs[j].imp_PR_end, cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000);

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

void generateVCFresults(CNV_Array **equal_bin_cnv_array, Chromosome_Tracking *chrom_tracking, Simple_Stats *equal_window_stats, Stats_Info *stats_info, User_Input *user_inputs, FILE *fp, FILE *sfh) {
    // According to Eric, we also need to produce one more file for loading CNV for the review in the TSV format with the following fields
    // chr, start, end, type, quality, var reads, total reads
    // while the 'var reads' is the normalized number of reads in the region (var reads) to what is expected (total reads).
    // we will use zScore for quality
    //
    fprintf(sfh, "chr\tstart\tend\ttype\tquality\tvar-reads\ttotal-reads\n");
    
    // walk through chromosome list
    //
    uint32_t i, j, cnv_start, cnv_end;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        CNV_Array *cnv_array = equal_bin_cnv_array[i];      // pointer assignment, don't need to be free-ed
        for (j=0; j<cnv_array->size; j++) {
            //if (cnv_array->cnvs[j].equal_bin_start == 77091000)
            //    printf("stop vcf\n");

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

            if (cnv_end - cnv_start < user_inputs->min_cnv_length)
                continue;

            char CNV[10];
            (cnv_array->cnvs[j].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "INS");

            char FILTER[30];
            strcpy(FILTER, "PASS");

            char GT[10];
            strcpy(GT, "./1");

            int supporting_evidences = 0;

            if (left_breakpoint > 0 || right_breakpoint > 0 || cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 > 0) {
                if (left_num_bpoint >= 3) supporting_evidences++;
                if (left_num_geTLEN >= 2) supporting_evidences++;
                if (right_num_bpoint >= 3) supporting_evidences++;
                if (right_num_geTLEN >= 2) supporting_evidences++;
                if (cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000 >= 2) supporting_evidences += 2;

                // For DUP, there will have no left_num_geTLEN and right_num_geTLEN
                // I will relax the checking for left_num_bpoint and right_num_bpoint
                // if left_num_bpoint or right_num_bpoint >= 25, I will PASS it
                //
                /*if (cnv_array->cnvs[j].cnv_type == 'P') {
                    if (left_num_bpoint >= 25 || right_num_bpoint >= 25)
                        supporting_evidences++;
                }*/

                if (supporting_evidences == 1) {
                    strcpy(FILTER, "littleBreakpointAndImpSupport");
                    strcpy(GT, "./.");
                }
            } else {
                strcpy(FILTER, "noBreakpointAndImpSupport");
                strcpy(GT, "./.");
            }

            if (supporting_evidences == 0) {
                strcpy(FILTER, "noBreakpointAndImpSupport");
                strcpy(GT, "./.");
            }

            if (supporting_evidences >= 2)
                cnv_array->cnvs[j].inner_cnv.passed = true;

            int32_t svLen = cnv_end-cnv_start;      // must be signed
            if (svLen < user_inputs->min_cnv_length) continue;

            if (cnv_array->cnvs[j].cnv_type == 'L')
                svLen *= -1;

            // now calculate QUAL score using zscore for simple CNV output
            //
            double qual = (cnv_array->cnvs[j].ave_coverage - equal_window_stats->average_coverage) / equal_window_stats->stdev;

            if (strcmp(FILTER, "PASS") == 0) {
                // For the sample CNV output in TSV format
                // chr, start, end, type, quality, var reads, total reads
                // calculate the total number of reads based on Atlas CNV approach
                //
                int32_t total_reads = cnv_array->cnvs[j].ave_coverage * (cnv_end-cnv_start) / stats_info->read_cov_stats->read_length;

                fprintf(sfh, "%s\t%"PRIu32"\t%"PRIu32"\t%s\t%.2f\t.\t%"PRIu32"\n", chrom_tracking->chromosome_ids[i], cnv_start, cnv_end, CNV, qual, (uint32_t) total_reads);
            }

            fprintf(fp, "%s\t%"PRIu32"\t.\tN\t%s\t%.2f\t%s\tEND=%"PRIu32";SVLEN=%"PRId32";SVTYPE=%s;AVGCOV=%.2f;BPTL=%"PRIu32";BPTLCOUNT=%"PRIu8";BPTLTLEN=%"PRIu8";BPTR=%"PRIu32";BPTRCOUNT=%"PRIu8";BPTRTLEN=%"PRIu8";IMPPRLEN=%"PRIu8"\tGT\t%s\n", \
                    chrom_tracking->chromosome_ids[i], cnv_start, CNV, qual, FILTER, cnv_end, svLen, CNV, cnv_array->cnvs[j].ave_coverage, left_breakpoint, left_num_bpoint, left_num_geTLEN, right_breakpoint, right_num_bpoint, right_num_geTLEN, cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000, GT);

            // update the inner_CNV variable for CNV segmentation
            //
            cnv_array->cnvs[j].inner_cnv.start = cnv_start;
            cnv_array->cnvs[j].inner_cnv.end   = cnv_end;
            cnv_array->cnvs[j].inner_cnv.qual  = qual;
            cnv_array->cnvs[j].inner_cnv.ave_coverage = cnv_array->cnvs[j].ave_coverage;
            cnv_array->cnvs[j].inner_cnv.cnv_type = cnv_array->cnvs[j].cnv_type;
            cnv_array->cnvs[j].inner_cnv.ave_coverage = cnv_array->cnvs[j].ave_coverage;
            cnv_array->cnvs[j].inner_cnv.left_breakpoint = left_breakpoint;
            cnv_array->cnvs[j].inner_cnv.left_breakpoint_count = left_num_bpoint;
            cnv_array->cnvs[j].inner_cnv.num_larger_TLEN_left  = left_num_geTLEN;
            cnv_array->cnvs[j].inner_cnv.right_breakpoint = right_breakpoint;
            cnv_array->cnvs[j].inner_cnv.right_breakpoint_count = right_num_bpoint;
            cnv_array->cnvs[j].inner_cnv.last_right_breakpoint_count = right_num_bpoint;
            cnv_array->cnvs[j].inner_cnv.num_larger_TLEN_right  = right_num_geTLEN;
            cnv_array->cnvs[j].inner_cnv.last_num_larger_TLEN_right  = right_num_geTLEN;
            cnv_array->cnvs[j].inner_cnv.imp_PR_start = cnv_array->cnvs[j].imp_PR_start;
            cnv_array->cnvs[j].inner_cnv.imp_PR_end = cnv_array->cnvs[j].imp_PR_end;
            cnv_array->cnvs[j].inner_cnv.num_larger_imp_PR_TLEN = cnv_array->cnvs[j].num_of_imp_PR_TLEN_1000;
            cnv_array->cnvs[j].inner_cnv.evidence_count = supporting_evidences;
            if (strcmp(FILTER, "PASS") == 0) {
                cnv_array->cnvs[j].inner_cnv.passed = true;
            } else {
                cnv_array->cnvs[j].inner_cnv.passed = false;
            }
            if (supporting_evidences >= 2)
                cnv_array->cnvs[j].inner_cnv.valid_cnv = true;

        } // end equal_bin_cnv_array
    } // end chromosome list
}

void outputLog2Ratio(Binned_Data_Wrapper **binned_data_wrapper, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs) {
    FILE *log2r_fh = fopen(user_inputs->log2ratio_output_file, "w");

    uint32_t i=0, j=0;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        for (j=0; j<binned_data_wrapper[i]->size; j++) {
            if (binned_data_wrapper[i]->data[j].length > 0) {
                fprintf(log2r_fh, "%s\t%"PRIu32"\t%.2f\t%.2f\n", chrom_tracking->chromosome_ids[i], binned_data_wrapper[i]->data[j].start, binned_data_wrapper[i]->data[j].log2ratio, binned_data_wrapper[i]->data[j].ave_coverage);
            }
        }
    }

    fclose(log2r_fh);
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
        cnv_array[i]->chrom_length = chrom_tracking->chromosome_lengths[i];
        cnv_array[i]->capacity = 5*PR_INIT_SIZE;
        cnv_array[i]->cnvs = calloc(cnv_array[i]->capacity, sizeof(CNV));
        cnv_array[i]->size = 0;
    }
}

void cnvArrayDestroy(CNV_Array **cnv_array, uint32_t number_of_chromosomes) {
    uint32_t i, j, k;
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
                for (k=0; k<cnv_array[i]->cnvs[j].cnv_breakpoints_size; k++) {
                    if (cnv_array[i]->cnvs[j].cnv_breakpoints[k].paired_read_starts)
                        free(cnv_array[i]->cnvs[j].cnv_breakpoints[k].paired_read_starts);

                    if (cnv_array[i]->cnvs[j].cnv_breakpoints[k].paired_read_ends)
                        free(cnv_array[i]->cnvs[j].cnv_breakpoints[k].paired_read_ends);
                }
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
