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

void generateCNVs(CNV_Array **equal_bin_cnv_array, Binned_Data_Wrapper **equal_size_window_wrappers, Binned_Data_Wrapper **raw_bin_data_wrappers, uint32_t number_of_chromosomes, Simple_Stats *the_stats, User_Input *user_inputs) {
    uint32_t i;
#pragma omp parallel shared(the_stats) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        for (i=0; i<number_of_chromosomes; ++i) {
            if (strcmp(equal_size_window_wrappers[i]->chromosome_id, raw_bin_data_wrappers[i]->chromosome_id) != 0) {
                fprintf(stderr, "ERROR: chromosome ids for equal window bins (%s) and raw bins (%s) are not the same\n", \
                            equal_size_window_wrappers[i]->chromosome_id, raw_bin_data_wrappers[i]->chromosome_id);
                exit(EXIT_FAILURE);
            }
#pragma omp task
            {
                int thread_id = omp_get_thread_num();
                printf("Current thread id in generating CNV calls: %d\n", thread_id);

                mergeNeighboringBinsBasedOnZscore(equal_bin_cnv_array[i], equal_size_window_wrappers[i], the_stats, 2);

                expandMergedCNVWithRawBins(raw_bin_data_wrappers[i], equal_bin_cnv_array[i], the_stats);

            } // end task 
        } // end for loop
#pragma omp taskwait

      } // end single
    } // end parallel
}

// Type 1: raw varying sized bins;  Type 2: Equal window bins
//
void mergeNeighboringBinsBasedOnZscore(CNV_Array *cnv_array, Binned_Data_Wrapper *equal_size_window_wrapper, Simple_Stats *equal_window_stats, int type) {
    uint32_t j;
    double hap_cutoff = equal_window_stats->average_coverage - equal_window_stats->zScore;
    double dup_cutoff = equal_window_stats->average_coverage + equal_window_stats->zScore;

    // flag used to indicate which category current bin belongs to: 1 -> hap; 2 -> dip; 3 -> dup;
    //
    uint8_t prev_flag=0, cur_flag=0;     
    uint32_t cnv_counter=0, bin_counter=0, bin_capacity=0;
    uint32_t prev_start=0, prev_end=0, prev_len;
    double total_coverage=0;
    uint32_t conseutive_0_length_bins = 0;      // bins with length 0 means the bins are Ns-regions or repeatmasks etc.
    Equal_Window_Bin *merged_equal_bin_array = NULL;

    for (j=0; j<equal_size_window_wrapper->size; j++) {
        if (equal_size_window_wrapper->data[j].length == 0) {
            if (prev_flag == 0) {
                conseutive_0_length_bins = 0;
                continue;
            } else {
                conseutive_0_length_bins += equal_size_window_wrapper->data[j].end - equal_size_window_wrapper->data[j].start;
            }
        }

        // starts a new CNV
        //
        if (prev_flag == 0) {
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
                cur_flag = prev_flag;
            } else if (equal_size_window_wrapper->data[j].ave_coverage <= hap_cutoff) {
                cur_flag = 1;
            } else if (equal_size_window_wrapper->data[j].ave_coverage >= dup_cutoff) {
                cur_flag = 3;
            } else {
                cur_flag = 2;
            }

            // combine neighboring bins together
            // For raw bins, if the distance (or gap) between neighboring bins is <= 200, merge them
            // This 200 value is calculated by averaging all Ns regions, repeatmask regions, low mappability regions etc.
            //
            if ((cur_flag == prev_flag) && (conseutive_0_length_bins <= 500) && 
                    ((type == 2 && equal_size_window_wrapper->data[j].start == prev_end) || 
                     (type == 1 && equal_size_window_wrapper->data[j].start - prev_end <= 200))) {
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

                if (bin_counter == bin_capacity - 2) {
                    bin_capacity = bin_capacity * 2;
                    dynamicIncreaseBinArraySize(&merged_equal_bin_array, bin_capacity);
                }

            } else {
                // the previous CNV differs from current CNV. 
                // Let's store the previous CNV first if its total length > user specified length
                // NOTE: we won't count this as CNV if the length used for CNV call is < 500bp
                //   OR: if the ratio between prev_len / (prev_end - prev_start) < 0.5, skip it
                //
                if (prev_end - prev_start >= 1000) {
                    if (prev_len >= 500 && ((double)prev_len/((double)prev_end - (double)prev_start)) >= 0.5) {
                        double coverage = total_coverage / prev_len;
                        storeCurrentCNVtoArray(&cnv_array->cnvs[cnv_counter], \
                            prev_start, prev_end, prev_len, coverage, merged_equal_bin_array, bin_counter);
                        cnv_counter++;
                    }
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

                if (cur_flag == 2) {
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
}

void storeCurrentCNVtoArray(CNV *cnv, uint32_t start, uint32_t end, uint32_t length, double coverage, Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_size) {
    cnv->start  = start;
    cnv->end    = end; 
    cnv->length = length;
    cnv->ave_coverage = coverage;

    // store all equal window bin info to the current CNV
    //
    cnv->size = bin_size;
    cnv->capacity = bin_size;
    cnv->equal_bin_array = calloc(cnv->size, sizeof(Equal_Window_Bin));
    uint32_t k;
    for (k=0; k<cnv->size; k++) {
        cnv->equal_bin_array[k].start  = merged_equal_bin_array[k].start;
        cnv->equal_bin_array[k].end    = merged_equal_bin_array[k].end;
        cnv->equal_bin_array[k].length = merged_equal_bin_array[k].length;
        cnv->equal_bin_array[k].ave_coverage = merged_equal_bin_array[k].ave_coverage;
    }
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
        for (j=raw_bin_start; j<binned_data_wrapper->size; j++) {
            if (cnv_array->cnvs[i].start > binned_data_wrapper->data[j].end) { 
                // no intersect, record restart position and then skip
                //
                raw_bin_start = j;
                continue;
            } else if (cnv_array->cnvs[i].end < binned_data_wrapper->data[j].start) {
                // no intersect. The current raw bin has pass the current CNV bin
                // no need to continue, just break out
                //
                break;
            } else {
                raw_bin_start = j;

                // they intersect. Now elongate the current CNV
                //
                if (cnv_array->cnvs[i].start > binned_data_wrapper->data[j].start) {
                    cnv_array->cnvs[i].start = binned_data_wrapper->data[j].start;
                    addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                    
                    // now let's extend it further at the left end
                    //
                    j--;
                    while ((binned_data_wrapper->data[j].ave_coverage <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                           (binned_data_wrapper->data[j].ave_coverage >= dup_cutoff 
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        // check if the distance is within 200bp away
                        // 200 bp is the average length of all excluded regions (such as Ns regions, repeatmasks)
                        //
                        if (cnv_array->cnvs[i].start - binned_data_wrapper->data[j].end <= 200) {
                            cnv_array->cnvs[i].start = binned_data_wrapper->data[j].start;
                            addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);
                            j--;
                        } else {
                            break;
                        }
                    }
                }

                j = raw_bin_start;
                if (cnv_array->cnvs[i].end < binned_data_wrapper->data[j].end) {
                    cnv_array->cnvs[i].end = binned_data_wrapper->data[j].end;
                    addRawBinToCNV(binned_data_wrapper, j, cnv_array->cnvs, i);

                    // now let's extend the CNV further at the right end
                    //
                    j++;
                    while ((binned_data_wrapper->data[j].ave_coverage <= hap_cutoff
                                        && cnv_array->cnvs[i].ave_coverage <= hap_cutoff) || 
                            (binned_data_wrapper->data[j].ave_coverage >= dup_cutoff
                                        && cnv_array->cnvs[i].ave_coverage >= dup_cutoff)) {
                        // check if the distance is within 200bp away
                        // 200 bp is the average length of all excluded regions (such as Ns regions, repeatmasks)
                        //
                        if (binned_data_wrapper->data[j].start - cnv_array->cnvs[i].end <= 200) {
                            cnv_array->cnvs[i].end = binned_data_wrapper->data[j].end;
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
    cnv[cnv_index].size++;
}

void outputCNVArray(CNV_Array **cnv_array, uint32_t number_of_chromosomes, int type) {
    FILE *fp;
    if (type == 1) {        // For raw varying size bin output
        fp = fopen("CNV_Array_raw_bin.txt", "w");
    } else {                // For equal size window output
        fp = fopen("CNV_Array_equal_bin.txt", "w");
    }
    uint32_t i,j;
    for (i=0; i<number_of_chromosomes; i++) {
        for (j=0; j<cnv_array[i]->size; j++) {
            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array[i]->chromosome_id, \
                    cnv_array[i]->cnvs[j].start, cnv_array[i]->cnvs[j].end, \
                    cnv_array[i]->cnvs[j].end - cnv_array[i]->cnvs[j].start,\
                    cnv_array[i]->cnvs[j].length, cnv_array[i]->cnvs[j].ave_coverage);

            uint32_t k;
            for (k=0; k<cnv_array[i]->cnvs[j].size; k++) {
                fprintf(fp, "\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array[i]->cnvs[j].equal_bin_array[k].start, \
                        cnv_array[i]->cnvs[j].equal_bin_array[k].end, \
                        cnv_array[i]->cnvs[j].equal_bin_array[k].length, \
                        cnv_array[i]->cnvs[j].equal_bin_array[k].ave_coverage);
            }
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
