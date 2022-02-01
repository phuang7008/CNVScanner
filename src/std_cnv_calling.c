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

double calculate99Percentile(Binned_Data_Wrapper **equal_size_window_wrapper, uint32_t number_of_chromosomes) {
    uint32_t i, j;
    uint32_t total_bins=0;
    for (i=0; i<number_of_chromosomes; i++) {
        for (j=0; j<equal_size_window_wrapper[i]->size; j++) {
            // Need to skip those windows with repeat maskers, mappbility <0.2, Ns regions etc
            //
            //if (equal_size_window_wrapper[i]->data[j].length < 250)
            if (equal_size_window_wrapper[i]->data[j].length == 0)
                continue;

            total_bins++;
        }
    }

    DoubleArray *equal_window_coverages = calloc(1, sizeof(DoubleArray));
    equal_window_coverages->size = total_bins;
    equal_window_coverages->array = calloc(equal_window_coverages->size, sizeof(double));
    uint32_t count=0;
    for (i=0; i<number_of_chromosomes; i++) { 
        for (j=0; j<equal_size_window_wrapper[i]->size; j++) {
            //if (equal_size_window_wrapper[i]->data[j].length < 250) 
            if (equal_size_window_wrapper[i]->data[j].length == 0) 
                continue;
            equal_window_coverages->array[count] = equal_size_window_wrapper[i]->data[j].ave_coverage;
            count++;
        }
    }

    double percentile =  CalcualtePercentile(equal_window_coverages, 99);
    free(equal_window_coverages);
    return percentile;
}

void generateStatsForNormalizedData(Binned_Data_Wrapper **equal_size_window_wrapper, uint32_t number_of_chromosomes, Simple_Stats *equal_window_stats) {
    // setup one pass stdev data-structure
    //
    double total_cov_sum  = 0;
    uint64_t total_bins   = 0;
    double total_sum_of_bin_cov_square = 0;

    uint32_t i, j;
    for (i=0; i<number_of_chromosomes; i++) {
        for (j=0; j<equal_size_window_wrapper[i]->size; j++) {
            // // Need to skip those windows with repeat maskers, mappbility <0.2, Ns regions etc
            //
            //if (equal_size_window_wrapper[i]->data[j].length < 250 || equal_size_window_wrapper[i]->data[j].ave_coverage >= cutoff_99p)
            //if (equal_size_window_wrapper[i]->data[j].length == 0 || equal_size_window_wrapper[i]->data[j].ave_coverage >= cutoff_99p)
            if (equal_size_window_wrapper[i]->data[j].length == 0)
                continue;

            total_bins++;
            total_cov_sum += equal_size_window_wrapper[i]->data[j].ave_coverage;
            total_sum_of_bin_cov_square += equal_size_window_wrapper[i]->data[j].ave_coverage * equal_size_window_wrapper[i]->data[j].ave_coverage;
        }
    }

    // get the stats here and store them into equal_window_stats
    //
    equal_window_stats->average_coverage  = total_cov_sum / total_bins;
    equal_window_stats->total_bases_used  = total_bins;
    equal_window_stats->stdev = sqrt((total_sum_of_bin_cov_square / total_bins) - pow(equal_window_stats->average_coverage, 2));
    //equal_window_stats->outlier_cutoff = 1.645 * equal_window_stats->stdev;     // 90% confident inverval for z score
    equal_window_stats->outlier_cutoff = 1.96 * equal_window_stats->stdev;    // 95% confident inverval for z score
    //equal_window_stats->outlier_cutoff = 3 * equal_window_stats->stdev;
}

void mergeNeighboringBinsBasedOnZscore(CNV_Array **cnv_array, Binned_Data_Wrapper **equal_size_window_wrapper, uint32_t number_of_chromosomes, Simple_Stats *equal_window_stats) {
    uint32_t i,j;
    double hap_cutoff = equal_window_stats->average_coverage - equal_window_stats->outlier_cutoff;
    double dup_cutoff = equal_window_stats->average_coverage + equal_window_stats->outlier_cutoff;

    for (i=0; i<number_of_chromosomes; i++) {
        // flag used to indicate which category current bin belongs to: 1 -> hap; 2 -> dip; 3 -> dup;
        //
        uint8_t prev_flag=0, cur_flag=0;     
        uint32_t cnv_counter=0, bin_counter=0, bin_capacity=0;
        uint32_t prev_start=0, prev_end=0, prev_len;
        double total_coverage=0;
        Equal_Window_Bin *merged_equal_bin_array = NULL;

        for (j=0; j<equal_size_window_wrapper[i]->size; j++) {
            if (equal_size_window_wrapper[i]->data[j].length == 0) {
                // store current CNV if pass the length criteria
                //
                if ( prev_end == 7053500) {
                    printf("stop here\n");
                }
                if (prev_end - prev_start >= 1000) { 
                    double coverage = total_coverage / prev_len;
                    storeCurrentCNVtoArray(&cnv_array[i]->CNVs_per_chromosome->cnvs[cnv_counter], \
                            prev_start, prev_end, prev_len, coverage, merged_equal_bin_array, bin_counter);
                    cnv_counter++;

                    // clean-up and reset variables
                    //
                    /*if (merged_equal_bin_array != NULL) {
                        free(merged_equal_bin_array);
                        merged_equal_bin_array=NULL;
                    }*/
                    prev_start = 0;
                    prev_end   = 0;
                    prev_len   = 0;
                    bin_counter  = 0;
                    prev_flag = 0;
                }
                continue;
            }

            if (prev_flag == 0) {
                if (equal_size_window_wrapper[i]->data[j].ave_coverage <= hap_cutoff) {
                    prev_flag = 1;
                } else if (equal_size_window_wrapper[i]->data[j].ave_coverage >= dup_cutoff) {
                    prev_flag = 3;
                } else {
                    //prev_flag = 2;
                    continue;
                }

                prev_start = equal_size_window_wrapper[i]->data[j].start;
                prev_end   = equal_size_window_wrapper[i]->data[j].end;
                prev_len   = equal_size_window_wrapper[i]->data[j].length;
                total_coverage = equal_size_window_wrapper[i]->data[j].ave_coverage * prev_len;

                // store current bin info (clean it before using it)
                //
                /*if (merged_equal_bin_array) {
                    free(merged_equal_bin_array);
                    merged_equal_bin_array=NULL;
                    bin_counter = 0;
                    bin_capacity = 0;
                }

                bin_capacity = EQUAL_BIN_SIZE;
                merged_equal_bin_array = calloc(bin_capacity, sizeof(Equal_Window_Bin));
                merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper[i]->data[j].start;
                merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper[i]->data[j].end;
                merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper[i]->data[j].length;
                merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper[i]->data[j].ave_coverage;
                bin_counter++;*/
            } else {
                if (equal_size_window_wrapper[i]->data[j].ave_coverage <= hap_cutoff) {
                    cur_flag = 1;
                } else if (equal_size_window_wrapper[i]->data[j].ave_coverage >= dup_cutoff) {
                    cur_flag = 3;
                } else {
                    cur_flag = 2;
                }

                if (cur_flag == prev_flag && prev_end == equal_size_window_wrapper[i]->data[j].start) {
                    // elongate bin size by resetting the end position and length
                    //
                    prev_end  = equal_size_window_wrapper[i]->data[j].end;
                    prev_len += equal_size_window_wrapper[i]->data[j].length;
                    total_coverage += equal_size_window_wrapper[i]->data[j].ave_coverage * equal_size_window_wrapper[i]->data[j].length;

                    // add current window bin info
                    //
                    /*merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper[i]->data[j].start;
                    merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper[i]->data[j].end;
                    merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper[i]->data[j].length;
                    merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper[i]->data[j].ave_coverage;
                    bin_counter++;
                    */

                    if (bin_counter == bin_capacity - 2) {
                        bin_capacity = bin_capacity * 2;
                        dynamicIncreaseBinArraySize(merged_equal_bin_array, bin_capacity);
                    }

                } else {
                    // the previous CNV differs from current CNV. 
                    // Let's store the previous CNV first if its total length > user specified length
                    //
                    if (prev_end - prev_start >= 1000) {
                        double coverage = total_coverage / prev_len;
                        storeCurrentCNVtoArray(&cnv_array[i]->CNVs_per_chromosome->cnvs[cnv_counter], \
                                prev_start, prev_end, prev_len, coverage, merged_equal_bin_array, bin_counter);
                        cnv_counter++;
                    }

                    // clean-up
                    //
                    fprintf(stderr, "bin counter is %"PRIu32"\n", bin_counter);
                    fprintf(stderr, "cnv counter is %"PRIu32"\n", cnv_counter);
                    /*if (merged_equal_bin_array != NULL) {
                        free(merged_equal_bin_array);
                        merged_equal_bin_array=NULL;
                    }*/

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
                    prev_start = equal_size_window_wrapper[i]->data[j].start;
                    prev_end   = equal_size_window_wrapper[i]->data[j].end;
                    prev_len   = equal_size_window_wrapper[i]->data[j].length;
                    total_coverage = equal_size_window_wrapper[i]->data[j].ave_coverage * prev_len;

                    bin_capacity = EQUAL_BIN_SIZE;
                    merged_equal_bin_array = calloc(bin_capacity, sizeof(Equal_Window_Bin));
                    merged_equal_bin_array[bin_counter].start  = equal_size_window_wrapper[i]->data[j].start;
                    merged_equal_bin_array[bin_counter].end    = equal_size_window_wrapper[i]->data[j].end;
                    merged_equal_bin_array[bin_counter].length = equal_size_window_wrapper[i]->data[j].length;
                    merged_equal_bin_array[bin_counter].ave_coverage =  equal_size_window_wrapper[i]->data[j].ave_coverage;
                    bin_counter++;
                }
            }
        }
        cnv_array[i]->CNVs_per_chromosome->size = cnv_counter;
        cnv_counter = 0;
        if (merged_equal_bin_array) {
            free(merged_equal_bin_array);
            merged_equal_bin_array=NULL;
        }
    }
}

void storeCurrentCNVtoArray(CNV *cnv, uint32_t start, uint32_t end, uint32_t length, double coverage, Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_size) {
    cnv->start  = start;
    cnv->end    = end; 
    cnv->length = length;
    cnv->ave_coverage = coverage;

    // store all equal window bin info to the current CNV
    //
    /*cnv->equal_bin_array = calloc(bin_size, sizeof(Equal_Window_Bin));
    cnv->size = bin_size;
    uint32_t k;
    for (k=0; k<bin_size; k++) {
        cnv->equal_bin_array[k].start  = merged_equal_bin_array[k].start;
        cnv->equal_bin_array[k].end    = merged_equal_bin_array[k].end;
        cnv->equal_bin_array[k].length = merged_equal_bin_array[k].length;
        cnv->equal_bin_array[k].ave_coverage = merged_equal_bin_array[k].ave_coverage;
    }*/
}

void outputCNVArray(CNV_Array **cnv_array, uint32_t number_of_chromosomes) {
    FILE *fp = fopen("CNV_Array.txt", "w");
    uint32_t i,j;
    for (i=0; i<number_of_chromosomes; i++) {
        for (j=0; j<cnv_array[i]->CNVs_per_chromosome->size; j++) {
            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\n", cnv_array[i]->chromosome_id, \
                    cnv_array[i]->CNVs_per_chromosome->cnvs[j].start, cnv_array[i]->CNVs_per_chromosome->cnvs[j].end, \
                    cnv_array[i]->CNVs_per_chromosome->cnvs[j].end - cnv_array[i]->CNVs_per_chromosome->cnvs[j].start,\
                    cnv_array[i]->CNVs_per_chromosome->cnvs[j].length, cnv_array[i]->CNVs_per_chromosome->cnvs[j].ave_coverage);
        }
    }
    fclose(fp);
}

void dynamicIncreaseBinArraySize(Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_capacity) {
    if (merged_equal_bin_array) {
        merged_equal_bin_array = realloc(merged_equal_bin_array, bin_capacity * sizeof(Equal_Window_Bin));
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
        cnv_array[i]->CNVs_per_chromosome = calloc(1, sizeof(CNVs_Per_Chromosome));
        cnv_array[i]->CNVs_per_chromosome->capacity = 5*PR_INIT_SIZE;
        cnv_array[i]->CNVs_per_chromosome->cnvs = calloc(cnv_array[i]->CNVs_per_chromosome->capacity, sizeof(CNV));
        cnv_array[i]->CNVs_per_chromosome->size = 0;
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

        for (j=0; j<cnv_array[i]->CNVs_per_chromosome->size;j++) {
            if (cnv_array[i]->CNVs_per_chromosome->cnvs[j].equal_bin_array) {
                free(cnv_array[i]->CNVs_per_chromosome->cnvs[j].equal_bin_array);
                cnv_array[i]->CNVs_per_chromosome->cnvs[j].equal_bin_array = NULL;
            }
        }

        if (cnv_array[i]->CNVs_per_chromosome->cnvs) {
            free(cnv_array[i]->CNVs_per_chromosome->cnvs);
            cnv_array[i]->CNVs_per_chromosome->cnvs = NULL;
        }
    }
    if (cnv_array)
        free(cnv_array);
}
