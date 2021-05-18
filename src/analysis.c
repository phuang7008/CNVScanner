/*
 * =====================================================================================
 *
 *		Filename:		analysis.c
 *
 *		Description:	The implementation file for the analysis functions
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang, phuang@bcm.edu
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access() and getopt()
#include <dirent.h>		// for checking output directory
#include <libgen.h>		// for function basename()

#include "analysis.h"
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"
#include "fileProcessing.h"

void coverageBinningWrapper(Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Binned_Data_Wrapper *binned_data_wrapper, int32_t chrom_idx) {
    FILE *wgs_binned_coverage_fp = fopen(user_inputs->wgs_binning_file, "a");

    writeCoverageBins(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_idx, user_inputs, stats_info, wgs_binned_coverage_fp, binned_data_wrapper);

    fclose(wgs_binned_coverage_fp);
}

void writeCoverageBins(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, int32_t chrom_idx, User_Input *user_inputs, Stats_Info *stats_info, FILE *fh_binned_coverage, Binned_Data_Wrapper *binned_data_wraper) {
    // NOTE: for the bed format, the end position is included!
    // next, we need to provide various boundaries for binning
    //       haploid             diploid               duplication               10x of haploid         anything > 10x of haploid
    //  [0x –– Dip*0.66)  [Dip*0.66 – Dip*1.32)  [Dip*1.32 – Dip + Dip *0.6) [Dip + Dip*0.6 – 10xHap]            > 10xHap
    //
    double hap_upper  = user_inputs->average_coverage * 0.66;   // from paper
    double dip_upper  = user_inputs->average_coverage * 1.32;   // from paper
    double dup_upper = user_inputs->average_coverage * 10 / 2;

    fprintf(stderr, "#haploid upper bound is:\t%f\n", hap_upper);
    fprintf(stderr, "#diploid upper bound is:\t%f\n", dip_upper);
    fprintf(stderr, "#dup upper bound is:\t%f\n", dup_upper);

    uint32_t i=0;
    for (i=begin; i<begin+length; i++) {
        uint32_t start=0, end=0;
        uint64_t cov_total=0;

        // check if the size is approach the capacity
        //
        if (binned_data_wraper->size + 20 > binned_data_wraper->capacity) {
            fprintf(stderr, "before dynamic size increase %d\n", binned_data_wraper->capacity);
            dynamicIncreaseBinSize(binned_data_wraper);
            fprintf(stderr, "after dynamic size increase %d\n", binned_data_wraper->capacity);
        }

        // For the excluded bases
        //
        if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] == -1)) {
            start = i;
            while ( i < begin+length && (chrom_tracking->coverage[chrom_idx][i] == -1)) {
                i++;
            }
            end = i;
            if (user_inputs->debug_ON)
                fprintf(fh_binned_coverage, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t-1\n",
                        chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start);
        } else if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] < hap_upper)) {
            start = i;
            while( i < begin+length && (0 <= chrom_tracking->coverage[chrom_idx][i]) && (chrom_tracking->coverage[chrom_idx][i] < hap_upper)) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;

            stats_info->wgs_cov_stats->total_genome_coverage += cov_total;

            if (start <= end) {
                processBinnedData(start, end, cov_total, binned_data_wraper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
            }
        } else if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] < dip_upper)) {
            start = i;
            while( i < begin+length && (hap_upper <= chrom_tracking->coverage[chrom_idx][i]) 
                    && (chrom_tracking->coverage[chrom_idx][i] < dip_upper)) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;

            stats_info->wgs_cov_stats->total_genome_coverage += cov_total;

            if (start <= end) {
                processBinnedData(start, end, cov_total, binned_data_wraper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
            }
        } else if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] < dup_upper)) {
            start = i;
            while( i < begin+length && (dip_upper <= chrom_tracking->coverage[chrom_idx][i]) &&
                    (chrom_tracking->coverage[chrom_idx][i] < dup_upper)) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;

            stats_info->wgs_cov_stats->total_genome_coverage += cov_total;

            if (start <= end) {
                processBinnedData(start, end, cov_total, binned_data_wraper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
            }
        } else {
            start = i;
            while( i < begin+length && (chrom_tracking->coverage[chrom_idx][i] >= dup_upper)) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;

            stats_info->wgs_cov_stats->total_genome_coverage += cov_total;

            if (start <= end) {
                processBinnedData(start, end, cov_total, binned_data_wraper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
            }
        }

        i--;
    }
}

/* 
 * During the insertion, we are also conducting edge smoothing
 * Case 1: we will check the neighboring bins to see if their coverage average is closer (such as < 5x)
 * In this case, we will merge them together
 *
 * Case 2: After binning, the neighboring bins' average coverage are less than 5x for the previous two bins
 * In this case, we will also merge them together
 *
 */
void processBinnedData(uint32_t start, uint32_t end, uint32_t coverage, Binned_Data_Wrapper *binned_data_wraper, Chromosome_Tracking *chrom_tracking, FILE *fh_binned_coverage, int32_t chrom_idx, User_Input *user_inputs) {
    double ave_coverage = (double)coverage / (double)(end - start);
    if (user_inputs->debug_ON)
        fprintf(fh_binned_coverage, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n",
                chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);

    // here we need to check with the previous bin block,
    // if they are only 5x away, combined them together
    //
    double cov_total = 0.0;
    if (binned_data_wraper->size > 1) {
        uint32_t prev_index = binned_data_wraper->size - 1;
        if ( binned_data_wraper->data[prev_index].end == start &&
                abs(binned_data_wraper->data[prev_index].ave_coverage - ave_coverage) <= DIFF_COV_TO_MERGE) {
            cov_total = (double) coverage + 
                binned_data_wraper->data[prev_index].ave_coverage * binned_data_wraper->data[prev_index].length;
            binned_data_wraper->data[prev_index].end    = end;
            binned_data_wraper->data[prev_index].length += end - start;
            binned_data_wraper->data[prev_index].ave_coverage = cov_total / (double)binned_data_wraper->data[prev_index].length;
        } else {
            // check the previous two bins and see if they need to be merged!
            //
            if (binned_data_wraper->size > 2) {
                uint32_t pre_prev_index = binned_data_wraper->size - 2;
                if (abs(binned_data_wraper->data[pre_prev_index].ave_coverage - 
                            binned_data_wraper->data[prev_index].ave_coverage) <= DIFF_COV_TO_MERGE) {
                    cov_total = binned_data_wraper->data[prev_index].ave_coverage * binned_data_wraper->data[prev_index].length +
                        binned_data_wraper->data[pre_prev_index].ave_coverage * binned_data_wraper->data[pre_prev_index].length;
                    binned_data_wraper->data[pre_prev_index].end = binned_data_wraper->data[prev_index].end;
                    binned_data_wraper->data[pre_prev_index].length += binned_data_wraper->data[prev_index].length;
                    binned_data_wraper->data[pre_prev_index].ave_coverage = cov_total / binned_data_wraper->data[pre_prev_index].length;
                    binned_data_wraper->size--;     // need to reduce the size by one

                    // now add the new data in
                    //
                    insertBinData(start, end, end-start, ave_coverage, binned_data_wraper);
                } else {
                    if (binned_data_wraper->data[pre_prev_index].length < SMALL_LENGTH_CUTOFF) {
                        // if the length of the pre_prev_index entry is < 50 (SMALL_LENGTH_CUTOFF), remove it
                        //
                        /*binned_data_wraper->data[pre_prev_index].start  = binned_data_wraper->data[prev_index].start;
                        binned_data_wraper->data[pre_prev_index].end    = binned_data_wraper->data[prev_index].end;
                        binned_data_wraper->data[pre_prev_index].length = binned_data_wraper->data[prev_index].length;
                        binned_data_wraper->data[pre_prev_index].ave_coverage = binned_data_wraper->data[prev_index].ave_coverage;
                        binned_data_wraper->size--;*/
                    }
                        
                    insertBinData(start, end, end-start, ave_coverage, binned_data_wraper);
                }
            } else {
                insertBinData(start, end, end-start, ave_coverage, binned_data_wraper);
            }
        }
    } else {
        insertBinData(start, end, end-start, ave_coverage, binned_data_wraper);
    }
}

void insertBinData(uint32_t start, uint32_t end, uint32_t length, double ave_coverage, Binned_Data_Wrapper *binned_data_wrapper) {
    binned_data_wrapper->data[binned_data_wrapper->size].start  = start;
    binned_data_wrapper->data[binned_data_wrapper->size].end    = end;
    binned_data_wrapper->data[binned_data_wrapper->size].length = length;
    binned_data_wrapper->data[binned_data_wrapper->size].ave_coverage = ave_coverage;
    binned_data_wrapper->starts[binned_data_wrapper->size] = start;
    binned_data_wrapper->ends[binned_data_wrapper->size]   = end;

    binned_data_wrapper->size++;
    //fprintf(stderr, "Current size is %d\n", binned_data_wrapper->size);

}

void outputBinnedData(Binned_Data_Wrapper *binned_data_wrapper, char* chrom_id) {
    FILE *binned_coverage_fp = fopen("Binned_Data.txt", "w");

    uint32_t i;
    for (i=0; i<binned_data_wrapper->size; i++) {
        fprintf(binned_coverage_fp, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n", 
                chrom_id, binned_data_wrapper->data[i].start, binned_data_wrapper->data[i].end,
                binned_data_wrapper->data[i].length, binned_data_wrapper->data[i].ave_coverage);
    }

    if (binned_coverage_fp) fclose(binned_coverage_fp);
}

void reportStatsForDebugging(Stats_Info *stats_info, User_Input *user_inputs) {
    if (stats_info->read_cov_stats->total_reads_aligned == 0) {
        fprintf(stderr, "ERROR: No reads aligned. Aborting.\n");
        exit(EXIT_FAILURE);
    }

    uint32_t non_duplicate_reads = stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads;
    if (non_duplicate_reads == 0) {
        fprintf(stderr, "ERROR: All reads are duplicates. Aborting.\n");
        exit(EXIT_FAILURE);
    }

    uint64_t sum=0;
    int32_t i=0;
    double average_coverage=0.0;
    khiter_t k_iter;
    uint32_t coverage_bins[20000] = {0};    // initialize array contents to 0

    //Do not consider the Ns for Median calculation.
    //
    uint64_t total_genome_non_Ns_bases = stats_info->wgs_cov_stats->total_genome_bases - stats_info->wgs_cov_stats->total_excluded_bases;

    for (k_iter=0; k_iter!=kh_end(stats_info->wgs_cov_stats->genome_coverage_for_median); k_iter++) {
        if (kh_exist(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter)) {
            if (kh_key(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter) < 20000)
                coverage_bins[kh_key(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter)] =
                    kh_value(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter);
        }
    }

    for (i=0; i<20000; i++) {

        if (sum >= (total_genome_non_Ns_bases/2)) {
            stats_info->wgs_cov_stats->median_genome_coverage = i--;
            break;
        }else{
            sum += coverage_bins[i];
            if (i == 0)
                sum -= stats_info->wgs_cov_stats->total_excluded_bases;
        }
    }

    // open WGS coverage summary report file handle
    //
    FILE *out_fp = fopen(user_inputs->wgs_cov_report, "a");

    average_coverage = (double) stats_info->wgs_cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;
    outputGeneralInfo(out_fp, stats_info, average_coverage);

    //fprintf(out_fp, "WGS Max_Coverage\t%"PRIu32"\n", stats_info->wgs_cov_stats->wgs_max_coverage);
    //fprintf(out_fp, "Number_of_Bases_with_WGS_Max_Coverage\t%"PRIu32"\n", stats_info->wgs_cov_stats->base_with_wgs_max_coverage);

    fprintf(out_fp, "#Base_Stats\n");
    fprintf(out_fp, "Total_Genome_Base_Targeted_w/o_Ns\t%"PRIu64"\n", total_genome_non_Ns_bases);
    fprintf(out_fp, "Total_Ns_Bases_in_Ns_Regions\t%"PRIu32"\n", stats_info->wgs_cov_stats->total_excluded_bases);
    fprintf(out_fp, "Total_Mapped_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->total_mapped_bases);
    fprintf(out_fp, "Total_Uniquely_Aligned_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->total_uniquely_aligned_bases);

    float percent = calculatePercentage64(stats_info->wgs_cov_stats->total_overlapped_bases, stats_info->wgs_cov_stats->total_mapped_bases);
    fprintf(out_fp, "Total_Overlapped_Aligned_Bases\t%"PRIu32"\n", stats_info->wgs_cov_stats->total_overlapped_bases);
    fprintf(out_fp, "PCT_Overlapped_Aligned_Bases\t%0.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->wgs_cov_stats->base_quality_20, stats_info->wgs_cov_stats->total_mapped_bases);
    fprintf(out_fp, "Aligned_Q20_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->base_quality_20);
    fprintf(out_fp, "PCT_Aligned_Q20_Bases\t%0.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->wgs_cov_stats->base_quality_30, stats_info->wgs_cov_stats->total_mapped_bases);
    fprintf(out_fp, "Aligned_Q30_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->base_quality_30);
    fprintf(out_fp, "PCT_Aligned_Q30_Bases\t%0.2f%%\n", percent);

    //fprintf(out_fp, "Median_Coverage\t%d\n", stats_info->wgs_cov_stats->median_genome_coverage);

    fprintf(out_fp, "\n");
    fclose(out_fp);
}

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage) {
    fprintf(fp, "#Read_Stats\n");

    fprintf(fp, "Total_Reads(TR)\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_produced);

    uint64_t yield = stats_info->read_cov_stats->read_length * (uint64_t) stats_info->read_cov_stats->total_reads_produced;
    fprintf(fp, "Sequenced_Read_Length\t%d\n", stats_info->read_cov_stats->read_length);
    fprintf(fp, "Total_Yield\t%"PRIu64"\n", yield);

    float percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_aligned,stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Aligned_Reads(AR)\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_Reads_Aligned\t%.2f%%\n", percent);

    yield = stats_info->read_cov_stats->read_length * (uint64_t) (stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads);
    uint64_t uniquely_aligned = stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads;
    percent = calculatePercentage64(uniquely_aligned, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Uniquely_Aligned_Yield\t%"PRIu64"\n", yield);
    fprintf(fp, "Unique_Aligned_Reads\t%"PRIu64"\n", uniquely_aligned);
    fprintf(fp, "PCT_of_Unique_Aligned_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(uniquely_aligned, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Unique_Aligned_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_duplicate_reads, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Duplicate_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_duplicate_reads);
    fprintf(fp, "PCT_of_Duplicate_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_duplicate_reads,stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Duplicate_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_supplementary_reads,stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Supplementary_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_supplementary_reads);
    fprintf(fp, "PCT_of_Supplementary_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_supplementary_reads,stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Supplementary_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_paired, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Paired_READ\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_paired);
    fprintf(fp, "PCT_of_Paired_READ_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_paired, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Paired_READ_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_paired_reads_with_mapped_mates, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Paired_Reads_With_Mapped_Mates\t%"PRIu64"\n", stats_info->read_cov_stats->total_paired_reads_with_mapped_mates);
    fprintf(fp, "PCT_of_Paired_Reads_With_Mapped_Mates_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_paired_reads_with_mapped_mates, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Paired_Reads_With_Mapped_Mates_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_proper_paired, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Properly_Paired_Reads\t%"PRIu64"\n",stats_info->read_cov_stats->total_reads_proper_paired);
    fprintf(fp, "PCT_of_Properly_Paired_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_proper_paired, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Properly_Paired_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_chimeric_reads, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Chimeric_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_chimeric_reads);
    fprintf(fp, "PCT_of_Chimeric_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_chimeric_reads, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Chimeric_Reads_(agst_AR)\t%.2f%%\n", percent);

    fprintf(fp, "Average_Coverage\t%.2f\n", average_coverage);
}

void binnedDataWrapperInit(Binned_Data_Wrapper** binned_data_wrapper, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // here we will initialize each binned data array to be 100,000 defined as INIT_BIN_SIZE
        // the size will dynamically expanded when there are more needed
        //
        binned_data_wrapper[i] = calloc(1, sizeof(Binned_Data_Wrapper));
        binned_data_wrapper[i]->chromosome_id = calloc(strlen(chrom_tracking->chromosome_ids[i])+1, sizeof(char*));
        strcpy(binned_data_wrapper[i]->chromosome_id, chrom_tracking->chromosome_ids[i]);
        binned_data_wrapper[i]->size = 0;
        binned_data_wrapper[i]->capacity = INIT_BIN_SIZE;
        binned_data_wrapper[i]->data   = calloc(INIT_BIN_SIZE, sizeof(Binned_Data));
        binned_data_wrapper[i]->starts = calloc(INIT_BIN_SIZE, sizeof(uint32_t));
        binned_data_wrapper[i]->ends   = calloc(INIT_BIN_SIZE, sizeof(uint32_t));
    }
}

void dynamicIncreaseBinSize(Binned_Data_Wrapper* binned_data_wrapper) {
    if (!binned_data_wrapper) return;

    binned_data_wrapper->capacity += INIT_BIN_SIZE;

    binned_data_wrapper->data = 
        realloc(binned_data_wrapper->data, (binned_data_wrapper->capacity)*sizeof(Binned_Data));

    exitWithFailure(binned_data_wrapper->data);

    binned_data_wrapper->starts = 
        realloc(binned_data_wrapper->starts, (binned_data_wrapper->capacity)*sizeof(uint32_t));

    exitWithFailure(binned_data_wrapper->starts);

    binned_data_wrapper->ends = 
        realloc(binned_data_wrapper->ends, (binned_data_wrapper->capacity)*sizeof(uint32_t));

    exitWithFailure(binned_data_wrapper->ends);

}

void exitWithFailure(void * data_point_in) {
    if (data_point_in == NULL) {
        fprintf(stderr, "ERROR: Dynamic Memory allocation for start positions array\n");
        exit(EXIT_FAILURE);
    }
}

void binnedDataWrapperDestroy(Binned_Data_Wrapper** binned_data_wrapper, Chromosome_Tracking *chrom_tracking) {

    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        free(binned_data_wrapper[i]->chromosome_id);
        free(binned_data_wrapper[i]->starts);
        free(binned_data_wrapper[i]->data);
        free(binned_data_wrapper[i]->ends);
        free(binned_data_wrapper[i]);
    }
    free(binned_data_wrapper);
}

void mappabilityNormalization(Binned_Data_Wrapper *binned_data_wraper, khash_t(khIntStr) *map_starts, khash_t(khIntStr) *map_ends, uint32_t total_map_lines) {
    // create an all starts and ends array, the size will be dynamically increased later
    //
    AllStartsEndsArray *all_starts_ends_array = calloc(1, sizeof(AllStartsEndsArray));
    all_starts_ends_array->capacity = binned_data_wraper->size*2 + total_map_lines;
    all_starts_ends_array->array = calloc(all_starts_ends_array->capacity, sizeof(uint32_t));
    all_starts_ends_array->size = 0;

    khash_t(khIntStr) *binned_starts  = kh_init(khIntStr);      // key: start, value: "start end length ave_cov"
    khash_t(khIntStr) *binned_ends    = kh_init(khIntStr);      // key: end,   value: "start end length ave_cov"

    generateHashFromDynamicBins(binned_data_wraper, binned_starts, binned_ends, all_starts_ends_array, 1);
    combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, map_starts);
    combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, map_ends);

    // sort the all_starts_ends_array
    //
    qsort(all_starts_ends_array->array, all_starts_ends_array->size, sizeof(uint32_t), compare);

    // do intersect
    //
    uint32_t i=0, counter=0;
    uint32_t prev_start0=0, prev_start1=0, map_position=0;

    for (i=0; i<all_starts_ends_array->size; i++) {
        if (checkKhashKey(map_ends, all_starts_ends_array->array[i]) || 
                checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {
            counter--;

            if (counter == 1) {
                // the first part of the annotation should come from binned_starts or binned_ends
                // the second part of the annotation should come from map_starts or map_ends
                //
                char *binned_string=NULL, *map_string=NULL;

                if (checkKhashKey(binned_starts, prev_start0))
                    binned_string = getKhashValue(binned_starts, prev_start0);

                if (checkKhashKey(map_starts, map_position)) {
                    map_string = getKhashValue(map_starts, map_position);
                } else if (checkKhashKey(map_ends, map_position)) {
                    map_string = getKhashValue(map_ends, map_position);
                }

                generateNormalizedMappabilityForCurrentBin(binned_data_wraper, 
                        binned_string, map_string, all_starts_ends_array->array[i], prev_start1);

                if (binned_string) free(binned_string);
                if (map_string) free(map_string);
            }

            // because bed file starts with 0 (or 0-indxed), so one value will appear in both starts and ends hash-tables
            // we have to remove those appeas in end position, so the next round, it will be the start position only
            //
            khiter_t iter;
            if (checkKhashKey(map_ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, map_ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(map_ends))
                    kh_del(khIntStr, map_ends, iter);
            } else if (checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, binned_ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(binned_ends))
                    kh_del(khIntStr, binned_ends, iter);
            }
        } else {
            // it should be in the start position
            // need to record the start position of current intersect
            //
            if (checkKhashKey(binned_starts, all_starts_ends_array->array[i])) {
                prev_start0 = all_starts_ends_array->array[i];      // this is needed for the dynamic bin annotation
                counter++;
                prev_start1 = all_starts_ends_array->array[i];      // this could be either from the dynamic binned or from map
            } else if (checkKhashKey(map_starts, all_starts_ends_array->array[i])) {
                counter++;
                prev_start1 = all_starts_ends_array->array[i];
            }
        }

        // need to record the map position of current intersect for the map annotation
        //
        if (checkKhashKey(map_ends, all_starts_ends_array->array[i]) || 
                checkKhashKey(map_starts, all_starts_ends_array->array[i])) {
            map_position = all_starts_ends_array->array[i];
        }
    }

    // clean-up
    //
    cleanAllStartsEndsArray(all_starts_ends_array);
    cleanKhashIntStr(binned_starts);
    cleanKhashIntStr(binned_ends);
}

void generateNormalizedMappabilityForCurrentBin(Binned_Data_Wrapper *binned_data_wraper, char *bin_string, char* map_string, uint32_t current_position, uint32_t prev_start) {
    StringArray *binned_array = calloc(1, sizeof(StringArray));
    stringArrayInit(binned_array, 10);
    splitStringToArray(bin_string, binned_array);

    StringArray *mapped_array = calloc(1, sizeof(StringArray));
    stringArrayInit(mapped_array, 10);
    splitStringToArray(map_string, mapped_array);

    // Note: the binned_array->theArray[0] is the index to the binned_array_wrapper->data
    //
    int length = current_position - prev_start;
    int orig_len = strtoul(binned_array->theArray[2], NULL, 10) - strtoul(binned_array->theArray[1], NULL, 10);
    double ave = strtod(binned_array->theArray[3], NULL);
    double map_ratio = strtod(mapped_array->theArray[3], NULL);
    double weighted_mappability = ((double)length * ave) / (orig_len * map_ratio);
    binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_map_normalized += weighted_mappability;

    // output for debugging
    //
    fprintf(stderr, "%"PRIu32"\t%"PRIu32"\t%.2f\t%s\t%s\t%.2f\n", prev_start, current_position, weighted_mappability, bin_string, map_string, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_map_normalized);

    // clean up
    //
    stringArrayDestroy(binned_array);
    stringArrayDestroy(mapped_array);
}

// type 1 is for mappability normalization, while type 2 is for gc normalization
//
void generateHashFromDynamicBins(Binned_Data_Wrapper *binned_data_wrapper, khash_t(khIntStr) *binned_starts, khash_t(khIntStr) *binned_ends, AllStartsEndsArray *all_starts_ends_array, int type) {
    // create a string pointer to store values
    //
    char * insert_value = calloc(100, sizeof(char));

    uint32_t i=0;
    for (i = 0; i < binned_data_wrapper->size; i++) {
        // Note here we need to store index i to the insert_value for later usage
        //
        if (type == 1) {
            sprintf(insert_value, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f", i, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i], binned_data_wrapper->data[i].ave_coverage);
        } else {
            sprintf(insert_value, "%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f", i, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i], binned_data_wrapper->data[i].ave_cov_map_normalized);
        }

        khashInsertion(binned_starts, binned_data_wrapper->starts[i], insert_value);
        khashInsertion(binned_ends, binned_data_wrapper->ends[i], insert_value);

        all_starts_ends_array->array[all_starts_ends_array->size] = binned_data_wrapper->starts[i];
        all_starts_ends_array->size++;
        all_starts_ends_array->array[all_starts_ends_array->size] = binned_data_wrapper->ends[i];
        all_starts_ends_array->size++;
    }
}

void combineAllStartsAndEndsFromOtherSource(AllStartsEndsArray *all_starts_ends_array, khash_t(khIntStr) *hash_in) {
    khiter_t iter;
    for (iter=kh_begin(hash_in); iter!=kh_end(hash_in); ++iter) {
        if (kh_exist(hash_in, iter)) {
            all_starts_ends_array->array[all_starts_ends_array->size] = kh_key(hash_in, iter);
            all_starts_ends_array->size++;
        }
    }
}
