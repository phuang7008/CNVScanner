/*
 * =====================================================================================
 *
 *		Filename:		reports.c
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

#include "reports.h"

void coverageBinningWrapper(Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Binned_Data_Wrapper *binned_data_wrapper, int32_t chrom_idx, Simple_Stats *wgs_simple_stats, int thread_id) {
    FILE *wgs_binned_coverage_fp = fopen(user_inputs->wgs_binning_file, "a");
    fileOpenError(wgs_binned_coverage_fp, user_inputs->wgs_binning_file);

    writeCoverageBins(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_idx, user_inputs, stats_info, wgs_binned_coverage_fp, binned_data_wrapper, wgs_simple_stats, thread_id);

    fclose(wgs_binned_coverage_fp);
}

void writeCoverageBins(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, int32_t chrom_idx, User_Input *user_inputs, Stats_Info *stats_info, FILE *fh_binned_coverage, Binned_Data_Wrapper *binned_data_wrapper, Simple_Stats *wgs_simple_stats, int thread_id) {
    // NOTE: for the bed format, the end position is included!
    // next, we need to provide various boundaries for binning
    //       haploid             diploid               duplication               10x of haploid         anything > 10x of haploid
    //  [0x –– Dip*0.66)  [Dip*0.66 – Dip*1.32)  [Dip*1.32 – Dip + Dip *0.6) [Dip + Dip*0.6 – 10xHap]            > 10xHap
    //
    double hap_upper  = user_inputs->average_coverage * 0.66;   // from paper
    double dip_upper  = user_inputs->average_coverage * 1.32;   // from paper
    double dup_upper = user_inputs->average_coverage * 10 / 2;

    fprintf(stderr, "For chromosome %s using thread %d\n", chrom_tracking->chromosome_ids[chrom_idx], thread_id);
    /*fprintf(stderr, "\t#haploid upper bound is:\t%f \n", hap_upper);
    fprintf(stderr, "\t#diploid upper bound is:\t%f \n", dip_upper);
    fprintf(stderr, "\t#dup upper bound is:\t%f \n", dup_upper);
    */

    uint32_t i=0;
    for (i=begin; i<begin+length; i++) {
        uint32_t start=0, end=0;
        uint64_t cov_total=0;
        //if (i==6818612) {
        //    printf("6818612 stopped!");
        //}

        // check if the size is approach the capacity
        //
        if (binned_data_wrapper->size + 20 > binned_data_wrapper->capacity) {
            fprintf(stderr, "before dynamic size increase %d\n", binned_data_wrapper->capacity);
            dynamicIncreaseBinSize(binned_data_wrapper);
            fprintf(stderr, "after dynamic size increase %d\n", binned_data_wrapper->capacity);
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
        } else if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] >= wgs_simple_stats->outlier_cutoff)) {
           start = i;
           while (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] >= wgs_simple_stats->outlier_cutoff)) {
               i++;
           }
           end = i;
           if (user_inputs->debug_ON)
               fprintf(fh_binned_coverage, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t-2\n",
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
                processBinnedData(start, end, cov_total, binned_data_wrapper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
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
                processBinnedData(start, end, cov_total, binned_data_wrapper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
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
                processBinnedData(start, end, cov_total, binned_data_wrapper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
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
                processBinnedData(start, end, cov_total, binned_data_wrapper, chrom_tracking, fh_binned_coverage, chrom_idx, user_inputs);
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
void processBinnedData(uint32_t start, uint32_t end, uint32_t coverage, Binned_Data_Wrapper *binned_data_wrapper, Chromosome_Tracking *chrom_tracking, FILE *fh_binned_coverage, int32_t chrom_idx, User_Input *user_inputs) {
    double ave_coverage = (double)coverage / (double)(end - start);
    if (user_inputs->debug_ON)
        fprintf(fh_binned_coverage, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n",
                chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);

    // here we need to check with the previous bin block,
    // if they are only 5x away, combined them together
    //
    double cov_total = 0.0;
    if (binned_data_wrapper->size > 1) {
        uint32_t prev_index = binned_data_wrapper->size - 1;
        if ( binned_data_wrapper->data[prev_index].end == start && 
                abs(binned_data_wrapper->data[prev_index].ave_coverage - ave_coverage) <= DIFF_COV_TO_MERGE) { 
            cov_total = (double) coverage + 
                binned_data_wrapper->data[prev_index].ave_coverage * binned_data_wrapper->data[prev_index].length;
            binned_data_wrapper->data[prev_index].end    = end;
            binned_data_wrapper->data[prev_index].length += end - start;
            binned_data_wrapper->data[prev_index].ave_coverage = cov_total / (double)binned_data_wrapper->data[prev_index].length;
            binned_data_wrapper->ends[prev_index] = end;        // we don't have to update the start position
        } else {
            // check the previous two bins and see if they need to be merged!
            //
            if (binned_data_wrapper->size > 2) {
                uint32_t pre_prev_index = binned_data_wrapper->size - 2;
                if (abs(binned_data_wrapper->data[pre_prev_index].ave_coverage - 
                            binned_data_wrapper->data[prev_index].ave_coverage) <= DIFF_COV_TO_MERGE &&
                            binned_data_wrapper->data[pre_prev_index].end == binned_data_wrapper->data[prev_index].start) {
                    cov_total = binned_data_wrapper->data[prev_index].ave_coverage * binned_data_wrapper->data[prev_index].length +
                        binned_data_wrapper->data[pre_prev_index].ave_coverage * binned_data_wrapper->data[pre_prev_index].length;
                    binned_data_wrapper->data[pre_prev_index].end = binned_data_wrapper->data[prev_index].end;
                    binned_data_wrapper->data[pre_prev_index].length += binned_data_wrapper->data[prev_index].length;
                    binned_data_wrapper->data[pre_prev_index].ave_coverage = cov_total / binned_data_wrapper->data[pre_prev_index].length;
                    binned_data_wrapper->ends[pre_prev_index] = binned_data_wrapper->data[prev_index].end;
                    binned_data_wrapper->size--;     // need to reduce the size by one

                    // now add the new data in
                    //
                    insertBinData(start, end, end-start, ave_coverage, binned_data_wrapper);
                } else {
                    //if (binned_data_wrapper->data[pre_prev_index].length < SMALL_LENGTH_CUTOFF) {
                        // if the length of the pre_prev_index entry is < 50 (SMALL_LENGTH_CUTOFF), remove it
                        //
                        /*binned_data_wrapper->data[pre_prev_index].start  = binned_data_wrapper->data[prev_index].start;
                        binned_data_wrapper->data[pre_prev_index].end    = binned_data_wrapper->data[prev_index].end;
                        binned_data_wrapper->data[pre_prev_index].length = binned_data_wrapper->data[prev_index].length;
                        binned_data_wrapper->data[pre_prev_index].ave_coverage = binned_data_wrapper->data[prev_index].ave_coverage;
                        binned_data_wrapper->size--;*/
                    //}
                        
                    insertBinData(start, end, end-start, ave_coverage, binned_data_wrapper);
                }
            } else {
                insertBinData(start, end, end-start, ave_coverage, binned_data_wrapper);
            }
        }
    } else {
        insertBinData(start, end, end-start, ave_coverage, binned_data_wrapper);
    }
}

void insertBinData(uint32_t start, uint32_t end, uint32_t length, double ave_coverage, Binned_Data_Wrapper *binned_data_wrapper) {
    binned_data_wrapper->data[binned_data_wrapper->size].start  = start;
    binned_data_wrapper->data[binned_data_wrapper->size].end    = end;
    binned_data_wrapper->data[binned_data_wrapper->size].length = length;
    binned_data_wrapper->data[binned_data_wrapper->size].index  = binned_data_wrapper->size;
    binned_data_wrapper->data[binned_data_wrapper->size].ave_coverage = ave_coverage;
    binned_data_wrapper->data[binned_data_wrapper->size].weighted_mappability = 0;
    binned_data_wrapper->data[binned_data_wrapper->size].ave_cov_map_normalized = 0;
    binned_data_wrapper->data[binned_data_wrapper->size].weighted_gc_scale = 0;
    binned_data_wrapper->data[binned_data_wrapper->size].ave_cov_map_gc_normalized = 0;

    binned_data_wrapper->starts[binned_data_wrapper->size] = start;
    binned_data_wrapper->ends[binned_data_wrapper->size]   = end;

    binned_data_wrapper->size++;
    //fprintf(stderr, "Current size is %d\n", binned_data_wrapper->size);

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
    fileOpenError(out_fp, user_inputs->wgs_cov_report);

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

void binnedDataWrapperInit(Binned_Data_Wrapper** binned_data_wrappers, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // here we will initialize each binned data array to be 100,000 defined as INIT_BIN_SIZE
        // the size will dynamically expanded when there are more needed
        //
        binned_data_wrappers[i] = calloc(1, sizeof(Binned_Data_Wrapper));
        binned_data_wrappers[i]->chromosome_id = calloc(strlen(chrom_tracking->chromosome_ids[i])+1, sizeof(char*));
        strcpy(binned_data_wrappers[i]->chromosome_id, chrom_tracking->chromosome_ids[i]);

        binned_data_wrappers[i]->size = 0;
        binned_data_wrappers[i]->capacity = INIT_BIN_SIZE;
        binned_data_wrappers[i]->data   = calloc(INIT_BIN_SIZE, sizeof(Binned_Data));
        binned_data_wrappers[i]->starts = calloc(INIT_BIN_SIZE, sizeof(uint32_t));
        binned_data_wrappers[i]->ends   = calloc(INIT_BIN_SIZE, sizeof(uint32_t));
    }
}

void dynamicIncreaseBinSize(Binned_Data_Wrapper* binned_data_wrapper) {
    if (!binned_data_wrapper) return;

    binned_data_wrapper->capacity += INIT_BIN_SIZE;

    if (binned_data_wrapper->data) {
        binned_data_wrapper->data = 
            realloc(binned_data_wrapper->data, (binned_data_wrapper->capacity)*sizeof(Binned_Data));
        exitWithFailure(binned_data_wrapper->data, "binned_data_wrapper->data");
    } else {
        fprintf(stderr, "Error: The binned_data_wrapper->data is NULL\n");
        exit(EXIT_FAILURE);
    }


    binned_data_wrapper->starts = 
        realloc(binned_data_wrapper->starts, (binned_data_wrapper->capacity)*sizeof(uint32_t));

    exitWithFailure(binned_data_wrapper->starts, "binned_data_wrapper->starts");

    binned_data_wrapper->ends = 
        realloc(binned_data_wrapper->ends, (binned_data_wrapper->capacity)*sizeof(uint32_t));

    exitWithFailure(binned_data_wrapper->ends, "binned_data_wrapper->ends");

}

void exitWithFailure(void * data_point_in, char* message) {
    if (data_point_in == NULL) {
        fprintf(stderr, "ERROR: Dynamic Memory allocation for %s failed!\n", message);
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
