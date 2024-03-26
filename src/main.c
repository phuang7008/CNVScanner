/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  The main() of the WGS CNV program
 *
 *        Version:  1.0
 *        Created:  03/31/2021 10:09:05 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include "htslib/include/htslib/sam.h"

#include "data_structure.h"
#include "coverage_tracking.h"
#include "utility.h"

#include "terms.h"
#include "storage.h"

#include "analysis.h"
#include "breakpoints.h"
#include "cnv_segmentation.h"
#include "improperly_paired_reads.h"
#include "calculate_stdev.h"
#include "excluded_regions.h"
#include "fileProcessing.h"
#include "reports.h"
#include "stats.h"
#include "std_cnv_calling.h"
#include "user_inputs.h"
#include "utils.h"

// Need to declaration a global variable defined in data_structure.h
// For khash: key -> string (char*);    value -> int
//
int khStrInt = 34;

// used for KHASH_MAP_INIT_INT(kh32, char*)
//
int khIntStr = 35;

// used for KHASH_MAP_INIT_STR(khIntPrArray, Paired_Reads_Across_A_Breakpoint*)
//
int khIntPrArray = 36;

int main(int argc, char *argv[]) {
    // get user input options and then processing it accordingly
    //
    User_Input *user_inputs = userInputInit();
    processUserOptions(user_inputs, argc, argv);
    //outputUserInputOptions(user_inputs);

    // now for the bam/cram file open it for read for multi-threading
    //
    int t;
    samFile **sfh = calloc(user_inputs->num_of_threads, sizeof(samFile*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        sfh[t] = sam_open(user_inputs->bam_file, "r");

        if (sfh[t] == 0) {
            fprintf(stderr, "ERROR: Cannot open file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // since we are going to handle one chromosome per thread, we need to get the index file
    //
    hts_idx_t **sfh_idx = calloc(user_inputs->num_of_threads, sizeof(hts_idx_t*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        sfh_idx[t] = sam_index_load(sfh[t], user_inputs->bam_file);
        if (sfh_idx[t] == NULL) {
            fprintf(stderr, "ERROR: Can't locate the index file\n");
            return -1;
        }
    }

    // Set the reference if it is the cram file
    //
    char * fn_ref = 0;
    if (user_inputs->reference_file) {
        fn_ref = getReferenceFaiPath(user_inputs->reference_file);

        for (t=0; t<user_inputs->num_of_threads; t++) {
            if (hts_set_fai_filename(sfh[t], fn_ref) != 0) {
                fprintf(stderr, "ERROR: Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
                if (fn_ref) free(fn_ref);
                return -1;
            }
        }
    } else {
        fprintf(stderr, "ERROR: Please provide the reference sequences for the input BAM/CRAM/SAM file \n%s\n", user_inputs->bam_file);
        return -1;
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs->bam_file, sfh[0]);

    // use sam_hdr_read to process both bam and cram headers
    //
    bam_hdr_t **headers = calloc(user_inputs->num_of_threads, sizeof(bam_hdr_t*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        if ((headers[t] = sam_hdr_read(sfh[t])) == 0) return -1;
    }

    // for the overall stats_info
    //
    Stats_Info *stats_info = calloc(1, sizeof(Stats_Info));
    statsInfoInit(stats_info);

    // setup a tracking variable to track chromosome working status
    //
    Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));

    // setup a variable to store chromosomes that specified by the user
    //
    khash_t(khStrInt) *wanted_chromosome_hash = kh_init(khStrInt);

    // because hash keys are not in order, therefore, I need to store the chromosome ids in an array
    // to make them the same order as those in bam/cram file
    //
    if (user_inputs->chromosome_bed_file != NULL) {
        stats_info->wgs_cov_stats->total_genome_bases = loadWantedChromosomes(wanted_chromosome_hash, 
                user_inputs->reference_version, user_inputs->chromosome_bed_file);
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, headers[0]);

        // here we need to verify if the chromosome naming convention matches 
        // between the bam/cram file and the chromosome bed file specified by the end user
        //
        checkNamingConvention(headers[0], wanted_chromosome_hash);

    } else {
        stats_info->wgs_cov_stats->total_genome_bases = 
            loadGenomeInfoFromBamHeader(wanted_chromosome_hash, headers[0], user_inputs->reference_version);
        chrom_tracking->number_of_chromosomes = headers[0]->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, headers[0]);
    }

    fprintf(stderr, "The total genome bases is %"PRIu64"\n", stats_info->wgs_cov_stats->total_genome_bases);

    // setup the stats_info for each individual chromosome for each threads and initialize them
    //
    Stats_Info **stats_info_per_chr = calloc(chrom_tracking->number_of_chromosomes, sizeof(Stats_Info*));
    uint32_t chrom_index = 0;
    for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
        stats_info_per_chr[chrom_index] = calloc(1, sizeof(Stats_Info));
        statsInfoInit(stats_info_per_chr[chrom_index]);
    }

    // For the excluded region bed file
    //
    Bed_Info *excluded_bed_info=NULL;

    if (user_inputs->excluded_region_file) {
        excluded_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, excluded_bed_info, stats_info,
                wanted_chromosome_hash, user_inputs->excluded_region_file); 
        fprintf(stderr, "The Ns base is %"PRIu32"\n", stats_info->wgs_cov_stats->total_excluded_bases);
    }

    setupOutputReportFiles(user_inputs);

    // need to setup data struture to store the raw binned regions
    //
    Binned_Data_Wrapper **binned_data_wrappers = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(binned_data_wrappers, "Binned_Data_Wrapper **binned_data_wrappers");
    binnedDataWrapperInit(binned_data_wrappers, chrom_tracking);

    // setup data structure to store the equal-sized window bins
    //
    Binned_Data_Wrapper **equal_size_window_wrappers = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(equal_size_window_wrappers, "Binned_Data_Wrapper **equal_size_window_wrappers");
    binnedDataWrapperInit(equal_size_window_wrappers, chrom_tracking);

    // setup the Breakpoint_Array and Paired_Reads_Across_Breakpoints_Array
    //
    Breakpoint_Array **breakpoint_array = calloc(chrom_tracking->number_of_chromosomes, sizeof(Breakpoint_Array));
    BreakpointArrayInit(breakpoint_array, chrom_tracking);

    // not properly paired/aligned reads
    //
    Not_Properly_Paired_Reads_Array** improperly_PR_array = calloc(chrom_tracking->number_of_chromosomes, sizeof (Not_Properly_Paired_Reads_Array*));
    NotProperlyPairedReadsInit(improperly_PR_array, chrom_tracking);

    // store anchor breakpoints into a hash array (array is for each chromosome)
    //
    khash_t(m32) **anchor_breakpoints_hash_array = calloc(chrom_tracking->number_of_chromosomes, sizeof(khash_t(m32)));
    AnchorBreakpointsHashArrayInit(anchor_breakpoints_hash_array, chrom_tracking);

    // calculate the whole genome base coverage mean and standard deviation
    //
    Simple_Stats *wgs_simple_stats = calloc(1, sizeof(Simple_Stats));
    SimpleStatsInit(wgs_simple_stats);
    OnePassCalculateSedev(user_inputs, headers, sfh_idx, sfh, excluded_bed_info, wgs_simple_stats, breakpoint_array, anchor_breakpoints_hash_array, improperly_PR_array);

    // The following is for debugging purpose
    //
    if (user_inputs->debug_ON) {
        findDebugPoint();
        forDebug();

        //removeDebugFiles(user_inputs);
    }

    fflush(stdout);

    // set random seed (should only be called ONCE)
    //
    if(user_inputs->percentage < 1.0)
        srand((uint32_t)time(NULL));    // set random seed

#pragma omp parallel shared(chrom_tracking) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index)
        {
            // need to check if we need to process current chromosome
            //
            if (wanted_chromosome_hash != NULL) {
              khiter_t iter = kh_get(khStrInt, wanted_chromosome_hash, chrom_tracking->chromosome_ids[chrom_index]);
              if (iter == kh_end(wanted_chromosome_hash)) {
                // chrom_id is not one of the primary chromosomes, so skip it!
                //
                continue;
              }
            }
#pragma omp task
          {
            //int num_of_threads = omp_get_num_threads();
            int thread_id = omp_get_thread_num();
            printf("In Main: current thread id: %d for chr %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);

            // get the iterator for the current chromosome
            //
            hts_itr_t *iter = sam_itr_querys(sfh_idx[thread_id], headers[thread_id], chrom_tracking->chromosome_ids[chrom_index]);
            if (iter == NULL) {
                fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", chrom_tracking->chromosome_ids[chrom_index]);
                exit(EXIT_FAILURE);
            }

            chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[chrom_index], chrom_index);

            int result=0;
            bam1_t *b = bam_init1();
            while ((result = sam_itr_next(sfh[thread_id], iter, b)) >= 0) {
              processCurrentRecord(user_inputs, b, headers[thread_id], stats_info_per_chr[chrom_index], chrom_tracking, chrom_index, NULL, NULL, NULL);
            }

            if (result < -1) {
                fprintf(stderr, "ERROR: bam1_t read has error with less than -1 result: %d\n", result);
                exit(EXIT_FAILURE);
            } else if (result <0 && result > -1) {
                fprintf(stderr, "ERROR: bam1_r read has 0-1 error %.2d\n", result);
            }

            bam_destroy1(b);
            hts_itr_destroy(iter);

            //printf("Set all excluded region bases to -1\n");
            if (user_inputs->excluded_region_file)
                zeroAllExcludedRegions(chrom_tracking, chrom_index, excluded_bed_info);

            // now need to do the binning
            //
            printf("Thread %d is conducting binning for chr %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);

            coverageBinningWrapper(chrom_tracking, user_inputs, stats_info, binned_data_wrappers[chrom_index], chrom_index, wgs_simple_stats, thread_id);
            //if (user_inputs->debug_ON)
            //    outputBinnedData(binned_data_wrappers[chrom_index], user_inputs, 1, chrom_tracking->chromosome_ids[chrom_index]);

            // clean-up array
            //
            if (chrom_tracking->coverage[chrom_index]) {
              free(chrom_tracking->coverage[chrom_index]);
              chrom_tracking->coverage[chrom_index] = NULL;
            }

            // without normalization section
            //
            uint32_t i;
            for (i = 0; i < binned_data_wrappers[chrom_index]->size; i++) {
                binned_data_wrappers[chrom_index]->data[i].ave_cov_gc_normalized = binned_data_wrappers[chrom_index]->data[i].ave_coverage;
                binned_data_wrappers[chrom_index]->data[i].ave_cov_map_gc_normalized = binned_data_wrappers[chrom_index]->data[i].ave_coverage;

                binned_data_wrappers[chrom_index]->data[i].weighted_mappability = 1;
                binned_data_wrappers[chrom_index]->data[i].weighted_gc_scale = 1;
            }

            // create equal-sized-bins
            //
            int total_lines = processFile(chrom_tracking->chromosome_ids[chrom_index], \
                                            user_inputs->equal_size_window_file, equal_size_window_wrappers[chrom_index]);

            generateEqualSizedBins(user_inputs, binned_data_wrappers[chrom_index], \
                    equal_size_window_wrappers[chrom_index],  total_lines, chrom_tracking->chromosome_ids[chrom_index]);

          } // omp task
        } // for loop
#pragma omp taskwait

        fflush(stdout);
      } // omp single
    } // omp parallel

    for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
        combineCoverageStats(stats_info, stats_info_per_chr[chrom_index]);

        if (stats_info_per_chr[chrom_index])
            statsInfoDestroy(stats_info_per_chr[chrom_index]);
    }

    // output report for debugging
    //
    reportStatsForDebugging(stats_info, user_inputs);

    if (user_inputs->debug_ON) {
        // output final normalized binned results
        //
        outputFinalBinnedData(binned_data_wrappers, user_inputs, chrom_tracking, 1);

        // output final equal sized window 
        //
        outputFinalBinnedData(equal_size_window_wrappers, user_inputs, chrom_tracking, 2);
    }

    // calculate the statistics here for WGS
    //
    Simple_Stats *the_stats = calloc(1, sizeof(Simple_Stats));
    calculateMeanAndStdev(equal_size_window_wrappers, the_stats, chrom_tracking);
    fprintf(stderr, "Mean:  %.2f\n", the_stats->average_coverage);
    fprintf(stderr, "Median:  %.2f\n", the_stats->median);
    fprintf(stderr, "Stdev: %.2f\n", the_stats->stdev);
    fprintf(stderr, "99_percentile: %.2f\n", the_stats->ninty_nine_percentile);
    fprintf(stderr, "98_percentile: %.2f\n", the_stats->ninty_eight_percentile);
    fprintf(stderr, "Z Score: %.2f\n", the_stats->zScore);
    fprintf(stderr, "Haploid cutoff: %.2f\n", the_stats->average_coverage - the_stats->zScore);
    fprintf(stderr, "Duplicate cutoff: %.2f\n\n", the_stats->average_coverage + the_stats->zScore);

    // once the median is calculated, we need to obtain the log2ratio for each equal-window bin
    // and output to a log2R file
    //
    calculateLog2Ratio(equal_size_window_wrappers, the_stats, chrom_tracking, user_inputs);

    // for Equal bin window CNV Calls
    //
    CNV_Array **equal_bin_cnv_array = calloc(chrom_tracking->number_of_chromosomes, sizeof(CNV_Array*));
    checkMemoryAllocation(equal_size_window_wrappers, "CNV_Array **equal_bin_cnv_array");
    cnvArrayInit(equal_bin_cnv_array, chrom_tracking);

    // merge and expand the CNV calls using raw bin data
    //
    generateCNVs(equal_bin_cnv_array, equal_size_window_wrappers, binned_data_wrappers, anchor_breakpoints_hash_array, improperly_PR_array, chrom_tracking, the_stats, stats_info,  user_inputs, headers, sfh_idx, sfh);

    // Create segmentation storage
    //
    Segment_Array **segment_arrays = calloc(chrom_tracking->number_of_chromosomes, sizeof (Segment_Array*));
    checkMemoryAllocation(segment_arrays, "Binned_Data_Wrapper **binned_data_wrappers");
    segmentArrysInit(segment_arrays, chrom_tracking);

    // Do segmentation
    //
    cnvSegmentation(chrom_tracking, segment_arrays, user_inputs);

    Segmented_CNV_Array **seg_cnv_array = calloc(chrom_tracking->number_of_chromosomes, sizeof(Segmented_CNV_Array*));
    checkMemoryAllocation(seg_cnv_array, "CNV_Array **final_cnv_array");
    SegmentedCNVArrayInit(seg_cnv_array, chrom_tracking);
    storeSegmentsLocallyAndInit(segment_arrays, seg_cnv_array, chrom_tracking);

    // find final CNVs after segmentation
    //
    processSegmentationData(equal_bin_cnv_array, seg_cnv_array, chrom_tracking, anchor_breakpoints_hash_array, user_inputs, the_stats, stats_info);

    fprintf(stderr, "After all processed. Now clean up everything\n");

    // clean up
    //
    if (the_stats) free(the_stats);
    printf("step 00\n");
    if (wgs_simple_stats) free(wgs_simple_stats);
    printf("step 0\n");
    cnvArrayDestroy(equal_bin_cnv_array, chrom_tracking->number_of_chromosomes);
    printf("step 1\n");
    segmentArrysDestroy(segment_arrays, chrom_tracking);
    printf("step 2\n");
    SegmentedCNVArrayDestroy(seg_cnv_array, chrom_tracking);
    printf("step 3\n");
    
    binnedDataWrapperDestroy(binned_data_wrappers, chrom_tracking);
    printf("step 4\n");
    binnedDataWrapperDestroy(equal_size_window_wrappers, chrom_tracking);
    printf("step 5\n");

    BreakpointArrayDestroy(breakpoint_array, chrom_tracking);
    //BreakpointStatsArrayDestroy(bpt_stats_array);
    printf("step 6\n");
    NotProperlyPairedReadsDestroy(improperly_PR_array, chrom_tracking);
    printf("step 7\n");
    AnchorBreakpointsHashArrayDestroy(anchor_breakpoints_hash_array, chrom_tracking);
    printf("step 8\n");

    if (excluded_bed_info != NULL)
        cleanBedInfo(excluded_bed_info);

    printf("step 9\n");
    chromosomeTrackingDestroy(chrom_tracking);
    printf("step 10\n");
    if (chrom_tracking)
        free(chrom_tracking);
    printf("step 11\n");

    if (stats_info)
        statsInfoDestroy(stats_info);
    printf("step 12\n");

    for (t=0; t<user_inputs->num_of_threads; t++) {
        bam_hdr_destroy(headers[t]);
        hts_idx_destroy(sfh_idx[t]);
        sam_close(sfh[t]);
    }
    printf("step 13\n");

    if (fn_ref) free(fn_ref);
    printf("step 14\n");

    if (wanted_chromosome_hash != NULL)
        cleanKhashStrInt(wanted_chromosome_hash);
    printf("step 15\n");

    userInputDestroy(user_inputs);
    printf("step 16\n");

    return 0;
}
