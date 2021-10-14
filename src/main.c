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
#include "htslib/sam.h"

#include "data_structure.h"
#include "coverage_tracking.h"
#include "utility.h"

#include "analysis.h"
#include "breakpoints.h"
#include "calculate_stdev.h"
#include "excluded_regions.h"
#include "fileProcessing.h"
#include "reports.h"
#include "stats.h"
#include "terms.h"
#include "user_inputs.h"
#include "utils.h"

// Need to declaration a global variable defined in data_structure.h
// For khash: key -> string (char*);    value -> int
//
int khStrInt = 34;

// used for KHASH_MAP_INIT_INT(kh32, char*)
//
int khIntStr = 35;

// used for KHASH_MAP_INIT_STR(khIntPrArray, Paired_Reads_Cross_A_Breakpoint*)
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

        if (hts_set_fai_filename(sfh[0], fn_ref) != 0) {
            fprintf(stderr, "ERROR: Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
            if (fn_ref) free(fn_ref);
            return -1;
        }
    } else {
        if ( sfh[0]->is_cram || sfh[0]->format.format == cram ) {
            fprintf(stderr, "ERROR: Please provide the reference sequences for the input CRAM file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs->bam_file, sfh[0]);

    // use sam_hdr_read to process both bam and cram header
    //
    bam_hdr_t **header = calloc(user_inputs->num_of_threads, sizeof(bam_hdr_t*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        if ((header[t] = sam_hdr_read(sfh[t])) == 0) return -1;
    }

    // for the overall stats_info
    //
    Stats_Info *stats_info = calloc(1, sizeof(Stats_Info));
    statsInfoInit(stats_info);

    // setup a tracking variable to track chromosome working status
    //
    Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));

    // Target_Buffer_Status need to be set even though there is no target or Ns region specified
    // It is because one of the method processRecord() need chromosome lengths information to be set
    //
    Target_Buffer_Status *target_buffer_status = NULL;

    // setup a variable to store chromosomes that specified by the user
    //
    khash_t(khStrInt) *wanted_chromosome_hash = kh_init(khStrInt);

    // because hash keys are not in order, therefore, I need to store the chromosome ids in an array
    // to make them the same order as those in bam/cram file
    //
    if (user_inputs->chromosome_bed_file != NULL) {
        stats_info->wgs_cov_stats->total_genome_bases = loadWantedChromosomes(wanted_chromosome_hash, 
                user_inputs->reference_version, user_inputs->chromosome_bed_file);
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, header[0]);

        // here we need to verify if the chromosome naming convention matches 
        // between the bam/cram file and the chromosome bed file specified by the end user
        //
        checkNamingConvention(header[0], wanted_chromosome_hash);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit2(target_buffer_status, wanted_chromosome_hash);
    } else {
        stats_info->wgs_cov_stats->total_genome_bases = 
            loadGenomeInfoFromBamHeader(wanted_chromosome_hash, header[0], user_inputs->reference_version);
        chrom_tracking->number_of_chromosomes = header[0]->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, header[0]);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit(target_buffer_status, header[0]);
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
        processBedFiles(user_inputs, excluded_bed_info, stats_info, target_buffer_status,
                wanted_chromosome_hash, user_inputs->excluded_region_file, chrom_tracking->number_of_chromosomes); 
        fprintf(stderr, "The Ns base is %"PRIu32"\n", stats_info->wgs_cov_stats->total_excluded_bases);
    }

    setupOutputReportFiles(user_inputs);

    // need to setup data struture to store the binned regions
    //
    Binned_Data_Wrapper **binned_data_wrappers = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(binned_data_wrappers, "Binned_Data_Wrapper **binned_data_wrappers");
    binnedDataWrapperInit(binned_data_wrappers, chrom_tracking);

    // setup data structure to store the equal-sized window bins
    //
    Binned_Data_Wrapper **equal_size_window_wrappers = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(equal_size_window_wrappers, "Binned_Data_Wrapper **equal_size_window_wrappers");
    binnedDataWrapperInit(equal_size_window_wrappers, chrom_tracking);
    
    // setup the Breakpoint_Array, Paired_Reads_Cross_Breakpoints_Array and Breakpoint_Stats_Array
    //
    Breakpoint_Array *breakpoint_array = calloc(1, sizeof(Breakpoint_Array));
    BreakpointArrayInit(breakpoint_array, chrom_tracking);

    Paired_Reads_Cross_Breakpoints_Array *preads_x_bpts_array = calloc(1, sizeof(Paired_Reads_Cross_Breakpoints_Array));
    PairedReadsCrossBreakpointsArrayInit(preads_x_bpts_array, chrom_tracking);

    Breakpoint_Stats_Array *bpt_stats_array = calloc(1, sizeof(Breakpoint_Stats_Array));
    BreakpointStatsArrayInit(bpt_stats_array, chrom_tracking);

    // calculate the whole genome base coverage mean and standard deviation
    //
    Simple_Stats *wgs_simple_stats = calloc(1, sizeof(Simple_Stats));
    SimpleStatsInit(wgs_simple_stats);
    OnePassCalculateSedev(user_inputs, header, sfh_idx, sfh, excluded_bed_info, wgs_simple_stats, target_buffer_status, breakpoint_array, preads_x_bpts_array);

    // The following is for debugging purpose
    //
    if (user_inputs->debug_ON) {
        findDebugPoint();
        forDebug();

        removeDebugFiles(user_inputs);
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
            printf("Current thread id: %d\n", thread_id);

            // get the iterator for the current chromosome
            //
            hts_itr_t *iter = sam_itr_querys(sfh_idx[thread_id], header[thread_id], chrom_tracking->chromosome_ids[chrom_index]);
            if (iter == NULL) {
                fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", chrom_tracking->chromosome_ids[chrom_index]);
                exit(EXIT_FAILURE);
            }

            chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[chrom_index], chrom_index);

            bam1_t *b = bam_init1();
            while (sam_itr_next(sfh[thread_id], iter, b) >= 0) {
              processCurrentRecord(user_inputs, b, header[thread_id], stats_info_per_chr[chrom_index], chrom_tracking, chrom_index, NULL, 0, NULL);
            }
            bam_destroy1(b);
            hts_itr_destroy(iter);

            //printf("Set all excluded region bases to -1\n");
            if (user_inputs->excluded_region_file)
            zeroAllNsRegions(chrom_tracking->chromosome_ids[chrom_index], excluded_bed_info, chrom_tracking, target_buffer_status, -1);

            // now need to do the binning
            //
            printf("Thread %d is conducting binning for chr %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);

            coverageBinningWrapper(chrom_tracking, user_inputs, stats_info, binned_data_wrappers[chrom_index], chrom_index, wgs_simple_stats, thread_id);
            if (user_inputs->debug_ON)
                outputBinnedData(binned_data_wrappers[chrom_index], user_inputs, 1);

            // clean-up array
            //
            if (chrom_tracking->coverage[chrom_index]) {
              free(chrom_tracking->coverage[chrom_index]);
              chrom_tracking->coverage[chrom_index] = NULL;
            }

            // now handling mappability normalization
            //
            uint32_t total_lines;
            if (user_inputs->mappability_file) {
              khash_t(khIntStr) * map_starts = kh_init(khIntStr);
              khash_t(khIntStr) * map_ends   = kh_init(khIntStr);

              total_lines = processFile(chrom_tracking->chromosome_ids[chrom_index],
                                        user_inputs->mappability_file, map_starts, map_ends, NULL);
              printf("The total number of mappability lines is %i\n", total_lines);
              outputHashTable(map_starts, 1, user_inputs);

              mappabilityGcNormalization(binned_data_wrappers[chrom_index], 
                                         user_inputs, map_starts, map_ends, total_lines, 1);

              // clean-up
              //
              cleanKhashIntStr(map_starts);
              cleanKhashIntStr(map_ends);
            }

            // GC normalization section
            //
            if (user_inputs->gc_content_file) {
              khash_t(khIntStr) *gc_starts = kh_init(khIntStr);
              khash_t(khIntStr) *gc_ends   = kh_init(khIntStr);

              total_lines = processFile(chrom_tracking->chromosome_ids[chrom_index], 
                                        user_inputs->gc_content_file, gc_starts, gc_ends, NULL);

              printf("The total number of GC%% lines is %i\n", total_lines);
              outputHashTable(gc_starts, 2, user_inputs);

              mappabilityGcNormalization(binned_data_wrappers[chrom_index], user_inputs, gc_starts, gc_ends, total_lines, 2);

              // clean-up
              //
              cleanKhashIntStr(gc_starts);
              cleanKhashIntStr(gc_ends);
            }

            // create equal-sized-bins
            //
            khash_t(khIntStr) *window_starts = kh_init(khIntStr);
            khash_t(khIntStr) *window_ends   = kh_init(khIntStr);
            total_lines = processFile(chrom_tracking->chromosome_ids[chrom_index], user_inputs->equal_size_window, 
                                        window_starts, window_ends, equal_size_window_wrappers[chrom_index]);

            generateEqualSizedBins(user_inputs, binned_data_wrappers[chrom_index],
                                        equal_size_window_wrappers[chrom_index],  total_lines);

            // clean-up
            //
            cleanKhashIntStr(window_starts);
            cleanKhashIntStr(window_ends);

          }
        }
#pragma omp taskwait

        // now need to combine all the stats_info for the final results
        //
#pragma omp critical
        {
            for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
                combineCoverageStats(stats_info, stats_info_per_chr[chrom_index]);

                if (stats_info_per_chr[chrom_index])
                    statsInfoDestroy(stats_info_per_chr[chrom_index]);
            }
        }

        //printf("\n");
        fflush(stdout);
      }
    }

    // output report for debugging
    //
    reportStatsForDebugging(stats_info, user_inputs);

    // output final normalized binned results
    //
    outputFinalBinnedData(binned_data_wrappers, user_inputs, chrom_tracking, 1);

    // output final equal sized window 
    //
    outputFinalBinnedData(equal_size_window_wrappers, user_inputs, chrom_tracking, 2);

    // calculate the statistics here
    //
    Stats *the_stats = calloc(1, sizeof(Stats));
    calculateMeanAndStdev(user_inputs, equal_size_window_wrappers, the_stats, chrom_tracking);
    //calculateMeanAndStdev(user_inputs, binned_data_wrappers, the_stats, chrom_tracking);
    fprintf(stderr, "Mean:  %.2f\n", the_stats->mean);
    fprintf(stderr, "Stdev: %.2f\n", the_stats->stdev);
    fprintf(stderr, "99_percentile: %.2f\n", the_stats->ninty_nine_percentile);
    fprintf(stderr, "98_percentile: %.2f\n", the_stats->ninty_eight_percentile);
    if (the_stats) free(the_stats);

    // clean up
    //
    TargetBufferStatusDestroy(target_buffer_status, chrom_tracking->number_of_chromosomes);

    binnedDataWrapperDestroy(binned_data_wrappers, chrom_tracking);

    BreakpointArrayDestroy(breakpoint_array);
    PairedReadsCrossBreakpointsArrayDestroy(preads_x_bpts_array);
    BreakpointStatsArrayDestroy(bpt_stats_array);

    if (excluded_bed_info != NULL)
        cleanBedInfo(excluded_bed_info);

    chromosomeTrackingDestroy(chrom_tracking);
    if (chrom_tracking)
        free(chrom_tracking);

    if (stats_info)
        statsInfoDestroy(stats_info);

    for (t=0; t<user_inputs->num_of_threads; t++) {
        sam_close(sfh[t]);
        bam_hdr_destroy(header[t]);
        hts_idx_destroy(sfh_idx[t]);
    }

    if (fn_ref) free(fn_ref);

    if (wanted_chromosome_hash != NULL)
        cleanKhashStrInt(wanted_chromosome_hash);

    userInputDestroy(user_inputs);

    return 0;
}
