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
#include "excluded_regions.h"
#include "fileProcessing.h"
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


int main(int argc, char *argv[]) {
    // get user input options and then processing it accordingly
    //
    User_Input *user_inputs = userInputInit();
    processUserOptions(user_inputs, argc, argv);
    //outputUserInputOptions(user_inputs);

    // now for the bam/cram file open it for read
    //
    samFile *sfh = sam_open(user_inputs->bam_file, "r");
    if (sfh == 0) {
        fprintf(stderr, "ERROR: Cannot open file \n%s\n", user_inputs->bam_file);
        return -1;
    }

    // since we are going to handle one chromosome per thread, we need to get the index file
    //
    hts_idx_t *sfh_idx = sam_index_load(sfh, user_inputs->bam_file);
    if (sfh_idx == NULL) {
        fprintf(stderr, "ERROR: Can't locate the index file\n");
        return -1;
    }

    // Set the reference if it is the cram file
    //
    char * fn_ref = 0;
    if (user_inputs->reference_file) {
        fn_ref = getReferenceFaiPath(user_inputs->reference_file);

        if (hts_set_fai_filename(sfh, fn_ref) != 0) {
            fprintf(stderr, "ERROR: Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
            if (fn_ref) free(fn_ref);
            return -1;
        }
    } else {
        if ( sfh->is_cram || sfh->format.format == cram ) {
            fprintf(stderr, "ERROR: Please provide the reference sequences for the input CRAM file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs->bam_file, sfh);

    // use sam_hdr_read to process both bam and cram header
    //
    bam_hdr_t *header=NULL;
    if ((header = sam_hdr_read(sfh)) == 0) return -1;

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
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, header);

        // here we need to verify if the chromosome naming convention matches 
        // between the bam/cram file and the chromosome bed file specified by the end user
        //
        checkNamingConvention(header, wanted_chromosome_hash);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit2(target_buffer_status, wanted_chromosome_hash);
    } else {
        stats_info->wgs_cov_stats->total_genome_bases = 
            loadGenomeInfoFromBamHeader(wanted_chromosome_hash, header, user_inputs->reference_version);
        chrom_tracking->number_of_chromosomes = header->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, header);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit(target_buffer_status, header);
    }

    fprintf(stderr, "The total genome bases is %"PRIu64"\n", stats_info->wgs_cov_stats->total_genome_bases);

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
    Binned_Data_Wrapper **binned_data_wrapper = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(binned_data_wrapper, "Binned_Data_Wrapper **binned_data_wrapper");
    binnedDataWrapperInit(binned_data_wrapper, chrom_tracking);

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
        uint32_t chrom_index = 0;
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
            hts_itr_t *iter = sam_itr_querys(sfh_idx, header, chrom_tracking->chromosome_ids[chrom_index]);
            if (iter == NULL) {
                fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", chrom_tracking->chromosome_ids[chrom_index]);
                exit(EXIT_FAILURE);
            }

            chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[chrom_index], chrom_index);

            bam1_t *b = bam_init1();
            while (sam_itr_next(sfh, iter, b) >= 0) {
              processCurrentRecord(user_inputs, b, header, stats_info, chrom_tracking, chrom_index);
            }
            bam_destroy1(b);
            hts_itr_destroy(iter);

            //printf("Set all excluded region bases to -1\n");
            if (user_inputs->excluded_region_file)
            zeroAllNsRegions(chrom_tracking->chromosome_ids[chrom_index], excluded_bed_info, chrom_tracking, target_buffer_status, -1);

            // now need to do the binning
            //
            printf("Thread %d is conducting binning for chr %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);

            //AllStartsEndsArray *all_starts_ends_array = calloc(1, sizeof(AllStartsEndsArray));
            //all_starts_ends_array->capacity = 100000;
            //all_starts_ends_array->array = calloc(all_starts_ends_array->capacity, sizeof(uint32_t));
            //all_starts_ends_array->size = 0;

            //khash_t(khIntStr) *binned_starts  = kh_init(khIntStr);      // key: start, value: "start end length ave_cov"
            //khash_t(khIntStr) *binned_ends    = kh_init(khIntStr);      // key: end,   value: "start end length ave_cov"

            coverageBinningWrapper(chrom_tracking, user_inputs, stats_info, binned_data_wrapper[chrom_index], chrom_index, thread_id);
            if (user_inputs->debug_ON)
                outputBinnedData(binned_data_wrapper[chrom_index], chrom_tracking->chromosome_ids[chrom_index], user_inputs);

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

              total_lines =
                    processFile(chrom_tracking->chromosome_ids[chrom_index], user_inputs->mappability_file, map_starts, map_ends);
              printf("The total number of mappability lines is %i\n", total_lines);
              outputHashTable(map_starts, 1, user_inputs);

              mappabilityGcNormalization(binned_data_wrapper[chrom_index], user_inputs, map_starts, map_ends, total_lines, 1);

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

              total_lines =
                  processFile(chrom_tracking->chromosome_ids[chrom_index], user_inputs->gc_content_file, gc_starts, gc_ends);

              printf("The total number of GC%% lines is %i\n", total_lines);
              outputHashTable(gc_starts, 2, user_inputs);

              mappabilityGcNormalization(binned_data_wrapper[chrom_index], user_inputs, gc_starts, gc_ends, total_lines, 2);

              // clean-up
              //
              cleanKhashIntStr(gc_starts);
              cleanKhashIntStr(gc_ends);
            }
          }
        }
        //printf("\n");
        fflush(stdout);
#pragma omp taskwait
      }
    }

    sam_close(sfh);
    hts_idx_destroy(sfh_idx);

    // output report for debugging
    //
    reportStatsForDebugging(stats_info, user_inputs);

    // output final normalized binned results
    //
    outputFinalBinnedData(binned_data_wrapper, user_inputs, chrom_tracking);

    // clean up
    //
    TargetBufferStatusDestroy(target_buffer_status, chrom_tracking->number_of_chromosomes);

    binnedDataWrapperDestroy(binned_data_wrapper, chrom_tracking);

    if (excluded_bed_info != NULL)
        cleanBedInfo(excluded_bed_info);

    chromosomeTrackingDestroy(chrom_tracking);
    if (chrom_tracking)
        free(chrom_tracking);

    if (stats_info)
        statsInfoDestroy(stats_info);

    bam_hdr_destroy(header);

    if (fn_ref) free(fn_ref);

    if (wanted_chromosome_hash != NULL)
        cleanKhashStrInt(wanted_chromosome_hash);

    userInputDestroy(user_inputs);

    return 0;
}
