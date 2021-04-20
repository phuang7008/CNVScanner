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
#include "stats.h"
#include "terms.h"
#include "user_inputs.h"
#include "utils.h"

// Need to declaration a global variable defined in data_structure.h
// For khash: key -> string (char*);    value -> int
//
int khStrInt = 34;


int main(int argc, char *argv[]) {
    // get user input options and then processing it accordingly
    //
    User_Input *user_inputs = userInputInit();
    processUserOptions(user_inputs, argc, argv);
    //outputUserInputOptions(user_inputs);

    // now for the bam/cram file open it for read
    //
    samFile *sfd = sam_open(user_inputs->bam_file, "r");
    if (sfd == 0) {
        fprintf(stderr, "ERROR: Cannot open file \n%s\n", user_inputs->bam_file);
        return -1;
    }

    // Set the reference if it is the cram file
    //
    char * fn_ref = 0;
    if (user_inputs->reference_file) {
        fn_ref = getReferenceFaiPath(user_inputs->reference_file);

        if (hts_set_fai_filename(sfd, fn_ref) != 0) {
            fprintf(stderr, "ERROR: Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
            if (fn_ref) free(fn_ref);
            return -1;
        }
    } else {
        if ( sfd->is_cram || sfd->format.format == cram ) {
            fprintf(stderr, "ERROR: Please provide the reference sequences for the input CRAM file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs->bam_file, sfd);

    // use sam_hdr_read to process both bam and cram header
    //
    bam_hdr_t *header=NULL;
    if ((header = sam_hdr_read(sfd)) == 0) return -1;

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
    uint32_t i, j;
    Binned_Data_Wrapper **binned_data_wrapper = calloc(chrom_tracking->number_of_chromosomes, sizeof(Binned_Data_Wrapper*));
    checkMemoryAllocation(binned_data_wrapper, "Binned_Data_Wrapper **binned_data_wrapper");
    binnedDataWrapperInit(binned_data_wrapper, chrom_tracking);
    //binnedDataWrapperInit(gc_binned_data_wrapper, chrom_tracking);
    //binnedDataWrapperInit(mappability_gc_binned_data_wrapper, chrom_tracking);

    // can't set to be static as openmp won't be able to handle it
    // check the bam/cram file size first
    //
    uint64_t input_bam_file_size = check_file_size(user_inputs->bam_file);

    uint32_t total_chunk_of_reads = 200000;     // for small bam/cram file
    if (input_bam_file_size > 5000000000)       // anything > 5Gb
        total_chunk_of_reads = 1000000;         // Good for 3 threads with 9gb  of memory

    // try to allocate the bam1_t array here for each thread, so they don't have to create and delete the array at each loop
    // Here there are two layers of read_buff are created through calloc(), which need to be freed at the end for two levels
    //
    Read_Buffer *read_buff = calloc(user_inputs->num_of_threads, sizeof(Read_Buffer));
    for (i=0; i<user_inputs->num_of_threads; i++) {
        read_buff[i].chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t*));
        for (j=0; j<total_chunk_of_reads; j++)
            read_buff[i].chunk_of_reads[j]=NULL;

        read_buff[i].capacity = total_chunk_of_reads;
        read_buff[i].size = 0;
    }

    // The following is for debugging purpose
    //
    findDebugPoint();

    fflush(stdout);

    // set random seed (should only be called ONCE)
    //
    if(user_inputs->percentage < 1.0)
        srand((uint32_t)time(NULL));    // set random seed

    // now let's do the parallelism
    //
    uint64_t total_reads=0;
    while(chrom_tracking->more_to_read) {
#pragma omp parallel shared(read_buff, chrom_tracking) num_threads(user_inputs->num_of_threads)
      {
        //int num_of_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
        readBufferInit(&read_buff[thread_id]);      // third level of memory at created through bam_init1()
        uint32_t num_records = 0;
        //printf("Before File Reading: number of threads is %d and current thread id is %d\n", num_of_threads, thread_id);

#pragma omp critical 
        {
          // this part of the code need to run atomically, that is only one thread should allow access to read
          num_records = readBam(sfd, header, chrom_tracking, &read_buff[thread_id]);
          read_buff[thread_id].size = num_records;
          total_reads += num_records;
        }

        Stats_Info *tmp_stats_info = calloc(1, sizeof(Stats_Info));
        statsInfoInit(tmp_stats_info);
        //Coverage_Stats *cov_stats = calloc(1, sizeof(Coverage_Stats)); 
        //coverageStatsInit(cov_stats);
        khash_t(str) *coverage_hash = kh_init(str);     // hash_table using string as key

        // can not use the else {} with the previous if {} block, otherwise the barrier will wait forever!
        //
        if (num_records > 0) {
          printf("Reading: %"PRIu32" records\t\tTotal: %"PRIu64"\t\tThread id: %d.\n", num_records, total_reads, thread_id);

          processBamChunk(user_inputs, tmp_stats_info, coverage_hash, header, &read_buff[thread_id],
                  target_buffer_status, thread_id, wanted_chromosome_hash, chrom_tracking->number_of_chromosomes);

        }

        // release the allocated chunk of buffer for aligned reads after they have been processed!
        //
        //printf("cleaning the read buffer hash for thread %d...\n\n", thread_id);
        readBufferDestroy(&read_buff[thread_id]);       // third level of memory allocation is destroyed and freed here

        if (num_records == 0) {
          printf("No more to read for thread %d !!!!!!!!!!!!\n", thread_id);
          if (chrom_tracking->more_to_read) chrom_tracking->more_to_read = false;
        }

#pragma omp critical
        {
          if (num_records > 0) {
            combineThreadResults(chrom_tracking, coverage_hash);
            combineCoverageStats(stats_info, tmp_stats_info);

            // if all reads have been processed for the entire file, we need to set the status to 2 for all
            // 
            if (!chrom_tracking->more_to_read) {
              for(i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++) {
                if (chrom_tracking->chromosome_status[i] == 1)
                  chrom_tracking->chromosome_status[i] = 2;
              }
            }
          }
        }

        cleanKhashStr(coverage_hash, 1);
        statsInfoDestroy(tmp_stats_info);
        tmp_stats_info = NULL;   // to eleminate dangling pointer

// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier 

#pragma omp single
        {
          if (num_records > 0) {
            for (i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++) {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                // check to see if any of the chromosomes has finished. If so, write the results out
                // for the whole genome, we need to use the file that contains regions of all Ns in the reference
                // As they will be not used, so we are going to set the count info in these regions to 0
                //
                printf("Zero all excluded regions\n");
                if (user_inputs->excluded_region_file)
                  zeroAllNsRegions(chrom_tracking->chromosome_ids[i], excluded_bed_info, chrom_tracking, target_buffer_status, -1);
              }
            }
          }
        }

#pragma omp barrier

        i = 0;
        while (i<(uint32_t)chrom_tracking->number_of_chromosomes) {

#pragma omp sections
          {
            if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2)
              printf("\n");

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("Thread %d is now producing coverage bins for chromosome %s\n", 
                        thread_id, chrom_tracking->chromosome_ids[i]);
                        
                // First, we need to find the index that is used to track current chromosome chrom_id
                //
                int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_tracking->chromosome_ids[i], chrom_tracking);
                coverageBinningWrapper(chrom_tracking, user_inputs, stats_info, binned_data_wrapper[chrom_idx], chrom_idx);
                if (user_inputs->debug_ON) outputBinnedData(binned_data_wrapper[chrom_idx], chrom_tracking->chromosome_ids[i]);
              }
            }

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                // For Whole Genome Annotations (use MySQL for the annotation)
                // if the annotation is not on, it will just output . . . . . . . )
                //
                /*if (user_inputs->wgs_coverage) {
                  printf("Thread %d is now writing WGS annotation for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
                  writeAnnotations(chrom_tracking->chromosome_ids[i], chrom_tracking, user_inputs, intronic_regions, db_exon_regions);

                  // if user specifies the range information (usually for graphing purpose), need to handle it here
                  //
                  printf("Thread %d is now writing coverage uniformity data for chromosome %s\n", 
                          thread_id, chrom_tracking->chromosome_ids[i]);
                  coverageRangeInfoForGraphing(chrom_tracking->chromosome_ids[i], chrom_tracking, user_inputs);
                }*/
              }
            }

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {

                /*if (TARGET_FILE_PROVIDED) {
                  printf("Thread %d is now calculating gene/transcript/cds coverage percentage for chromosome %s\n", 
                          thread_id, chrom_tracking->chromosome_ids[i]);

                  // process one annotation file at a time
                  //
                  int x;
                  for (x=0; x<user_inputs->num_of_target_files; x++) {
                    khash_t(khStrLCG) *transcript_hash = kh_init(khStrLCG);
                    khash_t(khStrStrArray) *gene_transcripts = kh_init(khStrStrArray);

                    if (USER_DEFINED_DATABASE) {
                        khash_t(khStrLCG) *user_defined_cds_gene_hash = kh_init(khStrLCG);

                        if (raw_user_defined_databases[x]) {
                          userDefinedGeneCoverageInit(user_defined_cds_gene_hash, 
                                  chrom_tracking->chromosome_ids[i], raw_user_defined_databases[x], gene_transcripts);
                        } else {
                          genePercentageCoverageInit(user_defined_cds_gene_hash, 
                                  chrom_tracking->chromosome_ids[i], dbs, gene_transcripts);
                          intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], 
                                  target_bed_info[x], chrom_tracking, user_defined_cds_gene_hash);
                        }

                        calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, user_inputs, user_defined_cds_gene_hash);
                        transcriptPercentageCoverageInit(transcript_hash, user_defined_cds_gene_hash);
                        storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], user_inputs, 
                                transcript_hash, gene_transcripts, hgmd_transcripts, gene_transcript_percentage_hash, x);

                        // clean-up
                        genePercentageCoverageDestroy(user_defined_cds_gene_hash);
                        user_defined_cds_gene_hash=NULL;
                    } else {
                        // Allocate memories for the low_cov_gene_hash, which get all its content from the MySQL database (official RefSeq DB)
                        // For calculating the percentage of gene bases with low coverge for capture only
                        // and use the intersect regions between refseq_cds_genes for official annotation and low_cov_genes for targets
                        //
                        khash_t(khStrLCG) *low_cov_gene_hash = kh_init(khStrLCG);

                        genePercentageCoverageInit(low_cov_gene_hash, chrom_tracking->chromosome_ids[i], dbs, gene_transcripts);
                        intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, low_cov_gene_hash);

                        calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, user_inputs, low_cov_gene_hash);
                        transcriptPercentageCoverageInit(transcript_hash, low_cov_gene_hash);

                        storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], user_inputs, 
                                transcript_hash, gene_transcripts, hgmd_transcripts, gene_transcript_percentage_hash, x);

                        // clean-up the memory space
                        genePercentageCoverageDestroy(low_cov_gene_hash);
                        low_cov_gene_hash=NULL;
                    }

                    // More clean-ups
                    genePercentageCoverageDestroy(transcript_hash);
                    cleanKhashStrStrArray(gene_transcripts);
                    transcript_hash=NULL;
                    gene_transcripts=NULL;
                  }
                }*/
              }
            }
          }

#pragma omp barrier 

#pragma omp single
          {
            if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("\n");

              // clean up the array allocated
              //
              if (chrom_tracking->coverage[i]) {
                free(chrom_tracking->coverage[i]);
                chrom_tracking->coverage[i] = NULL;
              }
              chrom_tracking->chromosome_status[i] = 3;
            }
            i++;
          }
        }

#pragma omp barrier 
      }     // End of while loop before flush for thread

      //printf("\n");
      fflush(stdout);
    }

    sam_close(sfd);

    // output report for debugging
    //
    reportStatsForDebugging(stats_info, user_inputs);
}
