/*
 * =====================================================================================
 *
 *        Filename:        targets.c
 *
 *        Description:    The implementation of the target.h file
 *
 *      Version:        1.0
 *      Created:        02/06/2017 04:45:04 PM
 *      Revision:        none
 *      Compiler:        gcc
 *
 *      Author:            Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */

#include "excluded_regions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>        // for file access()
#include <time.h>
#include "utils.h"

void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, Stats_Info *stats_info, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name, char *ref_version) {
    // First, let's get the total number of lines(items or count) within the target file
    //
    bed_info->size = getLineCount(bedfile_name);

    // Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
    //
    bed_info->coords = calloc(bed_info->size, sizeof(Bed_Coords));
    
    // load target file or Ns bed file again and store the information (starts, stops and chromosome ids)
    //
    uint32_t total_size = loadBedFiles(user_inputs->reference_version, bedfile_name, bed_info, wanted_chromosome_hash, ref_version);
    //printf("total size is %"PRIu32"\n", total_size);

    // Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
    //
    if (stats_info) {
        generateBedBufferStats(bed_info, stats_info, wanted_chromosome_hash);

        // Here we need to check if the bed file is merged and uniqued by comparing two different ways of addition of bases
        //
        if (total_size != stats_info->wgs_cov_stats->total_excluded_bases) {
            printf("\nERROR: The excluded region bed file \n%s\n needs to be bedtools sorted, merged and uniqued!\n", bedfile_name);
            printf("\ttotal size: %"PRIu32"\n", total_size);
            printf("\ttotal excluded bases: %"PRIu32"\n", stats_info->wgs_cov_stats->total_excluded_bases);
            //printf("\tThe chromosome ids listed in the Ns-region file MUST appear in the chromosome input file (--chr_list option)!\n");
            exit(EXIT_FAILURE);
        }
    }
}

void outputForDebugging(Bed_Info *bed_info) {
    // Output stored coordinates for verification
    uint32_t i=0;
    for (i=0; i<bed_info->size; i++) {
        printf("%s\t%"PRIu32"\t%"PRIu32"\n", bed_info->coords[i].chrom_id, bed_info->coords[i].start, bed_info->coords[i].end);
    }
}

void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, khash_t(khStrInt)* wanted_chromosome_hash) {
    uint32_t i=0, j=0, chrom_len=0;
    char cur_chrom_id[50];
    strcpy(cur_chrom_id, "something");

    for (i = 0; i < bed_info->size; i++) {

        // skip if the chromosome is not going to be processed
        //
        if (wanted_chromosome_hash != NULL) {
            khiter_t iter_p = kh_get(khStrInt, wanted_chromosome_hash, bed_info->coords[i].chrom_id);
            if (iter_p == kh_end(wanted_chromosome_hash)) {
                continue;
            } else {
                chrom_len = kh_value(wanted_chromosome_hash, iter_p);
            }
        }

        for (j=bed_info->coords[i].start; j<bed_info->coords[i].end; j++) {
            if (j >= chrom_len) continue;

            if (j < bed_info->coords[i].end) {
                if (j >= chrom_len) continue;

                stats_info->wgs_cov_stats->total_excluded_bases += 1;

                if ((strcmp(bed_info->coords[i].chrom_id, "chrX") == 0) || (strcmp(bed_info->coords[i].chrom_id, "X") == 0))
                    stats_info->wgs_cov_stats->total_excluded_bases_on_chrX += 1;

                if ((strcmp(bed_info->coords[i].chrom_id, "chrY") == 0) || (strcmp(bed_info->coords[i].chrom_id, "Y") == 0))
                    stats_info->wgs_cov_stats->total_excluded_bases_on_chrY += 1;
            }
        }
    }
}

void zeroAllExcludedRegions(Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Bed_Info *excluded_bed_info) {
    uint32_t i=0, j=0;

    for (i=0; i<excluded_bed_info->size; i++) {
        if (strcmp(excluded_bed_info->coords[i].chrom_id, chrom_tracking->chromosome_ids[chrom_index]) == 0) {
            for (j=excluded_bed_info->coords[i].start; j<excluded_bed_info->coords[i].end; j++) {
                if (j>=chrom_tracking->chromosome_lengths[chrom_index]) continue;

                // set all excluded bases to -1
                //
                chrom_tracking->coverage[chrom_index][j] = -1;
            }
        }
    }
}
