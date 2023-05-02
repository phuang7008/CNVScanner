/*
 * =====================================================================================
 *
 *       Filename:  calculate_stdev.c
 *
 *    Description:  the detailed implementation of calculate_stdev.h
 *
 *        Version:  1.0
 *        Created:  08/30/2021 02:21:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, Tx
 *
 * =====================================================================================
 */

#include "calculate_stdev.h"

void OnePassCalculateSedev(User_Input *user_inputs, bam_hdr_t **header, hts_idx_t **sfh_idx, samFile **sfh, Bed_Info *excluded_bed_info,  Simple_Stats *simple_stats, Breakpoint_Array **breakpoint_array, khash_t(m32) **anchor_breakpoints_hash_array, Not_Properly_Paired_Reads_Array** improperly_paired_reads_array, khash_t(khStrInt) *unmapped_read_hash) {

    // for tmp chromosome tracking
    //
    Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));

    // setup a variable to store chromosome info that user provided
    //
    khash_t(khStrInt) *wanted_chromosome_hash = kh_init(khStrInt);

    if (user_inputs->chromosome_bed_file != NULL) {
        loadWantedChromosomes(wanted_chromosome_hash, user_inputs->reference_version, user_inputs->chromosome_bed_file);
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, header[0]);

        checkNamingConvention(header[0], wanted_chromosome_hash);
    } else {
        loadGenomeInfoFromBamHeader(wanted_chromosome_hash, header[0], user_inputs->reference_version);
        chrom_tracking->number_of_chromosomes = header[0]->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, header[0]);
    }

    // setup the stats_info for each individual chromosome for each threads and initialize them
    //
    Stats_Info **stats_info_per_chr = calloc(chrom_tracking->number_of_chromosomes, sizeof(Stats_Info*));
    uint32_t chrom_index = 0;
    for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
        stats_info_per_chr[chrom_index] = calloc(1, sizeof(Stats_Info));
        statsInfoInit(stats_info_per_chr[chrom_index]);
    }

    // setup one pass stdev data-structure
    //
    OnePassStdev **one_pass_stdev = calloc(chrom_tracking->number_of_chromosomes, sizeof(OnePassStdev*));
    OnePassStdevInit(one_pass_stdev, chrom_tracking);

#pragma omp parallel shared(chrom_tracking) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
        {
            for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index)
            {
                if (wanted_chromosome_hash != NULL) {
                    khiter_t iter = kh_get(khStrInt, wanted_chromosome_hash, chrom_tracking->chromosome_ids[chrom_index]);
                    if (iter == kh_end(wanted_chromosome_hash)) {
                        continue;
                    }
                }
#pragma omp task
                {
                    int thread_id = omp_get_thread_num();
                    printf("In One PASS: current thread id: %d for chr %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);

                    // get the iterator for the current chromosome
                    //
                    hts_itr_t *iter = sam_itr_querys(sfh_idx[thread_id], header[thread_id], chrom_tracking->chromosome_ids[chrom_index]);
                    if (iter == NULL) {
                        fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", chrom_tracking->chromosome_ids[chrom_index]);
                        exit(EXIT_FAILURE);
                    }

                    // fetch the breakpoint array chromosome index here
                    //
                    uint32_t bpt_chr_idx = fetchBreakpointArrayChrIndex(breakpoint_array, chrom_tracking, chrom_index);
                    uint32_t imp_chr_idx = fetchImproperPRArrayChrIndex(improperly_paired_reads_array, chrom_tracking, chrom_index);

                    // create a lookup table for paired reads of each breakpoint
                    // key: read name (str), while value: the index in the breakpoint array for this chromosome
                    //
                    //khash_t(khStrInt) *breakpoint_pairs_hash = kh_init(khStrInt);

                    chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[chrom_index], chrom_index);

                    int result = 0;
                    bam1_t *b = bam_init1();
                    while ((result = sam_itr_next(sfh[thread_id], iter, b)) >= 0) {
                        processCurrentRecord(user_inputs, b, header[thread_id], stats_info_per_chr[chrom_index], chrom_tracking, chrom_index, breakpoint_array[bpt_chr_idx], one_pass_stdev[chrom_index], improperly_paired_reads_array[imp_chr_idx], unmapped_read_hash);
                    }

                    if (result < -1) {
                        fprintf(stderr, "ERROR: bam1_t read has error with less than -1 result: %d\n", result);
                        exit(EXIT_FAILURE);
                    } else if (result <0 && result > -1) {
                        fprintf(stderr, "ERROR: bam1_r read has 0-1 error %.2d\n", result);
                    }

                    bam_destroy1(b);
                    hts_itr_destroy(iter);

                    // process all breakpoints and generate anchor breakpoint hash
                    //
                    processBreakpoints(breakpoint_array[bpt_chr_idx], anchor_breakpoints_hash_array[chrom_index], user_inputs);

                    // output improperly paired reads for debugging
                    //
                    if (user_inputs->debug_ON)
                        outputGroupedImproperlyPairedReads(improperly_paired_reads_array, chrom_tracking);

                    if (user_inputs->excluded_region_file)
                        zeroAllExcludedRegions(chrom_tracking, chrom_index, excluded_bed_info);

                    // walk through each base for coverage info
                    //
                    get_coverage_info(chrom_tracking, chrom_index,  one_pass_stdev[chrom_index]);

                    // clean-up array 
                    //
                    if (chrom_tracking->coverage[chrom_index]) {
                        free(chrom_tracking->coverage[chrom_index]);
                        chrom_tracking->coverage[chrom_index] = NULL;
                    }

                    // output breakpoint array and paired_reads across breakpointsfor debugging
                    //
                    if (user_inputs->debug_ON) {
                        outputBreakpointArray(breakpoint_array[bpt_chr_idx]);
                    }
                } // omp task
            } // for loop
#pragma omp taskwait
        } // omp single
    } // omp parallel

    StdevCalculation(one_pass_stdev, chrom_tracking, simple_stats);

    for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
        // clean up the stats info per chromosome
        //
        if (stats_info_per_chr[chrom_index])
            statsInfoDestroy(stats_info_per_chr[chrom_index]);
    }
    if (stats_info_per_chr)
        free(stats_info_per_chr);

    // clean-up
    //
    OnePassStdevDestroy(one_pass_stdev, chrom_tracking);

    chromosomeTrackingDestroy(chrom_tracking);
    if (chrom_tracking)
        free(chrom_tracking);

    if (wanted_chromosome_hash != NULL)
        cleanKhashStrInt(wanted_chromosome_hash);
}

void OnePassStdevInit(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        one_pass_stdev[i] = calloc(1, sizeof(OnePassStdev));
        one_pass_stdev[i]->total_sum   = 0.0;
        one_pass_stdev[i]->total_bases = 0.0;
        one_pass_stdev[i]->sum_of_base_cov_square = 0.0;

        one_pass_stdev[i]->total_tlen_sum = 0;
        one_pass_stdev[i]->total_reads = 0;
        one_pass_stdev[i]->total_reads_0_tlen = 0;
        one_pass_stdev[i]->sum_of_tlen_square = 0;
    }
}

void OnePassStdevDestroy(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (one_pass_stdev[i] != NULL)
            free(one_pass_stdev[i]);
    }

    if (one_pass_stdev)
        free(one_pass_stdev);
}

void StdevCalculation(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking, Simple_Stats *simple_stats) {
    uint32_t i;
    uint64_t total_tlen_sum=0, total_tlen_square_sum=0, total_tlen_reads=0, total_reads_0_tlen=0;
    double total_sum=0, total_square_sum=0, total_bases=0;

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        total_sum   += one_pass_stdev[i]->total_sum;
        total_bases += one_pass_stdev[i]->total_bases;
        total_square_sum += one_pass_stdev[i]->sum_of_base_cov_square;

        total_tlen_sum += one_pass_stdev[i]->total_tlen_sum;
        total_tlen_square_sum += one_pass_stdev[i]->sum_of_tlen_square;
        total_tlen_reads += one_pass_stdev[i]->total_reads;
        total_reads_0_tlen += one_pass_stdev[i]->total_reads_0_tlen;
    }

    // For overall base stats
    //
    simple_stats->average_coverage  = total_sum / total_bases;
    simple_stats->total_bases_used  = total_bases;
    simple_stats->stdev = sqrt(total_square_sum / total_bases - simple_stats->average_coverage * simple_stats->average_coverage);
    simple_stats->outlier_cutoff    = simple_stats->average_coverage + 3 * simple_stats->stdev;

    // For TLEN stats
    // but the reads with TLEN 0 is counted twice, so we need to adjust it here
    //
    total_tlen_reads = total_tlen_reads + (total_reads_0_tlen / 2);
    simple_stats->tlen_mean = (double)total_tlen_sum / (double)total_tlen_reads;
    simple_stats->tlen_total_reads_used = total_tlen_reads;
    simple_stats->tlen_stdev = sqrt((double)total_tlen_square_sum / (double)total_tlen_reads - \
                                            simple_stats->tlen_mean * simple_stats->tlen_mean);
    simple_stats->total_reads_w_0_tlen = total_reads_0_tlen;
    simple_stats->tlen_outlier_cutoff = simple_stats->tlen_mean + 3 * simple_stats->tlen_stdev;

    fprintf(stderr, "After one pass calculation of mean and stdev\n");
    fprintf(stderr, "Average Coverage: %.2f\n", simple_stats->average_coverage);
    fprintf(stderr, "Standard Deviation: %.2f\n", simple_stats->stdev);
    fprintf(stderr, "Outlier Cutoff: %.2f\n", simple_stats->outlier_cutoff);
    fprintf(stderr, "Total bases Used: %"PRIu32"\n\n", simple_stats->total_bases_used);

    fprintf(stderr, "After one pass calculation of TLEN mean and TLEN stdev\n");
    fprintf(stderr, "Total reads Used: %"PRIu32"\n", simple_stats->tlen_total_reads_used);
    fprintf(stderr, "TLEN Mean: %.2f\n", simple_stats->tlen_mean);
    fprintf(stderr, "TLEN StDev: %.2f\n", simple_stats->tlen_stdev);
    fprintf(stderr, "TLEN TLEN outlier cutoff (3 * StDev + TLEN_Mean): %.2f\n", simple_stats->tlen_outlier_cutoff);
    fprintf(stderr, "Total reads with TLEN = 0: %"PRIu32"\n\n", simple_stats->total_reads_w_0_tlen);
}

void SimpleStatsInit(Simple_Stats *simple_stats) {
    simple_stats->average_coverage = 0.0;
    simple_stats->stdev = 0.0;
    simple_stats->total_bases_used = 0;

    simple_stats->tlen_mean = 0.0;
    simple_stats->tlen_stdev = 0.0;
    simple_stats->tlen_total_reads_used = 0;
    simple_stats->total_reads_w_0_tlen = 0;
    simple_stats->tlen_outlier_cutoff = 0.0;
}

void get_coverage_info(Chromosome_Tracking *chrom_tracking, int32_t chrom_idx, OnePassStdev *one_pass_stdev) {
    uint32_t i=0;
    for (i=0; i<chrom_tracking->chromosome_lengths[chrom_idx]; i++) {

        if (chrom_tracking->coverage[chrom_idx][i] != -1) {
            // remove bases with coverage extremely high
            //
            //if (chrom_tracking->coverage[chrom_idx][i] >= 1000)
            if (chrom_tracking->coverage[chrom_idx][i] >= 500)
                //chrom_tracking->coverage[chrom_idx][i] = 100;
                continue;

            one_pass_stdev->total_sum += (double) chrom_tracking->coverage[chrom_idx][i];
            one_pass_stdev->total_bases++;
            one_pass_stdev->sum_of_base_cov_square += (double) chrom_tracking->coverage[chrom_idx][i] * chrom_tracking->coverage[chrom_idx][i];
        }
    }
}
