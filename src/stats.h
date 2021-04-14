/*
 * =====================================================================================
 *
 *       Filename:  stats.h
 *
 *    Description:  The header file for sequence stats analysis
 *
 *        Version:  1.0
 *        Created:  02/24/2017 03:47:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef STATS_H
#define STATS_H

//#include <my_global.h>
//#include <mysql.h>
//#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"

// The function is for debugging purpose
//
void findDebugPoint();

/**
 * This function is used to process one chunk of aligned reads from bam file and process it using current thread
 * @param user_inputs: variable that contains all the user input info
 * @param cov_stats: variable used to store all the statistical information regarding bases and reads
 * @param coverage_hash: the hash table that is used to store temp calculation results
 * @param header: bam/cram file header information that has chromosome id and length info 
 * @param read_buff_in: the array of Read_Buffer with each read that will be processed by current function
 * @param target_buffer_status: it contains the target and buffer info for capture sequencing
 * @param thread_id: the current thread id
 * @param primary_chromosome_hash: handle primary chromosomes only if it is not NULL
 */
void processBamChunk(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in, Target_Buffer_Status *target_buffer_status, int thread_id, khash_t(khStrInt) * primary_chromosome_hash, int number_of_chromosomes);

/**
 * This function is used to process individual aligned read and put the results into a hash table where the key is the chromosome location
 * while the value is the count
 * @param user_inputs: variable that contains all the user input info
 * @param cov_stats: variable used to store all the statistical information regarding bases and reads
 * @param coverage_hash: the hash table used to store the coverage count
 * @param chrom_id: current chromosome id to be handled
 * @param rec: the individual alignment record to be processed
 * @param target_buffer_status: it contains the target and buffer info for capture sequencing
 * @param same_chr: a boolean value to indicate if chr id has been changed or not
 * @param iter_in_out: an khiter_t variable to store the point value for current chromosome id key
 */
void processRecord(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, char *chrom_id, bam1_t *rec, Target_Buffer_Status * target_buffer_status, bool same_chr, khiter_t *iter_in_out, int number_of_chromosomes);

/**
 * This function is used to combine each thread coverage results and put them onto a big array with length of each chromosome
 * This function should be critical as no more than one thread should write to the same array at the same time
 * @param chrom_tracking: the variable that is used to tracking the status of all chromosomes
 * @param coverage_hash: the coverage hash table whose contents will write into the chrom_coverage
 */
void combineThreadResults(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash);

bool getOverlapInfo(User_Input *user_inputs, Stats_Info *stats_info, bam1_t *rec, uint32_t *m_pos_r_end);

#endif // STATS_H
