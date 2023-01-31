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
#include "breakpoints.h"
#include "htslib/sam.h"
#include "terms.h"
#include "improperly_paired_reads.h"

// The function is for debugging purpose
//
void findDebugPoint();

/**
 * This function is used to process one chunk of aligned reads from bam file and process it using current thread
 * @param user_inputs: variable that contains all the user input info
 * @param rec: the individual alignment record to be processed
 * @param header: bam/cram file header information that has chromosome id and length info 
 * @param tmp_stats_info: variable used to store all the statistical information regarding bases and reads
 * @param primary_chromosome_hash: handle primary chromosomes only if it is not NULL
 * @param chrom_tracking, a variable used to track count information at each chromosome position
 * @param chrom_index, a index points to the current chromosome index at the chromo_tracking array variable
 */
void processCurrentRecord(User_Input *user_inputs, bam1_t *rec, bam_hdr_t *header, Stats_Info *tmp_stats_info, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Breakpoint_Array *breakpoint_array, Not_Properly_Paired_Reads_Array* improperly_paired_reads_array, khash_t(khStrInt) *unmapped_read_hash);

/**
 * This function is used to process individual aligned read and put the results into a chrom_tracking variable
 * @param user_inputs: variable that contains all the user input info
 * @param tmp_stats_info: variable used to store all the statistical information regarding bases and reads
 * @param rec: the individual alignment record to be processed
 * @param chrom_tracking, a variable used to track count information at each chromosome position
 * @param chrom_index, a index points to the current chromosome index at the chromo_tracking array variable
 */
void processRecord(User_Input *user_inputs, Stats_Info *tmp_stats_info, bam1_t *rec, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Breakpoint_Array *breakpoint_array);

bool getOverlapInfo(User_Input *user_inputs, Stats_Info *stats_info, bam1_t *rec, uint32_t *m_pos_r_end);

void calculateMeanAndStdev(Binned_Data_Wrapper **binned_data_wraper, Simple_Stats *the_stats, Chromosome_Tracking *chrom_tracking);

int compareDouble(const void * val1, const void * val2);

double CalcualtePercentile(DoubleArray *array_in, int percentile);

#endif // STATS_H
