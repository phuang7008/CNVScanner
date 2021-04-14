/*
 * =====================================================================================
 *
 *      Filename:       analysis.h
 *
 *      Description:    The header file for analysis methods
 *
 *      Version:        1.0
 *      Created:        02/06/2017 04:45:04 PM
 *      Revision:       none
 *      Compiler:       gcc
 *
 *      Author:         Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdbool.h>	// for bool return type
#include <ctype.h>		// for isdigit()
#include <errno.h>
#include "terms.h"

/*
 * This is a wrapper function to help generate the coverage range info using writeCoverageRanges()
 * This will not generate annotation information for the speed reason
 * @param chrom_id: the chromosome id to be handed
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param user_inputs: contains all the user_inputs options
 * @param stats_info: a variable that contains all reads and base information
 */
void coverageBinningWrapper(char *chrom_id, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info);


/*
 * This is used to generate average coverage information for a range of position based on different binning strategies
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_idx: the chromosome index 
 * @param user_inputs: contains all the user_inputs options
 * @param wgs_binned_coverage_fp: the file handle for depositing the binned data
 * */
void writeCoverageBins(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, uint16_t chrom_idx, User_Input *user_inputs, Stats_Info *stats_info, FILE *fh_binned_coverage);

/*
 * The following method is for debugging only
 */
void reportStatsForDebugging(Stats_Info *stats_info, User_Input *user_inputs);

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage);

#endif //UTILS_H
