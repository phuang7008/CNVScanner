/*
 * =====================================================================================
 *
 *      Filename:       reports.h
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
#ifndef REPORTS_H
#define REPORTS_H

#include <stdbool.h>	// for bool return type
#include <ctype.h>		// for isdigit()
#include <errno.h>
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"

/*
 * This is a wrapper function to help generate the coverage range info using writeCoverageRanges()
 * This will not generate annotation information for the speed reason
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param user_inputs: contains all the user_inputs options
 * @param stats_info: a variable that contains all reads and base information
 */
void coverageBinningWrapper(Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Binned_Data_Wrapper *binned_data_wrapper, int32_t chrom_idx, int thread_id);

/*
 * This is used to generate average coverage information for a range of position based on different binning strategies
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_idx: the chromosome index 
 * @param user_inputs: contains all the user_inputs options
 * @param wgs_binned_coverage_fp: the file handle for depositing the binned data
 * */
void writeCoverageBins(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, int32_t chrom_idx, User_Input *user_inputs, Stats_Info *stats_info, FILE *fh_binned_coverage, Binned_Data_Wrapper *binned_data_wraper, int thread_id);

/*
 * This is the function used for insertion. But before doing any insertion, 
 * it will try to eliminate the edge effect by checking the binned data in front of it
 * if it is 5x apart, merge them
 * @param start: the start position of the bin
 * @param end: the end position of the bin
 * @param coverage: the total coverage of the bin to be inserted
 * @param binned_data_wraper: the viable that holding binned data sequentially
 */
void insertBinData(uint32_t start, uint32_t end, uint32_t length, double ave_coverage, Binned_Data_Wrapper *binned_data_wrapper);

/*
 * This function is used to process binned data and also insert data into bin
 * It will handle edge effects by combining neighboring bins
 */
void processBinnedData(uint32_t start, uint32_t end, uint32_t coverage, Binned_Data_Wrapper *binned_data_wraper, Chromosome_Tracking *chrom_tracking, FILE *fh_binned_coverage, int32_t chrom_idx, User_Input *user_inputs, uint32_t prev_negative_one);

/*
 * The following method is for debugging only
 */
void reportStatsForDebugging(Stats_Info *stats_info, User_Input *user_inputs);

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage);

void binnedDataWrapperInit(Binned_Data_Wrapper** binned_data_wrapper, Chromosome_Tracking *chrom_tracking);

void binnedDataWrapperDestroy(Binned_Data_Wrapper** binned_data_wrapper, Chromosome_Tracking *chrom_tracking);

void dynamicIncreaseBinSize(Binned_Data_Wrapper* binned_data_wrapper);

void exitWithFailure(void * data_point_in, char* message);


#endif //UTILS_H
