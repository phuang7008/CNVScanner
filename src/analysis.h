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
#include "utils.h"

/*
 * This function performs mappability or GC% normalization
 * @param binned_data_wraper: a variable stores the binned results including the normalized data
 * @param user_inputs, a variable to store all user's inputs
 * @param starts, a hash table stores start as key, the mappability/GC% scale line as value
 * @param ends,   a hash table stores end as key, the mappability/GC% scale line as value
 * @param total_lines, the total number of entries in the mappability or GC% file for this chromosome
 * @param type, the normalization type: type 1 for gc while type 2 for mappability
 * @chrom_id, the chromosome id it is currently handling
 */ 
void mappabilityGcNormalization(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, khash_t(khIntStr) *starts, khash_t(khIntStr) *ends, uint32_t total_lines, Simple_Stats *wgs_simple_stats, int type, char *chrom_id);

/*
 * This function will generate all starts and ends array from dynamic bins.
 * It will also also create two hash tables, one for all starts and one for all ends
 * @param binned_data_wraper: a variable stores the binned results including the normalized data
 * @param user_inputs, a variable to store all user's inputs
 */
void performNormalizationForCurrentBin(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, char *bin_string, char* map_gc_string, uint32_t current_position, uint32_t prev_start, Simple_Stats *wgs_simple_stats, int type, FILE *map_gc_detail_fp);

/*
 * This function is used to obtain certain variable data from Dynamic Bins stored at 'binned_data_wrapper'
 * if type == 1, we will fetch the raw average coverage so that we can perform GC% normalization
 * if type == 2, we will fetch the GC% normalized average coverage so that we can perform mappability normalization
 * if type == 3, we will fetch the GC% & mappability normalized average coverage, so that we can calculate weighted equal window average coverage for the equal bins
 * if type == 4, we simply fetch final equal bin data for the CNV calling
 */
void generateHashFromDynamicBins(Binned_Data_Wrapper *binned_data_wrapper, khash_t(khIntStr) *binned_starts, khash_t(khIntStr) *binned_ends, AllStartsEndsArray *all_starts_ends_array, int type);

void generateEqualSizedBins(User_Input *user_inputs, Binned_Data_Wrapper *binned_data_wrapper, Binned_Data_Wrapper *equal_size_window_wrapper, uint32_t total_lines, char *chrom_id);

void storeWindowResults(Binned_Data_Wrapper *binned_data_wraper, Binned_Data_Wrapper *equal_size_window_wrappers, User_Input *user_inputs, char *binned_string, char *interval_string, uint32_t current_position, uint32_t prev_start, FILE *equal_window_fp);

void combineAllStartsAndEndsFromOtherSource(AllStartsEndsArray *all_starts_ends_array, khash_t(khIntStr) *hash_in);


#endif //UTILS_H
