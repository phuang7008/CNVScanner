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
 * @param starts, a hash table stores start as key, the mappability/GC% scale value as value
 * @param user_inputs, a variable to store all user's inputs
 * @param ends,   a hash table stores end as key and the mappability/GC% as value
 * @param totalline, the total number of entries in the mappability or GC% file for this chromosome
 */ 
void mappabilityGcNormalization(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, khash_t(khIntStr) *starts, khash_t(khIntStr) *ends, uint32_t total_lines, int type);

/*
 * This function will generate all starts and ends array from dynamic bins.
 * It will also also create two hash tables, one for all starts and one for all ends
 * @param binned_data_wraper: a variable stores the binned results including the normalized data
 * @param user_inputs, a variable to store all user's inputs
 */
void performNormalizationForCurrentBin(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, char *bin_string, char* map_gc_string, uint32_t current_position, uint32_t prev_start, int type);

void generateHashFromDynamicBins(Binned_Data_Wrapper *binned_data_wrapper, khash_t(khIntStr) *binned_starts, khash_t(khIntStr) *binned_ends, AllStartsEndsArray *all_starts_ends_array, int type);

void generateEqualSizedBins(User_Input *user_inputs, Binned_Data_Wrapper *binned_data_wrapper, Binned_Data_Wrapper *equal_size_window_wrapper, uint32_t total_lines);

void store_window_results(Binned_Data_Wrapper *binned_data_wraper, Binned_Data_Wrapper *equal_size_window_wrappers, User_Input *user_inputs, char *binned_string, char *interval_string, uint32_t current_position, uint32_t prev_start);

void combineAllStartsAndEndsFromOtherSource(AllStartsEndsArray *all_starts_ends_array, khash_t(khIntStr) *hash_in);


#endif //UTILS_H
