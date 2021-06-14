/*
 * =====================================================================================
 *
 *      Filename:       utilss.h
 *
 *      Description:    For the general utility functionalities
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
#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>	// for bool return type
#include <ctype.h>		// for isdigit()
#include <errno.h>
#include "terms.h"

/*
 * it is used to calculation the size of low coverage regions from a StrInt Hash table
 * In addition, it will generate the low coverage regions in sorted order in string format
 * @param low_cov_regions: a StrInt hash table contain the start and end position of low coverage region
 * @param output: output string array
 * @return total low coverage region size
 */
uint32_t processLowCovRegionFromKhash(khash_t(khStrInt) *low_cov_regions, char *output);

/* this is used to process low coverage regions for output
 * @param low_cov_regions: a string array that contains all low coverage regions
 * @param output: combined all low coverage regions
 * @return number of low coverage regions
 */
uint32_t processLowCovRegionFromStrArray(StringArray *low_cov_regions, char *output);

/*
 * This function is used to clean the khash_t (uint32_t key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void cleanKhashInt(khash_t(m32) *hash_to_clean);

/* 
 * It is used to clean the kash_t (char* as key, but string array as value) hash table
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories including keys and char* array
 */
void cleanKhashStrStrArray(khash_t(khStrStrArray) * hash_to_clean);

/*
 * This function is used to clean the khash_t (string key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 * @param type: if it is set to 1, it will clean Temp_Coverage_Array struct variable as well
 */
void cleanKhashStr(khash_t(str) *hash_to_clean, uint8_t type);

/*
 * This is used to free the memory used by the hash_to_clean
 * @param hash_to_clean: a khash_t variable to be cleaned
 */
void cleanKhashStrStr(khash_t(khStrStr) * hash_to_clean);

/*
 * Initialize the member of the Stats_Info variable
 * @param stats_info: an instance of Stats_Info to store the whole genome coverage stat info
 */
void statsInfoInit(Stats_Info *stats_info);

/*
 * for XXX_Coverage_Stats variable initialization
 * @param xxx_cov_stats: an instance of XXX_Coverage_Stats to store the coverage statistics for current chrom
 */
void WGSCoverageStatsInit(WGS_Coverage_Stats * wgs_cov_stats);
void readCoverageStatsInit(Read_Coverage_Stats * read_cov_stats);

/*
 * to destroy everything allocated for stats_info
 * @param stats_info
 */
void statsInfoDestroy(Stats_Info *stats_info);

/*
 * To add value into a hash table by the key
 * @param hash_in: the hash table to be modified
 * @param pos_key: the key for a specific position
 * @param val: value to be added
 */
void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val);
void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint32_t val);
void addValueToKhashBucketStrStr(khash_t(khStrStr) *hash_in, char *key, char * val);

/*
 * Get value from the hash table by the key
 * @param hash_in
 * @param key
 * @return the value pointed by the key
 */
uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key);
uint16_t getValueFromKhash16(khash_t(m16) *hash16, uint32_t pos_key);
char * getValueFromKhashStrStr(khash_t(khStrStr) *hash_in, char* key);

/* converts fractions into percentages with 2 decimal positions
 * @param num: the numeric value 
 * @param dom: the donominator value
 * @return a float value with 2 decimal points
 */
float calculatePercentage64(uint64_t num, uint64_t dom);
float calculatePercentage32(uint32_t num, uint32_t dom);
float calculatePercentage32_64(uint32_t num, uint64_t dom);

/* 
 * combine all the coverage stats from each individual thread and store them in the stats_info
 * @param stats_info: a storage place for all summarized sequencing stats
 * @param cov_stats: detailed coverage stats
 */
void combineCoverageStats(Stats_Info *stats_info, Stats_Info *tmp_stats_info);

void copyReadCoverageStats(Read_Coverage_Stats *read_cov_stats, Read_Coverage_Stats *tmp_read_cov_stats);
void copyWGSCoverageStats(WGS_Coverage_Stats *wgs_cov_stats, WGS_Coverage_Stats *tmp_wgs_cov_stats);

/* obtain the hash key for the integer value passed in                                                                
 * the keys are defined every 1000 position for quick lookup
 * For example: 3241645 (pos) -> 3241000 (key)
 * @param position_in, the position for the key                                                               
 * @return hash key                                                                                           
 */
uint32_t getHashKey(uint32_t position_in);

void outputFreqDistribution(User_Input *user_inputs, khash_t(m32) *cov_freq_dist);

void splitStringToArray(char* string_to_split, StringArray *string_array);

/*
 * It is used to print a string array before (OR after sorting) for viewing and comparison.
 * @param strings_in: the string array to be printed!
 * @param length_in: specifies the size of input string array
 */
void print_string_array(char** strings_in, size_t length_in);

void checkMemoryAllocation(void* newly_created_object, char* message);

bool checkKhashKey(khash_t(khIntStr) *hash_in, uint32_t key);

char* getKhashValue(khash_t(khIntStr) *hash_in, uint32_t key);

/**
 * This is used to free the memory used by the hash_to_clean
 * @param hash_to_clean: a khash_t variable to be cleaned
 */
void cleanKhashIntStr(khash_t(khIntStr) * hash_to_clean);

void cleanAllStartsEndsArray(AllStartsEndsArray *all_starts_ends_array);

/*
 * This function will output final normalized data to the file
 * @param binned_data_wrapper: a variable that contains the final normalized binned results
 * @param user_inputs: a variable that contains the user provided options
 * @param chrom_tracking: a variable that contains the chromosome informaiton to be processed
 */
void outputFinalBinnedData(Binned_Data_Wrapper **binned_data_wrapper, User_Input *user_inputs, Chromosome_Tracking *chrom_tracking);

void outputBinnedData(Binned_Data_Wrapper *binned_data_wrapper, char* chrom_id, User_Input *user_inputs);

/*
 * this function clears the debugging output file produced in early runs
 * because most of using append mode
 */
void removeDebugFiles(User_Input *user_inputs);

void outputAllPositionArray(AllStartsEndsArray *all_starts_ends_array, int status);


#endif //UTILS_H
