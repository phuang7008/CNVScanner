/*
 * =====================================================================================
 *
 *		Filename:		utils.c
 *
 *		Description:	The implementation file for the utility functions
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang, phuang@bcm.edu
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access() and getopt()
#include <dirent.h>		// for checking output directory
#include <libgen.h>		// for function basename()
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"
#include "utility.h"


void cleanKhashStr(khash_t(str) *hash_to_clean, uint8_t type) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean key if the key exist
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));

			if (type == 1) {
				// clean Temp_Coverage_Array
				free(kh_value(hash_to_clean, k)->cov_array);
				free(kh_value(hash_to_clean, k));
				//kh_del(str, hash_to_clean, k);
			}
		}
	}
	//printf("before clean hash string\n");

	if (hash_to_clean) kh_destroy(str, hash_to_clean);
	//printf("after clean hash string\n");
}

void cleanKhashStrStrArray(khash_t(khStrStrArray) * hash_to_clean) {
	khint_t k;
	int i;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean the value first
			//
			for (i=0; i<kh_value(hash_to_clean, k)->size; i++) {
				if (kh_value(hash_to_clean, k)->theArray[i] != NULL)
					free(kh_value(hash_to_clean, k)->theArray[i]);
			}

			// clean theArray pointer
			//
			if (kh_value(hash_to_clean, k)->theArray)
				free(kh_value(hash_to_clean, k)->theArray);

			// clean key if the key exist
			//
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));

			// clean value if it exists
			//
			if (kh_value(hash_to_clean, k)) free (kh_value(hash_to_clean, k));
		}
	}

	if (hash_to_clean) kh_destroy(khStrStrArray, hash_to_clean);
}

void cleanKhashStrStr(khash_t(khStrStr) * hash_to_clean) {
    khint_t k;
    for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
        if (kh_exist(hash_to_clean, k)) {
            // clean key if the key exist
            //
            if (kh_key(hash_to_clean, k)) 
                free((char *) kh_key(hash_to_clean, k));

            // clean value if it exists
            //
            if (kh_value(hash_to_clean, k)) 
                free ((char *) kh_value(hash_to_clean, k));
        }
    }

    if (hash_to_clean) kh_destroy(khStrStr, hash_to_clean);
}

void cleanKhashIntPrArray(khash_t(khIntPrArray) *hash_to_clean) {
    khint_t k;
    for (k=kh_begin(hash_to_clean); k!=kh_end(hash_to_clean); ++k) {
        if (kh_exist(hash_to_clean, k)) { 
            uint32_t i;
            for (i=0; i<kh_value(hash_to_clean, k)->size; i++) {
                if (kh_value(hash_to_clean, k)->pread_x_a_bpt[i].read_name != NULL) {
                    free(kh_value(hash_to_clean, k)->pread_x_a_bpt[i].read_name);
                    kh_value(hash_to_clean, k)->pread_x_a_bpt[i].read_name = NULL;
                }
            }

            if (kh_value(hash_to_clean, k)->pread_x_a_bpt != NULL) {
                free(kh_value(hash_to_clean, k)->pread_x_a_bpt);
                kh_value(hash_to_clean, k)->pread_x_a_bpt = NULL;
            }

            if (kh_value(hash_to_clean, k) != NULL) {
                free(kh_value(hash_to_clean, k));
                kh_value(hash_to_clean, k) = NULL;
            }
        }
    }
    if (hash_to_clean)
        free(hash_to_clean);
}

bool checkKhashKey(khash_t(khIntStr) *hash_in, uint32_t key) {
    khiter_t iter = kh_get(khIntStr, hash_in, key);

    if (iter == kh_end(hash_in))
        return false;

    return true;
}

bool checkm32KhashKey(khash_t(m32) *hash_in, uint32_t key) {
    khiter_t iter = kh_get(m32, hash_in, key);

    if (iter == kh_end(hash_in))
        return false;

    return true;
}

char* getKhashValue(khash_t(khIntStr) *hash_in, uint32_t key) {
    khiter_t iter = kh_get(khIntStr, hash_in, key);

    if (iter == kh_end(hash_in))
        return NULL;

    //return strdup(kh_value(hash_in, key));      // need to remember to free the char* in the caller's method
    char *tmp_str = calloc(strlen(kh_value(hash_in, iter))+1, sizeof(char));
    strcpy(tmp_str, kh_value(hash_in, iter));
    return tmp_str;
}

void statsInfoInit(Stats_Info *stats_info) {
	if (!stats_info) {
		fprintf(stderr, "ERROR: Memory allocation failed for Stats_Info\n");
		exit(EXIT_FAILURE);
	}

    stats_info->read_cov_stats = calloc(1, sizeof(Read_Coverage_Stats));
    stats_info->wgs_cov_stats  = calloc(1, sizeof(WGS_Coverage_Stats));

    int i;
    for (i=0; i<=1000; i++) {
        stats_info->wgs_cov_stats->genome_cov_histogram[i]=0;
    }

    stats_info->wgs_cov_stats->genome_base_with_N_coverage = kh_init(m32);
    stats_info->wgs_cov_stats->genome_coverage_for_median  = kh_init(m32);

    readCoverageStatsInit(stats_info->read_cov_stats);
    WGSCoverageStatsInit(stats_info->wgs_cov_stats);   
}

void readCoverageStatsInit(Read_Coverage_Stats * read_cov_stats) {
    if (!read_cov_stats) {
        fprintf(stderr, "ERROR: Memory allocation failed for Read_Coverage_Stats\n");
        exit(EXIT_FAILURE);
    }

    read_cov_stats->total_reads_produced = 0;
    read_cov_stats->total_reads_aligned = 0;
    read_cov_stats->total_reads_paired = 0;
    read_cov_stats->total_reads_proper_paired = 0;
    read_cov_stats->total_duplicate_reads = 0;
    read_cov_stats->total_chimeric_reads = 0;
    read_cov_stats->total_supplementary_reads = 0;
    read_cov_stats->total_paired_reads_with_mapped_mates = 0;
    read_cov_stats->read_length = 0;
}

void WGSCoverageStatsInit(WGS_Coverage_Stats * wgs_cov_stats) {
    if (!wgs_cov_stats) {
        fprintf(stderr, "ERROR: Memory allocation failed for WGS Coverage Stats\n");
        exit(EXIT_FAILURE);
    }

    wgs_cov_stats->total_genome_bases = 0;
    wgs_cov_stats->total_excluded_bases = 0;
    wgs_cov_stats->total_excluded_bases_on_chrX = 0;
    wgs_cov_stats->total_excluded_bases_on_chrY = 0;
    wgs_cov_stats->total_mapped_bases = 0;
    wgs_cov_stats->total_uniquely_aligned_bases = 0;
    wgs_cov_stats->total_genome_coverage = 0;
    wgs_cov_stats->base_quality_20 = 0;
    wgs_cov_stats->base_quality_30 = 0;
    wgs_cov_stats->total_overlapped_bases = 0;

    wgs_cov_stats->wgs_max_coverage = 0;
    wgs_cov_stats->base_with_wgs_max_coverage = 0;
    wgs_cov_stats->median_genome_coverage = 0;

}

void statsInfoDestroy(Stats_Info *stats_info) {
    if (stats_info->read_cov_stats) {
        free(stats_info->read_cov_stats);
        stats_info->read_cov_stats = NULL;
    }

    if (stats_info->wgs_cov_stats) {
        //kh_destroy(m32, stats_info->wgs_cov_stats->genome_cov_histogram);
        kh_destroy(m32, stats_info->wgs_cov_stats->genome_base_with_N_coverage);
        kh_destroy(m32, stats_info->wgs_cov_stats->genome_coverage_for_median);
        free(stats_info->wgs_cov_stats);
        stats_info->wgs_cov_stats = NULL;
    }

    if (stats_info != NULL) free(stats_info);

}

void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val) {
    int ret;
    khiter_t k_iter = kh_put(m16, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(EXIT_FAILURE);
    }

    kh_value(hash_in, k_iter) += val;

    return;
}

void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint32_t val) {
    int ret;
    khiter_t k_iter = kh_put(m32, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
		//printf("add value is 0 %d ret, with key %"PRIu32"\n", ret, pos_key);
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(EXIT_FAILURE);
    }

    kh_value(hash_in, k_iter) += val;
	//if (pos_key == 1 && kh_value(hash_in, k_iter)%1000000 == 0) printf("key 1 value is %"PRIu32"\n", kh_value(hash_in, k_iter));
	if (kh_value(hash_in, k_iter) > 4294967290) printf("larger value %"PRIu32"\n", kh_value(hash_in, k_iter));

    return;
}

uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key) {
    khiter_t k_iter;
	if (hash32 != NULL) {
		k_iter = kh_get(m32, hash32, pos_key);

		if (k_iter == kh_end(hash32))
			// this is needed as if the bucket is never used, the value will be whatever left there, and the value is undefined!
			//
			return 0;

        return kh_value(hash32, k_iter);
    }

	return 0;
}

void addValueToKhashBucketStrStr(khash_t(khStrStr) *hash_in, char *key, char * val) {
    int ret;
    khiter_t k_iter = kh_put(khStrStr, hash_in, key, &ret);
    if (ret == 1) {
        kh_key(hash_in, k_iter) = strdup(key);
        kh_value(hash_in, k_iter) = strdup(val);
        //strcpy(kh_value(hash_in, k_iter), val);
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key  %s\n", key);
        exit(EXIT_FAILURE);
    }
}

char * getValueFromKhashStrStr(khash_t(khStrStr) *hash_in, char* key) {
    khiter_t k_iter;
    if (hash_in != NULL) {
        k_iter = kh_get(khStrStr, hash_in, key);
        
        if (k_iter == kh_end(hash_in))
            // this is needed as if the bucket is never used, the value will be whatever left there, and the value is undefined!
            //
            return NULL;
        
        return kh_value(hash_in, k_iter);
    }
    
    return NULL;
}

uint16_t getValueFromKhash(khash_t(m16) *hash16, uint32_t pos_key) {
    khiter_t k_iter;
    if (hash16 != NULL) {
		k_iter = kh_get(m16, hash16, pos_key);

		if (k_iter == kh_end(hash16))
			return 0;

        return kh_value(hash16, k_iter);
    }

    return 0;
}

float calculatePercentage32(uint32_t num, uint32_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

float calculatePercentage32_64(uint32_t num, uint64_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

float calculatePercentage64(uint64_t num, uint64_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

void combineCoverageStats(Stats_Info *stats_info, Stats_Info *tmp_stats_info) {
	// For read length
    //
    if (stats_info->read_cov_stats->read_length == 0)
        stats_info->read_cov_stats->read_length = tmp_stats_info->read_cov_stats->read_length;

    copyReadCoverageStats(stats_info->read_cov_stats, tmp_stats_info->read_cov_stats);
    copyWGSCoverageStats(stats_info->wgs_cov_stats, tmp_stats_info->wgs_cov_stats);
}

void copyReadCoverageStats(Read_Coverage_Stats *read_cov_stats, Read_Coverage_Stats *tmp_read_cov_stats) {
	read_cov_stats->total_reads_produced  += tmp_read_cov_stats->total_reads_produced;
	read_cov_stats->total_reads_aligned   += tmp_read_cov_stats->total_reads_aligned;
	read_cov_stats->total_chimeric_reads  += tmp_read_cov_stats->total_chimeric_reads;
	read_cov_stats->total_duplicate_reads += tmp_read_cov_stats->total_duplicate_reads;
	read_cov_stats->total_reads_paired    += tmp_read_cov_stats->total_reads_paired;
	read_cov_stats->total_reads_proper_paired += tmp_read_cov_stats->total_reads_proper_paired;
	read_cov_stats->total_supplementary_reads += tmp_read_cov_stats->total_supplementary_reads;
	read_cov_stats->total_paired_reads_with_mapped_mates += tmp_read_cov_stats->total_paired_reads_with_mapped_mates;
}

void copyWGSCoverageStats(WGS_Coverage_Stats *wgs_cov_stats, WGS_Coverage_Stats *tmp_wgs_cov_stats) {
	// base related stats
	//
	wgs_cov_stats->base_quality_20        += tmp_wgs_cov_stats->base_quality_20;
	wgs_cov_stats->base_quality_30        += tmp_wgs_cov_stats->base_quality_30;
	wgs_cov_stats->total_mapped_bases     += tmp_wgs_cov_stats->total_mapped_bases;
	wgs_cov_stats->total_overlapped_bases += tmp_wgs_cov_stats->total_overlapped_bases;
	wgs_cov_stats->total_uniquely_aligned_bases += tmp_wgs_cov_stats->total_uniquely_aligned_bases;
}

// get the key for hash table for everything 1000 bases
//
uint32_t getHashKey(uint32_t position_in) {
	uint32_t tmp_key = (uint32_t) position_in / 1000;
	return tmp_key * 1000;
}

void outputFreqDistribution(User_Input *user_inputs, khash_t(m32) *cov_freq_dist) {
	// open WGS coverage summary report file handle
	//
	FILE *out_fp = fopen(user_inputs->wgs_cov_report, "a");
    fileOpenError(out_fp, user_inputs->wgs_cov_report);

	fprintf(out_fp, "\n#Smoothed_Coverage_Frequency_Distribution_for_Whole_Genome\n");
	fprintf(out_fp, "==");
	khiter_t iter;
	for (iter=kh_begin(cov_freq_dist); iter!=kh_end(cov_freq_dist); iter++) {
		if (kh_exist(cov_freq_dist, iter))
			fprintf(out_fp, "%"PRIu32",", kh_value(cov_freq_dist, iter));
	}
	fprintf(out_fp, "\n");
	fclose(out_fp);
}

void splitStringToArray(char* string_to_split, StringArray *string_array) {
    // since the strtok_r is destructive, so I have to use a copy for this
    // the bin_string has 5 entries like the follow order
    //      index    chrom    start    end    value
    // the map_gc_string has 4 entries as the following:
    //      chrom    start    end    mappability/gc%-scale
    //
    char *tmp_string = calloc(strlen(string_to_split)+1, sizeof(char));
    strcpy(tmp_string, string_to_split);

    char *tokPtr;
    char *savePtr = tmp_string;

    while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
        string_array->theArray[string_array->size] = strdup(tokPtr);
        string_array->size++;
    }

    if (tmp_string!=NULL)
        free(tmp_string);
}

// To view/print the content of string array before OR after sorting
// Note:  Any arrays will decay to a pointer to the first element when passing to a function.
// Therefore, we will have to pass the size info into the function to make it work!
//
void print_string_array(char** strings_in, size_t length_in) {
	size_t i;
	
	for (i=0; i<length_in; i++) {
		printf("%s\t", strings_in[i]);
	}

	printf("\n");
}

void checkMemoryAllocation(void* newly_created_object, char* message) {
    if (!newly_created_object) {
        fprintf(stderr, "ERROR: Memory allocation for %s failed\n", message);
        exit(EXIT_FAILURE);
    }
}

void cleanKhashIntStr(khash_t(khIntStr) * hash_to_clean) {
    khint_t k;

    for (k=kh_begin(hash_to_clean); k!=kh_end(hash_to_clean); ++k) {
        if (kh_exist(hash_to_clean, k)) {
            if (kh_value(hash_to_clean, k))
                free(kh_value(hash_to_clean,k));
        }
    }

    if (hash_to_clean) kh_destroy(khIntStr, hash_to_clean);
}

void cleanAllStartsEndsArray(AllStartsEndsArray *all_starts_ends_array) {
    if (all_starts_ends_array->array) free(all_starts_ends_array->array);
    if (all_starts_ends_array) free(all_starts_ends_array);
}

void outputFinalBinnedData(Binned_Data_Wrapper **binned_data_wrapper, User_Input *user_inputs, Chromosome_Tracking *chrom_tracking, int type) {
    FILE *fp=NULL;
    if (type == 1) {
        fp = fopen(user_inputs->normalized_result_file, "w");
    } else {
        fp = fopen("Equal_window_bins_norm_final.txt", "w");
    }
    fileOpenError(fp, "Equal_window_bins_norm_final.txt or user_inputs->normalized_result_file");

    uint32_t i=0, j=0;

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        for (j=0; j<binned_data_wrapper[i]->size; j++) {

            if (type == 2) {
                // Need to skip those windows with repeat maskers, mappbility <0.2, Ns regions etc
                //
                if (binned_data_wrapper[i]->data[j].length == 0)
                    continue;
            }

            // need to round the final normalized data for stats analysis
            //
            uint32_t round_normalized_value = (uint32_t) (binned_data_wrapper[i]->data[j].ave_cov_map_gc_normalized + 0.5);

            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%"PRIu32"\n", 
                binned_data_wrapper[i]->chromosome_id, binned_data_wrapper[i]->data[j].start, 
                binned_data_wrapper[i]->data[j].end, binned_data_wrapper[i]->data[j].length,
                binned_data_wrapper[i]->data[j].weighted_mappability,
                binned_data_wrapper[i]->data[j].weighted_gc_scale,
                binned_data_wrapper[i]->data[j].ave_coverage,
                binned_data_wrapper[i]->data[j].ave_cov_gc_normalized,
                binned_data_wrapper[i]->data[j].ave_cov_map_gc_normalized, round_normalized_value);
        }
    }
    if (fp) fclose(fp);
}

void outputBinnedData(Binned_Data_Wrapper *binned_data_wrapper, User_Input *user_inputs, int type) {
    FILE *binned_coverage_fp;
    if (type == 1) {
        binned_coverage_fp = fopen(user_inputs->merged_bin_file, "a");
    } else {
        binned_coverage_fp = fopen("Raw_equal_window_bins.txt", "a");
    }
    fileOpenError(binned_coverage_fp, "Raw_window_bins_w_index.txt or user_inputs->merged_bin_file");
                    
    uint32_t i;
    for (i=0; i<binned_data_wrapper->size; i++) {
        fprintf(binned_coverage_fp, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n",
                binned_data_wrapper->chromosome_id, binned_data_wrapper->data[i].start, binned_data_wrapper->data[i].end,
                binned_data_wrapper->data[i].length, binned_data_wrapper->data[i].ave_coverage);
    }

    if (binned_coverage_fp) fclose(binned_coverage_fp);
}

void removeDebugFiles(User_Input *user_inputs) {
    if (user_inputs->debug_ON) {
        if (user_inputs->merged_bin_file && checkFile(user_inputs->merged_bin_file))
            remove(user_inputs->merged_bin_file);
        
        if (user_inputs->map_gc_details_file && checkFile(user_inputs->map_gc_details_file))
            remove(user_inputs->map_gc_details_file);
                
        if (user_inputs->normalized_result_file && checkFile(user_inputs->normalized_result_file))
            remove(user_inputs->normalized_result_file);

        if (user_inputs->mappability_outfile && checkFile(user_inputs->mappability_outfile))
            remove(user_inputs->mappability_outfile);

        if (user_inputs->gc_content_outfile && checkFile(user_inputs->gc_content_outfile))
            remove(user_inputs->gc_content_outfile);
    }
}

// status: 1 is before sorting, while 2 is after sorting
void outputAllPositionArray(AllStartsEndsArray *all_starts_ends_array, int status) {
    FILE *fp = NULL;

    if (status == 1) {
        fp = fopen("array_data_before_sorting.txt", "w");
    } else {
        fp = fopen("array_data_after_sorting.txt", "w");
    }
    fileOpenError(fp, "Array data before/after sorting");

    uint32_t i;
    for (i=0; i<all_starts_ends_array->size; i++) {
        fprintf(fp, "%"PRIu32"\t%"PRIu32"\n", i, all_starts_ends_array->array[i]);
    }

    fclose(fp);
}

void fileOpenError(FILE *fp, char *message) {
    if (fp == NULL) {
        fprintf(stderr, "File %s open failed\n", message);
        exit(EXIT_FAILURE);
    }
}

void failureExit(void * data_point_in, char* message) {
    if (data_point_in == NULL) {
        fprintf(stderr, "ERROR: Dynamic Memory allocation for %s failed!\n", message);
        exit(EXIT_FAILURE);
    }
}

