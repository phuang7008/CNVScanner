/*
 * =====================================================================================
 *
 *		Filename:		analysis.c
 *
 *		Description:	The implementation file for the analysis functions
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

#include "analysis.h"
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"
#include "fileProcessing.h"

void mappabilityGcNormalization(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, khash_t(khIntStr) *starts, khash_t(khIntStr) *ends, uint32_t total_lines, int type) {
    // create an all starts and ends array, the size will be dynamically increased later
    //
    AllStartsEndsArray *all_starts_ends_array = calloc(1, sizeof(AllStartsEndsArray));
    //all_starts_ends_array->capacity = binned_data_wraper->size * 4 + total_lines * 4;
    all_starts_ends_array->capacity = binned_data_wraper->size * 2 + total_lines * 2 + 10;
    all_starts_ends_array->array = calloc(all_starts_ends_array->capacity, sizeof(uint32_t));
    all_starts_ends_array->size = 0;

    // key: start (or end), value: "start end length ave_cov (type==1), first_normalized ave (type==2)"
    //
    khash_t(khIntStr) *binned_starts  = kh_init(khIntStr);      
    khash_t(khIntStr) *binned_ends    = kh_init(khIntStr);

    generateHashFromDynamicBins(binned_data_wraper, binned_starts, binned_ends, all_starts_ends_array, type);
    combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, starts);
    combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, ends);

    // before sorting
    //
    if (user_inputs->debug_ON)
        //outputAllPositionArray(all_starts_ends_array, 1);
        fprintf(stderr, "Before sorting, all_starts_ends_array size is %"PRIu32"\n", all_starts_ends_array->size);

    // sort the all_starts_ends_array
    //
    qsort(all_starts_ends_array->array, all_starts_ends_array->size, sizeof(uint32_t), compare);

    // after sorting
    //
    if (user_inputs->debug_ON)
        //outputAllPositionArray(all_starts_ends_array, 2);
        fprintf(stderr, "After sorting, all_starts_ends_array size is %"PRIu32"\n", all_starts_ends_array->size);

    // do intersect
    //
    uint32_t i=0;
    int16_t counter=0;      // the counter needs to be defined as signed int, otherwise, it will never be negative
    uint32_t prev_start0=0, prev_start1=0, map_gc_position=0;

    for (i=0; i<all_starts_ends_array->size; i++) {
        //if (user_inputs->debug_ON)
        //    fprintf(stderr, "%"PRIu32"\t%"PRIu32"\n", i, all_starts_ends_array->array[i]);

        if (checkKhashKey(ends, all_starts_ends_array->array[i]) || 
                checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {

            // always decrease counter if it is the end position
            //
            counter--;

            if (counter < 0) {
                fprintf(stderr, "Error: At type %d, the counter %"PRId16" should NOT be negative", type, counter);
                exit(EXIT_FAILURE);
            }

            if (counter == 1) {
                // the first part of the annotation should come from binned_starts or binned_ends
                // the second part of the annotation should come from starts or ends
                //
                char *binned_string=NULL, *map_gc_string=NULL;

                if (checkKhashKey(binned_starts, prev_start0))
                    binned_string = getKhashValue(binned_starts, prev_start0);

                if (checkKhashKey(starts, map_gc_position)) {
                    map_gc_string = getKhashValue(starts, map_gc_position);
                } else if (checkKhashKey(ends, map_gc_position)) {
                    map_gc_string = getKhashValue(ends, map_gc_position);
                }

                if (binned_string && map_gc_string) {
                    performNormalizationForCurrentBin(binned_data_wraper, user_inputs,
                        binned_string, map_gc_string, all_starts_ends_array->array[i], prev_start1, type);
                } else {
                    fprintf(stderr, "Something is wrong: Binned_string=%s; mapping_gc_string=%s with index %"PRIu32"\n",
                            binned_string, map_gc_string, i);
                    continue;
                }

                if (binned_string) { free(binned_string); binned_string=NULL; }
                if (map_gc_string) { free(map_gc_string); map_gc_string=NULL; }
            }

            // because bed file starts with 0 (or 0-indxed), so one value will appear in both starts and ends hash-tables
            // we have to remove those appeas in end position, so the next round, it will be a new start position only
            //
            khiter_t iter;
            if (checkKhashKey(ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(ends)) {
                    if (kh_value(ends, iter))
                        free(kh_value(ends, iter));
                    kh_del(khIntStr, ends, iter);
                }
            } else if (checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, binned_ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(binned_ends)) {
                    if (kh_value(binned_ends, iter))
                        free(kh_value(binned_ends, iter));      // this deletes the value first
                    kh_del(khIntStr, binned_ends, iter);        // this deletes the key second
                }
            }
        } else {
            // it should be in the 'start' position, always increment counter
            //
            counter++;

            // need to record the start position of current intersect
            // 
            // prev_start0-------------------------------------end0     => this is from the binned_array
            //              prev_start1-----------end1                  => this is from the map_gc array
            //
            // we see that we need to record both prev_start0 and prev_start1
            // it is because the intersect is from the map or gc array
            // But we need annotations from both binned array and map/gc array
            // I will use prev_start0 for the binned array annotation
            // I will use the next block for the map/gc annotation
            //
            if (checkKhashKey(binned_starts, all_starts_ends_array->array[i])) {
                prev_start0 = all_starts_ends_array->array[i];  // this is needed for the dynamic bin annotation
            } else if (!checkKhashKey(starts, all_starts_ends_array->array[i])) {
                fprintf(stderr, "Warning: index %"PRIu32" value %"PRIu32" is not in the starts hash, for type %d\n", 
                        i, all_starts_ends_array->array[i], type);
            }
            prev_start1 = all_starts_ends_array->array[i];      // this could be either from the dynamic binned or from map
        }

        // need to remember the map/gc% position of current intersect for the map/gc annotation
        //
        if (checkKhashKey(ends, all_starts_ends_array->array[i]) || 
                checkKhashKey(starts, all_starts_ends_array->array[i])) {
            map_gc_position = all_starts_ends_array->array[i];
        }
    }

    // clean-up
    //
    cleanAllStartsEndsArray(all_starts_ends_array);
    cleanKhashIntStr(binned_starts);
    cleanKhashIntStr(binned_ends);
}

void performNormalizationForCurrentBin(Binned_Data_Wrapper *binned_data_wraper, User_Input *user_inputs, char *bin_string, char* map_gc_string, uint32_t current_position, uint32_t prev_start, int type) {
    StringArray *binned_array = calloc(1, sizeof(StringArray));
    stringArrayInit(binned_array, 10);
    splitStringToArray(bin_string, binned_array);       // bin_string has 5 entries: index,chr,start,end,ave_cov

    StringArray *map_gc_array = calloc(1, sizeof(StringArray));
    stringArrayInit(map_gc_array, 10);
    splitStringToArray(map_gc_string, map_gc_array);    // map_gc_string has 4 entries: chr,start,end,mappability(or gc% scale)

    // Note: the binned_array->theArray[0] is the index to the binned_array_wrapper->data
    //
    uint32_t length = current_position - prev_start;        // it is the current window length after intersect
    uint32_t orig_len = strtoul(binned_array->theArray[3], NULL, 10) - strtoul(binned_array->theArray[2], NULL, 10); // the orig raw binned window
    double ave = strtod(binned_array->theArray[4], NULL);
    double scale_ratio = strtod(map_gc_array->theArray[3], NULL);

    // output for debugging
    //
    FILE *fp=NULL;
    if (user_inputs->debug_ON) {
        fp = fopen(user_inputs->map_gc_details_file, "a");
        fileOpenError(fp, user_inputs->map_gc_details_file);
    }

    if (type == 1) {
        binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_gc_scale += ((double) length / (double)orig_len) * scale_ratio;
        binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_gc_normalized = 
            ave * binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_gc_scale;

        // output for debugging
        //
        if (user_inputs->debug_ON)
            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%s\t%s\t%.2f\n", binned_array->theArray[1], prev_start, current_position, length, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_mappability, bin_string, map_gc_string, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_gc_normalized);

    } else {
        binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_mappability += ((double) length / (double)orig_len) * scale_ratio;
        binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_map_gc_normalized =
            ave / binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_mappability;

        // output for debugging
        //
        if (user_inputs->debug_ON)
            fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%s\t%s\t%.2f\n", binned_array->theArray[1], prev_start, current_position, length, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_gc_scale, bin_string, map_gc_string, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].ave_cov_map_gc_normalized);
    }

    if (fp != NULL) fclose(fp);

    // clean up
    //
    stringArrayDestroy(binned_array);
    stringArrayDestroy(map_gc_array);

    if (binned_array)
        free(binned_array);

    if (map_gc_array)
        free(map_gc_array);
}

// type 1 is for gc content normalization, while type 2 is for mappability normalization
// type 3 is for binned array used for equal window bins and type 4 is for equal window bins
// and store dynamically merged bin's starts and ends to all_starts_ends_array
//
void generateHashFromDynamicBins(Binned_Data_Wrapper *binned_data_wrapper, khash_t(khIntStr) *binned_starts, khash_t(khIntStr) *binned_ends, AllStartsEndsArray *all_starts_ends_array, int type) {
    // create a string pointer to store values
    //
    char * insert_value = calloc(250, sizeof(char));

    uint32_t i=0;
    for (i = 0; i < binned_data_wrapper->size; i++) {
        // Note here we need to store index i to the insert_value for later usage
        //
        if (type == 1) {
            sprintf(insert_value, "%"PRIu32"\t%s\t%"PRIu32"\t%"PRIu32"\t%.2f", i, binned_data_wrapper->chromosome_id, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i], binned_data_wrapper->data[i].ave_coverage);
        } else if (type == 2) {
            sprintf(insert_value, "%"PRIu32"\t%s\t%"PRIu32"\t%"PRIu32"\t%.2f", i, binned_data_wrapper->chromosome_id, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i], binned_data_wrapper->data[i].ave_cov_gc_normalized);
        } else if (type == 3) {
            sprintf(insert_value, "%"PRIu32"\t%s\t%"PRIu32"\t%"PRIu32"\t%.2f\t%.2f", i, binned_data_wrapper->chromosome_id, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i], binned_data_wrapper->data[i].weighted_mappability, binned_data_wrapper->data[i].ave_cov_map_gc_normalized);
        } else {
            sprintf(insert_value, "%"PRIu32"\t%s\t%"PRIu32"\t%"PRIu32, i, binned_data_wrapper->chromosome_id, binned_data_wrapper->starts[i], binned_data_wrapper->ends[i]);
        }

        khashInsertion(binned_starts, binned_data_wrapper->starts[i], insert_value);
        khashInsertion(binned_ends, binned_data_wrapper->ends[i], insert_value);

        all_starts_ends_array->array[all_starts_ends_array->size] = binned_data_wrapper->starts[i];
        all_starts_ends_array->size++;
        all_starts_ends_array->array[all_starts_ends_array->size] = binned_data_wrapper->ends[i];
        all_starts_ends_array->size++;
    }

    if (insert_value) free(insert_value);
}

void generateEqualSizedBins(User_Input *user_inputs, Binned_Data_Wrapper *binned_data_wrapper, Binned_Data_Wrapper *equal_size_window_wrapper, uint32_t total_lines) {
    // create an all starts and ends array, the size will be dynamically increased later
    //
    AllStartsEndsArray *all_starts_ends_array = calloc(1, sizeof(AllStartsEndsArray));
    all_starts_ends_array->capacity = binned_data_wrapper->size * 2 + total_lines * 2 + 10;
    all_starts_ends_array->array = calloc(all_starts_ends_array->capacity, sizeof(uint32_t));
    all_starts_ends_array->size = 0;

    khash_t(khIntStr) *binned_starts  = kh_init(khIntStr);      // key: start, value: "index chr start end mapp ave_cov_norm"
    khash_t(khIntStr) *binned_ends    = kh_init(khIntStr);      // key: end,   value: "index chr start end mapp ave_cov_norm" 

    generateHashFromDynamicBins(binned_data_wrapper, binned_starts, binned_ends, all_starts_ends_array, 3);

    khash_t(khIntStr) *window_starts  = kh_init(khIntStr);      // key: start, value: "index chr start end"
    khash_t(khIntStr) *window_ends    = kh_init(khIntStr);      // key: end,   value: "index chr start end" 

    generateHashFromDynamicBins(equal_size_window_wrapper, window_starts, window_ends, all_starts_ends_array, 4);
    
    // for debugging, output the equal_size_window_wrapper
    //
    if (user_inputs->debug_ON) {
        fprintf(stderr, "Before output raw window bins\n");
        outputBinnedData(equal_size_window_wrapper, user_inputs, 2);
        fprintf(stderr, "After output raw window bins\n");
    }

    //combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, window_starts);
    //combineAllStartsAndEndsFromOtherSource(all_starts_ends_array, window_ends);

    // sort the all_starts_ends_array
    //
    qsort(all_starts_ends_array->array, all_starts_ends_array->size, sizeof(uint32_t), compare);
    fprintf(stderr, "the sorted all starts and ends array has the size as %"PRIu32"\n", all_starts_ends_array->size);

    // do intersect
    //
    uint32_t i=0;
    int16_t counter=0;      // the counter needs to be defined as signed int, otherwise, it will never be negative
    uint32_t prev_start0=0, prev_start1=0, interval_position=0;

    for (i=0; i<all_starts_ends_array->size; i++) {
        //if (user_inputs->debug_ON)
        //    fprintf(stderr, "%"PRIu32"\t%"PRIu32"\n", i, all_starts_ends_array->array[i]);
        //
        if (checkKhashKey(window_ends, all_starts_ends_array->array[i]) ||
                checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {

            // always decrease counter if it is the end position
            //
            counter--;

            if (counter < 0) {
                fprintf(stderr, "Error: the counter %"PRId16" should NOT be negative", counter);
                exit(EXIT_FAILURE);
            }

            if (counter == 1) {
                // the first part of the annotation should come from binned_starts or binned_ends
                // the second part of the annotation should come from starts or ends
                //
                char *binned_string=NULL, *interval_string=NULL;

                if (checkKhashKey(binned_starts, prev_start0))
                    binned_string = getKhashValue(binned_starts, prev_start0);

                if (checkKhashKey(window_starts, interval_position)) {
                    interval_string = getKhashValue(window_starts, interval_position);
                } else if (checkKhashKey(window_ends, interval_position)) {
                    interval_string = getKhashValue(window_ends, interval_position);
                }

                if (binned_string && interval_string) {
                    storeWindowResults(binned_data_wrapper, equal_size_window_wrapper, user_inputs,
                            binned_string, interval_string, all_starts_ends_array->array[i], prev_start1);

                } else {
                    fprintf(stderr, "Something is wrong: Binned_string=%s; interval_string=%s with index %"PRIu32"\n",
                                binned_string, interval_string, i);
                    continue;
                }

                if (binned_string) { free(binned_string); binned_string=NULL; }
                if (interval_string) { free(interval_string); interval_string=NULL; }
            }

            // because bed file starts with 0 (or 0-indxed), so one value will appear in both starts and ends hash-tables
            // we have to remove those appeas in end position, so the next round, it will be a new start position only
            //
            // In addition:
            // NOTE: here we need to use if .. else if .. clause it is because of the following cases:
            // (base) [phuang@sug-login1 equal_interval_combining]$ grep 375000 1000_chr19_*
            //  1000_chr19_500bp_windows:19 374500  375000
            //  1000_chr19_500bp_windows:19 375000  375500
            //
            //  1000_chr19_binned_data:19   374643  375000  357 1.00    1.00    39.87   39.87   40.04   40
            //  1000_chr19_binned_data:19   375000  375275  275 1.00    0.95    31.00   31.00   29.35   29
            //
            //  both binned data and window interval contain the same coordinates
            //  Need to handle them one at a time, NOT BOTH
            //                                                                                                 # 
            khiter_t iter;
            if (checkKhashKey(binned_ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, binned_ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(binned_ends)) {
                    if (kh_value(binned_ends, iter))
                        free(kh_value(binned_ends, iter));
                    kh_del(khIntStr, binned_ends, iter);
                }
            } else if (checkKhashKey(window_ends, all_starts_ends_array->array[i])) {
                iter = kh_get(khIntStr, window_ends, all_starts_ends_array->array[i]);
                if (iter != kh_end(window_ends)) {
                    if (kh_value(window_ends, iter))
                        free(kh_value(window_ends, iter));      // this deletes the value first
                    kh_del(khIntStr, window_ends, iter);        // this deletes the key second
                }
            }
        } else {
            // it should be in the 'start' position, always increment counter
            //
            counter++;

            if (checkKhashKey(binned_starts, all_starts_ends_array->array[i])) {
                prev_start0 = all_starts_ends_array->array[i];  // this is needed for the dynamic bin annotation
            } else if (!checkKhashKey(window_starts, all_starts_ends_array->array[i])) {
                fprintf(stderr, "Warning: index %"PRIu32" value %"PRIu32" is not in the starts hash, in generateEqualBinSize()\n",
                            i, all_starts_ends_array->array[i]);
            }
                
            prev_start1 = all_starts_ends_array->array[i];  // this could be either from the dynamic bins or from window_interval
        }

        // need to record the map/gc% position of current intersect for the map/gc annotation
        //
        if (checkKhashKey(window_ends, all_starts_ends_array->array[i]) ||
                checkKhashKey(window_starts, all_starts_ends_array->array[i])) {
            interval_position = all_starts_ends_array->array[i];
        }
    }

    // clean-up
    //
    cleanAllStartsEndsArray(all_starts_ends_array);
    cleanKhashIntStr(binned_starts);
    cleanKhashIntStr(binned_ends);
    cleanKhashIntStr(window_starts);
    cleanKhashIntStr(window_ends);

    // now we need to divide the average coverage by total length to get the real average coverage
    //
    for (i=0; i<equal_size_window_wrapper->size; i++) {
        equal_size_window_wrapper->data[i].weighted_mappability /= equal_size_window_wrapper->data[i].length;
        equal_size_window_wrapper->data[i].ave_coverage /= equal_size_window_wrapper->data[i].length;
    }
}

void storeWindowResults(Binned_Data_Wrapper *binned_data_wraper, Binned_Data_Wrapper *equal_size_window_wrapper, User_Input *user_inputs, char *binned_string, char *interval_string, uint32_t current_position, uint32_t prev_start) {
    StringArray *binned_array = calloc(1, sizeof(StringArray));
    stringArrayInit(binned_array, 10);
    splitStringToArray(binned_string, binned_array);    // bin_string has 6 entries: index,chr,start,end,mapp,ave_cov_normalized

    StringArray *window_bin_array = calloc(1, sizeof(StringArray));
    stringArrayInit(window_bin_array, 10);
    splitStringToArray(interval_string, window_bin_array);  // interval_string has 3 entries: chr,start,end

    // Note: the binned_array->theArray[0] is the index to the binned_array_wrapper->data
    //
    uint32_t length = current_position - prev_start;
    //uint32_t orig_len = strtoul(window_bin_array->theArray[3], NULL, 10) - strtoul(window_bin_array->theArray[2], NULL, 10);
    double ave = strtod(binned_array->theArray[5], NULL);
    double mappability = strtod(binned_array->theArray[4], NULL);

    // output for debugging
    //
    FILE *fp=NULL;
    if (user_inputs->debug_ON) {
        fp = fopen(user_inputs->window_details_file, "a");
        fileOpenError(fp, user_inputs->window_details_file);
    }

    equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].length += length;
    //equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].ave_coverage += length * ave;
    equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].ave_coverage += length * ave;
    equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].weighted_mappability += length * mappability;

    if (user_inputs->debug_ON) {
        fprintf(fp, "%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%.2f\t%s\t%s\t%"PRIu32"\t%.2f\t%.2f\n", binned_array->theArray[1], prev_start, current_position, length, binned_data_wraper->data[strtoul(binned_array->theArray[0], NULL, 10)].weighted_mappability, binned_string, interval_string, equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].length, equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].ave_coverage, equal_size_window_wrapper->data[strtoul(window_bin_array->theArray[0], NULL, 10)].weighted_mappability);
    }

    if (fp) fclose(fp);

    //clean-up
    //
    stringArrayDestroy(binned_array);
    stringArrayDestroy(window_bin_array);

    if (binned_array)
        free(binned_array);

    if (window_bin_array)
        free(window_bin_array);
}

void combineAllStartsAndEndsFromOtherSource(AllStartsEndsArray *all_starts_ends_array, khash_t(khIntStr) *hash_in) {
    khiter_t iter;
    for (iter=kh_begin(hash_in); iter!=kh_end(hash_in); ++iter) {
        if (kh_exist(hash_in, iter)) {
            all_starts_ends_array->array[all_starts_ends_array->size] = kh_key(hash_in, iter);
            all_starts_ends_array->size++;
        }
    }
}
