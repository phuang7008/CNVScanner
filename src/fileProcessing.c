/*
 * =====================================================================================
 *
 *       Filename:  fileProcessing.c
 *
 *    Description:  The detailed implementation of the header file
 *
 *        Version:  1.0
 *        Created:  04/21/2021 05:07:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston TX
 *
 * =====================================================================================
 */

#include <errno.h>

#include "fileProcessing.h"

uint32_t processFile(char* chrom_id, char* file_name, khash_t(khIntStr) * starts, khash_t(khIntStr) * ends, Binned_Data_Wrapper *binned_data_wrapper) {
    FILE *fp = fopen(file_name, "r");
    fileOpenError(fp, file_name);

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *tokPtr;
    uint32_t total_items=0;

    while ((read = getline(&line, &len, fp)) != -1) {
        // skip if it is a new line only or comment line
        //
        if (*line == '\n' || line[0] == '#') continue;

        // remove the new line character
        //
        line[strcspn(line, "\n")] = 0;

        // since strtok_r is disruptive, need to use *savePtr instead
        //
        char *savePtr = line;
        //fprintf(stderr, "%s", line);
        char *dup_line = strdup(line);

        uint32_t i=0, start=0;
        while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
            if (i == 0) {
                i++;
                if (strcmp(chrom_id, tokPtr) != 0)
                    break;

                total_items++;

            } else if (i == 1) {
                start = (uint32_t) strtol(tokPtr, NULL, 10);

                khashInsertion(starts, start, dup_line);
                i++;
            } else if (i == 2) {
                uint32_t end = (uint32_t) strtol(tokPtr, NULL, 10);
                khashInsertion(ends, end, dup_line);
                i++;

                if (binned_data_wrapper) {
                    // check the capacity ofthe binned_data_wrapper
                    //
                    if (binned_data_wrapper->size + 20 > binned_data_wrapper->capacity) {
                        fprintf(stderr, "before dynamic size increase %d in fileProcessing\n", binned_data_wrapper->capacity);
                        dynamicIncreaseBinSize(binned_data_wrapper);
                        fprintf(stderr, "after dynamic size increase %d in fileProcessing\n", binned_data_wrapper->capacity);
                    }

                    // need to set length to 0 here so it will be added later
                    //
                    insertBinData(start, end, 0, 0, binned_data_wrapper);
                }
            }
        }
        if (dup_line) free(dup_line);
    }

    if (fclose(fp) == -1)
        printf("feof=%d ferror=%d: %s\n", feof(fp), ferror(fp), strerror(errno));
    if (line) free(line);

    return total_items;
}

void khashInsertion(khash_t(khIntStr) *khash_in, uint32_t key, char* value) {
    int absent=0;
    khiter_t iter = kh_put(khIntStr, khash_in, key, &absent);
    if (absent) {
        kh_key(khash_in, iter) = key;
        kh_value(khash_in, iter) = strdup(value);
    } else {
        fprintf(stderr, "The key %"PRIu32" exist with value %s, maybe the input data is not sorted or merged?\n", key, value);
    }
}

void outputHashTable(khash_t(khIntStr) * khash_in, int type, User_Input *user_inputs, char *chrom_id) {
    char *filename=NULL;
    FILE *out_file;
    if (type == 1) {
        filename = calloc(strlen(user_inputs->mappability_outfile) + strlen(chrom_id) + 10, sizeof(char));
        sprintf(filename, "%s%s.txt", user_inputs->mappability_outfile, chrom_id);
        out_file = fopen(filename, "w");
    } else {
        filename = calloc(strlen(user_inputs->gc_content_outfile) + strlen(chrom_id) + 10, sizeof(char));
        sprintf(filename, "%s%s.txt", user_inputs->gc_content_outfile, chrom_id);
        out_file = fopen(filename, "w");
    }
    fileOpenError(out_file, filename);


    //fprintf(stderr, "start writing\n");
    khiter_t iter;
    for (iter=kh_begin(khash_in); iter!=kh_end(khash_in); ++iter) {
        if (kh_exist(khash_in, iter)) {
            fprintf(out_file, "key: %"PRIu32"\tvalue: %s\n", kh_key(khash_in, iter), kh_value(khash_in, iter));
        }
    }
    fclose(out_file);
    if (filename) free(filename);
}

void forDebug() {
    fprintf(stderr, "debugging\n");
}
