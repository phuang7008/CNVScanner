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

uint32_t processFile(char* chrom_id, char* file_name, khash_t(khIntStr) * starts, khash_t(khIntStr) * ends) {
    FILE *fp = fopen(file_name, "r");
    if (!fp) {perror("popen failed:"); exit(EXIT_FAILURE);}

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *tokPtr;
    uint32_t total_items=0;

    while ((read = getline(&line, &len, fp)) != -1) {
        // skip if it is a new line only or comment line
        //
        if (*line == '\n' || line[0] == '#') continue;
        total_items++;

        // remove the new line character
        //
        line[strcspn(line, "\n")] = 0;

        char *savePtr = line;
        //fprintf(stderr, "%s", line);
        char *dup_line = strdup(line);

        uint32_t i=0;
        while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
            if (i == 0) {
                i++;
                if (strcmp(chrom_id, tokPtr) != 0)
                    break;
            } else if (i == 1) {
                uint32_t start = (uint32_t) strtol(tokPtr, NULL, 10);

                // since strtok_r is disruptive, need to use *savePtr instead
                khashInsertion(starts, start, dup_line, file_name);
                i++;
            } else if (i == 2) {
                uint32_t end = (uint32_t) strtol(tokPtr, NULL, 10);
                khashInsertion(ends, end, dup_line, file_name);
                i++;
            }
        }
        free(dup_line);
    }

    if (fclose(fp) == -1)
        printf("feof=%d ferror=%d: %s\n", feof(fp), ferror(fp), strerror(errno));
    if (line) free(line);

    return total_items;
}

void khashInsertion(khash_t(khIntStr) *khash_in, uint32_t key, char* value, char *file_name) {
    int absent=0;
    khiter_t iter = kh_put(khIntStr, khash_in, key, &absent);
    if (absent) {
        kh_key(khash_in, iter) = key;
        kh_value(khash_in, iter) = calloc(strlen(value)+1, sizeof(char));
        strcpy(kh_value(khash_in, iter), value);
    } else {
        fprintf(stderr, "The file %s is not sorted or merged\n", file_name);
    }
}

void outputHashTable(khash_t(khIntStr) * khash_in) {
    FILE *out_file = fopen("mappability_output_test.txt", "w");
    if (out_file == NULL) fprintf(stderr, "Open mappability file writing failed\n");

    fprintf(stderr, "start writing\n");
    khiter_t iter;
    for (iter=kh_begin(khash_in); iter!=kh_end(khash_in); ++iter) {
        if (kh_exist(khash_in, iter)) {
            fprintf(out_file, "key: %"PRIu32"\tvalue: %s\n", kh_key(khash_in, iter), kh_value(khash_in, iter));
        }
    }
    fclose(out_file);
}

void forDebug() {
    fprintf(stderr, "debugging\n");
}
