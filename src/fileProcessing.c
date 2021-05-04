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

void processFile(char* chrom_id, char* file_name, Binned_Data_Wrapper* binned_data_wrapper) {
    char* cmd = calloc(strlen(file_name) * 2, sizeof(char));
    if (!cmd) { fprintf(stderr, "ERROR: char* memory allocation for cmd failed!"); exit(1); }

    // check file extension
    //
    const char* file_extension = getFileExtension(file_name);
    if (strcmp(file_extension, ".gz") == 0) {
        sprintf(cmd, "zgrep '^%s' %s", chrom_id, file_name);
    } else {
        sprintf(cmd, "grep '^%s' %s", chrom_id, file_name);
    }

    FILE *fp = popen(cmd, "r");
    if (!fp) {perror("popen failed:"); exit(1);}

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *tokPtr;
    uint32_t i;

    // need 2 loops, one for the binned_data_wrapper, while another one is for the file reading
    //
    for (i=0; i<binned_data_wrapper->size; i++) {

        while ((read = getline(&line, &len, fp)) != -1) {
            char *savePtr = line;
            //char *value = NULL;         // this value could be either mappability or gc content
            //char *chrom_id  = NULL;
            //int start=0, end=0;

            while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
                printf("data: '%s'\t", tokPtr);
           }
           break;
        }
        break;
    }

    printf("feof=%d ferror=%d: %s\n", feof(fp), ferror(fp), strerror(errno));
    pclose(fp);
    free(cmd);
}
