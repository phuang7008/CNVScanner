/*
 * =====================================================================================
 *
 *       Filename:  utility.c
 *
 *    Description:  The implementation file for the utility header file
 *
 *        Version:  1.0
 *        Created:  03/31/2021 11:46:44 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>        // for file access() and getopt()
#include <dirent.h>        // for checking output directory
#include <libgen.h>        // for function basename()

#include "utility.h"

// Before open any file for processing, it is always a good idea to make sure that file exists
//
bool checkFile(char * fName) {
    if(access(fName, F_OK|R_OK) == -1) {
        fprintf(stderr, "ERROR: No such file as \n%s\n  File not found.\n", fName);
        return false;
        //exit(EXIT_FAILURE);
    }
    return true;
}

uint64_t check_file_size(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0)
        return st.st_size;

    fprintf(stderr, "ERROR: Something is wrong when check the file size for file: \n%s\n", filename);
    exit(EXIT_FAILURE);
}

// This is used to check if a string (ie char *) is an int number
//
bool isNumber(const char * inStr) {
    if (strlen(inStr) == 0) return false;
    
    while( *inStr != '\0') {
        if (!isdigit(inStr[0])) {
            return false;
        }
        inStr++;
    }
    return true;
}

// This is used to check if a string (ie char *) is an int number
//
bool isFloat(const char *str, float *dest) {
    if (str == NULL) return false;

    char *endptr;
    *dest = (float) strtod(str, &endptr);
    if (str == endptr) return false;    // no conversion

    // look at trailing text
    //
    while (isspace((unsigned char) *endptr))
        endptr++;

    return *endptr == '\0';
}

// Here I need to pass in file_in name string as reference, otherwise, it will be by value and will get segmentation fault
//
void createFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append, char* version) {
    *file_in = calloc(strlen(output_dir)+strlen(base_name)+strlen(string_to_append)+2,  sizeof(char));
    strcpy(*file_in, output_dir);
    strcat(*file_in, "/");
    strcat(*file_in, base_name);
    strcat(*file_in, string_to_append);

    // need to write the version number to every output file
    //
    FILE *out_fp = fopen(*file_in, "w");
    fprintf(out_fp, "%s\n", version);
    fclose(out_fp);
}

void checkFileExtension(char* in_bam_file, samFile* sfd) {
    const char* file_extension = getFileExtension(in_bam_file);
    switch (sfd->format.format) {
        case bam:
            if (stristr(".bam", file_extension) == NULL) {
                fprintf(stderr, "\nERROR: The input file  \n%s\n is in bam format\n\n", in_bam_file);
                exit(EXIT_FAILURE);
            }
            break;
        case cram:
            if (stristr(".cram", file_extension) == NULL) {
                fprintf(stderr, "\nERROR: The input file \n%s\n is in cram format\n\n", in_bam_file);
                exit(EXIT_FAILURE);
            }
            break;
        case sam:
            if (stristr(".sam", file_extension) == NULL) {
                fprintf(stderr, "\nERROR: The input file \n%s\n is in sam format\n\n", in_bam_file);
                exit(EXIT_FAILURE);
            }
            break;
        default:
            fprintf(stderr, "\nERROR: Unknown input file format \n%s\n\n", in_bam_file);
            exit(EXIT_FAILURE);
    }

    if (file_extension != NULL) free((char*)file_extension);
}

const char* getFileExtension(const char* filename) {
    const char* dot = strrchr(filename, '.');
    if (!dot || dot == filename) return "";
    return strdup(dot);
}

const char* baseFilename(char const *path) {
    char *bn = strrchr(path, '/');
    if (bn == NULL)
        return strdup(path);
    else
        return strdup(bn+1);
}

char* getReferenceFaiPath(const char *fn_ref) {
    char *fn_fai = 0;
    if (fn_ref == 0) return 0;

    fn_fai = calloc(strlen(fn_ref) + 5, sizeof(char));
    strcat(strcpy(fn_fai, fn_ref), ".fai");

    if (access(fn_fai, R_OK) == -1) { // fn_fai is unreadable
        if (access(fn_ref, R_OK) == -1) {
            fprintf(stderr, "fail to read file \n%s\n in getReferenceFaiPath() \n", fn_ref);
        }
        fprintf(stderr, "fail to read file \n%s\n in getReferenceFaiPath() \n", fn_fai);
        free(fn_fai); fn_fai = 0;
    }

    return fn_fai;
}

void checkNamingConvention(bam_hdr_t *header, khash_t(khStrInt)* wanted_chromosome_hash) {
    bool match=false;
    int32_t i=0;
    for(i=0; i<header->n_targets; i++) {
        khiter_t iter = kh_get(khStrInt, wanted_chromosome_hash, header->target_name[i]);
        if (iter != kh_end(wanted_chromosome_hash)) {
            match=true;
            break;
        }
    }

    if (!match) {
        fprintf(stderr, "ERROR: The naming conventions for chromosomes are different\n");
        fprintf(stderr, "between the chromosome bed file and the input bam/cram file\n\n");
        exit(EXIT_FAILURE);
    }
}

void checkChromosomeID(char* chrom_id1, char* chrom_id2) {
    // make a copy and convert CHR in tmp_chrom_id to lowercase
    //
    char * tmp_chrom_id = calloc(strlen(chrom_id2) + 1, sizeof(char));
    strcpy(tmp_chrom_id, chrom_id2);

    int i=0;
    for (i=0; tmp_chrom_id[i]; i++) {
        if (tmp_chrom_id[i] == 'C' || tmp_chrom_id[i] == 'H' || tmp_chrom_id[i] == 'R')
            tmp_chrom_id[i] = tolower(tmp_chrom_id[i]);
    }

    if (strcmp(chrom_id1, "hg37") == 0 || strcmp(chrom_id1, "hg19") == 0) {
        // ignore the presence of 'chr' if it is in the middle of chr name
        //
        char * find = strstr(tmp_chrom_id, "chr");
        if (find != NULL && find == tmp_chrom_id) {
            fprintf(stderr, "ERROR: The chromosome ID for Human Genome hg37/hg19 shouldn't contain 'chr'\n!");
            exit(EXIT_FAILURE);
        }
    }

    if (strcmp(chrom_id1, "hg38") == 0 && strstr(tmp_chrom_id, "chr") == NULL) {
        fprintf(stderr, "ERROR: The chromosome ID for Human Genome hg38 should begin with 'chr'!");
        exit(EXIT_FAILURE);
    }

    if (tmp_chrom_id != NULL) free(tmp_chrom_id);
}

uint64_t loadWantedChromosomes(khash_t(khStrInt) *wanted_chromosome_hash, char* in_version, char* in_bed_file) {
    FILE * fp = fopen(in_bed_file, "r");
    size_t len = 0;
    ssize_t read;
    char *line = NULL, *tokPtr=NULL;
    char *chrom_id = NULL;
    uint64_t total_genome_bases = 0;

    while ((read = getline(&line, &len, fp)) != -1) {
        if (*line == '\n') continue;                // handle empty lines
        if (strstr(line, "#") != NULL) continue;    // commented out lines

        char *savePtr = line;
        int absent=0;
        khiter_t iter=0;
        uint32_t start=0, end=0;

        // now need to process input information with the following format
        // chr1    0       248956422
        //
        uint16_t j=0;
        while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
            if (j == 0) {
                // check to see if the hashkey exists for current chromosome id
                //
                iter = kh_put(khStrInt, wanted_chromosome_hash, tokPtr, &absent);

                // here the string key is not kept elsewhere, we have to make a deep copy.
                // as the function kh_put only use the  pointer to the tmp variable tokPtr
                //
                if (absent) {
                    kh_key(wanted_chromosome_hash, iter) = strdup(tokPtr);
                }
                kh_value(wanted_chromosome_hash, iter) = 1;

                if (chrom_id == NULL) {
                    dynamicStringExpansion(tokPtr, &chrom_id);
                }
            }

            if (j == 1) {
                start = (uint32_t) strtol(tokPtr, NULL, 10);
            }

            if (j == 2) {
                end = (uint32_t) strtol(tokPtr, NULL, 10);
                kh_value(wanted_chromosome_hash, iter) = end-start;

                // update the total genome size info
                //
                total_genome_bases += end - start;
                //printf("%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", start, end, stats_info->wgs_cov_stats->total_genome_bases);
            }
            j++;
        }
    }

    if (line != NULL) free(line);

    // check version using chrom_id (for human genome only)
    // Since the MySQL database won't be used for many cases, we can't rely on users to use -D option
    // As the Default reference setting is hg37, we will only handle cases where it is hg38.
    // Here we will dynamically check the chromosome format and set the -D option
    //
    if (strcasecmp(in_version, "hg38") != 0) {
        char * find = strstr(chrom_id, "chr");
        if (find != NULL && find == chrom_id ) {
            fprintf(stderr, "The reference version isn't set correctly\n");
            fprintf(stderr, "The current reference version was \n%s\n", in_version);
            exit(EXIT_FAILURE);
        }
    }

    if (chrom_id != NULL) free(chrom_id);
    fclose(fp);

    return total_genome_bases;
}

uint64_t loadGenomeInfoFromBamHeader(khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header, char* in_version) {
    int i, absent=0;
    uint64_t total_genome_bases = 0;

    for ( i = 0; i < header->n_targets; i++) {
        khiter_t iter = kh_put(khStrInt, wanted_chromosome_hash, header->target_name[i], &absent);
        if (absent) {
            kh_key(wanted_chromosome_hash, iter) = strdup(header->target_name[i]);
            kh_value(wanted_chromosome_hash, iter) = header->target_len[i];
        }
        total_genome_bases += header->target_len[i];
    }

    // Now need to find out the reference version
    //
    for ( i = 0; i < header->n_targets; i++) {
        if (strcasecmp(in_version, "hg38") != 0) {
            // sometimes, 'chr' could be appeared in the middle of chromosome name
            // therefore, here we only check if 'chr' is at the very beginning of the chromosome name
            //
            char * find = strstr(header->target_name[i], "chr");
            if (find != NULL && find == header->target_name[i]) {
                fprintf(stderr, "The reference version %s isn't set correctly\n", in_version);
                exit(EXIT_FAILURE);
            }
        }
    }

    return total_genome_bases;
}

void stringArrayInit(StringArray *string_array, uint32_t size_in) {
    string_array->capacity = size_in;
    string_array->size = 0;
    string_array->theArray = calloc(size_in, sizeof(char*));
}

void stringArrayDestroy(StringArray *arrayIn) {
    uint16_t i;
    for(i=0; i<arrayIn->size; i++) {
        if (arrayIn->theArray[i]) {
            free(arrayIn->theArray[i]);
        }
    }

    if (arrayIn->theArray != NULL)
        free(arrayIn->theArray);
}

// I found this on the internet
//
char* stristr( const char* str1, const char* str2 ) {
    const char* p1 = str1 ;
    const char* p2 = str2 ;
    const char* r = *p2 == 0 ? str1 : 0 ;

    while ( *p1 != 0 && *p2 != 0 ) {
        if ( tolower( (unsigned char)*p1 ) == tolower( (unsigned char)*p2 ) ) {
            if ( r == 0 ) {
                r = p1 ;
            }

            p2++ ;
        } else {
            p2 = str2 ;
            if ( r != 0 ) {
                p1 = r + 1 ;
            }

            if ( tolower( (unsigned char)*p1 ) == tolower( (unsigned char)*p2 ) ) {
                r = p1 ;
                p2++ ;
            } else {
                r = 0 ;
            }
        }

        p1++;
    } 
    
    return *p2 == 0 ? (char*)r : 0 ;
}

void dynamicStringExpansion(char *str_in, char **storage_str) {
    if (str_in == NULL) return;     // nothing to add

    int original_str_is_null = 1;   // the original *storage_str is NULL? 1 yes, 0 no

    if (*storage_str) {
        *storage_str = realloc(*storage_str, (strlen(*storage_str) + strlen(str_in) + 2)*sizeof(char));
        original_str_is_null = 0;
    } else {
        *storage_str = calloc(strlen(str_in) + 1, sizeof(char));
    }

    if (*storage_str == NULL) {
        fprintf(stderr, "ERROR: Dynamic Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    if (original_str_is_null == 0) {
        strcat(*storage_str, str_in);
    } else {
        strcpy(*storage_str, str_in);
    }
}


uint32_t getLineCount(char *bed_file) {
    FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("ERROR: target file \n%s\n open failed!", bed_file);
        exit(EXIT_FAILURE);
    }

    // Extract characters from file and store in character c
    //
    long count = 0; // To store total line number in the target file

    /*char c;         // To store a character read from file
    for (c = getc(fp); c != EOF; c = getc(fp)) {
        if (c == '\n')      // Increment count if this character is newline
            count = count + 1;
    }*/

    ssize_t read;
    size_t len = 0;
    char *line=NULL;

    while((read = getline(&line, &len, fp)) != -1) {
        if (*line == '\n')   continue;
        if (line[0] == '\0') continue;        // skip if it is a blank line
        if (line[0] == '#')  continue;        // skip if it is a comment line

        count = count + 1;
    }

    if (line != NULL) free(line);

    fclose(fp);     // Close the file

    return count;
}

uint32_t loadBedFiles(char* ref_version, char *bed_file, Bed_Info * bed_info, khash_t(khStrInt)* wanted_chromosome_hash, char *db_version) {
    FILE *fp = fopen(bed_file, "r");;

    char *error_message = calloc(strlen(bed_file)+100, sizeof(char));    // full path needs enough spaces for it
    if (fp == NULL) {       // ensure the target file open correctly!
        sprintf(error_message, "target or N region bed file \n%s\n open failed!.\n", bed_file);
        exitWithFailure(error_message);
    }

    // setup a variable to store chromosomes that have been seen
    // this is used to check if the bed file is sorted!
    // if the file is not sorted, exit and give an error message!
    //
    khash_t(khStrInt) *seen_chromosome_hash = kh_init(khStrInt);

    //only temp setup the memory, will be adjust dynamically later
    //
    char *prev_chr_id = calloc(50, sizeof(char));
    strcpy(prev_chr_id, "nothing99999");
    uint32_t prev_start=0, prev_end=0;
    bool sorted=true;

    uint32_t count=0;           // index for each target line (item)
    uint32_t total_size=0;      // used for target bed input file merge and uniq checking
    ssize_t read;
    size_t len = 0;
    char *p_token=NULL, *line=NULL;    // here p_ means a point. so p_token is a point to a token
    char *savePtr;

    while((read = getline(&line, &len, fp)) != -1) {
        //printf("%s\n", line);
        if (line[0] == '\0') continue;        // skip if it is a blank line
        if (line[0] == '#')  continue;        // skip if it is a comment line

        savePtr = line;

        // to keep track the tokens from string split using strtok_r()
        //
        uint32_t i = 0;

        // loop through the rest of items, but here we are only interested in start and end position
        //
        while ((p_token = strtok_r(savePtr, "\t", &savePtr))) {
            // get the first item, which is chromosome id
            //
            if (i == 0) {
                bed_info->coords[count].chrom_id=NULL;        // initialization

                // check the version before continue
                //
                checkReferenceVersion(p_token, db_version, bed_file);

                // skip if the chromosome is not going to be processed
                //
                if (wanted_chromosome_hash != NULL) {
                    khiter_t iter_p = kh_get(khStrInt, wanted_chromosome_hash, p_token);
                    if (iter_p == kh_end(wanted_chromosome_hash))
                        break;
                }

                checkChromosomeID(ref_version, p_token);       // check chromosome id for versioning
                dynamicStringExpansion(p_token, &bed_info->coords[count].chrom_id);
        
                if (bed_info->coords[count].chrom_id && strlen(bed_info->coords[count].chrom_id) == 0)
                    break;
            }

            if (i == 1)
                bed_info->coords[count].start = (uint32_t) strtol(p_token, NULL, 10);

            if (i == 2)
                bed_info->coords[count].end = (uint32_t) strtol(p_token, NULL, 10);

            i++;
        }

        // when using a user-spcified chr-list, it is possible that the next line is not part of the list
        // so the chrom_id will be NULL. Therefore, we need to check it here
        //
        if (bed_info->coords[count].chrom_id==NULL) continue;

        // checking if the bed file is sorted
        //
        //fprintf(stderr, "%s\t%s\n", bed_info->coords[count].chrom_id, prev_chr_id);
        if (strcmp(bed_info->coords[count].chrom_id, prev_chr_id) == 0) {
            if (prev_start >= bed_info->coords[count].start)
                sorted = false;

            if (prev_end >= bed_info->coords[count].start)
                sorted = false;
        } else {
            // its a new chromosome
            // need to check if we have seen the chromosome id before
            //
            int absent;
            khiter_t iter = kh_put(khStrInt, seen_chromosome_hash, bed_info->coords[count].chrom_id, &absent);

            // here the string key is not kept elsewhere, we have to make a deep copy.
            // as the function kh_put only use the same pointer to the variable bed_info->coords[count].chrom_id
            // 
            if (absent) {
                // haven't seen it, just add to the seen_chromosome_hash
                //
                kh_key(seen_chromosome_hash, iter) = strdup(bed_info->coords[count].chrom_id);
                kh_value(seen_chromosome_hash, iter) = 1;

                // need to free the prev_chr_id, otherwise, it will append to it
                //
                if (prev_chr_id) {
                    free(prev_chr_id);
                    prev_chr_id=NULL;
                }
                dynamicStringExpansion(bed_info->coords[count].chrom_id, &prev_chr_id);
            } else {
                // seen the current chromosome before, which mean chr id is not in order
                //
                sorted = false;
            }
        }

        if (!sorted) {
            fprintf(stderr, "ERROR: The input bed file \n%s\n", bed_file);
            fprintf(stderr, "\tis not sorted or merged or unique between chrom %s and %s!\n", prev_chr_id, bed_info->coords[count].chrom_id);
            fprintf(stderr, "\tis not sorted or merged or unique between %"PRIu32" and %"PRIu32"!\n", prev_start, bed_info->coords[count].start);
            fprintf(stderr, "Please make sure your input bed file is sorted, merged and unique!\n");
            exit(EXIT_FAILURE);
        }

        // reset the prev_start position and sum the total bed file size
        //
        prev_start = bed_info->coords[count].start;
        prev_end   = bed_info->coords[count].end;
        total_size += bed_info->coords[count].end - bed_info->coords[count].start;

        count++;
    }

    bed_info->size = count;

    if (line != NULL) free(line);
    if (error_message != NULL) free(error_message);
    if (prev_chr_id != NULL) free(prev_chr_id);
    if (seen_chromosome_hash != NULL)
        cleanKhashStrInt(seen_chromosome_hash);

    fclose(fp);

    return total_size;
}

// the following comparison is used to compare the uint32_t array for qsort()
//
int compare(const void * val1, const void * val2) {
    uint32_t tmp_val1 = *((uint32_t*) val1);
    uint32_t tmp_val2 = *((uint32_t*) val2);

    if (tmp_val1 == tmp_val2) return 0;
    else if (tmp_val1 < tmp_val2) return -1;
    else return 1;
}

void checkReferenceVersion(char *chrom_id, char *db_version, char *file_in) {
    char *error_message = calloc(strlen(file_in)+150, sizeof(char));    // full path needs enough spaces for it
    sprintf(error_message, "The reference version for file %s doesn't match that from the user option.\n", file_in);

    if (strstr(chrom_id, "chr") == NULL) {
        if (strcmp(db_version, "hg38") == 0)
            exitWithFailure(error_message);
    } else {
        if ( (strcmp(db_version, "hg37") == 0) || (strcmp(db_version, "hg19") == 0) )
            exitWithFailure(error_message);
    }

    if (error_message) free(error_message);
}

void exitWithFailure(char* message) {
    fprintf(stderr, "ERROR: %s\n", message);
    if (message != NULL) free(message);
    exit(EXIT_FAILURE);
}
