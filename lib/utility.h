/*
 * =====================================================================================
 *
 *       Filename:  utility.h
 *
 *    Description:  Include many useful functions
 *
 *        Version:  1.0
 *        Created:  03/31/2021 10:49:08 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX
 *
 * =====================================================================================
 */

#ifndef UTILITY_H
#define UTILITY_H

#include <ctype.h>      // for isdigit()
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>    // for bool return type

#include "data_structure.h"
#include "coverage_tracking.h"

/*
 * makes sure a file exists
 * @param fName (file name to check)
 * @return bool => return true if file exist, otherwise, exit(EXIT_FAILURE)
 */
bool checkFile(char * fName);

/*
 * check the size of input bam/cram file
 */
uint64_t check_file_size(const char *filename);

/*
 * check if the input string is a number
 * @param inStr: a variable need to be checked
 * @return true or false
 */
bool isNumber(const char *inStr);

/*
 * check if the input string is a float number
 * @param str, the string to be coverted
 * @param dest, the value to be copied to
 * @return true or false
 */
bool isFloat(const char *str, float *dest);

/*
 * this is a helper function that is used to create a file name from the existing base file name
 * @param output_dir: the directory name to store all the output reports/files
 * @param base_name: the file name to be created based on the base_name
 * @param file_in: the file_name to be created into the output_directory
 * @param string_to_append: the string to be attached to the base name
 */
void createFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append, char *version);

/*
 * This function is used to check if the input bam file named as .cram; or the input cram file named as .bam
 * @param user_inputs: it store the file name information
 * @param std: a samFile open file hand pointer
 */
void checkFileExtension(char* in_bam_file, samFile* sfd);

/*
 * This function is used to obtain the input filename extension
 * @param filename: the input filename
 * @return char*: the file extension
 */
const char* getFileExtension(const char* filename);

/*
 * Get the base name of the file including the extension, but without the full path
 * @param filepath: a full file path
 * @return char*: the basename from the file path
 */
const char* baseFilename(char const *filepath);

/*
 * The method is added for cram file input as cram file processing needs reference sequence
 * @param fn_ref: the file path to the reference sequences
 * @return the file path to the reference sequence's fai file
 */
char* getReferenceFaiPath(const char *fn_ref);

/* need to check for the naming convention between the chromosome bed file and the input bam/cram files
 * @param header: the bam/cram header to store the bam/cram header info
 * @param wanted_chromosome_hash: a kh_hash table to store the chromosome ids that need to be processed
 */
void checkNamingConvention(bam_hdr_t *header, khash_t(khStrInt)* wanted_chromosome_hash);

/*
 * since hg37 and hg38 used different chromosome nameing convention, 
 * we need to ensure that the input files follow this convention
 * @param chrom_id1: chromosome id to be checked
 * @param chrom_id2: chromosome id to be checked against
 */
void checkChromosomeID(char* chrom_id1, char* chrom_id2);

/*
 * it is a help function that is used to load all chromosomes need to be processed
 * @param primary_chromosome_hash: a hash table that is used to store the primary chromosomes for quick lookup
 * @param in_version: reference version to be checked
 * @parm in_bed_file: the bed file to be processed
 */
uint64_t loadWantedChromosomes(khash_t(khStrInt) *wanted_chromosome_hash, char* in_version, char* in_bed_file);

/*
 * calculate the total number of genome bases and load chromosome info from the bam/cram header info
 * @param header: bam/cram/sam header info
 * @param stats_info: to store the final number of total genome bases
 */
uint64_t loadGenomeInfoFromBamHeader(khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header, char* in_version);

/*
 * open the input bed-formated file and then load the coordinate regions into memory and leave them there
 * Note: the input file should be in .bed format
 * @param ref_version: reference version used
 * @param bed_file: input name for the bed file
 * @param bed_info: the Bed_Info structure variable to store coordinates of each bed section
 * @param db_version: the reference version used
 * @return the total number of bases covered in the target/N-regions bed file
 */
uint32_t loadBedFiles(char *ref_version, char *bed_file, Bed_Info * bed_info, khash_t(khStrInt)* wanted_chromosome_hash, char *db_version);

/*
 * To initializate a StringArray variable
 * @param string_array: a pointer to a StringArray variable to be initialized
 * @param size_in: the capacity of the StringArray variable to be initialized
 */
void stringArrayInit(StringArray *string_array, uint32_t size_in);

/*
 * This function is used to clean-up the StringArray variable
 * @param arrayIn: a StringArray variable to be cleaned
 */
void stringArrayDestroy(StringArray *arrayIn);

/*
 * find if a string contains a substr (case insensative)
 * it follows the behaviors per standard strstr()
 * @param string1, the string to be checked
 * @param substr2, the substring to search for
 * @return a string starts from the substr substr2 in the original string1 to the end of original string1
 */
char* stristr( const char* string1, const char* substr2 );

/*
 * open the input file (bed-formatted) and count how many items within the bed (such as target) files
 * @param bed_file: file in bed-formatted such as (target file or the Ns regions of reference sequences in bed format)
 * @return count    =>  Total number of targets or Ns regions
 */
uint32_t getLineCount(char *bed_file);

/**
 * This function is used to dynamically allocate string and grow it accordingly through either copy or concatenation
 * @param str_in, the string to be added
 * @param storage_str, the string to be dynamically allocated/expanded
 */
void dynamicStringExpansion(char *str_in, char **storage_str);

/*
 * the following comparison is used to compare the int elements in an array used by qsort()
 * @param val1: int array 1 used for comparison
 * @param val2: int array 2 used for comparison
 * @return: return value means
 * <0 The element pointed by val1 goes before the element pointed by val2
 * 0  The element pointed by val1 is equivalent to the element pointed by val2
 * >0 The element pointed by val1 goes after the element pointed by val2
 *
 */
int compare(const void * val1, const void * val2);

void checkReferenceVersion(char *chrom_id, char *db_version, char *file_in);

void exitWithFailure(char* message);

#endif
