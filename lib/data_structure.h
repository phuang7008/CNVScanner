/*
 * =====================================================================================
 *
 *       Filename:  data_structure.h
 *
 *    Description:  To define various data types for the library
 *
 *        Version:  1.0
 *        Created:  03/25/2021 09:40:47 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <inttypes.h>   // for PRIu32 and PRIu64 
#include <stdbool.h>    // for bool definition
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>     // for getopt() and usleep()
#include <math.h>
#include <zlib.h>
#include <sys/stat.h>
#include <stdio.h>      // for file read and write

// Users defined header files
#include "htslib/include/htslib/sam.h"
#include "htslib/include/htslib/khash.h"

// The following are data stractures that are defined in the library
//

/* each application that is trying to use this library needs to 
   have this variable defined in the main() function
**/
extern int  khStrInt;
extern int  m32;
extern int  m16;
extern int  m8;

KHASH_MAP_INIT_STR(khStrInt, uint32_t)

// Instantiate a hash map containing integer keys
// m32 means the key is 32 bit integer, while the value is of unsigned int type (ie uint32_t)
//
KHASH_MAP_INIT_INT(m32, uint32_t)
KHASH_MAP_INIT_INT(m16, uint16_t)
KHASH_MAP_INIT_INT(m8, uint16_t)

/* each application that is trying to use this library needs to 
 *    have this variable defined in the main() function
 **/
// const int get_khStrInt() { return khStrInt; }
// void set_khStrInt(int val) { khStrInt = val; }

/**
 * define a structure that holds the target/excluded region's coordinates
 **/
typedef struct {
    char *chrom_id;    // some chromosome ID would be quite long
    uint32_t start;
    uint32_t end;
} Bed_Coords;

/**
 * define a structure that holds the size and coordinates info of a bed file
 **/
typedef struct {
    Bed_Coords *coords;
    uint32_t size;          // number of regions in bed format from the original file
} Bed_Info;

/**
 * define a structure to hold a chunk of read buff to be processed by each thread
 **/
typedef struct {
    bam1_t **chunk_of_reads;
    uint32_t capacity;
    uint32_t size;
} Read_Buffer;

/**
 * data structure to store temp coverage results
 **/
typedef struct {
    uint32_t size;
    uint32_t * cov_array;
} Temp_Coverage_Array;

/**
 * define a strcuture for quick lookup of target information
 **/
typedef struct {
    char *chrom_id;                 // which chromosome it is tracking
    uint32_t size;                  // size of current chromosome
    int32_t index;                  // -1 means this chromosome is not important/processed, and we should skip it!
    uint8_t *target_status_array;   // array to store the target status for each chrom position
    uint8_t *buffer_status_array;   // array to store the buffer status for each chrom position
} Target_Buffer_Status;

/**
 *  define a strcuture that hold the chromosome ids status in array 
 *  (that is all of the member variables are array, except number_tracked)
 **/
typedef struct {
    uint32_t number_of_chromosomes;     // number of chromosomes users are interested in
    uint32_t number_tracked;            // the number of chromosomes we are tracking so far!
    int32_t **coverage;                 // the coverage count info for each base on each chromosome will be stored here!
                                        // this needs to be signed as the excluded bases will be assigned to -1 in wgs_cnv, 0 in scandium

    char **chromosome_ids;
    uint32_t *chromosome_lengths;
    uint8_t  *chromosome_status;    // 0 pending, 1 working, 2 finish coverage calculation, 3. done summary-report/annotation processing!
    bool more_to_read;
} Chromosome_Tracking;

/** define stringArray structure to store the annotation information
 * one for RefSeq, one for CCDS, one for VEGA and one for Gencode, one for miRNA
 */
typedef struct {
    char **theArray;
    uint16_t capacity;
    uint16_t size;
} StringArray;

typedef struct {
    uint32_t * array;
    uint32_t size;
    uint32_t capacity;
} AllStartsEndsArray;

#endif
