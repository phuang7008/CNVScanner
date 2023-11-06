/*
 * =====================================================================================
 *
 *       Filename:  coverage_tracking.h
 *
 *    Description:  To track coverage info for each base of each chromosome
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

#ifndef COVERAGE_TRACKING_H
#define COVERAGE_TRACKING_H

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
#include "data_structure.h"

/* each application that is trying to use this library needs to 
   have this variable defined in the main() function
const int get_khStrInt() { return khStrInt; }
void set_khStrInt(int val) { khStrInt = val; }
**/

/*
 * It is used to clean the kash_t (char* as key, with uint32_t as value) hash table
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories including keys
 **/
void cleanKhashStrInt(khash_t(khStrInt) *hash_to_clean);

/**
 * This function is used to clean the khash_t (uint32_t key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 **/
void cleanKhashInt(khash_t(m32) *hash_to_clean);

/**
 * This is used to initialize an array of read buffers of type bam1_t (pointers)
 * @param read_buff_in: the array of read_buff_in to be initialized
 * @param size: the size of the read_buff_in array
 **/
void readBufferInit(Read_Buffer *read_buff_in);

/**
 * This function is used to delete all the memory allocated for read_buff_in of type bam1_t
 * @param read_buff_in: the array of read_buff_in to be destroyed
 * @param size: the size of the read_buff_in array
 **/
void readBufferDestroy(Read_Buffer *read_buff_in);

/**
 * This function is used to read one chunk (10,000 -- 1,000,000) of aligned reads from bam into the memory 
 * @param sffh: samFile file hander: this is actually refer to the opened bam file handler
 * @param header: bam/cram file header information that has chromosome id and length info
 * @param chrom_tracking: it has a boolean flag (more_to_read) that is used to signal is bam file reading is finished (shared among all threads)
 * @param read_buff_in: the array of Read_Buffer to be used to store chunk of alignment reads to be processed by individual thread
 **/
uint32_t readBam(samFile *sffh, bam_hdr_t *header, Chromosome_Tracking *chrom_tracking,  Read_Buffer *read_buff_in);

/**
 * initialize the Chromosome_Tracking variable, this approach will process all chromosomes
 * @param header: a bam header that contains all the chromosome information
 * @return an instance of Chromosome_Tracking upon successful
 **/
void chromosomeTrackingInit1(Chromosome_Tracking *chrom_tracking, khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header);

/**
 * Initialize the chromosome_tracking variable using user specified region file
 * This approach will process those chromosomes specified by user only
 * @param wanted_chromosome_hash: a kh_hash table stores chromosomes to be processed
 **/
void chromosomeTrackingInit2(khash_t(khStrInt) *wanted_chromosome_hash, Chromosome_Tracking *chrom_tracking, bam_hdr_t *header);

/**
 * This function is used to update all members for the Chromosome_Tracking variable
 * @param chrom_tracking: a Chromosome_Tracking used for tracking
 * @param chrom_id: the current chromosome id
 * @param chrom_len: the length of current chromosome
 * @param index: the members in Chromosome_Tracking variable are stored in arrays, using index will help locate the chromosome info
 **/
void chromosomeTrackingUpdate(Chromosome_Tracking *chrom_tracking, uint32_t chrom_len, int index);

/**
 * To locate the index of a chromosome id in an array give the chromosome id information
 * @param chrom_id: current chromosome id to be handled
 * @param chrom_tracking: the Chromosome_Tracking variable to track the status of chromosome processed
 * @return a index at the tracking array
 **/
int32_t locateChromosomeIndexForChromTracking(char *chrom_id, Chromosome_Tracking *chrom_tracking);

/**
 * To clean up the allocated memory for chromosome tracking
 * @param chrom_tracking: the tracking variable to be cleaned
 **/
void chromosomeTrackingDestroy(Chromosome_Tracking * chrom_tracking);

/**
 * it will set the coverage for all of the Ns regions in the genome to zero
 * @param chrom_id: current chromosome id to be handled
 * @param Ns_info: the detailed Ns info
 * @param chrom_tracking: a storage used to track each chromosome in details
 * @param target_buffer_status: it is used to tell which regions are targets and which regions are buffer. 
 *      Sometimes, Ns/excluded regions will overlap with the buffer regions
 * @param value: the value to set to for Ns or excluded region bases
 **/
void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking, Target_Buffer_Status *target_buffer_status, int value);

/**
 * this function is used to free all the memories that are allocated by the program
 * @param bed_info: the declared Bed_Info variable
 **/
void cleanBedInfo(Bed_Info *bed_info);

void TargetBufferStatusInit(Target_Buffer_Status *target_buffer_status, bam_hdr_t *header);

void TargetBufferStatusInit2(Target_Buffer_Status *target_buffer_status, khash_t(khStrInt)* wanted_chromosome_hash);

void TargetBufferStatusUpdate(Target_Buffer_Status *target_buffer_status, int32_t target_buffer_index);

void TargetBufferStatusDestroyCurrentChromosome(Target_Buffer_Status *target_buffer_status, uint32_t number_of_chromosomes, char* chrom_id);

void TargetBufferStatusDestroy(Target_Buffer_Status *target_buffer_status, uint32_t number_of_chromosomes);


#endif
