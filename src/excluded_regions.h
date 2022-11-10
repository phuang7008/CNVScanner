/*
 * =====================================================================================
 *
 *      Filename:		targets.h
 *
 *		Description:	For the Capture/Target related functionalities
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef TARGETS_H
#define TARGETS_H

#include <stdbool.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "terms.h"
#include "user_inputs.h"
#include "utility.h"
#include "coverage_tracking.h"

/**
 * generate target and buffer lookup table for quick access. The khash.h file is used for hash table setup
 * @param bed_info: the bed information that is stored for the future usage
 * @param stat_info, statistical information for the reads/bases
 * @param header, the bam/sam/cram header pointer that hold the length info of each chromosome
 */
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, khash_t(khStrInt)* wanted_chromosome_hash);

/**
 * process bed-formatted file and populate the coordinates and lookup hash table
 * @param user_inputs: contains all the user inputs including the target or Ns bed file names
 * @param bed_info: the storage of bed coordinates and the size of the bed file
 * @param stats_info: a variable that contains various statistical information
 */
void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, Stats_Info *stats_info, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name);

void zeroAllExcludedRegions(Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Bed_Info *excluded_bed_info);

/**
 * just to output some information for debugging
 */
void outputForDebugging(Bed_Info *bed_info);

/**
 * It is used to store both target and buffer information
 * @param chromosome_id: the current chromosome to load
 * @param size: the size of the chromosome
 * @param target_buffer_regions: the 2D char array that store the flags either targets('T') or their surrounding buffers ('B')
 * @param target_coords: the struct array that stores the coordinates of targets
 * @return
 */
void getTargetAndBufferPositions(char * chromosome_id, int size, char *target_buffer_regions, Bed_Coords *coords);

#endif //TARGETS_H
