/*
 * =====================================================================================
 *
 *       Filename:  improperly_paired_reads.h
 *
 *    Description:  record all improperly paired reads for CNV validation
 *
 *        Version:  1.0
 *        Created:  01/23/2023 08:00:54 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef IMPROPERLY_PAIRED_READS
#define IMPROPERLY_PAIRED_READS

#include "htslib/include/htslib/sam.h"

#include "data_structure.h"
#include "coverage_tracking.h"

#include "terms.h"
#include "utils.h"
#include "utility.h"

void NotProperlyPairedReadsInit(Not_Properly_Paired_Reads_Array** improperly_paired_reads_array, Chromosome_Tracking *chrom_tracking);

void NotProperlyPairedReadsDestroy(Not_Properly_Paired_Reads_Array** improperly_PR_array, Chromosome_Tracking *chrom_tracking);

uint32_t fetchImproperPRArrayChrIndex(Not_Properly_Paired_Reads_Array** improperly_paired_reads_array, Chromosome_Tracking* chrom_tracking, uint32_t chrom_index);

void processPairedReadsWithinTheSameGroup(Not_Properly_Paired_Reads_Array* improperly_PR_array);

void processImproperlyPairedReads(Not_Properly_Paired_Reads_Array* improperly_paired_reads_array, bam1_t *rec);

void organizeImproperlyPairedReadArray(Not_Properly_Paired_Reads_Array* improperly_PR_array);

int32_t getMateMatchLengthFromMCTag(char *mate_cigar);

void outputGroupedImproperlyPairedReads(Not_Properly_Paired_Reads_Array** improperly_paired_reads_array, Chromosome_Tracking *chrom_tracking);

#endif
