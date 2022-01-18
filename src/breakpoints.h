/*
 * =====================================================================================
 *
 *       Filename:  breakpoints.h
 *
 *    Description:  the header file for breakpoints
 *
 *        Version:  1.0
 *        Created:  10/05/2021 01:54:20 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX, USA
 *
 * =====================================================================================
 */

#ifndef BREAKPOINTS_H
#define BREAKPOINTS_H

#include "data_structure.h"
#include "coverage_tracking.h"

#include "terms.h"
#include "utils.h"
#include "utility.h"

void BreakpointArrayInit(Breakpoint_Array *breakpoint_array, Chromosome_Tracking *chrom_tracking);

void BreakpointArrayDestroy(Breakpoint_Array *breakpoint_array);

void PairedReadsAcrossBreakpointsArrayInit(Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array, Chromosome_Tracking *chrom_tracking);

void PairedReadsAcrossBreakpointsArrayDestroy(Paired_Reads_Across_Breakpoints_Array *pread_x_bpts_array);

void PairedReadsAcrossABreakpointPerAnchorArrayDestroy(khash_t(khIntPrArray) * preads_x_per_anchor_arr);

void BreakpointStatsArrayInit(Breakpoint_Stats_Array *bpt_stats_array, Chromosome_Tracking *chrom_tracking);

void BreakpointStatsArrayDestroy(Breakpoint_Stats_Array *bpt_stats_array);

uint32_t fetchBreakpointArrayChrIndex(Breakpoint_Array *breakpoint_array, char * chrom_id);

uint32_t fetchPReadsXBreakpointArrayChrIndex(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, char * chrom_id);

void storeCurrentReadBreakpointInfo(uint32_t current_ref_pos, bam1_t *rec, Breakpoint_Array *breakpoint_array, uint32_t breakpoint_chr_index, khash_t(khStrInt) *breakpoint_pairs_hash, int type);

void storePairedReadsAcrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, Paired_Reads_Across_Breakpoints_Array *preads_x_bpts_arr, uint32_t pr_chr_ind, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh);

void eliminateUnwantedBreakpoints(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t pr_chr_ind);

void outputBreakpointArray(Breakpoint_Array *bpt_arr);

void outputPairedReadsAcrossBreakpointsArray(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr);

void dynamicBreakpointPerChrArraySizeIncrease(Breakpoints_Per_Chromosome *bpts_per_chr);

void dynamicPairedReadsAcrossABreakpointArraySizeIncrease(Paired_Reads_Across_Per_Anchor_Breakpoint_Array *preads_x_bpts_arr);

void dynamicBreakpointStatsPerChrSizeIncrease(Breakpoint_Stats_Per_Chromosome *bp_stats_per_chr);

void getSortedBreakpointArray(uint32_t *sorted_breakpoints, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx);

void recordAnchorBreakpoints(uint32_t *sorted_breakpoints, khash_t(m32) *anchor_breakpoints_hash, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx);

#endif
