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

void BreakpointArrayInit(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking);

void BreakpointArrayDestroy(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking);

void PairedReadsAcrossBreakpointsArrayInit(Paired_Reads_Across_Breakpoints_Array **pread_x_bpts_array, Chromosome_Tracking *chrom_tracking);

void PairedReadsAcrossBreakpointsArrayDestroy(Paired_Reads_Across_Breakpoints_Array **pread_x_bpts_array, Chromosome_Tracking *chrom_tracking);

uint32_t fetchBreakpointArrayChrIndex(Breakpoint_Array **breakpoint_array, Chromosome_Tracking *chrom_tracking, uint32_t chrom_idx);

uint32_t fetchPReadsXBreakpointArrayChrIndex(Paired_Reads_Across_Breakpoints_Array **preads_x_bpt_arr, Chromosome_Tracking *chrom_tracking, uint32_t chrom_idx);

void storeCurrentReadBreakpointInfo(uint32_t current_ref_pos, bam1_t *rec, Breakpoint_Array *bpt_arr, int type);

void storePairedReadsAcrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, Paired_Reads_Across_Breakpoints_Array *preads_x_bpts_arr, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh, User_Input *user_inputs);

void eliminateUnwantedBreakpoints(char* chr_id, Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr, uint32_t num_of_anchors);

void outputBreakpointArray(Breakpoint_Array *bpt_arr);

void outputPairedReadsAcrossBreakpointsArray(Paired_Reads_Across_Breakpoints_Array *preads_x_bpt_arr);

void dynamicBreakpointPerChrArraySizeIncrease(Breakpoint_Array *bpts_per_chr);

void dynamicPairedReadsAcrossABreakpointArraySizeIncrease(Paired_Reads_Across_Per_Anchor_Breakpoint_Array *preads_x_bpts_arr);

void getSortedBreakpointArray(uint32_t *sorted_breakpoints, Breakpoint_Array *bpt_arr, User_Input *user_inputs);

uint32_t recordAnchorBreakpoints(uint32_t *sorted_breakpoints, khash_t(m32) *anchor_breakpoints_hash, Breakpoint_Array *bpt_arr, User_Input *user_inputs);

#endif
