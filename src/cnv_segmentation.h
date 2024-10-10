/*
 * =====================================================================================
 *
 *       Filename:  cnv_segmentation.h
 *
 *    Description:  header file for cnv segmentation
 *
 *        Version:  1.0
 *        Created:  02/02/2024 02:26:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, pemhuang@gmail.com
 *        Company:  Burnaby, British Columbia, Canada
 *
 * =====================================================================================
 */

#ifndef CNV_SEGMENTATION_H
#define CNV_SEGMENTATION_H

#include <libgen.h>        // for function basename()
#include <omp.h>
#include <stdlib.h>

#include "fileProcessing.h"
#include "terms.h"
#include "storage.h"
#include "utils.h"
#include "utility.h"

typedef struct {
    INNER_CNV *seg_inner_cnv;
    int16_t seg_inner_cnv_size;
    int16_t seg_inner_cnv_capacity;

    // segmentation details
    //
    Segment_Details segment;

    // store CNVs that included in this final segment
    //
    uint32_t seg_cnv_index_size;
    uint32_t seg_cnv_index_capacity;
    uint32_t *seg_cnv_indices;        // array of index pointing to the position of each CNV in the CNV_Array

    // store related nearby breakpoint info
    //
    CNV_Breakpints *seg_breakpoints;
    uint32_t seg_breakpoints_size;
    uint32_t seg_breakpoints_capacity;
    int16_t seg_left_start_index;       // the signed index is set when there is a left-hand breakpoint (0-index is valid)
    int16_t seg_right_end_index;        // the signed index is set when there is a right-hand breakpoint (0-index is valid)

    // store improperly paired reads with perfect mapping and TLEN > 1000 and count >= 2
    //
    uint32_t seg_imp_PR_start;
    uint32_t seg_imp_PR_end;
    uint16_t seg_num_of_imp_PR_TLEN_1000;   // number of improperly paired-reads with TLEN >= 1000

    char* excluded_regions;                 // it contains Ns-regions, segdup and tandom repeat >= 10kb
} Segmented_CNV;

typedef struct {
    char *chromosome_id;
    uint32_t chrom_length;
    uint32_t size;
    uint32_t capacity;
    Segmented_CNV* seg_cnvs;              // CNV array per chromosome
} Segmented_CNV_Array;

void segmentArrysInit(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking);

void segmentArrysDestroy(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking);

void SegmentedCNVArrayInit(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking);

void SegmentedCNVArrayDestroy(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking);

void cnvSegmentation(Chromosome_Tracking *chrom_tracking, Segment_Array** segment_array, User_Input *user_inputs);

void storeSegmentsLocallyAndInit(Segment_Array** segment_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking);

void processSegmentationData(CNV_Array **cnv_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, khash_t(m32) **anchor_breakpoints_hash_array, User_Input *user_inputs, Simple_Stats *equal_window_stats, Stats_Info *stats_info, Bed_Info *low_mappability_bed_info, Bed_Info *gc_lt25pct_bed_info, Bed_Info *gc_gt85pct_bed_info);

void findIntersectedCNVs(CNV_Array *cnv_array, Segmented_CNV_Array *seg_cnv_array, User_Input *user_inputs);

void addCNVindexToSegment(uint32_t cnv_index, Segmented_CNV_Array *seg_cnv_array, uint32_t seg_cnv_index);

void checkBreakpointsForEachSegment(Segmented_CNV_Array *seg_cnv_array, khash_t(m32) *anchor_breakpoints_hash_array);

void addBreakpoints(Segmented_CNV_Array *seg_cnv_array, uint32_t seg_index, khash_t(m32) *anchor_breakpoints_hash, uint32_t anchor_breakpoint);

void setBreakpoinsAtSegmentEnds(Segmented_CNV_Array *seg_cnv_array);

void removeExcludedRegionsFromSegments(Segmented_CNV_Array *seg_cnv_array, Binned_Data_Wrapper *excluded_regions, int total_size);

void saveExcludedRegion(Segmented_CNV_Array *seg_cnv_array, uint32_t seg_cnv_index, char *regions_to_append);

void segmentInnerCNVInit(Segmented_CNV *seg_cnv, uint32_t capacity);

void segmentInnerCNVInit2(Segmented_CNV *seg_cnv, uint32_t index);

void copyCNV(Segmented_CNV *seg_cnv, uint32_t seg_index, CNV *cnv, uint32_t cnv_index);

void dynamicMemAllocateInnerCNVbreakpoints(INNER_CNV *seg_inner_cnv, uint32_t index);

int obtainSupportingEvidences(INNER_CNV* inner_cnv, int* total_breakpoints, double log2ratio, Simple_Stats *equal_window_stats, int type);

void mergeCNVsFromSameSegment(Segmented_CNV_Array *seg_cnv_array, CNV_Array *cnv_array, User_Input *user_inputs, Simple_Stats *equal_window_stats);

void checkingFalsePositives(Segmented_CNV_Array *seg_cnv_array, char * chr_id, Bed_Info *low_complexity_bed_info, int type);

void generateSegmentedCNVs(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, Simple_Stats *equal_window_stats, Stats_Info *stats_info, User_Input *user_inputs);

#endif
