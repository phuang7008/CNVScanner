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

#include "terms.h"
#include "storage.h"

typedef struct {
    uint32_t start;
    uint32_t end;
    uint32_t length;
    double ave_coverage;
    char cnv_type;      // L: for deletion, while P for Dup

    // segmentation details
    //
    uint32_t seg_start;
    uint32_t seg_end;
    uint32_t seg_length;
    double seg_ave_coverage;

    // store CNVs that included in this final segment
    //
    uint32_t size;
    uint32_t capacity;
    CNV* cnvs;

    // Used breakpoints
    //
    uint32_t *left_breakpoints;
    uint32_t *right_breakpoints;

    // Used paired-reads across a breakpoint
    //
    uint16_t *num_of_left_paired_reads;
    uint16_t *num_of_right_paired_reads;

    // store related nearby breakpoint info
    //
    CNV_Breakpints *cnv_breakpoints;
    uint32_t cnv_breakpoints_size;
    uint32_t cnv_breakpoints_capacity;
    int16_t left_start_index;            // the signed index is set when there is a left-hand breakpoint (0-index is valid)
    int16_t right_end_index;             // the signed index is set when there is a right-hand breakpoint (0-index is valid)

    // store improperly paired reads with perfect mapping and TLEN > 1000 and count >= 2
    //
    uint32_t imp_PR_start;
    uint32_t imp_PR_end;
    uint16_t num_of_imp_RP_TLEN_1000;   // number of improperly paired-reads with TLEN >= 1000
} Segmented_CNV;

typedef struct {
    char *chromosome_id;
    uint32_t chrom_length;
    uint32_t size;
    uint32_t capacity;
    Segmented_CNV* seg_cnvs;              // CNV array per chromosome
} Segmented_CNV_Array;

void segmentArrysInit(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking);

void cnv_segmentation(Chromosome_Tracking *chrom_tracking, Segment_Array** segment_array, User_Input *user_inputs);

void process_segmentation_data(Segment_Array **segment_array, CNV_Array **cnv_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, khash_t(m32) **anchor_breakpoints_hash_array, User_Input *user_inputs);

#endif
