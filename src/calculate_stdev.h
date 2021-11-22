/*
 * =====================================================================================
 *
 *       Filename:  calculate_stdev.h
 *
 *    Description:  calculate one pass standard deviation
 *
 *        Version:  1.0
 *        Created:  08/30/2021 02:14:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang at phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#ifndef CALCULATE_STDEV_H
#define CALCULATE_STDEV_H

#include <omp.h>
#include "terms.h"
#include "stats.h"
#include "utils.h"

#include "breakpoints.h"
#include "data_structure.h"
#include "coverage_tracking.h"
#include "utility.h"

void OnePassCalculateSedev(User_Input *user_inputs, bam_hdr_t **header, hts_idx_t **sfh_idx, samFile **sfh, Bed_Info *excluded_bed_info,  Simple_Stats *simple_stats, Target_Buffer_Status *target_buffer_status, Breakpoint_Array *breakpoint_array, Paired_Reads_Across_Breakpoints_Array *preads_x_bpts_array);

void OnePassStdevInit(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking);

void OnePassStdevDestroy(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking);

void StdevCalculation(OnePassStdev **one_pass_stdev, Chromosome_Tracking *chrom_tracking, Simple_Stats *simple_stats);

void SimpleStatsInit(Simple_Stats *simple_stats);

void get_coverage_info(Chromosome_Tracking *chrom_tracking, int32_t chrom_idx, OnePassStdev *one_pass_stdev);


#endif
