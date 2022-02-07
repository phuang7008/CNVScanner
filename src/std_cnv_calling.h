/*
 * =====================================================================================
 *
 *       Filename:  std_cnv_calling.h
 *
 *    Description:  Using Standard Deviation approach for CNV calling
 *
 *        Version:  1.0
 *        Created:  01/20/2022 09:23:38 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#ifndef STD_CNV_CALLING_H
#define STD_CNV_CALLING_H

#include "stats.h"
#include "terms.h"
#include "utils.h"

void generateCNVs(CNV_Array **equal_bin_cnv_array, Binned_Data_Wrapper **equal_size_window_wrappers, Binned_Data_Wrapper **raw_bin_data_wrappers, uint32_t number_of_chromosomes, Simple_Stats *equal_window_stats, User_Input *user_inputs);

void mergeNeighboringBinsBasedOnZscore(CNV_Array *cnv_array, Binned_Data_Wrapper *equal_size_window_wrapper, Simple_Stats *the_stats, int type);

void storeCurrentCNVtoArray(CNV *cnv, uint32_t start, uint32_t end, uint32_t length, double coverage, Equal_Window_Bin *merged_equal_bin_array, uint32_t bin_size);

void expandMergedCNVWithRawBins(Binned_Data_Wrapper *binned_data_wrapper, CNV_Array *cnv_array, Simple_Stats *equal_window_stats);

void addRawBinToCNV(Binned_Data_Wrapper *binned_data_wrapper, uint32_t raw_bin_index, CNV *cnv, uint32_t cnv_index);

void outputCNVArray(CNV_Array **cnv_array, uint32_t number_of_chromosomes, int type);

void dynamicIncreaseBinArraySize(Equal_Window_Bin **merged_equal_bin_array, uint32_t bin_capacity);

void cnvArrayInit(CNV_Array **cnv_array, Chromosome_Tracking *chrom_tracking);

void cnvArrayDestroy(CNV_Array **cnv_array, uint32_t number_of_chromosomes);

#endif
