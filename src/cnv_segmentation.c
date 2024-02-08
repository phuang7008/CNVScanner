/*
 * =====================================================================================
 *
 *       Filename:  cnv_segmentation.c
 *
 *    Description:  the implementation for the CNV segmentation
 *
 *        Version:  1.0
 *        Created:  02/02/2024 02:29:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, pemhuang@gmail.com
 *        Company:  Houston, TX
 *
 * =====================================================================================
 */

#include "cnv_segmentation.h"
#include "slmseg_connector.h"

void segmentArrysInit(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        segment_arrays[i] = calloc(1, sizeof(Segment_Array*));
        segment_arrays[i]->chrom_id = strdup(chrom_tracking->chromosome_ids[i]);

        // the following will be done at the connector/slmseg_connector.cpp
        //
        //segment_arrays[i]->size = 0;
        //segment_arrays[i]->capacity = INIT_SIZE;
        //segment_arrays[i]->segments = calloc(segment_arrays[i]->capacity, sizeof(Segment_Array));
    }
}

void segmentArrysDestory(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking) {
    uint32_t i; 
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (segment_arrays[i]->chrom_id) {
            free(segment_arrays[i]->chrom_id);
            segment_arrays[i]->chrom_id=NULL;
        }

        if (segment_arrays[i]->segments) {
            free(segment_arrays[i]->segments);
            segment_arrays[i]->segments = NULL;
        }
    }
}

void cnv_segmentation(Chromosome_Tracking *chrom_tracking, Segment_Array** segment_array, User_Input *user_inputs) {
    char *tmp_basename = basename(user_inputs->bam_file);

    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // process one chromosome at a time
        //
        char *input_file = calloc(5000, sizeof(char));
        sprintf(input_file, "%s%s%s", user_inputs->log2ratio_output_file, chrom_tracking->chromosome_ids[i], ".txt");

        char *output_file  = calloc(5000, sizeof(char));
        sprintf(output_file, "%s/%s_%s_%s", user_inputs->output_dir, tmp_basename, chrom_tracking->chromosome_ids[i], "segments.txt");
        slmseg_call(chrom_tracking->chromosome_ids[i], input_file, output_file, segment_array[i], 0.3, 0.00004, 1000000, 0);
        
        free(input_file);
        free(output_file);
        input_file = NULL;
        output_file = NULL;
    }

    free(tmp_basename);
}

void process_segmentation_data(Segment_Array **segment_array, CNV_Array **cnv_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, khash_t(m32) **anchor_breakpoints_hash_array, User_Input *user_inputs) {
#pragma omp parallel num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        uint32_t segmentation_index;
        for (segmentation_index=0; segmentation_index<chrom_tracking->number_of_chromosomes; ++segmentation_index) {
#pragma omp task
          {
            int thread_id = omp_get_thread_num();
            printf("Current thread id in process_segmentation_datacalls: %d with chr %s\n", thread_id, chrom_tracking->chromosome_ids[segmentation_index]);
            // find the corresponding index in cnv_array and seg_cnv_array
            //
            uint32_t cnv_array_idx, seg_cnv_array_idx;
            for (cnv_array_idx=0; cnv_array_idx<chrom_tracking->number_of_chromosomes; ++cnv_array_idx) {
                if (strcmp(cnv_array[cnv_array_idx]->chromosome_id, segment_array[segmentation_index]->chrom_id) == 0)
                    break;
            }

            for (seg_cnv_array_idx=0; seg_cnv_array_idx<chrom_tracking->number_of_chromosomes; ++seg_cnv_array_idx) {
                if (strcmp(seg_cnv_array[seg_cnv_array_idx]->chromosome_id, segment_array[segmentation_index]->chrom_id) == 0)
                    break;
            }

            // use intersect to find out which cnv segments are true CNVs
            //
            //find_true_cnvs();

          } // end task 
        } // end for loop
#pragma omp taskwait
      } // end single
    } // end parallel
}
