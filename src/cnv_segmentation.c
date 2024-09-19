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
        segment_arrays[i] = calloc(1, sizeof(Segment_Array));
        segment_arrays[i]->chrom_id = strdup(chrom_tracking->chromosome_ids[i]);

        // the following will be done at the connector/slmseg_connector.cpp
        //
        //segment_arrays[i]->size = 0;
        //segment_arrays[i]->capacity = INIT_SIZE;
        //segment_arrays[i]->segments = calloc(segment_arrays[i]->capacity, sizeof(Segment_Array));
    }
}

void segmentArrysDestroy(Segment_Array **segment_arrays, Chromosome_Tracking *chrom_tracking) {
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

        if (segment_arrays[i])
            free(segment_arrays[i]);
    }

    if (segment_arrays)
        free(segment_arrays);
}

void SegmentedCNVArrayInit(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        seg_cnv_array[i] = calloc(1, sizeof(Segmented_CNV_Array));
        seg_cnv_array[i]->chromosome_id = strdup(chrom_tracking->chromosome_ids[i]);
        seg_cnv_array[i]->chrom_length  = chrom_tracking->chromosome_lengths[i];

        seg_cnv_array[i]->size = 0;
        seg_cnv_array[i]->capacity = PR_INIT_SIZE * 50;     // set size of 10,000
        seg_cnv_array[i]->seg_cnvs = calloc(seg_cnv_array[i]->capacity, sizeof(Segmented_CNV));
        
        // the rest initialization will be done at the storeSegmentsLocallyAndInit
    }
}

void SegmentedCNVArrayDestroy(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i, j, k;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        if (seg_cnv_array[i]->chromosome_id)
            free(seg_cnv_array[i]->chromosome_id);

        for (j=0; j<seg_cnv_array[i]->size; ++j) {
            if (seg_cnv_array[i]->seg_cnvs[j].seg_cnv_indices)
                free(seg_cnv_array[i]->seg_cnvs[j].seg_cnv_indices);

            if (seg_cnv_array[i]->seg_cnvs[j].excluded_regions)
                free(seg_cnv_array[i]->seg_cnvs[j].excluded_regions);

            if (seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints) {
                for (k=0; k<seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints_size; k++) {
                    if (seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints[k].paired_read_starts)
                        free(seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints[k].paired_read_starts);

                    if (seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints[k].paired_read_ends)
                        free(seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints[k].paired_read_ends);
                }

                free(seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints);
            }

            if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv) {
                if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv->cnv_breakpoints)
                    free(seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv->cnv_breakpoints);
                free(seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv);
            }
        }

        if (seg_cnv_array[i]->seg_cnvs)
            free(seg_cnv_array[i]->seg_cnvs);

        if (seg_cnv_array[i])
            free(seg_cnv_array[i]);
    }

    if (seg_cnv_array)
        free(seg_cnv_array);
}

void cnvSegmentation(Chromosome_Tracking *chrom_tracking, Segment_Array** segment_array, User_Input *user_inputs) {
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

    // don't free the following
    //free(tmp_basename);
}

void storeSegmentsLocallyAndInit(Segment_Array** segment_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking) {
    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // process one chromosome at a time
        //
        seg_cnv_array[i]->size = segment_array[i]->size;
        seg_cnv_array[i]->capacity = segment_array[i]->size;
        seg_cnv_array[i]->seg_cnvs = calloc(seg_cnv_array[i]->capacity, sizeof(Segmented_CNV));

        uint32_t j;
        for (j = 0; j < segment_array[i]->size; ++j) {
            // CNV segments on each chromosome
            //
            seg_cnv_array[i]->seg_cnvs[j].segment.start = segment_array[i]->segments[j].start;
            seg_cnv_array[i]->seg_cnvs[j].segment.end   = segment_array[i]->segments[j].end;
            seg_cnv_array[i]->seg_cnvs[j].segment.ave_coverage = segment_array[i]->segments[j].ave_coverage;
            seg_cnv_array[i]->seg_cnvs[j].segment.log2ratio_mean = segment_array[i]->segments[j].log2ratio_mean;

            /* The following will be initialized at the mergeCNVsFromSameSegment()
             *
             * seg_cnv_array[i]->seg_cnvs[j].inner_cnv_size = 0;
            seg_cnv_array[i]->seg_cnvs[j].inner_cnv_capacity = 100;
            seg_cnv_array[i]->seg_cnvs[j].inner_cnv = \
                        calloc(seg_cnv_array[i]->seg_cnvs[j].inner_cnv_capacity, sizeof(INNER_CNV));
            */

            seg_cnv_array[i]->seg_cnvs[j].seg_cnv_index_size = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_cnv_index_capacity = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_cnv_indices = NULL;

            seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints = NULL;
            seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints_size = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_breakpoints_capacity = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_left_start_index  = -1;   // 0 index is valid, only -1 works
            seg_cnv_array[i]->seg_cnvs[j].seg_right_end_index = -1;     // 0 index is valid, only -1 works

            //seg_cnv_array[i]->seg_cnvs[j].left_breakpoints  = NULL;
            //seg_cnv_array[i]->seg_cnvs[j].right_breakpoints = NULL;
            //seg_cnv_array[i]->seg_cnvs[j].num_of_left_paired_reads  = NULL;
            //seg_cnv_array[i]->seg_cnvs[j].num_of_right_paired_reads = NULL;

            seg_cnv_array[i]->seg_cnvs[j].seg_imp_PR_start = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_imp_PR_end   = 0;
            seg_cnv_array[i]->seg_cnvs[j].seg_num_of_imp_PR_TLEN_1000 = 0;

            seg_cnv_array[i]->seg_cnvs[j].excluded_regions = NULL;
        }
    }
}

void processSegmentationData(CNV_Array **cnv_array, Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, khash_t(m32) **anchor_breakpoints_hash_array, User_Input *user_inputs, Simple_Stats *equal_window_stats, Stats_Info *stats_info, Bed_Info *low_mappability_bed_info, Bed_Info *gc_lt25pct_bed_info, Bed_Info *gc_gt85pct_bed_info) {
#pragma omp parallel num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        uint32_t chr_index;
        for (chr_index=0; chr_index<chrom_tracking->number_of_chromosomes; ++chr_index) {
#pragma omp task
          {
            int thread_id = omp_get_thread_num();
            printf("Current thread id in process_segmentation_datacalls: %d with chr %s\n", thread_id, chrom_tracking->chromosome_ids[chr_index]);
            // find the corresponding index in cnv_array and seg_cnv_array
            //
            uint32_t cnv_array_idx;
            for (cnv_array_idx=0; cnv_array_idx<chrom_tracking->number_of_chromosomes; ++cnv_array_idx) {
                if (strcmp(cnv_array[cnv_array_idx]->chromosome_id, chrom_tracking->chromosome_ids[chr_index]) == 0)
                    break;
            }

            uint32_t seg_cnv_array_idx;
            for (seg_cnv_array_idx=0; seg_cnv_array_idx<chrom_tracking->number_of_chromosomes; ++seg_cnv_array_idx) {
                if (strcmp(seg_cnv_array[seg_cnv_array_idx]->chromosome_id, chrom_tracking->chromosome_ids[chr_index]) == 0)
                    break;
            }

            /*uint32_t breakpoint_array_idx;
            for (breakpoint_array_idx=0; breakpoint_array_idx<chrom_tracking->number_of_chromosomes; ++breakpoint_array_idx) {
                if (strcmp(anchor_breakpoints_hash_array[breakpoint_array_idx]->chromosome_id, chrom_tracking->chromosome_ids[chr_index]) == 0)
                    break;
            }*/

            // remove Ns regions if they are on the edge
            //
            Binned_Data_Wrapper *excluded_regions = calloc(1, sizeof(Binned_Data_Wrapper));
            excluded_regions->chromosome_id = calloc(strlen(chrom_tracking->chromosome_ids[chr_index])+1, sizeof(char*));
            strcpy(excluded_regions->chromosome_id, chrom_tracking->chromosome_ids[chr_index]);
            excluded_regions->size = 0;
            excluded_regions->capacity = INIT_SIZE;
            excluded_regions->data   = calloc(INIT_SIZE, sizeof(Binned_Data));
            excluded_regions->starts = calloc(INIT_SIZE, sizeof(uint32_t));
            excluded_regions->ends   = calloc(INIT_SIZE, sizeof(uint32_t));
            int total_lines = processFile(chrom_tracking->chromosome_ids[chr_index], user_inputs->excluded_region_file, excluded_regions);
            removeExcludedRegionsFromSegments(seg_cnv_array[seg_cnv_array_idx], excluded_regions, total_lines);

            // use intersect to find out which cnv segments are true CNVs
            //
            fprintf(stderr, "Before intersect\n");
            findIntersectedCNVs(cnv_array[cnv_array_idx], seg_cnv_array[seg_cnv_array_idx], user_inputs);
            fprintf(stderr, "After intersect\n");
            checkBreakpointsForEachSegment(seg_cnv_array[seg_cnv_array_idx], anchor_breakpoints_hash_array[seg_cnv_array_idx]);
            fprintf(stderr, "After breakpoints check\n");
            setBreakpoinsAtSegmentEnds(seg_cnv_array[seg_cnv_array_idx]);
            fprintf(stderr, "After setBreakpoinsAtSegmentEnds \n");

            mergeCNVsFromSameSegment(seg_cnv_array[seg_cnv_array_idx], cnv_array[cnv_array_idx], user_inputs, equal_window_stats);
            fprintf(stderr, "After merging\n");

            // now checking the low mappability and GC% < 25%
            //
            if (low_mappability_bed_info)
                checkingFalsePositives(seg_cnv_array[seg_cnv_array_idx], chrom_tracking->chromosome_ids[chr_index], low_mappability_bed_info, 1);

            if (gc_lt25pct_bed_info)
                checkingFalsePositives(seg_cnv_array[seg_cnv_array_idx], chrom_tracking->chromosome_ids[chr_index], gc_lt25pct_bed_info, 2);

            if (gc_gt85pct_bed_info)
                checkingFalsePositives(seg_cnv_array[seg_cnv_array_idx], chrom_tracking->chromosome_ids[chr_index], gc_gt85pct_bed_info, 3);

            if (excluded_regions->data) free(excluded_regions->data);
            if (excluded_regions->starts) free(excluded_regions->starts);
            if (excluded_regions->ends) free(excluded_regions->ends);
            if (excluded_regions) free(excluded_regions);

          } // end task 
        } // end for loop
#pragma omp taskwait
      } // end single
    } // end parallel

    generateSegmentedCNVs(seg_cnv_array, chrom_tracking, equal_window_stats, stats_info, user_inputs);
}

void findIntersectedCNVs(CNV_Array *cnv_array, Segmented_CNV_Array *seg_cnv_array, User_Input *user_inputs) {
    // since only CNVs with length >= 1000bps will be considered, we need to find out how many here
    //
    uint32_t i=0;
    int32_t counter=0;      // counter might go negative
    for (i=0; i<cnv_array->size; i++) {
        if (cnv_array->cnvs[i].inner_cnv.end - cnv_array->cnvs[i].inner_cnv.start >= user_inputs->min_cnv_length)
            counter++;
    }

    // an array to store all positions
    //
    uint32_t capacity = (counter + seg_cnv_array->size) * 2;
    uint32_t *all_starts_ends = calloc(capacity, sizeof(uint32_t));
    counter = 0;

    khash_t(m32) *cnv_start_hash = kh_init(m32);
    khash_t(m32) *cnv_end_hash   = kh_init(m32);

    for (i=0; i<cnv_array->size; i++) {
        // skip those with CNV length < 1000, which is set at the generateVCF()
        //
        if (cnv_array->cnvs[i].inner_cnv.end - cnv_array->cnvs[i].inner_cnv.start < user_inputs->min_cnv_length)
            continue;

        all_starts_ends[counter] = cnv_array->cnvs[i].inner_cnv.start+1;
        setValueToKhashBucket32(cnv_start_hash, all_starts_ends[counter], i);   // store cnv index
        counter++;

        all_starts_ends[counter] = cnv_array->cnvs[i].inner_cnv.end;
        setValueToKhashBucket32(cnv_end_hash, all_starts_ends[counter], i);     // store cnv index
        counter++;
    }

    khash_t(m32) *segment_start_hash = kh_init(m32);
    khash_t(m32) *segment_end_hash   = kh_init(m32);

    for (i=0; i<seg_cnv_array->size; i++) {
        // segmentation follows bed format, so for start pos, we need to add 1
        //
        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.start + 1;
        setValueToKhashBucket32(segment_start_hash, all_starts_ends[counter], i);   // store seg_cnvs index
        counter++;

        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.end;
        setValueToKhashBucket32(segment_end_hash, all_starts_ends[counter], i);     // store seg_cnvs index
        counter++;
    }

    // sorting
    //
    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    // do intersect
    //
    khash_t(m32) *seen_cnv_starts_hash = kh_init(m32);
    khash_t(m32) *seen_seg_starts_hash = kh_init(m32);
    khash_t(m32) *seen_cnv_index_hash  = kh_init(m32);

    counter = 0;
    int32_t cnv_index = -1;             // need to use signed value as negative value might return
    int32_t seg_index = -1;             // need to use signed value as negative value might return

    for (i=0; i<capacity; i++) {
        //if (all_starts_ends[i] == 41627000)
        //    printf("Stop intersect 29060000\n");

        // if a key is present in both start and end position, we have to process start first
        //
        bool start_first = false;
        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(segment_start_hash, all_starts_ends[i]) && checkm32KhashKey(segment_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i]) && checkm32KhashKey(segment_end_hash, all_starts_ends[i]) 
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(segment_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (!start_first && (checkm32KhashKey(cnv_end_hash, all_starts_ends[i]) ||
                                    checkm32KhashKey(segment_end_hash, all_starts_ends[i]))) {

            // always decrease count if it is the end position
            //
            counter--;

            /*if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "CNV end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            } else if (checkm32KhashKey(segment_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "SEG end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }*/

            if (counter < 0) {
                if (i == 0) {
                    counter = 0;
                } else {
                    fprintf(stderr, "Error: chr id is %s and index i is %"PRIu32"\n", cnv_array->chromosome_id, i);
                    fprintf(stderr, "Error: the counter %"PRId16" should NOT be negative in breakpoint/cnv intersection\n", counter);
                    fprintf(stderr, "current pos: %"PRIu32" and the prev pos %"PRIu32"\n", all_starts_ends[i], all_starts_ends[i-1]);
                    exit(EXIT_FAILURE);
                }
            }

            // when current position is an end and counter >= 1, there is an intersect
            // There might be multiple CNV associated with this segment and I will store all of them
            //
            if (counter >= 1) {
                if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {   // it is the cnv end
                    cnv_index = getValueFromKhash32(cnv_end_hash, all_starts_ends[i]);

                    if (cnv_index >= 0 && checkm32KhashKey(seen_cnv_index_hash, cnv_index))
                        continue;

                } else {    // it is the segment end
                    seg_index = getValueFromKhash32(segment_end_hash, all_starts_ends[i]);
                }

                if (cnv_index < 0 || seg_index < 0) {
                    fprintf(stderr, "Error: cnv_index %"PRIu32" and seg_index %"PRIu32" should always be >= 0\n", cnv_index, seg_index);
                    exit(EXIT_FAILURE);
                }

                // it is possible that one CNV spans 2 segments. It is because we need to set the 
                // breakpoints on both end to extend the CNV. 
                // Need to check if the overlaping part is > half of the CNV length
                //
                uint32_t cnv_start = cnv_array->cnvs[cnv_index].inner_cnv.start;
                uint32_t cnv_end   = cnv_array->cnvs[cnv_index].inner_cnv.end;
                uint32_t seg_start = seg_cnv_array->seg_cnvs[seg_index].segment.start;
                uint32_t seg_end   = seg_cnv_array->seg_cnvs[seg_index].segment.end;

                if (cnv_start < seg_start) {
                    if (!checkm32KhashKey(seen_cnv_index_hash, cnv_index) && (cnv_end - seg_start > (cnv_end - cnv_start)/2)) {
                        addCNVindexToSegment(cnv_index, seg_cnv_array, seg_index);
                        setValueToKhashBucket32(seen_cnv_index_hash, cnv_index, 1);
                    }
                } else {
                    if (!checkm32KhashKey(seen_cnv_index_hash, cnv_index) && (seg_end - cnv_start > (cnv_end - cnv_start)/2)) {
                        addCNVindexToSegment(cnv_index, seg_cnv_array, seg_index);
                        setValueToKhashBucket32(seen_cnv_index_hash, cnv_index, 1);
                    }
                }
            }

            // delete end hash once it is processed
            // it is because some start/end has the same value, 
            // sometimes, start/start or end/end has the same value.
            //
            khiter_t iter;
            if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, cnv_end_hash, all_starts_ends[i]);
                if (iter != kh_end(cnv_end_hash)) {
                    kh_del(m32, cnv_end_hash, iter);
                }
            } else if (checkm32KhashKey(segment_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, segment_end_hash, all_starts_ends[i]);
                if (iter != kh_end(segment_end_hash)) {
                    kh_del(m32, segment_end_hash, iter);        // this deletes the key
                }
            }

        } else {    // it's a start position
            counter++;

            // get current CNV start position in cnv_index
            //
            if (checkm32KhashKey(cnv_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_cnv_starts_hash, all_starts_ends[i], 1);

                cnv_index = getSignedValueFromKhash32(cnv_start_hash, all_starts_ends[i]);
                if (cnv_index == -1) {
                    fprintf(stderr, "Something went wrong with cnv index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                //fprintf(stderr, "CNV start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }

            // get current segment start position index
            //
            if (checkm32KhashKey(segment_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_seg_starts_hash, all_starts_ends[i], 1);

                seg_index = getSignedValueFromKhash32(segment_start_hash, all_starts_ends[i]);
                if (seg_index == -1) {
                    fprintf(stderr, "Something went wrong with seg index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                //fprintf(stderr, "SEG start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }
        }
    }

    // clean-up
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, cnv_start_hash);
    kh_destroy(m32, cnv_end_hash);
    kh_destroy(m32, segment_start_hash);
    kh_destroy(m32, segment_end_hash);
    kh_destroy(m32, seen_cnv_starts_hash);
    kh_destroy(m32, seen_seg_starts_hash);
    kh_destroy(m32, seen_cnv_index_hash);
}

void addCNVindexToSegment(uint32_t cnv_index, Segmented_CNV_Array *seg_cnv_array, uint32_t seg_cnv_index) {
    if (seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices==NULL) {
        seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_size = 0;
        seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity = 500;
        seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices = \
                    calloc(seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity, sizeof(uint32_t));
    }

    seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices[seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_size] = cnv_index;
    seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_size++;

    if (seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_size + 5 >= seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity) {
        // need to dynamically increase the memory
        //
        seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity = seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity + 200;
        seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices = realloc(seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices, \
                seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_index_capacity*sizeof(uint32_t));
        failureExit(seg_cnv_array->seg_cnvs[seg_cnv_index].seg_cnv_indices, "seg_cnv_array->cnv_indices");
    }
}

void checkBreakpointsForEachSegment(Segmented_CNV_Array *seg_cnv_array, khash_t(m32) *anchor_breakpoints_hash) {
    // For breakpoints
    //
    int32_t capacity = PR_INIT_SIZE*500;       // the initial size is set to 100,000
    uint32_t *all_starts_ends = calloc(capacity, sizeof(uint32_t));

    khash_t(m32) *breakpoint_start_hash = kh_init(m32);
    khash_t(m32) *breakpoint_end_hash   = kh_init(m32);
    khash_t(m32) *breakpoint_start_end_lookup = kh_init(m32);
    khash_t(m32) *breakpoint_end_start_lookup = kh_init(m32);

    int32_t counter=0;

    khint_t k;
    for (k=kh_begin(anchor_breakpoints_hash); k!=kh_end(anchor_breakpoints_hash); ++k) {
        if (kh_exist(anchor_breakpoints_hash, k)) {
            // Only process breakpoint with count >= 2 and number of TLEN >=2
            //
            if (kh_value(anchor_breakpoints_hash, k) >= 2) {
                uint32_t breakpoint_pos = kh_key(anchor_breakpoints_hash, k);

                if ((int32_t)breakpoint_pos - DISTANCE_CUTOFF < 0) {
                    // prevent negative value or value overflow
                    //
                    all_starts_ends[counter] = 1;
                } else {
                    all_starts_ends[counter] = breakpoint_pos - DISTANCE_CUTOFF;
                }
                setValueToKhashBucket32(breakpoint_start_hash, all_starts_ends[counter], breakpoint_pos);
                counter++;

                if (breakpoint_pos + DISTANCE_CUTOFF > seg_cnv_array->chrom_length) {
                    all_starts_ends[counter] = seg_cnv_array->chrom_length - 1;
                } else {
                    all_starts_ends[counter] = breakpoint_pos + DISTANCE_CUTOFF;
                }
                setValueToKhashBucket32(breakpoint_end_hash, all_starts_ends[counter], breakpoint_pos);
                counter++;

                setValueToKhashBucket32(breakpoint_start_end_lookup, all_starts_ends[counter-2], all_starts_ends[counter-1]);
                setValueToKhashBucket32(breakpoint_end_start_lookup, all_starts_ends[counter-1], all_starts_ends[counter-2]);

                // dynamically increase the capacity of anchor breakpoint array
                //
                if (counter + 5 >= capacity) {
                    capacity += PR_INIT_SIZE*25;
                    all_starts_ends = realloc(all_starts_ends, capacity * sizeof(uint32_t));
                    failureExit(all_starts_ends, "all_starts_ends in intersecting CNV w/ breakpoint memory realloc failed\n");
                }
            }
        }
    }

    // reset the capacity to the size of counter plus the seg_cnv_array->size
    //
    capacity = counter + seg_cnv_array->size*2;
    all_starts_ends = realloc(all_starts_ends, capacity * sizeof(uint32_t));
    failureExit(all_starts_ends, "all_starts_ends in intersecting CNV/breakpoint memory realloc failed\n");

    // walk through the segment CNV array and save all CNV coordinates
    //
    khash_t(m32) *seg_start_hash = kh_init(m32);
    khash_t(m32) *seg_end_hash   = kh_init(m32);
    khash_t(m32) *seg_start_end_lookup = kh_init(m32);
    khash_t(m32) *seg_end_start_lookup = kh_init(m32);

    uint32_t i;
    for (i=0; i<seg_cnv_array->size;i++) {
        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.start+1;      // bed format 0-based
        setValueToKhashBucket32(seg_start_hash, all_starts_ends[counter], i);       // store index here
        counter++;

        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.end;
        setValueToKhashBucket32(seg_end_hash, all_starts_ends[counter], i);
        counter++;

        setValueToKhashBucket32(seg_start_end_lookup, all_starts_ends[counter-2], all_starts_ends[counter-1]);
        setValueToKhashBucket32(seg_end_start_lookup, all_starts_ends[counter-1], all_starts_ends[counter-2]);
    }

    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    // do intersect
    //
    khash_t(m32) *live_seg_start_hash = kh_init(m32);
    khash_t(m32) *live_bpt_start_hash = kh_init(m32);

    khash_t(m32) *seen_seg_starts_hash = kh_init(m32);
    khash_t(m32) *seen_bpt_starts_hash = kh_init(m32);

    khiter_t iter;
    counter = 0;
    int32_t seg_index = -1;             // need to use signed value as sometimes, no value found
    for (i=0; i<(uint32_t)capacity; i++) {
        //if (all_starts_ends[i] == 29093235 || all_starts_ends[i] == 29095235)
        //    printf("here it is\n");

        // if a key is present in both start and end position, we have to process start first
        //
        bool start_first = false;
        if (checkm32KhashKey(seg_start_hash, all_starts_ends[i]) && checkm32KhashKey(seg_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]))
            start_first = true;

        if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i]) && checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_bpt_starts_hash, all_starts_ends[i]))
            start_first = true;

        if (checkm32KhashKey(seg_start_hash, all_starts_ends[i]) && checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]))
            start_first = true;

        if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i]) && checkm32KhashKey(seg_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_bpt_starts_hash, all_starts_ends[i]))
            start_first = true;

        if (!start_first && (checkm32KhashKey(seg_end_hash, all_starts_ends[i]) ||
                    checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i]))) {

            // always decrease count if it is the end position
            //
            counter--;

            if (counter < 0) {
                if (i == 0) {
                    counter = 0;    // this happens when breakpoint_pos - DISTANCE_CUTOFF
                } else {
                    fprintf(stderr, "Error: chr id is %s and index i is %"PRIu32"\n", seg_cnv_array->chromosome_id, i);
                    fprintf(stderr, "Error: the counter %"PRId16" should NOT be negative in breakpoint/cnv intersection\n", counter);
                    fprintf(stderr, "current pos: %"PRIu32" and the prev pos %"PRIu32"\n", all_starts_ends[i], all_starts_ends[i-1]);
                    exit(EXIT_FAILURE);
                }
            }

            
            /*if (checkm32KhashKey(seg_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "%s\tSEG_End\t%"PRIu32"\t%"PRId32"\n", seg_cnv_array->chromosome_id, all_starts_ends[i], counter);
            } else if (checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "%s\tBreakpoint_End\t%"PRIu32"\t%"PRId32"\n", seg_cnv_array->chromosome_id, all_starts_ends[i], counter);
            }*/

            // when current position is an end and counter >= 1, there is an intersect
            // There might be multiple breakpoints associated with this CNV
            // and the breakpoint intervals are not Merged. Will store all of them
            //
            if (counter >= 1) {
                uint32_t seg_start = 0;
                uint32_t seg_end = 0;

                // because the max length of breakpoint interval used for intersection is only 600
                // so we don't have to worry that breakpoint interval completely engulf a CNV
                //
                uint32_t cur_anchor_breakpoint = 0;
                if (checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i])) {
                    cur_anchor_breakpoint = getValueFromKhash32(breakpoint_end_hash, all_starts_ends[i]);

                    // loop through current live_seg_start_hash
                    //
                    for (iter=kh_begin(live_seg_start_hash); iter!=kh_end(live_seg_start_hash); ++iter) {
                        if (kh_exist(live_seg_start_hash, iter)) {
                            seg_index = kh_value(live_seg_start_hash, iter);

                            // check distance between segment boundary and anchor breakpoint
                            // if it is more than DISTANCE_CUTOFF (1000 bps) away, don't store it
                            //
                            seg_start = seg_cnv_array->seg_cnvs[seg_index].segment.start;
                            seg_end   = seg_cnv_array->seg_cnvs[seg_index].segment.end;
                            if ( (abs((signed)(seg_start - cur_anchor_breakpoint)) > DISTANCE_CUTOFF) 
                                    && (abs((signed)(seg_end - cur_anchor_breakpoint)) > DISTANCE_CUTOFF) )
                                continue;

                            addBreakpoints(seg_cnv_array, seg_index, anchor_breakpoints_hash, cur_anchor_breakpoint);
                        }
                    }

                    // delete the breakpoint start in the live_bpt_start_hash
                    //
                    iter = kh_get(m32, live_bpt_start_hash, cur_anchor_breakpoint - DISTANCE_CUTOFF);
                    kh_del(m32, live_bpt_start_hash, iter);
                } else {
                    seg_start = getValueFromKhash32(seg_end_start_lookup, all_starts_ends[i]);
                    seg_index = getValueFromKhash32(live_seg_start_hash, seg_start);
                    seg_end   = seg_cnv_array->seg_cnvs[seg_index].segment.end;

                    // loop through the live_bpt_start_hash
                    //
                    for (iter=kh_begin(live_bpt_start_hash); iter!=kh_end(live_bpt_start_hash); ++iter) {
                        if (kh_exist(live_bpt_start_hash, iter)) {
                            uint32_t breakpoint_start = kh_key(live_bpt_start_hash, iter);
                            cur_anchor_breakpoint = getValueFromKhash32(breakpoint_start_hash, breakpoint_start);

                            // skip if the distances between the anchor breakpoint and segment boundaries > DISTANCE_CUTOFF
                            //
                            if (abs((signed)(seg_start - cur_anchor_breakpoint)) > DISTANCE_CUTOFF 
                                    && abs((signed)(seg_end - cur_anchor_breakpoint)) > DISTANCE_CUTOFF)
                                continue;

                            addBreakpoints(seg_cnv_array, seg_index, anchor_breakpoints_hash, cur_anchor_breakpoint);
                        }
                    }

                    // now delete current cnv start at the live_seg_start_hash
                    //
                    iter = kh_get(m32, live_seg_start_hash, seg_start);
                    kh_del(m32, live_seg_start_hash, iter);
                }
            } else {
                // clean-up both live_seg_start_hash and live_bpt_start_hash
                //
                for (iter=kh_begin(live_seg_start_hash); iter!=kh_end(live_seg_start_hash); ++iter) {
                    if (kh_exist(live_seg_start_hash, iter))
                        kh_del(m32, live_seg_start_hash, iter);
                }

                for (iter=kh_begin(live_bpt_start_hash); iter!=kh_end(live_bpt_start_hash); ++iter) {
                    if (kh_exist(live_bpt_start_hash, iter))
                        kh_del(m32, live_bpt_start_hash, iter);
                }
            }

            // as some start and end positions might be the same, it is quite confusion
            // so we need to delete the end ones one at a time once we have processed it
            //
            if (checkm32KhashKey(seg_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, seg_end_hash, all_starts_ends[i]);
                if (iter != kh_end(seg_end_hash)) {
                    kh_del(m32, seg_end_hash, iter);
                }
            } else if (checkm32KhashKey(breakpoint_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, breakpoint_end_hash, all_starts_ends[i]);
                if (iter != kh_end(breakpoint_end_hash)) {
                    kh_del(m32, breakpoint_end_hash, iter);        // this deletes the key
                }
            }
        } else {
            // start position
            //
            counter++;

            // get current segment CNV start position in seg_index
            //
            if (checkm32KhashKey(seg_start_hash, all_starts_ends[i])) {
                seg_index = getSignedValueFromKhash32(seg_start_hash, all_starts_ends[i]);
                if (seg_index == -1) {
                    fprintf(stderr, "Something went wrong with SEG/Breakpoint start at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }
                setValueToKhashBucket32(live_seg_start_hash, all_starts_ends[i], seg_index);
                setValueToKhashBucket32(seen_seg_starts_hash, all_starts_ends[i], seg_index);

                //fprintf(stderr, "%s\tSEG_Start\t%"PRIu32"\t%"PRId32"\n", seg_cnv_array->chromosome_id, all_starts_ends[i], counter);
            }

            // get current breakpoint start position
            //
            if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(live_bpt_start_hash, all_starts_ends[i], i);
                setValueToKhashBucket32(seen_bpt_starts_hash, all_starts_ends[i], i);
                
                //fprintf(stderr, "%s\tBreakpoint_Start\t%"PRIu32"\t%"PRId32"\n", seg_cnv_array->chromosome_id, all_starts_ends[i], counter);
            }

            // keep this here for future document.
            // Note: we don't want to delete anything from the hash table,
            // Instead, we will use seen* hashtable to guard the reuse of the hash keys
            /*
            if (checkm32KhashKey(seg_start_hash, all_starts_ends[i])) {
                iter = kh_get(m32, seg_start_hash, all_starts_ends[i]);
                if (iter != kh_end(seg_start_hash)) {
                    kh_del(m32, seg_start_hash, iter);
                }   
            } else if (checkm32KhashKey(breakpoint_start_hash, all_starts_ends[i])) {
                iter = kh_get(m32, breakpoint_start_hash, all_starts_ends[i]);
                if (iter != kh_end(breakpoint_start_hash)) {
                    kh_del(m32, breakpoint_start_hash, iter);        // this deletes the key
                }
            }*/
        }
    }

    // clean-up
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, breakpoint_end_start_lookup);
    kh_destroy(m32, breakpoint_start_end_lookup);
    kh_destroy(m32, breakpoint_start_hash);
    kh_destroy(m32, breakpoint_end_hash);
    kh_destroy(m32, seg_start_end_lookup);
    kh_destroy(m32, seg_end_start_lookup);
    kh_destroy(m32, seg_start_hash);
    kh_destroy(m32, seg_end_hash);
    kh_destroy(m32, live_seg_start_hash);
    kh_destroy(m32, live_bpt_start_hash);
    kh_destroy(m32, seen_seg_starts_hash);
    kh_destroy(m32, seen_bpt_starts_hash);
}

void addBreakpoints(Segmented_CNV_Array *seg_cnv_array, uint32_t seg_index, khash_t(m32) *anchor_breakpoints_hash, uint32_t anchor_breakpoint) {
    khint_t k = kh_get(m32, anchor_breakpoints_hash, anchor_breakpoint);
    if (k != kh_end(anchor_breakpoints_hash)) {
        // add breakpoint info
        //
        if (seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints == NULL) {
            seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_capacity = 300;
            seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints =
                calloc(seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_capacity, sizeof(CNV_Breakpints));
            seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_size = 0;
            seg_cnv_array->seg_cnvs[seg_index].seg_left_start_index = -1;   // 0-index is valid, so need to use -1 for it
            seg_cnv_array->seg_cnvs[seg_index].seg_right_end_index  = -1;   // 0-index is valid, so need to use -1 for it
        }

        uint16_t index = seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_size;

        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].breakpoint = anchor_breakpoint;
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].num_of_breakpoints = kh_value(anchor_breakpoints_hash, k);
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].num_of_TLEN_ge_1000=0;
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].num_of_paired_reads = 0;
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].capacity = 10;
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].paired_read_starts = \
                            calloc(seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].capacity, sizeof(uint32_t));
        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].paired_read_ends = \
                            calloc(seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints[index].capacity, sizeof(uint32_t));

        //storePairedReadsAcrossABreakpoint(seg_cnv_array, seg_index, anchor_breakpoint, header, sfh_idx, sfh, user_inputs);

        seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_size++;
        
        // resize the cnv_breakpoints array if needed
        //
        if ((seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_size + 3) > seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_capacity) {
            seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_capacity += 200;
            seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints = realloc(seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints, \
                    (seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints_capacity)*sizeof(CNV_Breakpints));
            failureExit(seg_cnv_array->seg_cnvs[seg_index].seg_breakpoints, "cnv_array->cnvs[cnv_index].cnv_breakpoints");
        }
    }
}

void setBreakpoinsAtSegmentEnds(Segmented_CNV_Array *seg_cnv_array) {
    uint32_t i;
    for (i=0; i<seg_cnv_array->size; ++i) {
        uint32_t seg_start = seg_cnv_array->seg_cnvs[i].segment.start;
        uint32_t seg_end   = seg_cnv_array->seg_cnvs[i].segment.end;

        if (seg_cnv_array->seg_cnvs[i].seg_breakpoints == NULL)
            continue;

        bool left_most_set = false, right_most_set = false;
        uint32_t j;
        for (j=0; j<seg_cnv_array->seg_cnvs[i].seg_breakpoints_size; ++j) {
            // number of count for each breakpoint should be >= 2 as checked earlier at addBreakpoint()
            //
            bool left=false;
            uint32_t cur_breakpoint = seg_cnv_array->seg_cnvs[i].seg_breakpoints[j].breakpoint;

            if (abs((signed) (cur_breakpoint - seg_start)) < abs((signed) (cur_breakpoint - seg_end))) {
                if (abs((signed) (cur_breakpoint - seg_start) > DISTANCE_CUTOFF*2)) {
                    continue;
                } else {
                    left = true;
                }
            } else {
                if ((signed) (cur_breakpoint - seg_end) > DISTANCE_CUTOFF*2) {
                    // no need to go any further
                    //
                    break;
                }

                if ((signed) seg_end - cur_breakpoint > DISTANCE_CUTOFF*2)
                    continue;
            }

            uint32_t prev_num_bpts=0, cur_num_bpts=0;
            if (cur_breakpoint < seg_start) {
                left = true;
            } else if (cur_breakpoint > seg_end) {
                left = false;
            } else {
                if (abs((signed) (cur_breakpoint - seg_start)) < abs((signed) (cur_breakpoint - seg_end)))
                    left = true;
            }

            if (left && cur_breakpoint > seg_end)
                left = false;

            if (left) {
                // the current anchor breakpoint closer to the left-hand side
                //
                if (!left_most_set) {
                    seg_cnv_array->seg_cnvs[i].seg_left_start_index = j;
                    left_most_set = true;
                } else {
                    // left most already set, let's do the comparison
                    //
                    prev_num_bpts = seg_cnv_array->seg_cnvs[i].seg_breakpoints[seg_cnv_array->seg_cnvs[i].seg_left_start_index].num_of_breakpoints;
                    //prev_tlen = cnv_array->cnvs[i].seg_breakpoints[cnv_array->cnvs[i].seg_left_start_index].num_of_TLEN_ge_1000;
                    cur_num_bpts = seg_cnv_array->seg_cnvs[i].seg_breakpoints[j].num_of_breakpoints;
                    //cur_tlen = cnv_array->cnvs[i].seg_breakpoints[j].num_of_TLEN_ge_1000;

                    if (cur_num_bpts > prev_num_bpts)
                        seg_cnv_array->seg_cnvs[i].seg_left_start_index = j;
                }
            } else {
                // the current anchor breakpoint closer to the right-hand side
                //
                if (!right_most_set) {
                    seg_cnv_array->seg_cnvs[i].seg_right_end_index = j;
                    right_most_set = true;
                } else {
                    // right hand side already set, let's do the comparison
                    //
                    prev_num_bpts = seg_cnv_array->seg_cnvs[i].seg_breakpoints[seg_cnv_array->seg_cnvs[i].seg_right_end_index].num_of_breakpoints;
                    //prev_tlen = seg_cnv_array->seg_cnvs[i].seg_breakpoints[seg_cnv_array->seg_cnvs[i].seg_right_end_index].num_of_TLEN_ge_1000;
                    cur_num_bpts  = seg_cnv_array->seg_cnvs[i].seg_breakpoints[j].num_of_breakpoints;
                    //cur_tlen = seg_cnv_array->seg_cnvs[i].seg_breakpoints[j].num_of_TLEN_ge_1000;

                    if (cur_num_bpts > prev_num_bpts)
                        seg_cnv_array->seg_cnvs[i].seg_right_end_index = j;
                }
            }
        }
    }
}

void removeExcludedRegionsFromSegments(Segmented_CNV_Array *seg_cnv_array, Binned_Data_Wrapper *excluded_regions, int total_size) {
    uint32_t capacity = (seg_cnv_array->size + total_size) * 2;
    uint32_t *all_starts_ends = calloc(capacity, sizeof(uint32_t));

    khash_t(m32) *excluded_start_hash = kh_init(m32);
    khash_t(m32) *excluded_end_hash   = kh_init(m32);

    uint32_t i=0, counter=0;
    for (i=0; i<total_size; ++i) {
        all_starts_ends[counter] = excluded_regions->data[i].start;
        setValueToKhashBucket32(excluded_start_hash, all_starts_ends[counter], i);     // store index
        counter++;

        all_starts_ends[counter] = excluded_regions->data[i].end;
        setValueToKhashBucket32(excluded_end_hash, all_starts_ends[counter], i);
        counter++;
    }

    khash_t(m32) *segment_start_hash = kh_init(m32);
    khash_t(m32) *segment_end_hash   = kh_init(m32);

    for (i=0; i<seg_cnv_array->size; i++) {
        // segmentation follows bed format, so for start pos, we need to add 1
        //
        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.start + 1;
        setValueToKhashBucket32(segment_start_hash, all_starts_ends[counter], i);   // store seg_cnvs index
        counter++;

        all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].segment.end;
        setValueToKhashBucket32(segment_end_hash, all_starts_ends[counter], i);     // store seg_cnvs index
        counter++;
    }

    // sorting
    //
    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    // do intersect
    //
    khash_t(m32) *seen_exc_starts_hash = kh_init(m32);
    khash_t(m32) *seen_seg_starts_hash = kh_init(m32);

    counter = 0;
    int32_t exc_index = -1;             // need to use signed value as negative value might return
    int32_t seg_index = -1;             // need to use signed value as negative value might return

    for (i=0; i<capacity; i++) {
        //if (all_starts_ends[i] == 99667000)
        //    printf("Stopped at remove excluded regions\n");

        // if a key is present in both start and end position, we have to process start first
        //
        bool start_first = false;
        if (checkm32KhashKey(excluded_start_hash, all_starts_ends[i]) && checkm32KhashKey(excluded_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_exc_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(segment_start_hash, all_starts_ends[i]) && checkm32KhashKey(segment_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(excluded_start_hash, all_starts_ends[i]) && checkm32KhashKey(segment_end_hash, all_starts_ends[i]) 
                && !checkm32KhashKey(seen_exc_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(segment_start_hash, all_starts_ends[i]) && checkm32KhashKey(excluded_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_seg_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (!start_first && (checkm32KhashKey(excluded_end_hash, all_starts_ends[i]) ||
                    checkm32KhashKey(segment_end_hash, all_starts_ends[i]))) {

            // always decrease count if it is the end position
            //
            counter--;

            /*if (checkm32KhashKey(cnv_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "CNV end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            } else if (checkm32KhashKey(segment_end_hash, all_starts_ends[i])) {
                fprintf(stderr, "SEG end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }*/

            // when current position is an end and counter >= 1, there is an intersect
            // There might be multiple CNV associated with this segment and I will store all of them
            //
            if (counter >= 1) {
                if (checkm32KhashKey(excluded_end_hash, all_starts_ends[i])) {   // it is the excluded end
                    exc_index = getValueFromKhash32(excluded_end_hash, all_starts_ends[i]);
                } else {    // it is the segment end
                    seg_index = getValueFromKhash32(segment_end_hash, all_starts_ends[i]);
                }

                if (exc_index < 0 || seg_index < 0) {
                    fprintf(stderr, "Error: cnv_index %"PRIu32" and seg_index %"PRIu32" should always be >= 0\n", exc_index, seg_index);
                    exit(EXIT_FAILURE);
                }

                // need to remove excluded region from the intersected segment
                //
                uint32_t excluded_start = excluded_regions->data[exc_index].start;
                uint32_t excluded_end   = excluded_regions->data[exc_index].end;
                uint32_t segment_start  = seg_cnv_array->seg_cnvs[seg_index].segment.start;
                uint32_t segment_end    = seg_cnv_array->seg_cnvs[seg_index].segment.end;

                // just skip it if the excluded regions is less than 500
                //
                if (excluded_end - excluded_start < 500)
                    continue;

                char regions_to_append[75];
                sprintf(regions_to_append, "%"PRIu32"-%"PRIu32, excluded_start, excluded_end); 

                if (excluded_end - excluded_start >= 500 && excluded_end - excluded_start < 1000) {
                    // the min CNV length is 1000, so the gap between 2 CNVs can be joined if they are separated by < 1000bps
                    // then save them for the output purpose
                    //
                    saveExcludedRegion(seg_cnv_array, seg_index, regions_to_append);
                    continue;
                }

                bool used = false;
                if (segment_start < excluded_end && excluded_end < segment_end) {
                    if (abs((signed) (excluded_end - segment_start)) < 1000) {
                        // 1000 is min CNV length, might be changed to 10% fraction
                        // need to remove left hand side of the segment
                        //
                        seg_cnv_array->seg_cnvs[seg_index].segment.start = excluded_end + 1;
                        used = true;
                    }
                }

                if (segment_start < excluded_start && excluded_start < segment_end) {
                    if (abs((signed) (segment_end - excluded_start)) < 1000) {
                        // need to remove right hand side of the segment
                        //
                        seg_cnv_array->seg_cnvs[seg_index].segment.end = excluded_start + 1;
                        used = true;
                    }
                }

                if (!used)
                    saveExcludedRegion(seg_cnv_array, seg_index, regions_to_append);
            }

            // delete processed end position from the hash table
            //
            khiter_t iter;
            if (checkm32KhashKey(excluded_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, excluded_end_hash, all_starts_ends[i]);
                if (iter != kh_end(excluded_end_hash)) {
                    kh_del(m32, excluded_end_hash, iter);
                }
            } else if (checkm32KhashKey(segment_end_hash, all_starts_ends[i])) {
                iter = kh_get(m32, segment_end_hash, all_starts_ends[i]);
                if (iter != kh_end(segment_end_hash)) {
                    kh_del(m32, segment_end_hash, iter);        // this deletes the key
                }
            }
        } else {    // it's a start position
            counter++;

            // get current excluded start position in exc_index
            //
            if (checkm32KhashKey(excluded_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_exc_starts_hash, all_starts_ends[i], 1);

                exc_index = getSignedValueFromKhash32(excluded_start_hash, all_starts_ends[i]);
                if (exc_index == -1) {
                    fprintf(stderr, "Something went wrong with cnv index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                //fprintf(stderr, "Excluded start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }

            // get current segment start position index
            //
            if (checkm32KhashKey(segment_start_hash, all_starts_ends[i])) {
                setValueToKhashBucket32(seen_seg_starts_hash, all_starts_ends[i], 1);

                seg_index = getSignedValueFromKhash32(segment_start_hash, all_starts_ends[i]);
                if (seg_index == -1) {
                    fprintf(stderr, "Something went wrong with seg index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                //fprintf(stderr, "SEG start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }
        }
    }
        
    // clean-up
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, excluded_start_hash);
    kh_destroy(m32, excluded_end_hash);
    kh_destroy(m32, segment_start_hash);
    kh_destroy(m32, segment_end_hash);
    kh_destroy(m32, seen_exc_starts_hash);
    kh_destroy(m32, seen_seg_starts_hash);
}

void saveExcludedRegion(Segmented_CNV_Array *seg_cnv_array, uint32_t seg_cnv_index, char *regions_to_append) {
    if (seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions == NULL) {
        seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions = calloc(strlen(regions_to_append)+1, sizeof(char));
        strcpy(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions, regions_to_append);
    } else {
        uint32_t tmp_size = strlen(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions)
                                                            + strlen(regions_to_append) + 30;
        seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions =
            realloc(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions, tmp_size*sizeof(char));
        failureExit(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions, "seg_cnv_array->seg_cnvs[i].excluded_regions");

        strcat(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions, ",");
        strcat(seg_cnv_array->seg_cnvs[seg_cnv_index].excluded_regions, regions_to_append);
    }
}

void segmentInnerCNVInit(Segmented_CNV *seg_cnv, uint32_t capacity) {
    seg_cnv->seg_inner_cnv_size = 0;
    seg_cnv->seg_inner_cnv_capacity = capacity;
    seg_cnv->seg_inner_cnv = calloc(seg_cnv->seg_inner_cnv_capacity, sizeof(INNER_CNV));
}

void segmentInnerCNVInit2(Segmented_CNV *seg_cnv, uint32_t index) {
    seg_cnv->seg_inner_cnv[index].cnv_breakpoints = NULL;
    seg_cnv->seg_inner_cnv[index].start = 0;
    seg_cnv->seg_inner_cnv[index].end = 0;
    seg_cnv->seg_inner_cnv[index].passed = false;
    seg_cnv->seg_inner_cnv[index].qual = 0.0;
    seg_cnv->seg_inner_cnv[index].ave_coverage = 0.0;
    seg_cnv->seg_inner_cnv[index].cnv_type = ' ';
    seg_cnv->seg_inner_cnv[index].valid_cnv = false;
    seg_cnv->seg_inner_cnv[index].left_breakpoint = 0;
    seg_cnv->seg_inner_cnv[index].right_breakpoint = 0;
    seg_cnv->seg_inner_cnv[index].left_breakpoint_count = 0;
    seg_cnv->seg_inner_cnv[index].right_breakpoint_count = 0;
    seg_cnv->seg_inner_cnv[index].last_num_larger_TLEN_right = 0;
    seg_cnv->seg_inner_cnv[index].last_right_breakpoint_count = 0;
    seg_cnv->seg_inner_cnv[index].num_larger_TLEN_left = 0;
    seg_cnv->seg_inner_cnv[index].num_larger_TLEN_right = 0;
    seg_cnv->seg_inner_cnv[index].imp_PR_start = 0;
    seg_cnv->seg_inner_cnv[index].imp_PR_end = 0;
    seg_cnv->seg_inner_cnv[index].num_larger_imp_PR_TLEN = 0;
    seg_cnv->seg_inner_cnv[index].evidence_count = 0;
    seg_cnv->seg_inner_cnv[index].num_merged_CNVs = 0;
    seg_cnv->seg_inner_cnv[index].low_mapp_length = 0;
    seg_cnv->seg_inner_cnv[index].gc_lt25pct_length = 0;
    
}

void copyCNV(Segmented_CNV *seg_cnv, uint32_t seg_index, CNV *cnvs, uint32_t cnv_index) {
    //seg_cnv->seg_inner_cnv[seg_index].cnv_breakpoints = NULL;
    seg_cnv->seg_inner_cnv[seg_index].start = cnvs[cnv_index].inner_cnv.start;
    seg_cnv->seg_inner_cnv[seg_index].end = cnvs[cnv_index].inner_cnv.end;
    seg_cnv->seg_inner_cnv[seg_index].passed = cnvs[cnv_index].inner_cnv.passed;
    seg_cnv->seg_inner_cnv[seg_index].qual = cnvs[cnv_index].inner_cnv.qual;;
    seg_cnv->seg_inner_cnv[seg_index].ave_coverage = cnvs[cnv_index].inner_cnv.ave_coverage;
    seg_cnv->seg_inner_cnv[seg_index].cnv_type = cnvs[cnv_index].inner_cnv.cnv_type;
    seg_cnv->seg_inner_cnv[seg_index].left_breakpoint = cnvs[cnv_index].inner_cnv.left_breakpoint;
    seg_cnv->seg_inner_cnv[seg_index].right_breakpoint = cnvs[cnv_index].inner_cnv.right_breakpoint;
    seg_cnv->seg_inner_cnv[seg_index].left_breakpoint_count = cnvs[cnv_index].inner_cnv.left_breakpoint_count;
    seg_cnv->seg_inner_cnv[seg_index].right_breakpoint_count = cnvs[cnv_index].inner_cnv.right_breakpoint_count;
    seg_cnv->seg_inner_cnv[seg_index].last_right_breakpoint_count = cnvs[cnv_index].inner_cnv.last_right_breakpoint_count;
    seg_cnv->seg_inner_cnv[seg_index].num_larger_TLEN_left = cnvs[cnv_index].inner_cnv.num_larger_TLEN_left;
    seg_cnv->seg_inner_cnv[seg_index].num_larger_TLEN_right = cnvs[cnv_index].inner_cnv.num_larger_TLEN_right;
    seg_cnv->seg_inner_cnv[seg_index].last_num_larger_TLEN_right = cnvs[cnv_index].inner_cnv.last_num_larger_TLEN_right;
    seg_cnv->seg_inner_cnv[seg_index].imp_PR_start = cnvs[cnv_index].inner_cnv.imp_PR_start;
    seg_cnv->seg_inner_cnv[seg_index].imp_PR_end = cnvs[cnv_index].inner_cnv.imp_PR_end;
    seg_cnv->seg_inner_cnv[seg_index].num_larger_imp_PR_TLEN = cnvs[cnv_index].inner_cnv.num_larger_imp_PR_TLEN;
    seg_cnv->seg_inner_cnv[seg_index].evidence_count = cnvs[cnv_index].inner_cnv.evidence_count;
    seg_cnv->seg_inner_cnv[seg_index].num_merged_CNVs = cnvs[cnv_index].inner_cnv.num_merged_CNVs;

    // only reset this if it is true
    //
    if (cnvs[cnv_index].inner_cnv.valid_cnv)
        seg_cnv->seg_inner_cnv[seg_index].valid_cnv = cnvs[cnv_index].inner_cnv.valid_cnv;
}

void dynamicMemAllocateInnerCNVbreakpoints(INNER_CNV *seg_inner_cnv, uint32_t index) {
    if (seg_inner_cnv[index].cnv_breakpoints == NULL) {
        seg_inner_cnv[index].breakpoint_size = 0;
        seg_inner_cnv[index].breakpoint_capacity = 10;
        seg_inner_cnv[index].cnv_breakpoints = calloc(seg_inner_cnv[index].breakpoint_capacity, sizeof(CNV_Breakpints));
    } else {
        seg_inner_cnv[index].breakpoint_capacity += 10;
        seg_inner_cnv[index].cnv_breakpoints = realloc(seg_inner_cnv[index].cnv_breakpoints, \
                            seg_inner_cnv[index].breakpoint_capacity * sizeof(CNV_Breakpints));
        failureExit(seg_inner_cnv[index].cnv_breakpoints, "dynamicMemAllocateInnerCNVbreakpoints\n");
    }
}

void mergeCNVsFromSameSegment(Segmented_CNV_Array *seg_cnv_array, CNV_Array *cnv_array, User_Input *user_inputs, Simple_Stats *equal_window_stats) {
    uint32_t i;
    for (i=0; i<seg_cnv_array->size; ++i) {
        uint32_t seg_start = seg_cnv_array->seg_cnvs[i].segment.start;
        uint32_t seg_end   = seg_cnv_array->seg_cnvs[i].segment.end;

        //if (strcmp(seg_cnv_array->chromosome_id, "22") == 0)
        //    printf("stopped in mergeCNV for chr22 \n");

        if (seg_start == 33838000 || seg_end == 106455000)
            printf("stopped in mergeCNV\n");

        // initialize the inner_cnv variable based on the cnv_index_size
        //
        if (seg_cnv_array->seg_cnvs[i].seg_cnv_index_size > 0)
            segmentInnerCNVInit(&seg_cnv_array->seg_cnvs[i], seg_cnv_array->seg_cnvs[i].seg_cnv_index_size);

        // walk through all enclosed fragmented CNVs
        //
        if (seg_cnv_array->seg_cnvs[i].seg_inner_cnv == NULL)
            continue;

        // Now process all enclosed CNVs one at a time
        //  =====================================================================  SLM segment
        //    ---------- CNV1(Dup/Del)       -------------- CNV2 (Dup/Del) etc.
        //
        //  To merge or not depends on the type of CNVs
        //  if they are the same type: yes
        //  if they are different types: no, will leave it to the review to decide
        //
        uint32_t prev_seg_inner_cnv_index=0;

        uint32_t j;
        for (j=0; j<seg_cnv_array->seg_cnvs[i].seg_cnv_index_size; ++j) {
            // initialize inner cnv contents
            //
            segmentInnerCNVInit2(&seg_cnv_array->seg_cnvs[i], j);

            // get current cnv index associated with this segment and use its inner_cnv to get the breakpoint coordinates
            //
            uint32_t cnv_idx   = seg_cnv_array->seg_cnvs[i].seg_cnv_indices[j];
            uint32_t cnv_start = cnv_array->cnvs[cnv_idx].inner_cnv.start;
            uint32_t cnv_end   = cnv_array->cnvs[cnv_idx].inner_cnv.end;

            //if (cnv_start == 20099245 || cnv_start == 18003430 || cnv_start == 17821157)
            //    printf("stop 11\n");

            //if (cnv_end == 20100939 || cnv_end == 39387669 || cnv_end == 18060457)
            //    printf("stop 12\n");

            // check if current CNV is outside the current segment as Ns regions might void this CNV
            // the current CNV is here might be because of the none-Ns regions is short, but equal-size 
            // window makes it larger.
            //  (CNV1) start1 ------------- end1                                (CNV2) start2 ----------------- end2
            //                                    start =================== end
            //                                                segment
            //
            if (cnv_end < seg_start || seg_end < cnv_start) {
                cnv_array->cnvs[cnv_idx].inner_cnv.valid_cnv = false;
                continue;
            }

            // always handle the left most breakpoint of the segment with the first CNV
            //
            uint32_t seg_inner_cnv_index = seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size;

            copyCNV(&seg_cnv_array->seg_cnvs[i], seg_inner_cnv_index, cnv_array->cnvs, cnv_idx);
            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].num_merged_CNVs = 1;
            
            // only needs to do the merge if seg_cnv_array->seg_cnvs[i].size >= 1 and j >= 1
            //
            if (j >= 1) {
                // if equal-bin size is 1000, I will do the merge earlier if two CNVs are within 1000bp
                // Therefore, any low-mappability regions with size < 1000 should be ignored.
                // To calculate the average length of low-mappability regions with the length > 1000bp
                // less GRCh37_lowmappabilityall.sorted.merged.bed | awk 'var=$3-$2 {print var}' | 
                // awk '$1>=1000 {print}' | awk '{sum+=$1} END {print sum}' 
                // 111258008
                // less GRCh37_lowmappabilityall.sorted.merged.bed | awk 'var=$3-$2 {print var}' | awk '$1>=1000 {print}' | wc
                //   26361   26361  133305
                // the ave length of low-mappability with size >= 1000 is: 111258008 / 26361 = 4220
                //
                // In addition, according to the PennCNV web, they suggested to merge adjacent CNVs if the gap between
                // them are within 20% of total combined CNV length.
                //
                uint32_t gap_length = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].start - seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].end;
                uint32_t total_combined_CNV_length = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].end - seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].start;
                double gap_ratio = (double) gap_length / (double) total_combined_CNV_length;

                if (seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].cnv_type == \
                        seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].cnv_type) {
                    if (gap_length <= 4220 || gap_ratio <= 0.2) {
                        // check to see if one of them are fully supported by the evidence.
                        // if so, don't merge
                        //
                        int total_breakpoints = 0;
                        int evidences1 = obtainSupportingEvidences(&(seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index]), &total_breakpoints, seg_cnv_array->seg_cnvs[i].segment.log2ratio_mean, 1);
                        int evidences2 = obtainSupportingEvidences(&(seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index]), &total_breakpoints, seg_cnv_array->seg_cnvs[i].segment.log2ratio_mean, 1);

                        // update the valid_cnv variable
                        //
                        if (evidences1 >= 2)
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].valid_cnv = true;

                        if (evidences2 >= 2)
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].valid_cnv = true;

                        // need to merge them 
                        // However, after merging, all the supporting evidences regarding breakpoint will be gone for the middle CNVs
                        // In this case, we will take half of the breakpoint count as evidence
                        //
                        if (evidences1 < 6 && evidences2 < 6) {
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].end = \
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].end;

                        /*if (seg_start == 45466000) {
                                fprintf(stderr, "4th: %d\t%d\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", i, j, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].right_breakpoint, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].right_breakpoint, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].right_breakpoint_count, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].right_breakpoint_count );
                        }*/
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].right_breakpoint = \
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].right_breakpoint;
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].last_right_breakpoint_count = \
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].last_right_breakpoint_count;

                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].right_breakpoint_count += \
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].right_breakpoint_count;
                            //fprintf(stderr, "Right at %d\t%"PRIu32"\t%"PRIu32"\t%d\n", i, seg_start, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].right_breakpoint_count, j);

                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].last_num_larger_TLEN_right = \
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].last_num_larger_TLEN_right;
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_larger_TLEN_right += \
                                 seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].num_larger_TLEN_right;

                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].left_breakpoint_count = \
                                (int)((float)seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].left_breakpoint_count/2.0 + 0.5) + seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].left_breakpoint_count;
                            //fprintf(stderr, "Left at %d\t%"PRIu32"\t%"PRIu32"\t%d\n", i, seg_start, seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].left_breakpoint_count, j);
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_larger_TLEN_left = \
                                (int)((float)seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].num_larger_TLEN_left/2.0 + 0.5) + seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_larger_TLEN_left;

                            // need to re-calculate the average coverage
                            //
                            uint32_t prev_length = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].end - seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].start;
                            uint64_t prev_total_coverage = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].ave_coverage * prev_length;
                            uint32_t curr_length = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].end - seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].start;
                            uint64_t curr_total_coverage = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].ave_coverage * curr_length; 
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].ave_coverage = (float) (prev_total_coverage + curr_total_coverage) /(float) (prev_length + curr_length);

                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].qual = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].ave_coverage / equal_window_stats->stdev;

                            if (seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].valid_cnv 
                                    || seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].valid_cnv)
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].valid_cnv = true;

                            if (seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].num_larger_imp_PR_TLEN >
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_larger_imp_PR_TLEN) {
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].imp_PR_start = \
                                    seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].imp_PR_start;
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].imp_PR_end = \
                                     seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].imp_PR_end;
                                seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_larger_imp_PR_TLEN = \
                                     seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].num_larger_imp_PR_TLEN;

                                // after the merge, delete the current one
                                //
                                //if (seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].cnv_breakpoints)
                                //    free(seg_cnv_array->seg_cnvs[i].seg_inner_cnv[seg_inner_cnv_index].cnv_breakpoints);
                            }

                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv[prev_seg_inner_cnv_index].num_merged_CNVs++;
                            seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size--;
                        }
                    }
                }
            }

            prev_seg_inner_cnv_index = seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size;
            seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size++;
        }
    }
}

// type: 1 for low mappability; 2 for GC% <= 25%
//
void checkingFalsePositives(Segmented_CNV_Array *seg_cnv_array, char * chr_id, Bed_Info *low_complexity_bed_info, int type) {
    // get the total size need for all_starts_ends array
    // since we process one chromosome at a time, we don't need to use the full size
    //
    int32_t i=0, j=0, counter=0, capacity = 0;
    for (i=0; i<low_complexity_bed_info->size; i++) {
        if (strcmp(low_complexity_bed_info->coords[i].chrom_id, chr_id) == 0) {
            capacity += 2;
        }
    }

    for (i=0; i<seg_cnv_array->size; i++)
        capacity += seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size * 2;

    uint32_t *all_starts_ends = calloc(capacity, sizeof(uint32_t));

    // store length info of a low_complexity_bed_info region
    //
    khash_t(m32) *lc_start_hash = kh_init(m32);
    khash_t(m32) *lc_end_hash   = kh_init(m32);
    //khash_t(m32) *lc_start_end_lookup = kh_init(m32);
    //khash_t(m32) *lc_end_start_lookup = kh_init(m32);

    for (i=0; i<low_complexity_bed_info->size; i++) {
        if (strcmp(low_complexity_bed_info->coords[i].chrom_id, chr_id) == 0) {
            //uint32_t length = low_complexity_bed_info->coords[i].end - low_complexity_bed_info->coords[i].start;
            all_starts_ends[counter] = low_complexity_bed_info->coords[i].start; 
            setValueToKhashBucket32(lc_start_hash, all_starts_ends[counter], low_complexity_bed_info->coords[i].end);
            counter++;

            all_starts_ends[counter] = low_complexity_bed_info->coords[i].end;
            setValueToKhashBucket32(lc_end_hash, all_starts_ends[counter], low_complexity_bed_info->coords[i].start);
            counter++;
        }
    }

    // store segmented CNVs info
    //
    khash_t(m32) *cnv_start_i_hash = kh_init(m32);
    khash_t(m32) *cnv_start_j_hash = kh_init(m32);
    khash_t(m32) *cnv_end_i_hash   = kh_init(m32);
    khash_t(m32) *cnv_end_j_hash   = kh_init(m32);
    //khash_t(m32) *cnv_start_end_lookup = kh_init(m32);
    //khash_t(m32) *cnv_end_start_lookup = kh_init(m32);

    for (i=0; i<seg_cnv_array->size; i++) {
        for (j=0; j<seg_cnv_array->seg_cnvs[i].seg_inner_cnv_size; j++) {
            // + 1 is added to the start position so that it will be 1-based and also avoid 
            // start/end the same for the neighbouring CNVs
            // but need to remember add 1 when using ke_del for the start key
            //
            all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[j].start+1; 
            setValueToKhashBucket32(cnv_start_i_hash, all_starts_ends[counter], i);
            setValueToKhashBucket32(cnv_start_j_hash, all_starts_ends[counter], j);
            counter++;

            all_starts_ends[counter] = seg_cnv_array->seg_cnvs[i].seg_inner_cnv[j].end;
            setValueToKhashBucket32(cnv_end_i_hash, all_starts_ends[counter], i);
            setValueToKhashBucket32(cnv_end_j_hash, all_starts_ends[counter], j);
            counter++;
        }
    }

    if (counter != capacity)
        fprintf(stderr, "The number doesn't match at the checkingFalsePositives\n");

    // sorting
    //
    qsort(all_starts_ends, capacity, sizeof(uint32_t), compare);

    // intersect
    //
    khash_t(m32) *live_lc_starts_hash = kh_init(m32);
    khash_t(m32) *live_cnv_starts_i_hash = kh_init(m32);
    khash_t(m32) *live_cnv_starts_j_hash = kh_init(m32);

    khash_t(m32) *seen_lc_starts_hash = kh_init(m32);
    khash_t(m32) *seen_cnv_starts_hash = kh_init(m32);

    int32_t lc_start    = -1;             // need to use signed value as negative value might return
    int32_t lc_end      = -1;             // need to use signed value as negative value might return
    int32_t cnv_i_index = -1;             // need to use signed value as negative value might return
    int32_t cnv_j_index = -1;             // need to use signed value as negative value might return
    khiter_t iter;
    counter = 0;

    for (i=0; i<capacity; i++) {

        if (all_starts_ends[i] == 29059774) 
            printf("stopped A\n");

        // if a key is present in both start and end position, we have to process start first
        //
        bool start_first = false;
        if (checkm32KhashKey(lc_start_hash, all_starts_ends[i]) && checkm32KhashKey(lc_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_lc_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(cnv_start_i_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_i_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(lc_start_hash, all_starts_ends[i]) && checkm32KhashKey(cnv_end_i_hash, all_starts_ends[i]) 
                && !checkm32KhashKey(seen_lc_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (checkm32KhashKey(cnv_start_i_hash, all_starts_ends[i]) && checkm32KhashKey(lc_end_hash, all_starts_ends[i]) \
                && !checkm32KhashKey(seen_cnv_starts_hash, all_starts_ends[i]) )
            start_first = true;

        if (!start_first && (checkm32KhashKey(lc_end_hash, all_starts_ends[i]) ||
                                    checkm32KhashKey(cnv_end_i_hash, all_starts_ends[i]))) {

            // always decrease count if it is the end position
            //
            counter--;

            /*if (checkm32KhashKey(cnv_end_i_hash, all_starts_ends[i])) {
                  fprintf(stderr, "CNV end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            } else if (checkm32KhashKey(lc_end_hash, all_starts_ends[i])) {
                  fprintf(stderr, "low complexity end: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }*/

            // when current position is an end and counter >= 1, there is an intersect
            // There might be multiple low complexity regions associated with this CNV; need store all of them
            //
            if (counter >= 1) {
                uint32_t start = 0, end = 0;
                if (checkm32KhashKey(lc_end_hash, all_starts_ends[i])) {   // it is the low complexity end
                    lc_start = getValueFromKhash32(lc_end_hash, all_starts_ends[i]);

                    // now loop through live_cnv_hash tables
                    //
                    for (iter=kh_begin(live_cnv_start_i_hash); iter!=kh_end(live_cnv_starts_i_hash); ++iter) {
                        if (kh_exist(live_cnv_starts_i_hash, iter)) {
                            cnv_i_index = kh_value(live_cnv_starts_i_hash, iter);

                            cnv_j_index = getValueFromKhash32(live_cnv_starts_j_hash, kh_key(live_cnv_starts_i_hash, iter));

                            // calculate the intersected length
                            //
                            if (lc_start >= seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start) {
                                start = lc_start;
                            } else {
                                start = seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start;
                            }

                            if (all_starts_ends[i] < seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].end) {
                                end = all_starts_ends[i];
                            } else {
                                end = seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].end;
                            }

                            if (end < start) {
                                // now overlapping, just skip
                                //
                                continue;
                                //fprintf(stderr, "Order is wrong start %"PRIu32" and end %"PRIu32"\n", start, end);
                            }

                            if (type == 1) {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].low_mapp_length += end - start + 1;
                            } else if (type == 2) {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].gc_lt25pct_length += end -start + 1;
                            } else {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].gc_gt85pct_length += end -start + 1;
                            }
                        }
                    }

                    // delete the low complexity start entry. Can't use all_starts_ends[i] because it is the end position
                    //
                    iter = kh_get(m32, live_lc_starts_hash, lc_start);
                    kh_del(m32, live_lc_starts_hash, iter);

                } else {    // it is the CNV end
                    cnv_i_index = getValueFromKhash32(cnv_end_i_hash, all_starts_ends[i]);
                    cnv_j_index = getValueFromKhash32(cnv_end_j_hash, all_starts_ends[i]);

                    // walk through live_lc_starts_hash
                    //
                    for (iter=kh_begin(live_lc_starts_hash); iter!=kh_end(live_lc_starts_hash); ++iter) {
                        if (kh_exist(live_lc_starts_hash, iter)) {
                            lc_end = kh_value(live_lc_starts_hash, iter);
                            lc_start = kh_key(live_lc_starts_hash, iter);

                            // calculate the intersected length
                            //
                            if (lc_start >= seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start) {
                                start = lc_start;
                            } else {
                                start = seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start;
                            }

                            if (all_starts_ends[i] < seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].end) {
                                end = all_starts_ends[i];
                            } else {
                                end = seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].end;
                            }

                            if (end < start) continue;

                            if (type == 1) {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].low_mapp_length += end - start + 1;
                            } else if (type == 2) {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].gc_lt25pct_length += end -start + 1;
                            } else {
                                seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].gc_gt85pct_length += end -start + 1;
                            }
                        }
                    }

                    // delete the live_cnv_starts at i and j; Dont use all_starts_ends[i] as it is the end position
                    // as it is the end of current CNV, we don't need it any more. So delete it
                    // because I added 1 for CNV start for 1-based bedfile, so need to adjust accordingly
                    //
                    iter = kh_get(m32, live_cnv_starts_i_hash, seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start+1);
                    kh_del(m32, live_cnv_starts_i_hash, iter);
                    iter = kh_get(m32, live_cnv_starts_j_hash, seg_cnv_array->seg_cnvs[cnv_i_index].seg_inner_cnv[cnv_j_index].start+1);
                    kh_del(m32, live_cnv_starts_j_hash, iter);
                }
            } else {
                // clean-up live_lc_starts_hash and live_cnv_starts_i_hash and live_cnv_starts_j_hash
                //
                for (iter=kh_begin(live_lc_starts_hash); iter!=kh_end(live_lc_starts_hash); ++iter) {
                    if (kh_exist(live_lc_starts_hash, iter))
                        kh_del(m32, live_lc_starts_hash, iter);
                }

                for (iter=kh_begin(live_cnv_starts_i_hash); iter!=kh_end(live_cnv_starts_i_hash); ++iter) {
                    if (kh_exist(live_cnv_starts_i_hash, iter))
                        kh_del(m32, live_cnv_starts_i_hash, iter);
                }

                for (iter=kh_begin(live_cnv_starts_j_hash); iter!=kh_end(live_cnv_starts_j_hash); ++iter) {
                    if (kh_exist(live_cnv_starts_j_hash, iter))
                        kh_del(m32, live_cnv_starts_j_hash, iter);
                }
            }
        } else {    // start position
            counter++;

            // get current low complexity start position in lc_length
            //
            if (checkm32KhashKey(lc_start_hash, all_starts_ends[i])) {
                lc_end = getSignedValueFromKhash32(lc_start_hash, all_starts_ends[i]);

                if (lc_end == -1) {
                    fprintf(stderr, "Something went wrong with low complexity index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                setValueToKhashBucket32(live_lc_starts_hash, all_starts_ends[i], lc_end);
                setValueToKhashBucket32(seen_lc_starts_hash, all_starts_ends[i], lc_end);

                /*if (type == 1) {
                    fprintf(stderr, "Low mappability start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
                } else if (type == 2) {
                    fprintf(stderr, "GC%% less than or equal to 25%% start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
                } else {
                    fprintf(stderr, "GC%% greater than or equal to 85%% start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
                }*/
            }

            // get current cnv start position index
            //
            if (checkm32KhashKey(cnv_start_i_hash, all_starts_ends[i])) {
                cnv_i_index = getSignedValueFromKhash32(cnv_start_i_hash, all_starts_ends[i]);
                cnv_j_index = getSignedValueFromKhash32(cnv_start_j_hash, all_starts_ends[i]);

                if (cnv_i_index == -1 || cnv_j_index == -1) {
                    fprintf(stderr, "Something went wrong with seg index at %"PRIu32"\n", all_starts_ends[i]);
                    exit(EXIT_FAILURE);
                }

                setValueToKhashBucket32(live_cnv_starts_i_hash, all_starts_ends[i], cnv_i_index);
                setValueToKhashBucket32(live_cnv_starts_j_hash, all_starts_ends[i], cnv_j_index);
                setValueToKhashBucket32(seen_cnv_starts_hash, all_starts_ends[i], cnv_i_index);

                //fprintf(stderr, "CNV start: %"PRIu32" with index %"PRIu32"\n", all_starts_ends[i], counter);
            }
                
        }
    }

    // clean-up
    //
    if (all_starts_ends != NULL) {
        free(all_starts_ends);
        all_starts_ends = NULL;
    }

    kh_destroy(m32, lc_start_hash);
    kh_destroy(m32, lc_end_hash);
    kh_destroy(m32, cnv_start_i_hash);
    kh_destroy(m32, cnv_start_j_hash);
    kh_destroy(m32, cnv_end_i_hash);
    kh_destroy(m32, cnv_end_j_hash);
    kh_destroy(m32, live_lc_starts_hash);
    kh_destroy(m32, live_cnv_starts_i_hash);
    kh_destroy(m32, live_cnv_starts_j_hash);
    kh_destroy(m32, seen_cnv_starts_hash);
    kh_destroy(m32, seen_lc_starts_hash);
}

int obtainSupportingEvidences(INNER_CNV* inner_cnv, int* total_breakpoints, double log2ratio, int type) {

    int evidences=0;
    if (inner_cnv->left_breakpoint > 0) {
        if (inner_cnv->left_breakpoint_count >= 3 \
                || (inner_cnv->left_breakpoint_count >= 2 && inner_cnv->num_larger_TLEN_left >= 2) ) 
            evidences++;
    }

    if (inner_cnv->right_breakpoint > 0 ) {
        if (inner_cnv->right_breakpoint_count >=3 \
                || (inner_cnv->right_breakpoint_count >=2 && inner_cnv->num_larger_TLEN_right >=2) )
            evidences++;
    }

    if (type == 2 && (inner_cnv->right_breakpoint_count != inner_cnv->last_right_breakpoint_count)) {
        uint32_t mid_right_bpt_count = inner_cnv->right_breakpoint_count - inner_cnv->last_right_breakpoint_count;
        inner_cnv->right_breakpoint_count = (int)((float)mid_right_bpt_count/2.0 + 0.5) + inner_cnv->last_right_breakpoint_count;
    }

    if (type == 2 && (inner_cnv->num_larger_TLEN_right != inner_cnv->last_num_larger_TLEN_right)) { 
        uint32_t mid_large_TLEN_count = inner_cnv->num_larger_TLEN_right - inner_cnv->last_num_larger_TLEN_right;
        inner_cnv->num_larger_TLEN_right = (int)((float)mid_large_TLEN_count/2.0 + 0.5) + inner_cnv->last_num_larger_TLEN_right;
    }

    *total_breakpoints = inner_cnv->left_breakpoint_count + inner_cnv->right_breakpoint_count;
    //if (*total_breakpoints >= 25) evidences++;

    if (inner_cnv->num_larger_TLEN_left >= 2)  evidences++;
    if (inner_cnv->num_larger_TLEN_right >= 2) evidences++;

    if (inner_cnv->num_larger_imp_PR_TLEN >= 4) {
        evidences += 2;
    } else if (inner_cnv->num_larger_imp_PR_TLEN >= 2) {
        evidences++;
    }

    // the following is needed for the merged ones as many merged ones will erase the breakpoint into
    //
    if (log2ratio <= -0.5 || log2ratio >= 0.5)
        evidences++;

    if ( (log2ratio >= 0.99 && log2ratio <= 1.1) || (log2ratio >= -1.1 && log2ratio <= -0.99) )
        evidences++;

    return evidences;
}

void generateSegmentedCNVs(Segmented_CNV_Array **seg_cnv_array, Chromosome_Tracking *chrom_tracking, Simple_Stats *equal_window_stats, Stats_Info *stats_info, User_Input *user_inputs) {
    FILE *fp = fopen(user_inputs->segmented_vcf_output_file, "w");
    FILE *sfh = fopen(user_inputs->simple_segmented_vcf_file, "w");
    fprintf(sfh, "chr\tstart\tend\ttype\tquality\tvar-reads\ttotal-reads\n");

    generateVCF_MetaData(user_inputs, chrom_tracking, fp);

    uint32_t i, j, k;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        for (j=0; j<seg_cnv_array[i]->size; j++) {
            //if ( seg_cnv_array[i]->seg_cnvs[j].segment.end == 106455000)
            //    printf("segment output stopped!\n");
            
            double log2ratio = seg_cnv_array[i]->seg_cnvs[j].segment.log2ratio_mean;

            if (user_inputs->debug_ON) {
                fprintf(fp, "\nSegment: %"PRIu32"\t%"PRIu32"\n", seg_cnv_array[i]->seg_cnvs[j].segment.start, seg_cnv_array[i]->seg_cnvs[j].segment.end);
            } else {
                if (log2ratio > -0.5 && log2ratio < 0.5 )
                    continue;
            }

            for (k=0; k<seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv_size; ++k) {
                char VALID_CNV[10];
                strcpy(VALID_CNV, "No");

                if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_merged_CNVs >=2
                                    && seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].valid_cnv) 
                    strcpy(VALID_CNV, "Yes");

                char CNV[10];
                (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].cnv_type == 'L') ? strcpy(CNV, "DEL") : strcpy(CNV, "INS");

                char FILTER[50];
                strcpy(FILTER, "PASS");

                char GT[10];
                strcpy(GT, "./1");

                // now calculate QUAL score using zscore for simple CNV output
                //
                double qual = (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].ave_coverage \
                                - equal_window_stats->average_coverage) / equal_window_stats->stdev;

                // need to re-access the evidence and set the FILTER accordingly
                //
                int total_breakpoints = 0;
                int evidences = obtainSupportingEvidences(&(seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k]), &total_breakpoints, log2ratio, 2);

                // 2.576 is z-score for 99% confident interval cutoff
                //
                if (abs(qual) < 2.576 && ((evidences <= 2 && total_breakpoints < 25) || (evidences == 0 && strcmp(VALID_CNV, "Yes") == 0))) {
                    strcpy(FILTER, "qualLessThan99pctCI;littleBreakpointAndImpSupport");
                    strcpy(GT, "./.");
                }

                if (evidences == 1 && strcmp(VALID_CNV, "No") == 0) {
                    strcpy(FILTER, "littleBreakpointAndImpSupport");
                    strcpy(GT, "./.");
                }

                if (evidences == 0 && strcmp(VALID_CNV, "No") == 0) {
                    strcpy(FILTER, "noBreakpointAndImpSupport");
                    strcpy(GT, "./.");
                }

                // must be signed
                //
                int32_t svLen = seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].end - seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].start;
                if (svLen < user_inputs->min_cnv_length) continue;

                // checking low mappability size and fraction before svLen set to negative for DEL
                //
                double lowMapRatio = (double)seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].low_mapp_length / (double) svLen;
                if (lowMapRatio > 0.5 || (lowMapRatio == 0.5 && evidences <= 2)) {
                    strcpy(FILTER, "lowMappability");
                    strcpy(GT, "./.");
                }

                // checking GC% <= 25% size and fraction before svLen set to negative for DEL
                //
                double gc_lt25_pct_ratio = (double) seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].gc_lt25pct_length / (double) svLen;
                if (gc_lt25_pct_ratio > 0.5 || (gc_lt25_pct_ratio == 0.5 && evidences <= 2)) {
                    strcpy(FILTER, "gcContent_lt25pct");;
                    strcpy(GT, "./.");
                }

                // checking GC% >= 85% size and fraction before svLen set to negative for DEL
                //
                double gc_gt85_pct_ratio = (double) seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].gc_gt85pct_length / (double) svLen;
                if (gc_gt85_pct_ratio > 0.5 || (gc_gt85_pct_ratio == 0.5 && evidences <= 2)) {
                    strcpy(FILTER, "gcContent_gt85pct");;
                    strcpy(GT, "./.");
                }

                double combined_low_complexity_ratio = lowMapRatio + gc_lt25_pct_ratio + gc_gt85_pct_ratio;

                if (lowMapRatio <= 0.5 && gc_lt25_pct_ratio <= 0.5 && gc_gt85_pct_ratio <= 0.5) {
                    if (combined_low_complexity_ratio > 0.5 || (combined_low_complexity_ratio == 0.5 && evidences <= 2)) {
                        strcpy(FILTER, "lowComplexityPct");
                        strcpy(GT, "./.");
                    }
                }

                // during HG002 concordance analysis, many FNs were rescued by the following condition
                //
                if (lowMapRatio > 0.5) {
                    double haploid_cutoff = equal_window_stats->average_coverage - equal_window_stats->zScore;
                    if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].ave_coverage * lowMapRatio < haploid_cutoff) {
                        if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_larger_imp_PR_TLEN >= 4 && evidences >= 3) {
                            strcpy(FILTER, "PASS");
                            strcpy(GT, "./1");
                        }
                    }
                }

                if (log2ratio > -0.5 && log2ratio < 0.5 ) {
                    if (evidences < 3 || seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_larger_imp_PR_TLEN < 3) {
                        strcpy(FILTER, "failLog2ratio");
                        strcpy(GT, "./.");
                    }
                }

                if (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].cnv_type == 'L')
                    svLen *= -1;

                // now calculate QUAL score using zscore for simple CNV output
                //
                double qual2 = (seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].ave_coverage - equal_window_stats->average_coverage) / equal_window_stats->stdev;

                if (strcmp(FILTER, "PASS") == 0) {
                    // For the sample CNV output in TSV format
                    // chr, start, end, type, quality, var reads, total reads
                    // calculate the total number of reads based on Atlas CNV approach
                    //
                    uint32_t cnv_start = seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].start;
                    uint32_t cnv_end   = seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].end;
                    int32_t total_reads = seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].ave_coverage * (cnv_end-cnv_start) / stats_info->read_cov_stats->read_length;

                    fprintf(sfh, "%s\t%"PRIu32"\t%"PRIu32"\t%s\t%.2f\t.\t%"PRIu32"\n", chrom_tracking->chromosome_ids[i], cnv_start, cnv_end, CNV, qual2, (uint32_t) total_reads);
                }

                fprintf(fp, "%s\t%"PRIu32"\t.\tN\t%s\t%.2f\t%s\tEND=%"PRIu32";SVLEN=%"PRId32";SVTYPE=%s;AVGCOV=%.2f;BPTL=%"PRIu32";BPTLCOUNT=%"PRIu32";BPTLTLEN=%"PRIu32";BPTR=%"PRIu32";BPTRCOUNT=%"PRIu32";BPTRTLEN=%"PRIu32";IMPPRLEN=%"PRIu16";MergedCNVs=%d;EvidenceCount=%d;ValidCNVForOneOfMergedCNVs=%s;lowMappabilityRatio=%.2f;gcContent_lt25pct=%.2f;gcContent_gt85pct=%.2f;log2ratio=%.2f\tGT\t%s\n", \
                        chrom_tracking->chromosome_ids[i],
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].start, CNV, qual, FILTER,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].end, svLen, CNV,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].ave_coverage,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].left_breakpoint,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].left_breakpoint_count,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_larger_TLEN_left,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].right_breakpoint,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].right_breakpoint_count,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_larger_TLEN_right,
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_larger_imp_PR_TLEN, 
                        seg_cnv_array[i]->seg_cnvs[j].seg_inner_cnv[k].num_merged_CNVs, evidences, 
                        VALID_CNV, lowMapRatio, gc_lt25_pct_ratio, gc_gt85_pct_ratio, log2ratio, GT);

            }
        }
    }

    fclose(fp);
    fclose(sfh);
}
