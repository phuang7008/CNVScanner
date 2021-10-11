/*
 * =====================================================================================
 *
 *       Filename:  breakpoints.c
 *
 *    Description:  the detail implementation of breakpoints header file
 *
 *        Version:  1.0
 *        Created:  10/05/2021 02:04:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX, USA
 *
 * =====================================================================================
 */

#include "breakpoints.h"

void BreakpointArrayInit(Breakpoint_Array *breakpoint_array, Chromosome_Tracking *chrom_tracking) {
    breakpoint_array->size  = chrom_tracking->number_of_chromosomes;
    breakpoint_array->chrom_ids = calloc(breakpoint_array->size, sizeof(char*));
    breakpoint_array->bpts_per_chr = 
            calloc(breakpoint_array->size, sizeof(Breakpoints_Per_Chromosome));

    uint32_t i;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        breakpoint_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);

        breakpoint_array->bpts_per_chr[i].capacity = INIT_SIZE;
        breakpoint_array->bpts_per_chr[i].size = 0;
        breakpoint_array->bpts_per_chr[i].breakpoints = 
                calloc(breakpoint_array->bpts_per_chr[i].capacity, sizeof(Breakpoint));
    }
}

void BreakpointArrayDestroy(Breakpoint_Array *breakpoint_array) {
    uint32_t i;
    for (i=0; i<breakpoint_array->size; i++) {
        if (breakpoint_array->chrom_ids[i] != NULL) {
            free(breakpoint_array->chrom_ids[i]);
            breakpoint_array->chrom_ids[i]=NULL;
        }
    }

    if (breakpoint_array->chrom_ids != NULL) {
        free(breakpoint_array->chrom_ids);
        breakpoint_array->chrom_ids = NULL;
    }

    if (breakpoint_array->bpts_per_chr) {
        free(breakpoint_array->bpts_per_chr);
        breakpoint_array->bpts_per_chr = NULL;
    }
}

void PairedReadsCrossBreakpointsArrayInit(Paired_Reads_Cross_Breakpoints_Array *pread_x_bpts_array, Chromosome_Tracking *chrom_tracking) {
    pread_x_bpts_array->size = chrom_tracking->number_of_chromosomes;
    pread_x_bpts_array->chrom_ids = calloc(pread_x_bpts_array->size, sizeof(char*));
    pread_x_bpts_array->preads_x_bpts_per_chr = 
                    calloc(pread_x_bpts_array->size, sizeof(Paired_Reads_Cross_Breakpoints_Per_Chromosome));

    uint32_t i;
    for (i=0; i<pread_x_bpts_array->size; i++) {
        pread_x_bpts_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);

        pread_x_bpts_array->preads_x_bpts_per_chr[i].capacity = INIT_SIZE;
        pread_x_bpts_array->preads_x_bpts_per_chr[i].size = 0;
        pread_x_bpts_array->preads_x_bpts_per_chr[i].preads_x_bpts =
                calloc(pread_x_bpts_array->preads_x_bpts_per_chr[i].capacity, sizeof(Paired_Reads_Cross_Breakpoints));
    }
}

void PairedReadsCrossBreakpointsArrayDestroy(Paired_Reads_Cross_Breakpoints_Array *pread_x_bpts_array) {
    uint32_t i, j;
    for (i=0; i< pread_x_bpts_array->size; i++) {
        if (pread_x_bpts_array->chrom_ids[i] != NULL) {
            free(pread_x_bpts_array->chrom_ids[i]);
            pread_x_bpts_array->chrom_ids[i] = NULL;
        }

        for (j=0; j<pread_x_bpts_array->preads_x_bpts_per_chr[i].size; j++) {
            if (pread_x_bpts_array->preads_x_bpts_per_chr[i].preads_x_bpts->read_name != NULL) {
                free(pread_x_bpts_array->preads_x_bpts_per_chr[i].preads_x_bpts->read_name);
                pread_x_bpts_array->preads_x_bpts_per_chr[i].preads_x_bpts->read_name = NULL;
            }
        }
    }

    if (pread_x_bpts_array->chrom_ids != NULL) {
        free(pread_x_bpts_array->chrom_ids);
        pread_x_bpts_array->chrom_ids = NULL;
    }

    if (pread_x_bpts_array->preads_x_bpts_per_chr != NULL) {
        free(pread_x_bpts_array->preads_x_bpts_per_chr);
        pread_x_bpts_array->preads_x_bpts_per_chr = NULL;
    }
}

void BreakpointStatsArrayInit(Breakpoint_Stats_Array *bpt_stats_array, Chromosome_Tracking *chrom_tracking) {
    bpt_stats_array->size = chrom_tracking->number_of_chromosomes;
    bpt_stats_array->chrom_ids = calloc(bpt_stats_array->size, sizeof(char*));
    bpt_stats_array->bp_stats_per_chr = calloc(bpt_stats_array->size, sizeof(Breakpoint_Stats_Per_Chromosome));

    uint32_t i;
    for (i=0; i<bpt_stats_array->size; i++) {
        bpt_stats_array->chrom_ids[i] = strdup(chrom_tracking->chromosome_ids[i]);

        bpt_stats_array->bp_stats_per_chr[i].size = 0;
        bpt_stats_array->bp_stats_per_chr[i].capacity = INIT_SIZE;
        bpt_stats_array->bp_stats_per_chr[i].bpts_stats =
                calloc(bpt_stats_array->bp_stats_per_chr[i].capacity, sizeof(Breakpoint_Stats));
    }
}

void BreakpointStatsArrayDestroy(Breakpoint_Stats_Array *bpt_stats_array) {
    uint32_t i;
    for (i=0; i< bpt_stats_array->size; i++) {
        if (bpt_stats_array->chrom_ids[i] != NULL) {
            free(bpt_stats_array->chrom_ids[i]);
            bpt_stats_array->chrom_ids[i] = NULL;
        }

        if (bpt_stats_array->bp_stats_per_chr[i].bpts_stats != NULL) {
            free(bpt_stats_array->bp_stats_per_chr[i].bpts_stats);
            bpt_stats_array->bp_stats_per_chr[i].bpts_stats = NULL;
        }
    }

    if (bpt_stats_array->chrom_ids != NULL) {
        free(bpt_stats_array->chrom_ids);
        bpt_stats_array->chrom_ids = NULL;
    }

    if (bpt_stats_array->bp_stats_per_chr != NULL) {
        free(bpt_stats_array->bp_stats_per_chr);
        bpt_stats_array->bp_stats_per_chr = NULL;
    }
    //if (bpt_stats_array->bp_stats_per_chr[i].breakpoints !=NULL) { }  // should be handled by BreakpointArrayDestroy()
}

uint32_t fetchBreakpointArrayChrIndex(Breakpoint_Array *breakpoint_array, char * chrom_id) {
    uint32_t i;
    for (i=0; i<breakpoint_array->size; i++) {
        if (strcmp(chrom_id, breakpoint_array->chrom_ids[i]) == 0) {
            return i;
        }
    }
    
    return breakpoint_array->size + 10;
}

// for type: 1 -> soft clipping while 2 -> hard clipping
//
void storeCurrentReadBreakpointInfo(uint32_t current_ref_pos, bam1_t *rec, Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, khash_t(khStrInt) *breakpoint_pairs_hash, int type) {

    // if paired reads are not on the same chromosome, we won't process it as it is translocation
    //
    if (rec->core.tid != rec->core.mtid) return;

    if (bpt_arr == NULL) return;

    uint32_t bpt_arr_ind = bpt_arr->bpts_per_chr[bpt_chr_idx].size;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].type = type;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_index = bpt_arr_ind;
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].breakpoint_position = current_ref_pos;

    // according to sam.h,  @field  pos: 0-based leftmost coordinate of the current read
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_read_start = rec->core.pos;

    // according to sam.h, @field  mpos: 0-based leftmost coordinate of next read in template
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_read_start = rec->core.mpos;

    // according to sam.h, @field  isize: observed template length ("insert size")
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].gap_distance_TLEN = rec->core.isize;

    // according to sam.h, @field  tid: chromosome ID, defined by sam_hdr_t
    //
    //bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].current_chr_id = rec->core.tid;

    // according to sam.h, @field  mtid: chromosome ID of next read in template, defined by sam_hdr_t
    //
    //bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_chr_id = rec->mtid;

    // @field  data: all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    // the read name should be accessed using bam_get_qname
    //
    bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].read_name = strdup(bam_get_qname(rec));

    // now we need to record this read name in the hash table for quick lookup
    //
    khiter_t iter = kh_get(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec));
    if (iter == kh_end(breakpoint_pairs_hash)) {
        // first time
        //
        int absent;
        iter = kh_put(khStrInt, breakpoint_pairs_hash, bam_get_qname(rec), &absent);
        if (absent) {
            kh_key(breakpoint_pairs_hash, iter) = strdup(bam_get_qname(rec));
            kh_value(breakpoint_pairs_hash, iter) = bpt_arr_ind;
        }
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_index = -1;
    } else {
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[bpt_arr_ind].mate_index = kh_value(breakpoint_pairs_hash, iter);

        // now need to go the mate data and reset its mate index
        //
        bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[kh_value(breakpoint_pairs_hash, iter)].mate_index = bpt_arr_ind;
    }

    bpt_arr->bpts_per_chr[bpt_chr_idx].size++;

    if (bpt_arr->bpts_per_chr[bpt_chr_idx].size + 10 >
            bpt_arr->bpts_per_chr[bpt_chr_idx].capacity) {

        // need to dynamically increase the array size
        // TODO
        dynamicBreakpointPerChrArraySizeIncrease(&bpt_arr->bpts_per_chr[bpt_chr_idx]);
    }
}

void dynamicBreakpointPerChrArraySizeIncrease(Breakpoints_Per_Chromosome *bpts_per_chr) {
    if (bpts_per_chr == NULL) return;

    bpts_per_chr->capacity += INIT_SIZE;

    if (bpts_per_chr->breakpoints) {
        bpts_per_chr->breakpoints = 
            realloc(bpts_per_chr->breakpoints, bpts_per_chr->capacity * sizeof(Breakpoint));

        failureExit(bpts_per_chr->breakpoints, "Breakpoints_Per_Chromosome->breakpoints");
    } else {
        fprintf(stderr, "Error: The bpts_per_chr.breakpoints is NULL\n");
        exit(EXIT_FAILURE);
    }
}

void outputBreakpointArray(Breakpoint_Array *bpt_arr) {
    FILE *fp = fopen("Breakpoint_Details_for_Debugging.txt", "w");
    uint32_t i, j;
    for (i=0; i<bpt_arr->size; i++) {
        fprintf(fp, "For chromosome: %s\n", bpt_arr->chrom_ids[i]);

        for (j=0; j<bpt_arr->bpts_per_chr[i].size ; j++) {
            fprintf(fp, "  Type: %d\n", bpt_arr->bpts_per_chr[i].breakpoints[j].type);
            fprintf(fp, "  breakpoint_position: %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].breakpoint_position);
            fprintf(fp, "  read_name:     %s\n", bpt_arr->bpts_per_chr[i].breakpoints[j].read_name);
            fprintf(fp, "  current_index: %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].current_index);
            fprintf(fp, "  mate_index:    %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].mate_index);
            fprintf(fp, "  current_read_start: %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].current_read_start);
            fprintf(fp, "  mate_read_start:    %"PRIu32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].mate_read_start);
            fprintf(fp, "  gap_distance_TLEN:  %"PRId32"\n", bpt_arr->bpts_per_chr[i].breakpoints[j].gap_distance_TLEN);
        }
    }
}

void storePairedReadsCrossBreakpointsPerChr(Breakpoint_Array *bpt_arr, uint32_t bpt_chr_idx, Paired_Reads_Cross_Breakpoints_Array *preads_x_bpts_arr, uint32_t preads_x_bpts_idx, bam_hdr_t *header, hts_idx_t *sfh_idx, samFile *sfh) {
    // need to store processed breakpoints
    //
    khash_t(m32) *seen_breakpoints_hash = kh_init(m32);

    // walk through breakpoints of the chromosome with the bpt_chr_idx
    //
    uint32_t k;
    for (k=0; k<bpt_arr->bpts_per_chr[bpt_chr_idx].size; k++) {
        // breakpoint position is 0-based
        // now form a region for search
        //
        char region[200];
        uint32_t bpt_pos = bpt_arr->bpts_per_chr[bpt_chr_idx].breakpoints[k].breakpoint_position;
        // check to see if we have seen this breakpoint position
        //
        khiter_t iter_h = kh_get(m32, seen_breakpoints_hash, bpt_pos);
        int absent;
        if (iter_h == kh_end(seen_breakpoints_hash)) {
            iter_h = kh_put(m32, seen_breakpoints_hash, bpt_pos, &absent);
            if (absent) {
                kh_key(seen_breakpoints_hash, iter_h) = bpt_pos;
                kh_value(seen_breakpoints_hash, iter_h) = bpt_pos;
            }
        } else {
            continue;
        }

        sprintf(region, "%s:%"PRIu32"-%"PRIu32"", bpt_arr->chrom_ids[bpt_chr_idx], bpt_pos, bpt_pos+1);
        hts_itr_t *hts_itr = sam_itr_querys(sfh_idx, header, region);

        if (hts_itr == NULL) {
            fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", region);
            exit(EXIT_FAILURE);
        }

        // need to record the query read name for quick lookup
        //
        khash_t(khStrInt) *breakpoint_pairs_hash = kh_init(khStrInt);

        bam1_t *b = bam_init1();
        while (sam_itr_next(sfh, hts_itr, b) >= 0) {
            // if paired reads are not on the same chromosome, don't process it, just skip it
            //
            if (b->core.tid != b->core.mtid) return;

            // now add the current read to the Paired_Reads_Cross_Breakpoints_Per_Chromosome
            //
            uint32_t arr_ind = preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].size;
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].breakpoint_position = bpt_pos;
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].current_start_position  = b->core.pos;
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].mate_start_position = b->core.mpos;
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].current_index = arr_ind;
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].tlen = abs(b->core.isize);
            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].read_name = strdup(bam_get_qname(b));

            khiter_t iter = kh_get(khStrInt, breakpoint_pairs_hash, bam_get_qname(b));
            if (iter == kh_end(breakpoint_pairs_hash)) {
                // first time
                //
                iter = kh_put(khStrInt, breakpoint_pairs_hash, bam_get_qname(b), &absent);
                if (absent) {
                    kh_key(breakpoint_pairs_hash, iter) = strdup(bam_get_qname(b));
                    kh_value(breakpoint_pairs_hash, iter) = arr_ind;
                }
                preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].mate_index = -1;
            } else {
                preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[arr_ind].mate_index = kh_value(breakpoint_pairs_hash, iter);

                // now need to go the mate data and reset its mate index
                //
                preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].preads_x_bpts[kh_value(breakpoint_pairs_hash, iter)].mate_index = arr_ind;
            }

            preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].size++;

            if (preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].size + 10 >
                    preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx].capacity)
                dynamicPairedReadsCrossBreakpointsPerChrArraySizeIncrease(&preads_x_bpts_arr->preads_x_bpts_per_chr[preads_x_bpts_idx]);
        }
        bam_destroy1(b);
        hts_itr_destroy(hts_itr);

    }
}

void dynamicPairedReadsCrossBreakpointsPerChrArraySizeIncrease(Paired_Reads_Cross_Breakpoints_Per_Chromosome *preads_x_bpts_per_chr) {
    if (preads_x_bpts_per_chr == NULL) return;

    preads_x_bpts_per_chr->capacity += INIT_SIZE;

    if (preads_x_bpts_per_chr->preads_x_bpts) {
        preads_x_bpts_per_chr->preads_x_bpts = realloc(preads_x_bpts_per_chr->preads_x_bpts, 
                        preads_x_bpts_per_chr->capacity *sizeof(Paired_Reads_Cross_Breakpoints));
        failureExit(preads_x_bpts_per_chr->preads_x_bpts, \
                "Paired_Reads_Cross_Breakpoints_Per_Chromosome->preads_x_bpts_per_chr");
    } else {
        fprintf(stderr, "Error: The preads_x_bpts_per_chr.preads_x_bpts is NULL\n");
        exit(EXIT_FAILURE);
    }
}

void dynamicBreakpointStatsPerChrSizeIncrease(Breakpoint_Stats_Per_Chromosome *bp_stats_per_chr) {
    if (bp_stats_per_chr == NULL) return;

    bp_stats_per_chr->capacity += INIT_SIZE;

    if (bp_stats_per_chr->bpts_stats) {
        bp_stats_per_chr->bpts_stats = 
            realloc(bp_stats_per_chr->bpts_stats, bp_stats_per_chr->capacity * sizeof(Breakpoint_Stats));
        failureExit(bp_stats_per_chr->bpts_stats, "Breakpoint_Stats_Per_Chromosome->bpts_stats");
    } else {
        fprintf(stderr, "Error: The Breakpoint_Stats_Per_Chromosome->bpts_stats is NULL\n");
        exit(EXIT_FAILURE);
    }
}
