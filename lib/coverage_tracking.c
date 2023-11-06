/*
 * =====================================================================================
 *
 *       Filename:  coverage_tracking.c
 *
 *    Description:  The detailed implementation of the header file
 *
 *        Version:  1.0
 *        Created:  03/25/2021 11:05:36 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dr. Fritz Mehner (mn), mehner@fh-swf.de
 *        Company:  FH Südwestfalen, Iserlohn
 *
 * =====================================================================================
 */

#include "coverage_tracking.h"

void cleanKhashStrInt(khash_t(khStrInt) *hash_to_clean) {
    khint_t k;
    for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
        if (kh_exist(hash_to_clean, k)) {
            // free the key if it is exist
            //
            if (kh_exist(hash_to_clean, k)) free((char*)kh_key(hash_to_clean, k));
        }
    }

    if (hash_to_clean) kh_destroy(khStrInt, hash_to_clean);

}

void cleanKhashInt(khash_t(m32) *hash_to_clean) {
    khint_t k;
    for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k)
        if (kh_exist(hash_to_clean, k))
            kh_del(m32, hash_to_clean, k);

    if (hash_to_clean) kh_destroy(m32, hash_to_clean);
}


void readBufferInit(Read_Buffer *read_buff_in) {
    uint32_t i=0;
    for (i=0; i<read_buff_in->capacity; i++) {
        read_buff_in->chunk_of_reads[i] = bam_init1();
    }
}

void readBufferDestroy(Read_Buffer *read_buff_in) {
    uint32_t i=0;
    for (i=0; i<read_buff_in->capacity;i++) {
        if (read_buff_in->chunk_of_reads[i] != NULL) {
            bam_destroy1(read_buff_in->chunk_of_reads[i]);
            read_buff_in->chunk_of_reads[i]=NULL;
        }
    }

    read_buff_in->size = 0;
}

uint32_t readBam(samFile *sfin, bam_hdr_t *header, Chromosome_Tracking *chrom_tracking, Read_Buffer *read_buff_in) {
    uint32_t record_idx = 0;
    while (record_idx < read_buff_in->capacity && chrom_tracking->more_to_read) {
        if (sam_read1(sfin, header, read_buff_in->chunk_of_reads[record_idx]) < 0) {
            chrom_tracking->more_to_read = false;
            //fprintf(stderr, "Reading Bam has encountered some problem\n");
            break;
        }
        ++record_idx;
    }

    return record_idx;
}

// Note: this function should never be used to update any information regarding the Chromosome_Tracking variable
// It is only used for initialization and dynamically allocate memories!
//
void chromosomeTrackingInit1(Chromosome_Tracking *chrom_tracking, khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header) {
    uint32_t i=0;

    chrom_tracking->coverage = calloc(chrom_tracking->number_of_chromosomes, sizeof(uint32_t*));
    if (!chrom_tracking->coverage) {
        fprintf(stderr, "ERROR: Memory allocation for %s failed\n", "chrom_tracking->coverage");
        exit(EXIT_FAILURE);
    }

    // since the thread finishes unevenly, some chromosome might be processed before its predecessor,
    // hence we need to initialize them here
    //
    for(i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++)
        chrom_tracking->coverage[i] = NULL;

    // the order of the chromosome id need to be in the same order as the bam/cram
    // 
    chrom_tracking->chromosome_ids = calloc(chrom_tracking->number_of_chromosomes, sizeof(char*));
    if (!chrom_tracking->chromosome_ids) {
        fprintf(stderr, "ERROR: Memory allocation for %s failed\n", "chrom_tracking->chromosome_ids");
        exit(EXIT_FAILURE);
    }

    // Now the chromosome lengths
    //
    chrom_tracking->chromosome_lengths = calloc(chrom_tracking->number_of_chromosomes, sizeof(uint32_t));
    if (!chrom_tracking->chromosome_lengths) {
        fprintf(stderr, "ERROR: Memory allocation for %s failed\n", "chrom_tracking->chromosome_lengths");
        exit(EXIT_FAILURE);
    }

    uint32_t j=0;
    for(i=0; i<(uint32_t) header->n_targets; i++) {
        // initialize the id here based on the chromosome ids need to be processed
        //
        khiter_t iter = kh_get(khStrInt, wanted_chromosome_hash, header->target_name[i]);
        if (iter != kh_end(wanted_chromosome_hash)) {
            // now set chromosome ids here
            //
            chrom_tracking->chromosome_ids[j] = calloc(strlen(header->target_name[i])+1, sizeof(char));
            if (chrom_tracking->chromosome_ids[j] == NULL) {
                printf("ERROR: Allocation failed for chromosome %s\n", header->target_name[i]);
                exit(EXIT_FAILURE);
            }
            strcpy(chrom_tracking->chromosome_ids[j], header->target_name[i]);

            // now set chromosome length info
            //
            chrom_tracking->chromosome_lengths[j] = kh_value(wanted_chromosome_hash, iter);

            j++;
        }
    }

    if (j > chrom_tracking->number_of_chromosomes) {
        fprintf(stderr, "ERROR: Number of chromosomes needs to be processed %d is larger than required %d!\n", j, chrom_tracking->number_of_chromosomes);
        exit(EXIT_FAILURE);
    }

    /*
    chrom_tracking->chromosome_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(uint8_t));
    if (!chrom_tracking->chromosome_status) {
        fprintf(stderr, "ERROR: Memory allocation for %s failed\n", "chrom_tracking->chromosome_status");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++)
        chrom_tracking->chromosome_status[i] = 0;
    */

    chrom_tracking->number_tracked = chrom_tracking->number_of_chromosomes;
    chrom_tracking->more_to_read = true;
}

void chromosomeTrackingInit2(khash_t(khStrInt) *wanted_chromosome_hash, Chromosome_Tracking *chrom_tracking, bam_hdr_t *header) {
    // find out how many chromosome we are dealing with
    //
    khiter_t iter;

    for (iter = kh_begin(wanted_chromosome_hash); iter != kh_end(wanted_chromosome_hash); ++iter) {
        if (kh_exist(wanted_chromosome_hash, iter))
            chrom_tracking->number_of_chromosomes++;
    }

    chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, header);
}

void chromosomeTrackingUpdate(Chromosome_Tracking *chrom_tracking, uint32_t chrom_len, int index) {
    // As the 0 position will be empty as the position will be 1-based
    // So I used it to store the index information for quick access
    // Also, I need to add 1 to the size to align with the 1-based position
    //
    chrom_tracking->coverage[index] = calloc(chrom_len + 1, sizeof(uint32_t));
    if (!chrom_tracking->coverage[index]) {
        fprintf(stderr, "ERROR: Memory allocation failed for chromosome_tracking->coverage");
        exit(EXIT_FAILURE);
    }

    //if (chrom_tracking->chromosome_status[index] == 0)
    //    chrom_tracking->chromosome_status[index] = 1;
    chrom_tracking->chromosome_lengths[index] = chrom_len;
}

int32_t locateChromosomeIndexForChromTracking(char *chrom_id, Chromosome_Tracking *chrom_tracking) {
    uint32_t i=0;
    for (i = 0; i < chrom_tracking->number_tracked; i++) {
        if (chrom_tracking->chromosome_ids[i]) {
            if (strcmp(chrom_id, chrom_tracking->chromosome_ids[i]) == 0) {
                return i;
            }
        }
    }

    //fprintf(stderr, "Something is wrong because the chromosome %s couldn't be found\n", chrom_id);
    return -1;
}

void chromosomeTrackingDestroy(Chromosome_Tracking *chrom_tracking) {
    uint32_t i = 0;
    for (i=0; i<chrom_tracking->number_tracked; i++) {
        if (chrom_tracking->chromosome_ids[i])
            free(chrom_tracking->chromosome_ids[i]);

        if (chrom_tracking->coverage[i])
            free(chrom_tracking->coverage[i]);
    }

    free(chrom_tracking->chromosome_ids);
    free(chrom_tracking->chromosome_lengths);
    //free(chrom_tracking->chromosome_status);
    free(chrom_tracking->coverage);
}

void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking, Target_Buffer_Status *target_buffer_status, int value) {
    // First, we need to find the index that is used to track current chromosome chrom_id
    //
    uint32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    uint32_t i=0,j=0, chrom_len=0;

    // here we need to find out the length of the current chromsome to prevent out of bound segmentation error
    // get the index for the target_buffer_status
    //
    for(i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++) {
        if (strcmp(target_buffer_status[i].chrom_id, chrom_id) == 0) {
            chrom_len = target_buffer_status[i].size;

            break;
        }
    }

    for (i=0; i<Ns_info->size; i++) {
        if (strcmp(Ns_info->coords[i].chrom_id, chrom_id) == 0) {
            //printf("%s\t%"PRIu32"\t%"PRIu32"\n", Ns_info->coords[i].chr, Ns_info->coords[i].start, Ns_info->coords[i].end);
            //
            for (j=Ns_info->coords[i].start; j<Ns_info->coords[i].end; j++) {
                if (j>=chrom_len) continue;

                //printf("value of j is %d\n", j);
                chrom_tracking->coverage[idx][j] = value;
            }
        }
    }
    
    //printf("Finished for zero all N zeros\n");
    //printf("\n");
}

void cleanBedInfo(Bed_Info *bed_info) {
    if (bed_info) {
        if (bed_info->coords) {
            uint32_t i;
            for (i=0; i<bed_info->size; i++) {
                if (bed_info->coords[i].chrom_id)
                    free(bed_info->coords[i].chrom_id);
            }
        
            free(bed_info->coords);
        }
        free(bed_info);
    }
}

void TargetBufferStatusInit(Target_Buffer_Status *target_buffer_status, bam_hdr_t *header) {
    int32_t i=0;
    for(i=0; i<header->n_targets; i++) {
        target_buffer_status[i].chrom_id = calloc(strlen(header->target_name[i])+1, sizeof(char));
        strcpy(target_buffer_status[i].chrom_id, header->target_name[i]);
        target_buffer_status[i].size = header->target_len[i];
        target_buffer_status[i].index = -1;
        //target_buffer_status[i].target_status_array = calloc(header->target_len[i], sizeof(uint8_t));
        //target_buffer_status[i].buffer_status_array = calloc(header->target_len[i], sizeof(uint8_t));
    }
}

void TargetBufferStatusInit2(Target_Buffer_Status *target_buffer_status, khash_t(khStrInt)* wanted_chromosome_hash) {
    khiter_t iter;
    uint32_t idx=0;
    for (iter = kh_begin(wanted_chromosome_hash); iter != kh_end(wanted_chromosome_hash); ++iter) {
        if (kh_exist(wanted_chromosome_hash, iter)) {
            target_buffer_status[idx].chrom_id = calloc(strlen(kh_key(wanted_chromosome_hash, iter))+1, sizeof(char));
            strcpy(target_buffer_status[idx].chrom_id, kh_key(wanted_chromosome_hash, iter));
            target_buffer_status[idx].size = kh_value(wanted_chromosome_hash, iter);
            target_buffer_status[idx].index = -1;
            //target_buffer_status[idx].target_status_array = calloc(kh_value(wanted_chromosome_hash, iter), sizeof(uint8_t));
            //target_buffer_status[idx].buffer_status_array = calloc(kh_value(wanted_chromosome_hash, iter), sizeof(uint8_t));
            idx++;
        }
    }
}

void TargetBufferStatusUpdate(Target_Buffer_Status *target_buffer_status, int32_t target_buffer_index) {
    target_buffer_status[target_buffer_index].target_status_array = calloc(target_buffer_status[target_buffer_index].size, sizeof(uint8_t));
    target_buffer_status[target_buffer_index].buffer_status_array = calloc(target_buffer_status[target_buffer_index].size, sizeof(uint8_t));
}

void TargetBufferStatusDestroyCurrentChromosome(Target_Buffer_Status *target_buffer_status, uint32_t number_of_chromosomes, char* chrom_id) {
    uint32_t idx = 0;
    for (idx=0; idx<number_of_chromosomes; idx++) {
        if (strcmp(target_buffer_status[idx].chrom_id, chrom_id) == 0) {
            if (target_buffer_status[idx].target_status_array) {
                free(target_buffer_status[idx].target_status_array);
                target_buffer_status[idx].target_status_array = NULL;
            }

            if (target_buffer_status[idx].buffer_status_array) {
                free(target_buffer_status[idx].buffer_status_array);
                target_buffer_status[idx].buffer_status_array = NULL;
            }
        }
    }
}

void TargetBufferStatusDestroy(Target_Buffer_Status *target_buffer_status, uint32_t number_of_chromosomes) {
    if (target_buffer_status) {
        uint32_t i;
        for (i=0; i<(uint32_t)number_of_chromosomes; i++) {
            if (target_buffer_status[i].chrom_id)
                free(target_buffer_status[i].chrom_id);

            if (target_buffer_status[i].target_status_array)
                free(target_buffer_status[i].target_status_array);

            if (target_buffer_status[i].buffer_status_array)
                free(target_buffer_status[i].buffer_status_array);
        }
        free(target_buffer_status);
    }
}
