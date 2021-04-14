/*
 * =====================================================================================
 *
 *       Filename:  stats.c
 *
 *    Description:  the detailed implementation of sequencing statistics
 *
 *        Version:  1.0
 *        Created:  02/22/2017 01:55:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang (phuang@bcm.edu)
 *        Company:  Baylor College of medicine
 *
 * =====================================================================================
 */

#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "utils.h"

void findDebugPoint() {
    printf("for debugging\n");
}

void processBamChunk(User_Input *user_inputs, Stats_Info *tmp_stats_info, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in, Target_Buffer_Status *target_buffer_status, int thread_id, khash_t(khStrInt)* primary_chromosome_hash, int number_of_chromosomes) {
	// it is the flag that is used to indicate if we need to add the khash into the coverage_hash
	bool not_added = true;	
	uint32_t i = 0;
	char cur_chr[150];
	bool same_chr=true;
	khiter_t *iter_in_out = calloc(1, sizeof(khiter_t));
	strcpy(cur_chr, "NOTTHEREALONE");

	for (i=0; i<read_buff_in->size; i++) {
		if(user_inputs->percentage < 1.0) {
			// set random seed and only need to be set ONCE
            //srand((uint32_t)time(NULL));	
			// it is already set at the main.c
			//
            float random_num = (float)rand() / (float)RAND_MAX;
            if(random_num > user_inputs->percentage) continue;
			fprintf(stderr, "Random selection (i.e. downsampling) is ON\n");
        }

		//tmp_stats_info->read_cov_stats->total_reads_produced++;

		// Need to check various 'READ' flags regarding the current read before doing statistics analysis
		// But the order here is quite important,
		// mapped vs unmapped first
		//
        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FUNMAP)	{			// Read Unmapped
		    tmp_stats_info->read_cov_stats->total_reads_produced++;
            continue;
		}

		// among mapped reads, need to check if we are only interested in the primary chromosomes
		//
		if (primary_chromosome_hash != NULL) {
			khiter_t iter_p = kh_get(khStrInt, primary_chromosome_hash, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]);
			if (iter_p == kh_end(primary_chromosome_hash)) {
				// chrom_id is not one of the primary chromosomes, so skip it!
				//
				continue;
			}
		}
		tmp_stats_info->read_cov_stats->total_reads_produced++;

		if(read_buff_in->chunk_of_reads[i]->core.qual < user_inputs->min_map_quality) {
            continue;
        }

        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FQCFAIL){           // Read Fails Vendor Quality Check
			continue;
        }

        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FSECONDARY){        // Read Alignment is not Primary Alignment
            continue;
        }
        
		// the total_mapped_bases should contain everything including soft-clipped bases
		//
		tmp_stats_info->wgs_cov_stats->total_mapped_bases += read_buff_in->chunk_of_reads[i]->core.l_qseq;	
		tmp_stats_info->read_cov_stats->total_reads_aligned++;

		if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FPAIRED) {			// Read is paired
			tmp_stats_info->read_cov_stats->total_reads_paired++;
            if(!(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FMUNMAP)) {	// Read is Paird with Mapped Mate 
				tmp_stats_info->read_cov_stats->total_paired_reads_with_mapped_mates++;
            }
        }

		if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FDUP) {				// Read is a Duplicate (either optical or PCR)
			tmp_stats_info->read_cov_stats->total_duplicate_reads++;
            if(user_inputs->remove_duplicate) continue;
        }

		if (read_buff_in->chunk_of_reads[i]->core.flag & BAM_FPROPER_PAIR) {	// Read is properly paired
			tmp_stats_info->read_cov_stats->total_reads_proper_paired++;
		} else {
			tmp_stats_info->read_cov_stats->total_chimeric_reads++;
            //continue;       // in hg37, it is also has supplementary flag set
		}

		if (read_buff_in->chunk_of_reads[i]->core.flag & BAM_FSUPPLEMENTARY) {
			tmp_stats_info->read_cov_stats->total_supplementary_reads++;
			//if (user_inputs->remove_supplementary_alignments)
			continue;
		}

		// skip if RNAME is "*"
		//
		if (strcmp("*", header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]) == 0)
			continue;

		// check to see if we have changed the chromosome
		if ( strcmp(cur_chr, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]) != 0) {

			// update the last_chr value, as we are going to handle new chromosome 
			strcpy(cur_chr, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]);

			printf("Processing chromosome id %s for thread %d\n", cur_chr, thread_id);

			not_added = true;
			same_chr = false;
		} else {
			if (!same_chr) same_chr=true;
		}

		// Add the instance of khash_t(m32) to the khash_t(str) object (ie coverage_hash) based on the chrom_id
		if (not_added) {
			int absent;
			// Create an instance of khash_t(m32) and initialize it, Note each thread should have its own khash instance
			Temp_Coverage_Array *temp_cov_array = calloc(1, sizeof(Temp_Coverage_Array));
			temp_cov_array->cov_array = calloc(header->target_len[read_buff_in->chunk_of_reads[i]->core.tid]+1, sizeof(uint32_t));
			temp_cov_array->size = header->target_len[read_buff_in->chunk_of_reads[i]->core.tid]+1;

			khiter_t k_iter = kh_put(str, coverage_hash, cur_chr, &absent);    // get key iterator for chrom_id
			if (absent)
				kh_key(coverage_hash, k_iter) = strdup(cur_chr);
			kh_value(coverage_hash, k_iter) = temp_cov_array;
			not_added = false;

			/////////////////Added on Feb 11, 2015/////////////////////////////
        	if (tmp_stats_info->read_cov_stats->read_length <= 0)
				tmp_stats_info->read_cov_stats->read_length = read_buff_in->chunk_of_reads[i]->core.l_qseq;
        	///////////////////////////////////////////////////////////////////
		}

        processRecord(user_inputs, tmp_stats_info, coverage_hash, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid], read_buff_in->chunk_of_reads[i], target_buffer_status, same_chr, iter_in_out, number_of_chromosomes);
    }

	if (iter_in_out != NULL) free(iter_in_out);

    //printf("Done read bam for thread %d\n", thread_id);
}

void processRecord(User_Input *user_inputs, Stats_Info *tmp_stats_info, khash_t(str) *coverage_hash, char * chrom_id, bam1_t *rec, Target_Buffer_Status *target_buffer_status, bool same_chr, khiter_t *iter_in_out, int number_of_chromosomes) {
	uint32_t i=0, chrom_len=0;

	// get the target_buffer_status index,
	// need to do the following no matter target file or Ns regions specified or not
	// because we need chromosome length information here
	//
	int32_t idx = -1;
	for (i=0; i<(uint32_t)number_of_chromosomes; i++) {
		if (strcmp(chrom_id, target_buffer_status[i].chrom_id) == 0) {
			chrom_len = target_buffer_status[i].size;
			idx = i;
			break;
		}
	}

    if (idx == -1) {
        printf("Can't matching chromosome at stats.c with idx = -1 for read %s\n", rec->data);
        return;
    }

	// get key iterator for chrom_id (outer hash table: khash_t(str))
	//
	if (!same_chr) {
		int	ret=0;
		*iter_in_out = kh_put(str, coverage_hash, chrom_id, &ret);
		if (ret)
			fprintf(stderr, "Something went wrong, hash key for %s is not added\n", chrom_id);
	}

	// get the quality information
	uint8_t *qual = bam_get_qual(rec);

	// Need to take care of soft-clip here as it will be part of the length, especially the soft-clip after a match string
	//
    uint32_t * cigar = bam_get_cigar(rec);		// get cigar info
	uint32_t pos_r = rec->core.pos;				// position at the reference 0-based (my own debugging finding)
    uint32_t pos_r_end = bam_endpos(rec)-1;     // end position at the reference is 1-based, so change it to 0-based
	uint32_t pos_q = 0;							// position at the query

    // before proceeding, we need to find out if there is overlapping between pair-end reads
    // here we only handle the first left most reads as we have trouble tracking them through khash table
    // as khash is not thread safe
    // Get its mate info and make sure they are on the same chromosome
    //
    uint32_t m_pos_r      = rec->core.mpos;         // mate alignment start position at the reference 0-based
    uint32_t m_pos_r_end  = 0;                      // mate alignment end position at the reference 0-based
    int32_t  isize = rec->core.isize;               // insertion size (isize or TLEN)
    bool flag_overlap = false;

    if ( isize > 0 && user_inputs->excluding_overlapping_bases) {    // check if users turn on the flag to excluding overlapped bases
        flag_overlap = getOverlapInfo(user_inputs, tmp_stats_info, rec, &m_pos_r_end);

        if (m_pos_r <= pos_r && pos_r_end <= m_pos_r_end) {
            // complete engulfed by the Reverse read, will skip this Forward read
            //
            //      --------------------------------- forward
            //  -------------------------------------------------- reverse
            //
            return;
        }

        //if (flag_overlap) fprintf(stderr, "%s\n", rec->data);
    }

	for (i=0; i<rec->core.n_cigar; ++i) {
        int cop = cigar[i] & BAM_CIGAR_MASK;    // operation
        int cln = cigar[i] >> BAM_CIGAR_SHIFT;  // length
        int j=0;    // need to make it signed as the comparison with unsigned will generated warnings

		/*
			The “M” CIGAR op (BAM_CMATCH) is forbidden in PacBio BAM files. 
			PacBio BAM files use the more explicit ops “X” (BAM_CDIFF) and “=” (BAM_CEQUAL). 
			PacBio software will abort if BAM_CMATCH is found in a CIGAR field.
		*/
		if (cop == BAM_CMATCH || cop == BAM_CEQUAL || cop == BAM_CDIFF) {
            //tmp_stats_info->wgs_cov_stats->total_uniquely_aligned_bases += cln;

			// For matched/mis-matched bases only. Thus, this portion doesn't contain soft-clipped
			//
			for (j=0; j<cln; ++j) {

				//if (pos_r < 0) continue;      // removed! As it is always false
				if (pos_r >= chrom_len) break;

                if ( (!flag_overlap) || (flag_overlap && (pos_r < m_pos_r || (m_pos_r_end > 0 && pos_r > m_pos_r_end)) ) ) {
				    if (qual[pos_q] >= 20) tmp_stats_info->wgs_cov_stats->base_quality_20++;
				    if (qual[pos_q] >= 30) tmp_stats_info->wgs_cov_stats->base_quality_30++;
                }

				if ( (user_inputs->min_base_quality) > 0 && (qual[pos_q] < user_inputs->min_base_quality) ) {	
					// filter based on the MIN_BASE_SCORE
					//
					pos_r++;
					pos_q++;
					continue;
				}

                if ( (!flag_overlap) || (flag_overlap && (pos_r < m_pos_r || (m_pos_r_end > 0 && pos_r > m_pos_r_end)) ) ) {
				    kh_value(coverage_hash, *iter_in_out)->cov_array[pos_r]++;
                    tmp_stats_info->wgs_cov_stats->total_uniquely_aligned_bases++;
                }
				pos_r++;
				pos_q++;
			}

		} else if (cop == BAM_CSOFT_CLIP) {
			// according to the website: https://github.com/lh3/samtools/blob/master/bam_md.c
			// at the BAM_CSOFT_CLIP, we shouldn't increment the pos_r
			// yes, the soft-clip position is not the beginning of the matched reference position
			//
			pos_q += cln;

		} else if (cop == BAM_CINS) {
			// as they are not part of the reference, we will not advance the pos_r
			// this is confirmed by the web:  https://github.com/lh3/samtools/blob/master/bam_md.c
			//
			//tmp_stats_info->wgs_cov_stats->total_uniquely_aligned_bases += cln;
			for (j=0; j<cln; ++j) {
				//if (pos_r < 0) continue;      // Removed! As it is always false!
				if (pos_r >= chrom_len) break; 

				// insertion bases will be counted as aligned bases
				//
				if (qual[pos_q] >= 20) tmp_stats_info->wgs_cov_stats->base_quality_20++;
				if (qual[pos_q] >= 30) tmp_stats_info->wgs_cov_stats->base_quality_30++;
				pos_q++;
			}

		} else if (cop == BAM_CREF_SKIP || cop == BAM_CDEL) {
			// here we only advance the reference position, but not the query position
			//
			pos_r += cln;

		} else if (cop == BAM_CPAD) {
			// according to the https://www.biostars.org/p/4211/
			// BAM_CPAD shouldn't advance pos_r
			// but according to the samtools specs, it shouldn't consume query also
            //

		} else {
            fprintf(stderr, "bam case not handled: %s\n", rec->data);
        }
	}
}

//Note: this function needs to be run unde the critical condition, that is only one thread will run this function
//
void combineThreadResults(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash) {
	khiter_t outer_iter;		// Note: khint_t and khiter_t are the same thing!
	uint32_t i=0;

	/**
		First we need to loop through the coverage_hash to find out how many chromosomes it tracks
		As the keys in hash is not in order, I need to order them here!!!
		The chromosome id list in the chrom_tracking is ordered, so I am going to use it!
	*/

	for (i = 0; i < chrom_tracking->number_tracked; i++) {
		if (chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] > 2) continue;

		// Now go through the coverage_hash map to see if it is currently tracking it
		// if it tracks new chromosome, we need to allocate memory for the coverage array of the new chromosome
		//
		for (outer_iter = kh_begin(coverage_hash); outer_iter != kh_end(coverage_hash); ++outer_iter) {
			if (kh_exist(coverage_hash, outer_iter)) {
				// compare with ordered chromosome ID , if it is matched, then process it!
				//
				if (strcmp(chrom_tracking->chromosome_ids[i], kh_key(coverage_hash, outer_iter)) == 0) {
					if (chrom_tracking->chromosome_status[i] == 0) {
						// Process the current chromosome for the first time!
						//
						chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[i], i);

						// if it goes to the next chromosome, the previous chromosome should be done processing
				        // Thus, we need to update the previous chromosome status
						// Here we will have to make sure that j is signed. Otherwise, you won't get negative value
                        //
                        int32_t j=0;
						for(j=i-1; j >= 0; j--) {
							if (chrom_tracking->chromosome_status[j] == 1)
								chrom_tracking->chromosome_status[j] = 2;
						}

					}

					// update the coverage information for each position at current chromosome
					// used unsigned variable uj as size is unsigned
                    //
                    uint32_t uj = 0;
					for(uj=0; uj<kh_value(coverage_hash, outer_iter)->size; uj++) {
						chrom_tracking->coverage[i][uj] += kh_value(coverage_hash, outer_iter)->cov_array[uj];
					}
				}
			}
		}
	}
}

// For return value: 0 => No overlapping; 1 => Overlapping Yes; 2 => Completely Overlap (need to skip)
//
bool getOverlapInfo(User_Input *user_inputs, Stats_Info *stats_info, bam1_t *rec, uint32_t *m_pos_r_end) {
    if ( user_inputs->excluding_overlapping_bases &&    // check if users turn on the flag to excluding overlapped bases
         (rec->core.tid == rec->core.mtid)              // on the same chromosome
       ) {
        // obtain overlapping info
        // need to find its mate coords
        //
        //int32_t  isize = rec->core.isize;         // insertion size (isize or TLEN)
        uint32_t pos_r = rec->core.pos;             // position at the reference 0-based (my own debugging finding)
        uint32_t pos_r_end = bam_endpos(rec)-1;     // end position at the reference is 1-based, so change it to 0-based
        uint32_t m_pos_r   = rec->core.mpos;        // mate alignment start position at the reference 0-based

        // For overlapping reads
        //
        uint8_t *tag_i = bam_aux_get(rec, "MC");
        if (user_inputs->non_MC_tag_ON) tag_i = NULL;   // developer option for testing to turn non_MC_tag algorithm ON

        if (tag_i == NULL) {
            // the MC tag is not set
            //
            if (rec->core.tid == rec->core.mtid) {          // on the same chromosome
                if ( (pos_r <= m_pos_r) && (m_pos_r <= pos_r_end) ) {
                    // The normal overlapping situation
                    //  pos_r ----------------- pos_r_end
                    //          m_pos_r --------------------- don't know
                    //
                    stats_info->wgs_cov_stats->total_overlapped_bases += pos_r_end - m_pos_r + 1;   // both are 0-based after I changed pos_r_end from 1-based to 0-based
                    return true;
                } else if (pos_r > m_pos_r) {
                    // Special cases
                    //          pos_r ---------------- pos_r_end
                    //  m_pos_r ---------------............. don't know
                    //
                    stats_info->wgs_cov_stats->total_overlapped_bases += pos_r_end - pos_r + 1;   // both are 0-based after I changed pos_r_end from 1-based to 0-based
                    return true;
                } else {
                    // non overlapping on the same chromosome
                    //  pos_r ----------- pos_r_end
                    //                  m_pos_r ----------------- don't know
                    //
                    return false;
                }
            } else {
                return false;
            }
        }

        // the following handles the case where MC tag is set
        //
        char *cigar_str = bam_aux2Z(tag_i);
        char *cigar_end = cigar_str;
        int cigar_length = 0;       // only add cases where the cigar string consumes the reference

        while (*cigar_end) {
            int n = strtol(cigar_str, &cigar_end, 10);
            switch (*cigar_end) {
                case 'D':
                    cigar_length += n;
                    break;
                case 'M':
                    cigar_length += n;
                    break;
                case 'N':
                    cigar_length += n;
                    break;
                case 'H':
                case 'I':
                case 'P':
                case 'S':
                    break;
                case 'X':
                    cigar_length += n;
                    break;
                case '=':
                    cigar_length += n;
                    break;
                case '?':
                    fprintf(stderr, "unknown and value is %d%c\n", n, *cigar_end);
                    break;
                default:
                    fprintf(stderr, "===>Non-option argument %c\n", *cigar_end);
                    break;
            }
            cigar_end++;
            cigar_str = cigar_end;
        }

        *m_pos_r_end = m_pos_r + cigar_length - 1;      // because cigar_length is 1-based, while m_pos_r is 0-based

        // because the uint32_t subtraction will cause over-flow is it is negative
        //
        int32_t diff_starts = pos_r_end - m_pos_r;
        int32_t diff_ends   = *m_pos_r_end - pos_r;

        //if (isize == 2 && pos_r == 47213293) 
        //    printf("matched\n");

        if ( diff_starts >= 0 && diff_ends >= 0 ) {
            // the overlapping algorithm was taken from 
            // https://stackoverflow.com/questions/15726825/find-overlap-between-collinear-lines
            //
            uint32_t start = (pos_r >= m_pos_r)? pos_r : m_pos_r;
            uint32_t end = (pos_r_end <= *m_pos_r_end) ? pos_r_end : *m_pos_r_end;
            if ((end - start + 1) == 0) {
                printf("handle [(end - start + 1) == 0] as no overlapping %s\n", rec->data);

                return false;
            } else if (end < start) {
                printf ("end < start => Special case to be handled later %s\n", rec->data);
            }

            stats_info->wgs_cov_stats->total_overlapped_bases += end - start + 1;   // should add 1 because endpos is 0-based

            return true;
        } else {
            //fprintf(stderr, "No Overlap: %s\n", rec->data);
            return false;
        }
    }

    fprintf(stderr, "Not Handled: %s\n", rec->data);
    return false;
}
