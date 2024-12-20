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

void processCurrentRecord(User_Input *user_inputs, bam1_t *rec, bam_hdr_t *header, Stats_Info *tmp_stats_info, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Breakpoint_Array *breakpoint_array, OnePassStdev *one_pass_stdev, Not_Properly_Paired_Reads_Array* improperly_paired_reads_array) {
    if(user_inputs->percentage < 1.0) {
        // set random seed and only need to be set ONCE
        //srand((uint32_t)time(NULL));    
        // it is already set at the main.c
        //
        float random_num = (float)rand() / (float)RAND_MAX;
        if(random_num > user_inputs->percentage) return;
        fprintf(stderr, "Random selection (i.e. downsampling) is ON\n");
    }

    tmp_stats_info->read_cov_stats->total_reads_produced++;

    // collect TLEN info
    //
    int32_t tlen = rec->core.isize;         // the right-hand read will have tlen < 0

    // for properly paired reads
    //
    if (one_pass_stdev != NULL && tlen >= 0 && (rec->core.flag & BAM_FPROPER_PAIR)) {
        one_pass_stdev->total_tlen_sum += tlen;
        one_pass_stdev->sum_of_tlen_square += tlen * tlen;

        if (tlen == 0) {
            one_pass_stdev->total_reads_0_tlen++;
        } else {
            one_pass_stdev->total_reads++;
        }
        //fprintf(stderr, "TLEN %"PRId32"\n", tlen);
    }

    // For improperly paired reads on the same chromosome
    /*
    if (one_pass_stdev != NULL && tlen >= 0 && !(rec->core.flag & BAM_FPROPER_PAIR)) {
        if (rec->core.tid == rec->core.mtid) {
            fprintf(stderr, "TLEN %"PRId32"\n", tlen);
        }
    }*/

    // Need to check various 'READ' flags regarding the current read before doing statistics analysis
    // But the order here is quite important,
    // mapped vs unmapped first
    //
    if(rec->core.flag & BAM_FUNMAP)    {            // Read Unmapped
        return;
    }

    if(rec->core.qual < user_inputs->min_map_quality) {
        return;
    }

    if(rec->core.flag & BAM_FQCFAIL){           // Read Fails Vendor Quality Check
        return;
    }

    if(rec->core.flag & BAM_FSECONDARY){        // Read Alignment is not Primary Alignment
        return;
    }
        
    // the total_mapped_bases should contain everything including soft-clipped bases
    //
    tmp_stats_info->wgs_cov_stats->total_mapped_bases += rec->core.l_qseq;    
    tmp_stats_info->read_cov_stats->total_reads_aligned++;

    if(rec->core.flag & BAM_FPAIRED) {            // Read is paired
        tmp_stats_info->read_cov_stats->total_reads_paired++;
        if(!(rec->core.flag & BAM_FMUNMAP)) {    // Read is Paird with Mapped Mate 
            tmp_stats_info->read_cov_stats->total_paired_reads_with_mapped_mates++;
        }
    }

    if(rec->core.flag & BAM_FDUP) {                // Read is a Duplicate (either optical or PCR)
        tmp_stats_info->read_cov_stats->total_duplicate_reads++;
        if(user_inputs->remove_duplicate) return;
    }

    if (rec->core.flag & BAM_FPROPER_PAIR) {    // Read is properly paired
        tmp_stats_info->read_cov_stats->total_reads_proper_paired++;
    } else {
        // improperly paired reads
        //
        tmp_stats_info->read_cov_stats->total_chimeric_reads++;
        if (improperly_paired_reads_array != NULL)
            processImproperlyPairedReads(improperly_paired_reads_array, rec);

        //return;
    }

    /*if (rec->core.flag & BAM_FSUPPLEMENTARY) {
        tmp_stats_info->read_cov_stats->total_supplementary_reads++;
        //if (user_inputs->remove_supplementary_alignments)
        return;
    }*/

    // skip if RNAME is "*"
    //
    if (strcmp("*", header->target_name[rec->core.tid]) == 0)
        return;

    if (tmp_stats_info->read_cov_stats->read_length <= 0)
        tmp_stats_info->read_cov_stats->read_length = rec->core.l_qseq;

    processRecord(user_inputs, tmp_stats_info, rec, chrom_tracking, chrom_index, breakpoint_array);

    //printf("Done read bam for thread %d\n", thread_id);
}

void processRecord(User_Input *user_inputs, Stats_Info *tmp_stats_info, bam1_t *rec, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Breakpoint_Array *bpt_arr) {
    uint32_t i=0;
    uint32_t chrom_len = chrom_tracking->chromosome_lengths[chrom_index];

    // get the quality information
    uint8_t *qual = bam_get_qual(rec);

    // Need to take care of soft-clip here as it will be part of the length, especially the soft-clip after a match string
    //
    uint32_t * cigar = bam_get_cigar(rec);      // get cigar info
    uint32_t pos_r = rec->core.pos;             // position at the reference 0-based (my own debugging finding)
    uint32_t pos_r_end = bam_endpos(rec)-1;     // end position at the reference is 1-based, so change it to 0-based
    uint32_t pos_q = 0;                         // position at the query

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
                    chrom_tracking->coverage[chrom_index][pos_r]++;
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

            // need to store the soft clip breakpoint here
            //
            if (bpt_arr != NULL)
                storeCurrentReadBreakpointInfo(pos_r, rec, bpt_arr, 1);

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

        } else if (cop == BAM_CPAD || cop == BAM_CHARD_CLIP) {
            // according to the https://www.biostars.org/p/4211/
            // BAM_CPAD shouldn't advance pos_r
            // but according to the samtools specs, it shouldn't consume query also
            //
            // May 17th 2021, need to add BAM_CHARD_CLIP in for Supplementary Reads
            // based on the sam.c, the BAM_CHARD_CLIP is handled the same as BAM_CPAD
            // My own Note: hard clip means the part of the sequence that is hard clipped
            // from the query sequence. Therefore, it is not part of query anymore
            // So no need to advance on query sequence anymore
            //

            if (cop == BAM_CHARD_CLIP && bpt_arr != NULL)
                storeCurrentReadBreakpointInfo(pos_r, rec, bpt_arr, 2);

        } else {
            if (user_inputs->debug_ON)
                fprintf(stderr, "bam case not handled: %s\tcigar operation:\t%d\t%d\n", rec->data, cop, cln);
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
                //printf("handle [(end - start + 1) == 0] as no overlapping %s\n", rec->data);

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

void calculateMeanAndStdev(Binned_Data_Wrapper **binned_data_wrapper, Simple_Stats *the_stats, Chromosome_Tracking *chrom_tracking) {

    uint32_t i=0, j=0, total_num=0;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        for (j=0; j<binned_data_wrapper[i]->size; j++) {
            // Need to skip those windows with Ns regions, segdup or tandem repeats (>10kb) etc.
            //
            if (binned_data_wrapper[i]->data[j].length == 0)
                continue;
        
            total_num++;
        }
    }

    DoubleArray *average_coverage_array = calloc(1, sizeof(DoubleArray));
    average_coverage_array->array = calloc(total_num, sizeof(double));
    average_coverage_array->size=0;

    double sum_of_squares=0.0, sum_of_values=0.0;
    uint32_t num_of_filtered=0;

    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        for (j=0; j<binned_data_wrapper[i]->size; j++) {
            if (binned_data_wrapper[i]->data[j].length > 0) {
                sum_of_values += binned_data_wrapper[i]->data[j].ave_coverage;
                sum_of_squares += pow(binned_data_wrapper[i]->data[j].ave_coverage, 2);

                average_coverage_array->array[average_coverage_array->size] = binned_data_wrapper[i]->data[j].ave_coverage;
                average_coverage_array->size++;
            } else {
                num_of_filtered++;
            }
        }
    }

    fprintf(stderr, "total sum of values is %.4f with size of %"PRIu32" with sum of squares is %.4f\n", sum_of_values, average_coverage_array->size, sum_of_squares);
    fprintf(stderr, "number of bins being filtered: %"PRIu32"\n", num_of_filtered);

    // Note the formula from: https://sureshemre.wordpress.com/2012/04/21/how-to-compute-standard-deviation-in-one-pass/
    //
    the_stats->average_coverage = sum_of_values / average_coverage_array->size;
    the_stats->stdev = sqrt((sum_of_squares - average_coverage_array->size * pow(the_stats->average_coverage,2)) / (average_coverage_array->size - 1 ));

    // now calculate the percentile
    //
    double percentile_cutoff = CalcualtePercentile(average_coverage_array, 99);
    the_stats->ninty_nine_percentile = percentile_cutoff;

    percentile_cutoff = CalcualtePercentile(average_coverage_array, 98);
    the_stats->ninty_eight_percentile = percentile_cutoff;

    // now calculate the z-score
    //
    //the_stats->zScore = 1.645 * the_stats->stdev;             // 90% confident inverval for z score
    the_stats->zScore = 1.96 * the_stats->stdev;                // 95% confident inverval for z score
    //the_stats->zScore_99_pct = 2.576 * the_stats->stdev;      // 99% confident interval for z score (99%)
    //the_stats->zScore_99_p7_pct = 3.00 * the_stats->stdev;    // 99% confident interval for z score (99.7%)

    // now calculate the median
    //
    qsort(average_coverage_array->array, average_coverage_array->size, sizeof(double), compare_double);
    int mid_point = (int) average_coverage_array->size/2;
    the_stats->median = average_coverage_array->array[mid_point];
    if (the_stats->median == 0)         // too many 0s, causing median to be 0, which will be NAN by division
        the_stats->median = the_stats->average_coverage;

    // clean up
    //
    free(average_coverage_array->array);
    free(average_coverage_array);
}

void calculateLog2Ratio(Binned_Data_Wrapper **binned_data_wrapper, Simple_Stats *the_stats, Chromosome_Tracking *chrom_tracking, User_Input* user_inputs) {
    // Calculate log2ration stats from autosomes and non-Ns regions only
    // not needed as log2 will transform the dataset into a non-normal distribution
    //
    //double log2ration_sum = 0.0, sum_of_squares = 0.0;
    //int32_t total_equal_window_size = 0;

    uint32_t i=0;
    for (i=0; i<chrom_tracking->number_of_chromosomes; i++) {
        // process one chromosome at a time
        //
        char *filename = calloc(5000, sizeof(char));
        sprintf(filename, "%s%s%s", user_inputs->log2ratio_output_file, chrom_tracking->chromosome_ids[i], ".txt");
        FILE *log2r_fh = fopen(filename, "w");

        // obtain equal window template
        //
        Binned_Data_Wrapper *equal_size_window = calloc(1, sizeof(Binned_Data_Wrapper));
        checkMemoryAllocation(equal_size_window, "Binned_Data_Wrapper *equal_size_window");
        equal_size_window->size     = 0;
        equal_size_window->capacity = INIT_SIZE;
        equal_size_window->data     = calloc(INIT_SIZE, sizeof(Binned_Data));
        equal_size_window->starts   = calloc(INIT_SIZE, sizeof(uint32_t));
        equal_size_window->ends     = calloc(INIT_SIZE, sizeof(uint32_t));

        int total_lines = processFile(chrom_tracking->chromosome_ids[i], user_inputs->equal_size_window_file, equal_size_window);
        uint32_t j=0, k=0, new_start=0, nan_counter=0;
        double small_value = 0.001;

        // The test showed that it is safe when the consecutive 0x coverage run is below 95, 
        // It breaks when using 96, but at 95 all test samples pass! So the cutoff should be 95
        // But, it is possible that for untested samples, the overfitting at 0x might be different
        // so to avoid overfitting around 0x, use a value 50 to give enough cushions for any samples.
        // But for chrX, the cutoff need to be 10
        // 
        int nan_cutoff = 50;
        if (strcmp(chrom_tracking->chromosome_ids[i], "X") == 0 || strcmp(chrom_tracking->chromosome_ids[i], "chrX") == 0 
                || strcmp(chrom_tracking->chromosome_ids[i], "CHRX") == 0)
            nan_cutoff = 10;

        if (strcmp(chrom_tracking->chromosome_ids[i], "Y") == 0 || strcmp(chrom_tracking->chromosome_ids[i], "chrY") == 0 
                || strcmp(chrom_tracking->chromosome_ids[i], "CHRY") == 0)
            nan_cutoff = 10;

        for (j=0; j<total_lines; j++) {
            for (k=new_start; k<binned_data_wrapper[i]->size; k++) {
                // this is to prevent the 'na' in log2 ratio
                // log2 0 is not defined. so to overcome this, I add 0.0001 before doing log2 ratio
                // This is adjusted based on ViZCNV
                // somehow, for the Y chr, it will cause the error: 
                // Error in JointSeg(data_matrix, eta, omega, muk, mi, smu, sepsilon) : 
                //   NA/NaN/Inf in foreign function call (arg 2)
                // Not used: it looks like for equal bin 1000, at 45949th bin needs not be too low like -13.29
                // Not used: instead, if it changed to -1.29, it will pass the error
                // Not used: Later, I found that this situation happens quite often
                // Not used: The best way is to allow same value run through 12,000 times and then switch to a different value
                // However, according to online document, it is because of model overfitting around 0x coverage.
                // also many coverage is between 0.01x - 0.09x, need to reduce it (consecutive runs) as well
                //
                double log2ratio = 0.0;
                if (equal_size_window->data[j].start == binned_data_wrapper[i]->data[k].start) {
                    //    binned_data_wrapper[i]->data[k].ave_coverage = 40.34;

                    if ((int)(binned_data_wrapper[i]->data[k].ave_coverage * 100000) <= 0.05 
                            || (int)(binned_data_wrapper[i]->data[k].ave_coverage * 10) <= 1 ) {
                        // this is because values with 0 is not easy to check
                        //
                        nan_counter++;
                        if (nan_counter > nan_cutoff) continue;
                    }
                        
                    log2ratio = log2((binned_data_wrapper[i]->data[k].ave_coverage / the_stats->median) + small_value);
                    fprintf(log2r_fh, "%s\t%"PRIu32"\t%.2f\t%.2f\n", chrom_tracking->chromosome_ids[i], 
                            binned_data_wrapper[i]->data[k].start, log2ratio, binned_data_wrapper[i]->data[k].ave_coverage);
                    new_start = k+1;
                    if (nan_counter > nan_cutoff) nan_counter = 0;

                    /*if (strcmp(chrom_tracking->chromosome_ids[i], "chrX") == 0 \
                            || strcmp(chrom_tracking->chromosome_ids[i], "chrY") == 0 \
                            || strcmp(chrom_tracking->chromosome_ids[i], "X") == 0 \
                            || strcmp(chrom_tracking->chromosome_ids[i], "Y") == 0 ) {
                        //continue;
                    } else {
                        total_equal_window_size++;
                        log2ration_sum += log2ratio;
                        sum_of_squares += pow(log2ratio, 2);
                    }*/
                } else {
                    nan_counter++;
                    if (nan_counter > nan_cutoff) continue;
                    
                    log2ratio = log2(small_value);
                    fprintf(log2r_fh, "%s\t%"PRIu32"\t%.2f\t%.2f\n", chrom_tracking->chromosome_ids[i], 
                                                                equal_size_window->data[j].start, log2ratio, 0.0);
                }
                break;
            }
        }
        if (equal_size_window->data) free(equal_size_window->data);
        if (equal_size_window->starts) free(equal_size_window->starts);
        if (equal_size_window->ends) free(equal_size_window->ends);
        if (equal_size_window) free(equal_size_window);

        free(filename);
        filename=NULL;
        fclose(log2r_fh);
    }

    // calculate log2ratio stats, but not needed as the log2 function transforms the original data into a non-normal data
    // The formula and the results are correct!!!
    //
    //the_stats->ave_log2ratio = log2ration_sum / (double) total_equal_window_size;
    //the_stats->stdev_log2ratio =
    //    sqrt((sum_of_squares - pow(log2ration_sum,2)/total_equal_window_size) / (total_equal_window_size - 1 ));
    //the_stats->zScore_log2_ratio = 1.645 * the_stats->stdev_log2ratio;   // 90% confident inverval for z score
    //the_stats->del_log2ratio = the_stats->ave_log2ratio - the_stats->zScore_log2_ratio;
    //the_stats->dup_log2ratio = the_stats->ave_log2ratio + the_stats->zScore_log2_ratio;

    // Not needed at this stage, might be needed in the future, who knows
    // The log2ratio should be calculated after the SLGSeg step!!!
    /*
    fprintf(stderr, "Lof2Ratio stats:\n");
    fprintf(stderr, "ave_log2ratio: %.4f\n", the_stats->ave_log2ratio);
    fprintf(stderr, "stdev_log2ratio: %.4f\n", the_stats->stdev_log2ratio);
    fprintf(stderr, "zScore_log2_ratio: %.4f\n", the_stats->zScore_log2_ratio);
    */

    the_stats->del_log2ratio = log2((the_stats->average_coverage - 2.576 * the_stats->stdev) / the_stats->median);
    //the_stats->del_log2ratio = log2((the_stats->average_coverage - 3 * the_stats->stdev) / the_stats->median);
    the_stats->dup_log2ratio = log2((the_stats->average_coverage + 2.576 * the_stats->stdev) / the_stats->median);
    //the_stats->dup_log2ratio = log2((the_stats->average_coverage + 3 * the_stats->stdev) / the_stats->median);

    fprintf(stderr, "del_log2ratio: %.4f\n", the_stats->del_log2ratio);
    fprintf(stderr, "dup_log2ratio: %.4f\n", the_stats->dup_log2ratio);
}

// the following function is used to for double number array qsort()
//
int compare_double(const void * val1, const void * val2) {
    double tmp_val1 = *((double*) val1);
    double tmp_val2 = *((double*) val2);

    if (tmp_val1 == tmp_val2) return 0;
    else if (tmp_val1 < tmp_val2) return -1;
    else return 1;
}

double CalcualtePercentile(DoubleArray *array_in, int percentile) {
    // sort arry first
    //
    qsort(array_in->array, array_in->size, sizeof(double), compare_double);

    // get the index where the 
    //
    uint32_t percentile_rank = (uint32_t) (floor(array_in->size * percentile)/100);
    fprintf(stderr, "The %d percentile index/rank is at %"PRIu32"\n", percentile, percentile_rank);

    return array_in->array[percentile_rank];
}
