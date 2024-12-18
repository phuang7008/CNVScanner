/*
 *
 * ===================================================================================
 *
 *        Filename:        terms.h
 *
 *        Description:    define constants/struct to be used in Sequencing Statistics
 *
 *        Version:        1.0
 *        Created:        01/30/2017 04:45:04 PM
 *        Revision:        none 
 *        Compiler:        gcc
 *
 *        Author:            Peiming (Peter) Huang (phuang@bcm.edu)
 *        Company:        Baylor College of Medicine
 *
 *        =====================================================================================
 */

#ifndef TERMS_H
#define TERMS_H

#include <inttypes.h>   // for PRIu32 and PRIu64 
#include <stdbool.h>    // for bool definition
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>     // for getopt() and usleep()
#include <math.h>
#include <zlib.h>
#include <sys/stat.h>
#include <stdio.h>      // for file read and write

// Users defined header files
#include "htslib/include/htslib/sam.h"
#include "htslib/include/htslib/khash.h"

#include "data_structure.h"
#include "coverage_tracking.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION_ "##CNVScanner_v2.0.0"
#define SOURCE_ "CNVScanner_v2.0.0"
#define INIT_SIZE 500000
#define PR_INIT_SIZE 200        // the init size for number of paired reads across a breakpoint
#define DIFF_COV_TO_MERGE 5
#define SMALL_LENGTH_CUTOFF 50
#define EQUAL_BIN_SIZE 20
//#define BREAKPOINT_DISTANCE_CUTOFF 300     // Qiaoyan: use 2 x seq-length = 2 x 150 = 300 on either size
//                                              it is now a user option
#define DISTANCE_CUTOFF 1000     // Dragen uses 1000bp for the junction detection
#define BREAKPOINT_DISTANCE_TO_GROUP 5      // group neighboring breakpoint within 5bp together

// We need to declared the followings as glabal since the program will change these values!!!
// The naming convention for this type of data is CAPTICAL_WORD1_WORD2_WORD3...extern bool EXCLUDED_FILE_PROVIDED;
// this is the file that contains regions of all Ns, Mappability regions, RepeatMaskers etc. in the reference
extern bool EXCLUDED_FILE_PROVIDED;
extern int  khStrStr;
extern int  khStrFloat;
extern int  khIntPrArray;

/**
 * define a structure that holds the file strings from user inputs
 */
typedef struct {
    // General file info
    char * bam_file;
    char * excluded_region_file;    // provide the regions with N/Repeats/Mappability in the reference genome in bed format
    char * output_dir;              // output directory (mandatory)
    char * reference_file;          // reference file name for cram input file
    char * chromosome_bed_file;     // a file contains chromosome ids and regions need to be processed in bed format
    char * equal_size_window_file;  // a bed file contains equal sized bin windows for all chromosomes

    // For whole genome (WGS) related outputs
    char * wgs_cov_file;        // output the whole genome coverage count information
    char * wgs_cov_report;      // output the whole genome coverage summary report    
    char * wgs_binning_file;    // for whole genome raw binned data file

    //misc
    char * reference_version;
    char * low_mappability_file;        // input file for low mappability (GA4GH)
    char * mappability_outfile;         // output mappability file for debugging
    float  mappability_cutoff;          // the value used to filter out low mappability (default 0.0)
    char * gc_lt25pct_file;             // input file that contains regions with the GC% less than 25% (GA4GH)
    char * gc_gt85pct_file;             // input file that contains regions with the GC% greater than 85% (GA4GH)
    char * gc_content_outfile;          // output GC% scale file for debugging
    char * map_gc_details_file;         // output map and gc calculation details
    char * window_details_file;         // output equal size window intersect details
    char * merged_bin_file;             // output the merged binned data from raw binned data
    char * normalized_result_file;
    char * vcf_output_file;             // produce output CNV file in VCF format
    char * segmented_vcf_output_file;   // produce segmented output CNV file in VCF format
    char * simple_vcf_output_file;      // Simple CNV output in TSV format for easy and quick review
    char * simple_segmented_vcf_file;   // Simple segmentated CNV output in TSV format for easy and quick review
    char * log2ratio_output_file;       // for log2ratio output file for segmentation of all chromosomes
    char * sample_name;
    uint16_t breakpoint_distance;       // the distance to a breakpoint for a CNV, default 300
    int16_t min_cnv_length;             // minimal length to pass a CNV; default 1000
    int8_t min_map_quality;
    int8_t min_base_quality;
    unsigned short num_of_threads;
    float percentage;                   // percentage (fraction) of total bam reads will be used for analysis
    bool remove_duplicate;
    bool remove_supplementary_alignments;
    bool excluding_overlapping_bases;   // to exclude overlapping reads/bases for double counting
    bool Write_WGS_cov_fasta;           // need two different flags, one for Capture and one for WGS

    // developer testing options
    bool non_MC_tag_ON;                 // use non_MC_tag approach for overlap base removal even though the bam file has MC_tags
    bool debug_ON;                      // output debug information for developer
} User_Input;

/* here is the data structure for binned data
 */
typedef struct {
    uint32_t start;
    uint32_t end;
    uint32_t length;
    uint32_t index;
    double   ave_coverage;
    double   log2ratio;
    float    gc_pct;
    double   weighted_mappability;
    double   ave_cov_gc_normalized;
    double   weighted_gc_scale;
    double   ave_cov_map_gc_normalized;
} Binned_Data;

typedef struct {
    char * chromosome_id;
    uint32_t * starts;
    uint32_t * ends;
    uint32_t size;
    uint32_t capacity;
    Binned_Data * data;     // the order will be sequentially from the small positions to the large positions
} Binned_Data_Wrapper;

typedef struct {
    double * array;
    uint32_t size;
} DoubleArray;

// calculation mean, stdev, tlen_mean and tlen_stdev in one pass
typedef struct {
    double total_sum;
    double total_bases;
    double sum_of_base_cov_square;

    uint64_t total_tlen_sum;
    uint32_t total_reads;
    uint32_t total_reads_0_tlen;
    uint64_t sum_of_tlen_square;
} OnePassStdev;

// store discovered CNVs
//
typedef struct {
    uint32_t start;
    uint32_t end;
    uint32_t length;
    double   ave_coverage;
    char type;                  // E: equal bin; R: raw bin
} Equal_Window_Bin;

typedef struct {
    uint32_t breakpoint;            // anchor breakpoint
    uint16_t num_of_breakpoints;    // the occurance of this specific breakpoint
    uint16_t num_of_TLEN_ge_1000;   // number of paired reads span more than 1000 bp
    uint32_t *paired_read_starts;   // an array of starts associated with the paired reads across the breakpoint
    uint32_t *paired_read_ends;     // an array of ends associated with the paired reads across the breakpoint
    uint16_t num_of_paired_reads;
    uint16_t capacity;
    char orientation;               // left 'L' or right 'R' breakpoint associated w/ merged CNV per segment
} CNV_Breakpints;

typedef struct {
    uint32_t start;
    uint32_t end;
    bool passed;
    double qual;
    double ave_coverage;
    char cnv_type;                  // L: for deletion, while P for Dup
    bool valid_cnv;                 // if this is a valid cnv (a passed CNV), if it is not, don't output
    uint32_t left_breakpoint;
    uint32_t right_breakpoint;
    uint32_t left_breakpoint_count;
    uint32_t right_breakpoint_count;
    uint32_t last_right_breakpoint_count;
    uint32_t last_num_larger_TLEN_right;
    uint32_t num_larger_TLEN_left;
    uint32_t num_larger_TLEN_right;
    uint32_t imp_PR_start;
    uint32_t imp_PR_end;
    uint16_t num_larger_imp_PR_TLEN;
    uint8_t evidence_count; 
    uint8_t num_merged_CNVs;
    CNV_Breakpints *cnv_breakpoints;    // Used to store breakpoints associated with CNV merging per segment
    uint16_t breakpoint_size;
    uint16_t breakpoint_capacity;
    uint32_t low_mapp_length;
    uint32_t gc_lt25pct_length;
    uint32_t gc_gt85pct_length;
} INNER_CNV;

typedef struct {
    uint32_t equal_bin_start;
    uint32_t equal_bin_end;
    uint32_t raw_bin_start;
    uint32_t raw_bin_end;
    uint32_t length;    // note, the length is not necessarily = end - start, as some bases like Ns regions are removed
    double  ave_coverage;
    char cnv_type;      // L: for deletion, while P for Dup
    bool combined;      // a combined CNV between neighboring CNVs separated by one Equal-bin

    // store merged bins 
    //
    uint32_t size;
    uint32_t capacity;
    Equal_Window_Bin *equal_bin_array;  // list of all equal window bins that are combined

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
    uint16_t num_of_imp_PR_TLEN_1000;   // number of improperly paired-reads with TLEN >= 1000

    INNER_CNV inner_cnv;
} CNV;

typedef struct {
    char *chromosome_id;
    uint32_t chrom_length;
    uint32_t size;
    uint32_t capacity;
    CNV* cnvs;              // CNV array per chromosome
} CNV_Array;                // CNV array of array

// store breakpoint info
// Breakpoint_Array
//                  chr_id_1        chr_id_2        chr_id_3        ...     chr_id_n
//                  bp_per_chr1     bp_per_chr2     bp_per_chr3     ...     bp_per_chr_n
//                      chr_1                       ...                     chr_n
//                      bp_1                        ...                     bp_n
//                          bp_position
//                          chr_id
//                          mate_chr_id
//                          cigar
//                          left_read
//                          read_start
//                          mate_start
//                          tlen
//============================================================================================
//                  read1                                                   read2
//      start1 -------------------->                        start2 <------------------------
//      TLEN = start2 - start1 + 1 + read2_length
//
typedef struct {
    int type;                       //1: soft-clip;   2: hard-clip;
    uint32_t breakpoint_position;
    //char * read_name;                 // samtools' qname: the current query name
    //int32_t current_index;          // the breakpoint index of itself in the breakpoint array
    //int32_t mate_index;             // the breakpoint index of its mate in the breakpoint array
    uint32_t current_read_start;    // samtools' pos:   Alignment position (1-based)
    uint32_t mate_read_start;       // samtools' pnext: Mate's alignment position (1-based)
    int32_t gap_distance_TLEN;      // samtools' tlen:  Template length (insert size)
                                    //                  The TLEN field is positive for the leftmost
                                    //                  segment of the template, negative for the rightmost,
} Breakpoint;

typedef struct {
    uint32_t size;
    uint32_t capacity;
    char * chrom_id;
    Breakpoint *breakpoints;        // there is a bpt_chr_idx associated with this array
} Breakpoint_Array;                 // all breakpoints are assigned to their own chromosome id

/**
 * define Not_Properly_Paired_Reads structure here
 * Note, I only need to record everything once for paired reads
 * Also, I am going to group nearby reads together if they are separated <= 150 bps (sequencing length)
 * Here we are going to handle only 3 cases
 * Case 1: unmapped read while its mate is mapped
 * Case 2: Paired reads mapped to different chromosome
 * Case 3: Paired reads mapped to the same chromosome, but with insertion size >= 1000
 */

/**
 * The schema: how the improperly paired read array works
 *  Not_Properly_Paired_Reads_Array: improperly_paired_reads_array => An array with size=num_of_chromosomes
 *    For each array element (one chr): there are things to be tracked
 *      chrom_id
 *      num_of_groups (simply size of the groups) all groups within a chromosome
 *      seen_paired_read_hash
 *      grouped_improperly_PRs (an array of type Grouped_Not_Properly_Paired_Reads)
 *          For each grouped improperly paired reads:
 *              group_start
 *              group_end
 *              total_paired_reads (within the group)
 *              num_of_pairs_TLEN_ge_1000
 *              num_of_mapped_reads_on_diff_chrom
 *              num_of_mapped_reads_on_same_chrom
 *    
 */
typedef struct {
    //One_Not_Properly_Paired_Reads * one_improperly_paired_reads;
    uint32_t group_start;                           // the first start position of all current group reads
    uint32_t group_start_end;                       // the last start position of all current group reads
    uint32_t group_mate_end;                        // both paired reads on the same chromosome with TLEN >= 1000
                                                    // it is the end of the last mate in the group
    uint32_t total_paired_reads;
    uint32_t num_of_pairs_TLEN_ge_1000;             // paired reads should be on the same chrom with TLEN >= 1000
    khash_t(m32) * mate_ends_hash;                  // hash of mate end positions that have the paired-read TLEN >= 1000
                                                    // the difference of neighboring ends (sorted) should be <= 150 
                                                    // just like the starts (the hash will be converted to array)
    uint32_t num_of_mapped_reads_on_diff_chrom;
    uint32_t num_of_mapped_reads_on_same_chrom;
} Grouped_Not_Properly_Paired_Reads;

/*
 * For single chromosome
 */
typedef struct {
    char * chrom_id;
    Grouped_Not_Properly_Paired_Reads* grouped_improperly_PRs;
    int32_t num_of_groups;                          // signed value, the size of groups, initialize to -1
    int32_t capacity;

    //khash_t(khStrInt) *seen_paired_read_hash;       // names of paired reads which already encountered
} Not_Properly_Paired_Reads_Array;                  // One array element for one chromosome

// khIntStr: the key as 32 bit integer, while the value is the string
KHASH_MAP_INIT_INT(khIntStr, char*)

/**
 * define a khash like structure that has string as key and various structures as values
 * Note: name part of init must be unique for the key, value types.
 * In our case, 33/32/31 (defined in main.c) are arbitrary symbolic names for hashtables
 * that contains string keys and Temp_Coverage_Array* and int values.
 */
KHASH_MAP_INIT_STR(khStrStr, char*)

KHASH_MAP_INIT_STR(khStrStrArray, StringArray*)

KHASH_MAP_INIT_STR(str, Temp_Coverage_Array*)


/**
 * define a coverage statistics structure for general read stats
 */
typedef struct {
    //read stats
    uint64_t total_reads_produced;          // total reads contains in the bam
    uint64_t total_reads_aligned;           // total reads aligned to a genomic region
    uint32_t total_duplicate_reads;         // total number of duplicate reads
    uint64_t total_reads_paired;            // total number of reads with mate pairs (if any)
    uint64_t total_reads_proper_paired;     // total number of reads with mate pairs (if any)
    uint32_t total_chimeric_reads;          // total number of duplicate reads
    uint32_t total_supplementary_reads;     // total number of reads with supplementary flag set
    uint64_t total_paired_reads_with_mapped_mates; // total number of aligned reads which have mapped mates

    uint16_t read_length;                   // the sequenced READ Length, it is taken from => read_buff_in->chunk_of_reads[i]->core.l_qseq
} Read_Coverage_Stats;

typedef struct {
    //base stats
    uint64_t total_genome_bases;            // total number of bases in Genome
    uint32_t total_excluded_bases;                // total number of bases that are N (unknown bases)
    uint32_t total_excluded_bases_on_chrX;        // total number of bases that are N (unknown bases) on X chromosome
    uint32_t total_excluded_bases_on_chrY;        // total number of bases that are N (unknown bases) on Y chromosome
    uint64_t total_mapped_bases;            // total number of mapped bases
    uint64_t total_uniquely_aligned_bases;  // aka. Reads Usable - where "Usable" is uniquely aligned, non-duplicate, on-target reads
    uint64_t total_genome_coverage;         // total number of read bases aligned to the Genome (used to calculate average coverage)
    uint64_t base_quality_20;               // total number of aligned bases with quality >= 20
    uint64_t base_quality_30;               // total number of aligned bases with quality >=30
    uint32_t total_overlapped_bases;        // total number of overlapped bases from pair-end reads

    //misc
    uint32_t wgs_max_coverage;
    uint32_t base_with_wgs_max_coverage;
    uint16_t median_genome_coverage;

    uint32_t genome_cov_histogram[1001];            // coverage histogram array
    khash_t(m32) *genome_base_with_N_coverage;      // here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100
    khash_t(m32) *genome_coverage_for_median;       // Used for calculating the median coverage for Whole Genome.
} WGS_Coverage_Stats;

/**
 * define a structure to store various information, such as coverage histogram etc
 */
typedef struct {
    Read_Coverage_Stats     *read_cov_stats;
    WGS_Coverage_Stats      *wgs_cov_stats;
} Stats_Info;


typedef struct {
    double average_coverage;
    double stdev;
    double median;

    double ave_log2ratio;
    double stdev_log2ratio;
    double zScore_log2_ratio;
    double del_log2ratio;
    double dup_log2ratio;

    uint32_t total_bases_used;      // excluding N-regsions

    double tlen_mean;               // average TLEN from all data point
    double tlen_stdev;
    uint32_t tlen_total_reads_used; // for a read pair, this should be counted only once (use positive TLEN)
    uint32_t total_reads_w_0_tlen;  // this will be counted twice, so the final results should be divided by 2
    double tlen_outlier_cutoff;

    double ninty_nine_percentile;
    double ninty_eight_percentile;
    double zScore;                  // 95%, 97.5%, 99%, 99.95%, 99.99% => zscore: 1.645, 1.96, 2.576, 3.291, 4
    //double zScore_99_p7_pct;        // 99.7% = mean + 3 * stdev
    double outlier_cutoff;
} Simple_Stats;

#endif //TERMS_H
