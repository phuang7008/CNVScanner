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
#include "htslib/sam.h"
#include "htslib/khash.h"

#include "data_structure.h"
#include "coverage_tracking.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION_ "##WGS_CNV_v1.0.0"
#define SOURCE_ "WGS_CNV_v1.0.0"
#define INIT_SIZE 500000
#define PR_INIT_SIZE 200        // the init size for number of paired reads across a breakpoint
#define DIFF_COV_TO_MERGE 5
#define SMALL_LENGTH_CUTOFF 50
#define EQUAL_BIN_SIZE 20
#define DISTANCE_CUTOFF 300     // Qiaoyan: use 2 x seq-length = 2 x 150 = 300 on either size

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
    char * equal_size_window;       // a bed file contains equal sized bin windows for all chromosomes

    // For whole genome (WGS) related outputs
    char * wgs_cov_file;        // output the whole genome coverage count information
    char * wgs_cov_report;      // output the whole genome coverage summary report    
    char * wgs_binning_file;    // for whole genome raw binned data file

    //misc
    char * reference_version;
    char * mappability_file;            // input mappability file
    char * mappability_outfile;         // output mappability file for debugging
    float  mappability_cutoff;          // the value used to filter out low mappability (default 0.0)
    char * gc_content_file;             // input GC% scale file
    char * gc_content_outfile;          // output GC% scale file for debugging
    char * map_gc_details_file;         // output map and gc calculation details
    char * window_details_file;         // output equal size window intersect details
    char * merged_bin_file;             // output the merged binned data from raw binned data
    char * normalized_result_file;
    char * vcf_output_file;             // produce output CNV file in VCF format
    char * sample_name;
    int16_t equal_bin_size;
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

// calculation stdev in one pass
typedef struct {
    double total_sum;
    double total_bases;
    double sum_of_base_cov_square;
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
    uint16_t num_of_breakpoints;     // the occurance of this specific breakpoint
    uint16_t num_of_TLEN_ge_1000;    // number of paired reads span more than 1000 bp
    uint8_t orientation;            // 0: not picked, 1: left breakpoint, 2: right breakpoint
} CNV_Breakpints;

typedef struct {
    uint32_t equal_bin_start;
    uint32_t equal_bin_end;
    uint32_t raw_bin_start;
    uint32_t raw_bin_end;
    uint32_t length;    // note, the length is not necessarily = end - start, as some bases like Ns regions are removed
    double  ave_coverage;
    char cnv_type;      // L: for deletion, while P for Dup

    // store merged bins 
    //
    uint32_t size;
    uint32_t capacity;
    Equal_Window_Bin *equal_bin_array;  // list of all equal window bins that are combined

    // store related nearby breakpoint info
    //
    CNV_Breakpints *cnv_breakpoints;
    uint16_t cnv_breakpoints_size;
    uint16_t cnv_breakpoints_capacity;
    int16_t left_start_index;            // the signed index is set when there is a left-hand breakpoint (0-index is valid)
    int16_t right_end_index;             // the signed index is set when there is a right-hand breakpoint (0-index is valid)
} CNV;

typedef struct {
    char *chromosome_id;
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

// consolidate all paired reads associated with a breakpoint
//                      breakpoint
//         read1_f -----------/----->     ...    <---------------  read1_r
//  read2_f ------------->              ...        <------------------  read2_r
//                      cross-breakpoint
//
// Here is the detailed storage structure
// Paired_Reads_Across_Breakpoints_Array: An array with size=num_of_chromosomes
//      Each array element has two variable:
//          chrom_id
//          preads_x_per_anchor_bpt_arr_hash
//              preads_x_per_anchor_bpt_hash (chr1)    preads_x_per_anchor_bpt_hash (chr2)
//                  key1: anchor_bpt_pos1
//                  val1: preads_x_per_anchor_bpt_array (chr1): need to initialize for the key1
//                    For paired reads info
//                      size=0
//                      capacity=PR_INIT_SIZE
//                      num_TLEN_ge_1000=0
//                      initialize seen_paired_read_hash
//                      initialize pread_x_a_bpt array to the size of PR_INIT_SIZE
//                          -> pread_x_a_bpt[0]: for current read, need to store the following info:
//                                current_start_position
//                                mate_start_position
//                                tlen
//                                read_name
//                          -> pread_x_a_bpt[1]
//                             ...
//                          -> pread_x_a_bpt[n]
//                    
//                    For breakpoint info:
//                      my_group_size = 0
//                      my_breakpoint_group[6]
//                      current_breakpoint_count = 0
//                           including neighboring breakpoints (within 5 bp distance) -> group them together
//                      num_of_soft_clipping
//                      num_of_hard_clipping
//                                              
//                  key2: anchor_bpt_pos2
//                  val2: ... 
//                  etc.
//
typedef struct {
    uint32_t current_start_position;
    uint32_t mate_start_position;
    uint32_t tlen;
    char * read_name;
} Paired_Reads_Across_A_Breakpoint;             // for a single breakpoint

typedef struct {
    // paired reads info
    //
    uint32_t total_paired_reads;                // number of total paired reads for this anchor breakpoint
    uint32_t size;                              // number of paired reads with tlen >= 1000 for this anchor breakpoint
    uint32_t capacity;
    uint16_t  num_TLEN_ge_1000;                  // number of paired reads with insertion size >= 1000
    khash_t(khStrInt) *seen_paired_read_hash;           // names of paired reads which already encountered
    Paired_Reads_Across_A_Breakpoint *pread_x_a_bpt;    // an array of paired reads with tlen >= 1000 in this anchor breakpoint group

    // breakpoint info
    //
    int my_group_size;
    uint32_t my_breakpoint_group[6];            // array of unique breakpoints (within 5 bp distance) -> group them together
    uint16_t current_breakpoint_count;          // number of breakpoint at this specific position (include those grouped)
    uint16_t num_of_soft_clipping;
    uint16_t num_of_hard_clipping;
} Paired_Reads_Across_Per_Anchor_Breakpoint_Array;      // for breakpoints on one anchor breakpoint

// key is the breakpoint position as uint32_t, while value is the the array of Paired_Reads_Across_A_Breakpoint
// the KHASH_MAP_INIT_INT, the last INT means the key is INT
//
KHASH_MAP_INIT_INT(khIntPrArray, Paired_Reads_Across_Per_Anchor_Breakpoint_Array*)

typedef struct {
    char * chrom_id;
    khash_t(khIntPrArray) *preads_x_per_anchor_bpt_hash;    //  a hashtable (on each chromosome)
                                                            //      key: anchor breakpoint, 
                                                            //      value: Paired_Reads_Across_A_Breakpoint_Per_Anchor_Array
} Paired_Reads_Across_Breakpoints_Array;                    // For all breakpoints on all chromosomes (one array element for one chromosome)
                                                            // The size of the array will be the number of chromosomes

/**
 * define Not_Properly_Paired_Reads structure here
 * Note, I only need to record everything once for paired reads
 * Also, I am going to group nearby reads together if they are separated <= 150 bps (sequencing length)
 */
typedef struct {
    char * chr_id;
    char * mate_chr_id;
    uint32_t start;
    uint32_t mate_start;
    uint32_t num_of_neighboring_pairs;
    uint32_t num_of_pairs_TLEN_ge_1000;
} Not_Properly_Paired_Reads;

/*
 * For single chromosome
 */
typedef struct {
    char * chrom_id;
    Not_Properly_Paired_Reads* improperly_paired_reads;
    uint32_t size;
    uint32_t capacity;

    khash_t(khStrInt) *seen_paired_read_hash;       // names of paired reads which already encountered
    khash_t(khStrInt) *seen_unmapped_read_hash;     // names of unmapped reads which already encountered
} Not_Properly_Paired_Reads_Array;

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

    uint32_t total_bases_used;      // excluding N-regsions, Repeat Masks, Mappability Excluded Regions and Seg-Dup

    double ninty_nine_percentile;
    double ninty_eight_percentile;
    double zScore;                  // 95%, 97.5%, 99%, 99.95%, 99.99% => zscore: 1.645, 1.96, 2.576, 3.291, 4
    //double zScore_99_p7_pct;        // 99.7% = mean + 3 * stdev
    double outlier_cutoff;
} Simple_Stats;

#endif //TERMS_H
