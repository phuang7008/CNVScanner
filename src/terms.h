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
#include <mysql.h>
#include <sys/stat.h>
#include <stdio.h>      // for file read and write

// Users defined header files
#include "htslib/sam.h"

#include "data_structure.h"
#include "coverage_tracking.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION_ "##WGS CNV v1.0.0"
#define INIT_BIN_SIZE 100000
#define DIFF_COV_TO_MERGE 5
#define SMALL_LENGTH_CUTOFF 50

// We need to declared the followings as glabal since the program will change these values!!!
// The naming convention for this type of data is CAPTICAL_WORD1_WORD2_WORD3...extern bool EXCLUDED_FILE_PROVIDED;
// this is the file that contains regions of all Ns, Mappability regions, RepeatMaskers etc. in the reference
extern bool EXCLUDED_FILE_PROVIDED;
extern int  khStrStr;
extern int  khStrFloat;
extern int  khStrStrArray;
extern int  khStrLCG;
extern int  khStrGTP;

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
    int16_t average_coverage;       // the average coverage of current sample, need to be signed for comparison

    // For whole genome (WGS) related outputs
    char * wgs_cov_file;        // output the whole genome coverage count information
    char * wgs_cov_report;      // output the whole genome coverage summary report    
    char * wgs_binning_file;    // for whole genome smoothly binned data file

    //misc
    char * reference_version;
    char * mappability_file;
    char * gc_content_file;
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
    double   stdev;
    double   z_score;
    double   ave_coverage;
    double   mappability;
    double   ave_cov_mappability_normalized;
    double   gc_ratio;
    double   ave_cov_mappability_gc_normalized;
} Binned_Data;

typedef struct {
    char * chromosome_id;
    uint32_t * starts;
    uint32_t * ends;
    uint32_t size;
    uint32_t capacity;
    Binned_Data * data;     // the order will be sequentially from the small positions to the large positions
} Binned_Data_Wrapper;

/*typedef struct {
    char * chrom_id;
    uint32_t * start;       // array of starts
    uint32_t * end;         // array of ends (note: start and end should be paired)
    float * mappability;    // array of mappability
    uint32_t size;          // array size
} Raw_Mappability;
*/

/** define stringArray structure to store the annotation information
 * one for RefSeq, one for CCDS, one for VEGA and one for Gencode, one for miRNA
 */
typedef struct {
    char **theArray;
    uint16_t capacity;
    uint16_t size;
} stringArray;

#include "htslib/khash.h"

// Instantiate a hash map containing integer keys
// m32 means the key is 32 bit integer, while the value is of unsigned int type (ie uint16_t)
//
KHASH_MAP_INIT_INT(m32, uint32_t)
KHASH_MAP_INIT_INT(m16, uint16_t)
KHASH_MAP_INIT_INT(m8, uint16_t)

// khIntStr: the key as 32 bit integer, while the value is the string
KHASH_MAP_INIT_INT(khIntStr, char*)

/**
 * define a khash like structure that has string as key and various structures as values
 * Note: name part of init must be unique for the key, value types.
 * In our case, 33/32/31 (defined in main.c) are arbitrary symbolic names for hashtables
 * that contains string keys and Temp_Coverage_Array* and int values.
 */
KHASH_MAP_INIT_STR(khStrStr, char*)

KHASH_MAP_INIT_STR(khStrStrArray, stringArray*)

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

#endif //TERMS_H
