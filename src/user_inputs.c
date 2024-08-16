/*
 * =====================================================================================
 *
 *        Filename:        user_inputs.c
 *
 *        Description:    The implementation file for the user inputs
 *
 *      Version:        1.0
 *      Created:        03/16/2020 04:45:04 PM
 *      Revision:        none
 *      Compiler:        gcc
 *
 *      Author:            Peiming (Peter) Huang, phuang@bcm.edu
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>        // for file access() and getopt()
#include <dirent.h>        // for checking output directory
#include <libgen.h>        // for function basename()
#include "terms.h"
#include "user_inputs.h"

#include "utility.h"

// define (initialize) global variables declared at the terms.h file
//
bool EXCLUDED_FILE_PROVIDED = false;

// print out the help information to show the general usage of the package
//
void usage() {
    printf("Version %s\n\n", VERSION_ );
    printf("Usage:  scandium -input_bam bam/cram -output_dir output_directory [options ...]\n");
    printf("Note:   this is a multi-threading program. Each thread needs 4Gb of memory. So please allocate them accordingly!\n");
    printf("\tfor example: 3 threads would use 12Gb of memory, while 4 threads would need 16Gb of memory, etc.\n\n");
    printf("Mandatory:\n");
    printf("--input_bam            -i  BAM/CRAM alignment file (multiple files are not allowed!).\n");
    printf("                           It Is Mandatory\n");
    printf("--output_dir           -o  output directory. It Is Mandatory\n");
    printf("--sample_name          -N  the sample name to be processed. It Is Mandatory\n");
    printf("--breakpoint_distance  -B  the distance checking to a breakpoint for a CNV. Default: 300\n");
    printf("                           It is based on sequencing length: 2x150=300. The next checking is 3x150=450 or any other user choice\n");
    printf("--equal_size_window    -w  the equal size window bed file. It Is Mandatory\n");
    printf("--reference            -R  the file path of the reference sequence. \n");
    printf("                           It is Mandatory for CRAM files\n\n");

    printf("The Followings Are Optional:\n");
    printf("--low_mappability_file -M  a file contains genomic regions with low mappability - sorted/merged from GA4GH.\n");
    printf("--gc_lt25pct_file      -G  a file contains genomic regions with GC%% less than 25%% - sorted/merged from GA4GH.\n");
    printf("--gc_gt85pct_file      -G  a file contains genomic regions with GC%% greater than 85%% - sorted/merged from GA4GH.\n");
    printf("--min_base_qual        -b  minimal base quality\n");
    printf("                           to filter out any bases with base quality less than b. Default 0\n");
    printf("--min_map_qual         -m  minimal mapping quality\n");
    printf("                           to filter out any reads with mapping quality less than m. Default 0\n");
    printf("--excluded_regions     -e  file name that contains regions to be excluded in bed format\n");
    printf("--percentage           -p  the percentage (fraction) of reads used for this analysis. Default 1.0 (ie, 100%%)\n");
    printf("--chr_list             -r  file name that contains chromosomes and their regions \n");
    printf("                           need to be processed in bed format. Default: Not Provided\n");
    printf("--threads              -T  the number of threads \n");
    printf("                           (Note: when used with HPC's msub, make sure that the number of\n"); 
    printf("                           processors:ppn matches to number of threads). Default 2\n");
    printf("--min_cnv_length       -S  the minimal length of CNV required. Default 1000\n");
    printf("--mappability_cutoff   -c  the minimal mappability used to filter out low mappability regions. Default 0.0\n");
    printf("--ref_version          -V  the reference version used. Default hg38\n");

    printf("The Followings Are Flags\n");
    printf("--duplicate            -d  Specify this flag only when you want to keep Duplicates reads.\n");
    printf("                           Default: Remove Duplicate is ON\n");
    printf("--supplemental         -s  Remove Supplementary alignments and DO NOT use them for statistics. Default: off\n");
    printf("--overlap              -O  Turn off Remove Overlapping Bases to avoid double counting. Default: on\n");
    printf("--wgs_depth            -W  Write/Dump the WGS base coverage depth into Coverage.fasta file \n");
    printf("                           (both -w and -W needed). Default: off\n");
    printf("--help                 -h  Print this help/usage message\n");
}

// Get command line arguments in and check the sanity of user inputs 
//
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]) {
    int arg;
    bool input_error_flag=false;
    bool flag_float=true;

    // Flag set by '--verbose'
    //static int verbose_flag;

    while (1) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose",  no_argument,  &verbose_flag,  9},
            //{"brief",    no_argument,  &verbose_flag,  0},
            /* These options don't set a flag. We distinguish them by their indices. */
            {"chr_list",             required_argument,  0,  'r'},
            {"input_bam",            required_argument,  0,  'i'},
            {"min_base_qual",        required_argument,  0,  'b'},
            {"min_map_qual",         required_argument,  0,  'm'},
            {"output_dir",           required_argument,  0,  'o'},
            {"reference",            required_argument,  0,  'R'},
            {"sample_name",          required_argument,  0,  'N'},
            {"breakpoint_distance",  required_argument,  0,  'B'},
            {"low_mappability_file", required_argument,  0,  'M'},
            {"gc_lt25pct_file",      required_argument,  0,  'L'},
            {"gc_gt85pct_file",      required_argument,  0,  'G'},
            {"ref_version",          required_argument,  0,  'V'},
            {"percentage",           required_argument,  0,  'p'},
            {"excluded_regions",     required_argument,  0,  'e'},
            {"min_cnv_length",       required_argument,  0,  'S'},
            {"equal_size_window",    required_argument,  0,  'w'},
            //{"mappability_cutoff",  required_argument,  0,  'c'},
            {"threads",              required_argument,  0,  'T'},
            {"duplicate",            no_argument,  0,  'd'},
            {"help",                 no_argument,  0,  'h'},
            {"overlap",              no_argument,  0,  'O'},
            {"supplemental",         no_argument,  0,  's'},
            {"wgs_depth",            no_argument,  0,  'W'},
            {"debug",                no_argument,  0,  'g'},
            {0,  0,  0,  0},
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        arg = getopt_long_only (argc, argv, "b:B:de:gG:hi:L:m:M:N:o:p:Or:R:s:S:T:V:w:W", long_options, &option_index);
        //arg = getopt_long_only (argc, argv, "b:c:de:gG:hi:m:M:N:o:p:Or:R:s:S:T:V:w:W", long_options, &option_index);

        /* Detect the end of the options. */
        if (arg == -1) break;

        //printf("User options for %c is %s\n", arg, optarg);
        switch(arg) {
            /*case 'a':
                user_inputs->average_coverage = atoi(optarg);
                break;
            */
            case 'b':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered base quality filter score %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->min_base_quality = atoi(optarg);
                break;
            case 'B':
                user_inputs->breakpoint_distance = atoi(optarg);
                break;
            case 'd': 
                user_inputs->remove_duplicate = false;
                break;
            case 'e':
                EXCLUDED_FILE_PROVIDED = true;
                user_inputs->excluded_region_file = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->excluded_region_file, optarg);
                break;
            case 'g':
                user_inputs->debug_ON = true;
                break;
            case 'G':
                user_inputs->gc_gt85pct_file = calloc(strlen(optarg)+1, sizeof(char));
                strcpy(user_inputs->gc_gt85pct_file, optarg);
                break;
            case 'h': usage(); exit(EXIT_FAILURE);
            case 'i':
                user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->bam_file, optarg);
                break;
            case 'L':
                user_inputs->gc_lt25pct_file = calloc(strlen(optarg)+1, sizeof(char));
                strcpy(user_inputs->gc_lt25pct_file, optarg);
                break;
            case 'm':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered map quality filter score %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->min_map_quality = atoi(optarg);
                break;
            case 'M':
                user_inputs->low_mappability_file = calloc(strlen(optarg)+1, sizeof(char));
                strcpy(user_inputs->low_mappability_file, optarg);
                break;
            case 'N':
                user_inputs->sample_name = calloc(strlen(optarg)+1, sizeof(char));
                strcpy(user_inputs->sample_name, optarg);
                break;
            case 'o':
                user_inputs->output_dir = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->output_dir, optarg);
                break;
            case 'O':
                user_inputs->excluding_overlapping_bases = false;
                break;
            case 'p': 
                   flag_float = isFloat(optarg, &(user_inputs->percentage)); 
                if (!flag_float) {
                    fprintf(stderr, "ERROR: Entered percentage value %s is not a float decimal number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'r': 
                user_inputs->chromosome_bed_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->chromosome_bed_file, optarg);
                break;
            case 'R':
                user_inputs->reference_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->reference_file, optarg);
                break;
            case 's': user_inputs->remove_supplementary_alignments = true; break;
            case 'S': user_inputs->min_cnv_length = atoi(optarg); break;
            case 'T':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered number of threads %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->num_of_threads = atoi(optarg);
                break;
            case 'V':
                strcpy(user_inputs->reference_version, optarg);
                break;
            case 'w':
                user_inputs->equal_size_window_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->equal_size_window_file, optarg);
                break;
            case 'W': 
                user_inputs->Write_WGS_cov_fasta = true;
                break;
            case '?':       //"b:de:gG:hi:m:M:o:p:Or:R:s:T:W:"
                if ( optopt == 'b' || optopt == 'e' || optopt == 'i' || optopt == 'm' 
                        || optopt == 'o' || optopt == 'p' || optopt == 'r'
                        || optopt == 'R' || optopt == 's' || optopt == 'T' || optopt == 'V')
                    fprintf(stderr, "ERROR: Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "ERROR: You have entered an Unknown option. See the above error message\n");
                else
                    fprintf (stderr, "ERROR: You have entered an unknown option. See the above error message\n");
                exit(EXIT_FAILURE);
                break;
            default: 
                fprintf(stderr, "ERROR: Non-option argument %c\n", optopt); input_error_flag=true; break;
                exit(EXIT_FAILURE);
        }
    }

    outputUserInputOptions(user_inputs);

    // check the mandatory arguments (will turn this on for the final test/run)
    if (user_inputs->bam_file == NULL) {
        fprintf(stderr, "ERROR: --input_bam (or -i)\toption is mandatory!\n");
        input_error_flag=true;
    }

    if (user_inputs->sample_name == NULL) {
        fprintf(stderr, "ERROR: --sample_name (or -N)\toption is mandatory!\n");
        input_error_flag=true;
    }

    /*if (user_inputs->mappability_file == NULL) {
        fprintf(stderr, "ERROR: --mappability_file (or -M)\toption is mandatory!\n");
        input_error_flag=true;
    }

    if (user_inputs->gc_content_file == NULL) {
        fprintf(stderr, "ERROR: --gc_content_file (or -G)\toption is mandatory!\n");
        input_error_flag=true;
    }*/

    if (user_inputs->equal_size_window_file == NULL) {
        fprintf(stderr, "ERROR: --equal_size_window (or -w)\toption is mandatory!\n");
        input_error_flag=true;
    }

    if (user_inputs->output_dir == NULL) {
        fprintf(stderr, "ERROR: --output_dir (or -o)\toption is mandatory!\n");
        input_error_flag=true;
    } else {
        // check to see if the directory exist!
        DIR* dir = opendir(user_inputs->output_dir);
        if (dir) {
            /* Directory exists */
            closedir(dir);
        } else if (ENOENT == errno) {
            fprintf(stderr, "ERROR: The output directory \n%s\n doesn't exist! \n", user_inputs->output_dir);
            fprintf(stderr, "Please double check the output directory and try again. Thanks!!\n");
            input_error_flag=true;
        } else {
            /* opendir() failed for some other reason, such as permission */
            fprintf(stderr, "ERROR: Can't open the output directory \n%s\n", user_inputs->output_dir);
            fprintf(stderr, "Please check to see if the permission is set correctly. Thanks!\n");
            input_error_flag=true;
        }
    }

    // Need to check out that all files user provided exist before proceeding
    if (user_inputs->bam_file && !checkFile(user_inputs->bam_file)) input_error_flag=true;
    if (EXCLUDED_FILE_PROVIDED && !checkFile(user_inputs->excluded_region_file)) input_error_flag=true;

    if (input_error_flag) {
        //usage();
        fprintf(stderr, "Please use --help (or -h) for all Scandium options\n");
        exit(EXIT_FAILURE);
    }

}

void setupOutputReportFiles(User_Input *user_inputs) {
    // need to get the basename from BAM/CRAM filename
    char *tmp_basename = basename(user_inputs->bam_file);
    if (!tmp_basename || strlen(tmp_basename) == 0) {
        fprintf(stderr, "ERROR: Something went wrong for extracting the basename from the input BAM/CRAM file\n");
        fprintf(stderr, "Please use --help (or -h) for all Scandium options\n");
        exit(EXIT_FAILURE);
    }

    char string_to_add[1000];

    // output WGS coverage summary report
    //
    sprintf(string_to_add, ".WGS_Coverage_Summary_Report.txt");
    createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_report, string_to_add, VERSION_);

    // output the binned data to files, which needs to be multi-threaded
    //
    if (user_inputs->debug_ON) {
        sprintf(string_to_add, ".WGS_binned_data_REPORT_");
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_binning_file, string_to_add);

        // output normalized binned result file at the end of the program. so no needs multi-threading here
        //
        sprintf(string_to_add, ".WGS_normalized_binned_results.txt");
        createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->normalized_result_file, string_to_add, VERSION_);

        // the output info here needs to be multi-threaded. 
        //
        sprintf(string_to_add, ".WGS_equal_window_details_");
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->window_details_file, string_to_add);
    }

    // for whole genome (wgs) file name
    if (user_inputs->Write_WGS_cov_fasta) {
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_file, ".WGS_cov_");
        //printf("Create wgs file name %s\n", user_inputs->wgs_file);
    }

    // output the mappability and gc% detailed calculation results, needs to be multi-threaded!
    //
    if (user_inputs->debug_ON) {
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->map_gc_details_file, ".map_gc_calculation_details_");
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->merged_bin_file, ".merged_bins_");
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->mappability_outfile, ".mappability_details_");
        generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->gc_content_outfile, ".gc_details_");
    }

    // output the final CNVs in VCF file format
    //
    sprintf(string_to_add, ".not_segmented.cnv.vcf");
    generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->vcf_output_file, string_to_add);

    // output the simple CNV vcf file in TSV format
    //
    sprintf(string_to_add, ".not_segmented.cnv.txt");
    generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->simple_vcf_output_file, string_to_add);

    sprintf(string_to_add, ".cnv.txt");
    generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->simple_segmented_vcf_file, string_to_add);

    // output segmented CNV in VCF format
    //
    sprintf(string_to_add, ".cnv.vcf");
    generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->segmented_vcf_output_file, string_to_add);

    // output log2ratio for segmentation
    //
    sprintf(string_to_add, ".log2ratio.");
    generateFileName(user_inputs->output_dir, tmp_basename, &user_inputs->log2ratio_output_file, string_to_add);

    // KEEP the following Please!
    // free the memory
    // However, here is the quote from the basename() official site
    // Both dirname() and basename() return pointers to null-terminated strings. (Do not pass these pointers to free(3))
    //if (tmp_basename) {
    //    free(tmp_basename);
    //    tmp_basename=NULL;
    //}

    // string_to_add is declared at the stack, so no need to free it!
}

void generateFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append) {
    *file_in = calloc(strlen(output_dir)+strlen(base_name)+strlen(string_to_append)+2,  sizeof(char));
    strcpy(*file_in, output_dir);
    strcat(*file_in, "/");
    strcat(*file_in, base_name);
    strcat(*file_in, string_to_append);
}

void outputUserInputOptions(User_Input *user_inputs) {
    fprintf(stderr, "The following are the options you have chosen:\n");
    fprintf(stderr, "\tInput bam/cram file: %s\n", user_inputs->bam_file);
    fprintf(stderr, "\tOutput directory is: %s\n", user_inputs->output_dir);

    if (user_inputs->excluded_region_file)
        fprintf(stderr, "\tThe file that contains all Ns regions is: %s\n", user_inputs->excluded_region_file);

    if (user_inputs->chromosome_bed_file)
        fprintf(stderr, "\tThe file that contains the chromosome IDs and regions to be processed: %s\n", user_inputs->chromosome_bed_file);

    if (user_inputs->reference_file)
        fprintf(stderr, "\tThe reference sequence file is: %s\n", user_inputs->reference_file);

    fprintf(stderr, "\n");

    fprintf(stderr, "\tThe minimum mapping quality is: %d\n", user_inputs->min_map_quality);
    fprintf(stderr, "\tThe minimum base quality is: %d\n", user_inputs->min_base_quality);
    fprintf(stderr, "\tThe distance checking to a breakpoint of a CNV: %d\n", user_inputs->breakpoint_distance);
    fprintf(stderr, "\tThe number of thread used is: %d\n", user_inputs->num_of_threads);
    fprintf(stderr, "\tThe percentage of reads used for analysis is: %.1f%%\n", user_inputs->percentage*100);

    if (user_inputs->Write_WGS_cov_fasta) {
        fprintf(stderr, "\tThe WGS coverage dump at base level is ON (this will create a cov.fasta file\n");
    } else {
        fprintf(stderr, "\tThe WGS coverage dump is OFF\n");
    }

    if (user_inputs->remove_duplicate) {
        fprintf(stderr, "\tRemove duplicate reads is ON\n");
    } else {
        fprintf(stderr, "\tRemove duplicate reads is OFF\n");
    }

    if (user_inputs->excluding_overlapping_bases) {
        fprintf(stderr, "\tExcluding Overlapping bases is ON\n");
    } else {
        fprintf(stderr, "\tExcluding Overlapping bases is Off\n");
    }

    if (user_inputs->remove_supplementary_alignments) {
        fprintf(stderr, "\tRemove supplementary alignments is ON\n");
    } else {
        fprintf(stderr, "\tRemove supplementary alignments is OFF\n");
    }

    if (user_inputs->low_mappability_file)
        fprintf(stderr, "\tThe sorted and merged GA4GH low mappability file used is %s\n", user_inputs->low_mappability_file);
    
    if (user_inputs->gc_lt25pct_file)
        fprintf(stderr, "\tThe sorted and merged GA4GH GC%% less than 25%% file used is %s\n", user_inputs->gc_lt25pct_file);

    if (user_inputs->gc_gt85pct_file)
        fprintf(stderr, "\tThe sorted and merged GA4GH GC%% greater than 25%% file used is %s\n", user_inputs->gc_gt85pct_file);

    if (user_inputs->debug_ON) {
        fprintf(stderr, "\tThe developer debugging is ON\n");
    } else {
        fprintf(stderr, "\tThe developer debugging is OFF\n");
    }

    fprintf(stderr, "User Input Options ===> DONE!\n\n");
    //printf("The  is: %d\n", user_inputs->);

}

User_Input * userInputInit() {
    User_Input * user_inputs = calloc(1, sizeof(User_Input));
    if (!user_inputs) {
        fprintf(stderr, "ERROR: Memory allocation failed in line %d!\n", __LINE__);
        exit(EXIT_FAILURE);
    }

    user_inputs->percentage = 1.0;
    //user_inputs->average_coverage = -1;
    user_inputs->min_map_quality  = 0;
    user_inputs->min_base_quality = 0;
    user_inputs->num_of_threads   = 2;
    user_inputs->min_cnv_length   = 1000;
    user_inputs->mappability_cutoff  = 0.0;
    user_inputs->breakpoint_distance = 300;
    user_inputs->Write_WGS_cov_fasta = false;
    user_inputs->excluding_overlapping_bases = true;
    user_inputs->remove_duplicate = true;
    user_inputs->remove_supplementary_alignments = false;

    user_inputs->bam_file = NULL;
    user_inputs->output_dir = NULL;
    user_inputs->sample_name = NULL;
    user_inputs->reference_file = NULL;
    user_inputs->merged_bin_file  = NULL;
    user_inputs->gc_lt25pct_file  = NULL;
    user_inputs->gc_gt85pct_file  = NULL;
    user_inputs->vcf_output_file  = NULL;
    user_inputs->mappability_outfile  = NULL;
    user_inputs->gc_content_outfile   = NULL;
    user_inputs->chromosome_bed_file  = NULL;
    user_inputs->map_gc_details_file  = NULL;
    user_inputs->window_details_file  = NULL;
    user_inputs->excluded_region_file = NULL;
    user_inputs->low_mappability_file = NULL;
    user_inputs->log2ratio_output_file = NULL;
    user_inputs->simple_vcf_output_file = NULL;
    user_inputs->normalized_result_file = NULL;
    user_inputs->equal_size_window_file = NULL;
    user_inputs->segmented_vcf_output_file = NULL;
    user_inputs->simple_segmented_vcf_file = NULL;

    user_inputs->reference_version = calloc(10, sizeof(char));
    strcpy(user_inputs->reference_version, "hg38");

    // WGS output file
    //
    user_inputs->wgs_cov_file  = NULL;
    user_inputs->wgs_cov_report = NULL;
    user_inputs->wgs_binning_file = NULL;

    // developer option for testing
    user_inputs->non_MC_tag_ON = false;
    user_inputs->debug_ON = false;
    
    return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {

    if (user_inputs->bam_file)
        free(user_inputs->bam_file);

    if (user_inputs->output_dir)
        free(user_inputs->output_dir);

    if (user_inputs->sample_name)
        free(user_inputs->sample_name);

    if (user_inputs->excluded_region_file)
        free(user_inputs->excluded_region_file);

    // whole genome output files clean-up
    //
    if (user_inputs->wgs_cov_file)
        free(user_inputs->wgs_cov_file);

    if (user_inputs->wgs_cov_report)
        free(user_inputs->wgs_cov_report);

    if (user_inputs->wgs_binning_file)
        free(user_inputs->wgs_binning_file);

    if (user_inputs->chromosome_bed_file)
        free(user_inputs->chromosome_bed_file);

    if (user_inputs->reference_file)
        free(user_inputs->reference_file);

    if (user_inputs->reference_version)
        free(user_inputs->reference_version);

    if (user_inputs->merged_bin_file)
        free(user_inputs->merged_bin_file);

    if (user_inputs->low_mappability_file)
        free(user_inputs->low_mappability_file);

    if (user_inputs->mappability_outfile)
        free(user_inputs->mappability_outfile);

    if (user_inputs->gc_lt25pct_file)
        free(user_inputs->gc_lt25pct_file);

    if (user_inputs->gc_gt85pct_file)
        free(user_inputs->gc_gt85pct_file);

    if (user_inputs->gc_content_outfile)
        free(user_inputs->gc_content_outfile);

    if (user_inputs->map_gc_details_file)
        free(user_inputs->map_gc_details_file);

    if (user_inputs->vcf_output_file)
        free(user_inputs->vcf_output_file);

    if (user_inputs->normalized_result_file)
        free(user_inputs->normalized_result_file);

    if (user_inputs->equal_size_window_file)
        free(user_inputs->equal_size_window_file);

    if (user_inputs->window_details_file)
        free(user_inputs->window_details_file);

    if (user_inputs->simple_segmented_vcf_file)
        free(user_inputs->simple_segmented_vcf_file);

    if (user_inputs->log2ratio_output_file)
        free(user_inputs->log2ratio_output_file);

    if (user_inputs->segmented_vcf_output_file)
        free(user_inputs->segmented_vcf_output_file);

    if (user_inputs->simple_vcf_output_file)
        free(user_inputs->simple_vcf_output_file);

    if (user_inputs)
        free(user_inputs);
}
