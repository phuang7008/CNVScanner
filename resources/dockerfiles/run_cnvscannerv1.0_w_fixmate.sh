#!/bin/bash

# start
echo " - start"
date

# number of args = 5, if not exit
if [ "$#" -ne 5 ]
then
    echo "Usage: run_cnvscannerv1.0.sh sample_name s3_input_CRAM_BAM_url s3_output_dir_url reference_version equal_window_size"
    echo -e "\treference_version : [hg38 or hg37]"
    exit 2
fi

set -euo pipefail
set -x

# set input variables
sample=$1
s3_input_CRAM_BAM_url=$2
s3_output_dir_url=$3
reference_version=$4
equal_window_size=$5

# function that prints the error statement
error() {
  echo "ERROR - ${1}" >&2
  exit 1
}


# function that checks for empty string
check_empty_string()
{
    if [ -z "${1}" ]
    then
        error "${1} s3 url can't be empty - check your input!"
        exit 1
    fi
}

# function that checks an s3 input url which should start with s3://
check_input_s3()
{
    check_empty_string "${1}"

    file_type="${1:0:5}"  # get first 5 chars from the left
    if [ $file_type != "s3://" ]
    then
        error "${1} is not a valid s3 url, expecting url starting with s3://"
    fi
}

# check input
echo " - checking inputs"
check_input_s3 "${s3_input_CRAM_BAM_url}"
check_input_s3 "${s3_output_dir_url}"
check_empty_string "${reference_version}"
check_empty_string "${sample}"
check_empty_string "${equal_window_size}"

# check if s3_output_dir_url ends with '/'
suffix="${s3_output_dir_url:(-1)}"
if [ $suffix != "/" ]
then
    error "${s3_output_dir_url} does not end in '/'"
fi


# set reference
if [ $reference_version == "hg38" ]
then
    s3_reference_file_url="s3://hgsccl-op-data/bwa_references/h/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    s3_reference_index_url="s3://hgsccl-op-data/bwa_references/h/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
    chr_list="s3://hgsccl-op-data/alignstats/regions/GRCh38_full_analysis_set_plus_decoy_hla.bed"
    excluded_regions="s3://hgsccl-op-data/cnv-scanner/hg38_resources/hg38_Ns_regions_seg_dups_gt10k_tandem_rep_gt10k_sorted_merged.bed"
    if [ $equal_window_size -eq 500 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg38_resources/hg38.500.windows.bed"
    elif [ $equal_window_size -eq 1000 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg38_resources/hg38.1000.windows.bed"
    elif [ $equal_window_size -eq 2000 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg38_resources/hg38.2000.windows.bed"
    else
        error "Error - there is no equal window for the size you listed: $equal_window_size"
        exit 1
    fi
elif [ $reference_version == "hg37" ]
then
    s3_reference_file_url="s3://hgsccl-op-data/bwa_references/h/hs37d5/hs37d5.fa"
    s3_reference_index_url="s3://hgsccl-op-data/bwa_references/h/hs37d5/hs37d5.fa.fai"
    chr_list="s3://hgsccl-op-data/alignstats/regions/hg19.bed"
    excluded_regions="s3://hgsccl-op-data/cnv-scanner/hg19_resources/hs37d5_N_regions_w_segdups_gt_10k_tandem_rep_gt_10k_sorted_merged.bed"

    if [ $equal_window_size -eq 500 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg19_resources/hg37.500.windows.bed"
    elif [ $equal_window_size -eq 1000 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg19_resources/hg37.1000.windows.bed"
    elif [ $equal_window_size -eq 2000 ]; then
        equal_window="s3://hgsccl-op-data/cnv-scanner/hg19_resources/hg37.2000.windows.bed"
    else
        error "Error - there is no equal window for the size you listed: $equal_window_size"
        exit 1
    fi
else
    error "ERROR - $reference_version - incorrect reference specified, it can only be hg38 or hg37"
    exit 1
fi

cd /scratch

echo " - copying files from s3 to container"
echo " - copying ${s3_input_CRAM_BAM_url} ..."
aws s3 cp ${s3_input_CRAM_BAM_url} /scratch/ || error "failed to download ${s3_input_CRAM_BAM_url}"
input_file_basename=`basename "${s3_input_CRAM_BAM_url}"`
input_filepath="/scratch/${input_file_basename}"

echo " - copying ${s3_input_CRAM_BAM_url} index file..."
aws s3 cp ${s3_input_CRAM_BAM_url}.crai /scratch/ || error "failed to download ${s3_input_CRAM_BAM_url}.crai"
input_index_basename=`basename "${s3_input_CRAM_BAM_url}.crai"`
input_index_filepath="/scratch/${input_index_basename}"

echo " - copying ${s3_reference_file_url} ..."
aws s3 cp ${s3_reference_file_url} /scratch/ || error "failed to download ${s3_reference_file_url}"
reference_file_basename=`basename "${s3_reference_file_url}"`
reference_filepath="/scratch/${reference_file_basename}"

echo " - copying ${s3_reference_index_url} ..."
aws s3 cp ${s3_reference_index_url} /scratch/ || error "failed to download ${s3_reference_index_url}"
reference_index_basename=`basename "${s3_reference_index_url}"`
reference_index_filepath="/scratch/${reference_index_basename}"

echo " - copying ${chr_list} ..."
aws s3 cp ${chr_list} /scratch/ || error "failed to download ${chr_list}"
chr_list_basename=`basename "${chr_list}"`
echo " - copying ${equal_window} ..."
aws s3 cp ${equal_window} /scratch/ || error "failed to download ${equal_window}"
equal_window_basename=`basename "${equal_window}"`

echo " - copying ${excluded_regions} ..."
aws s3 cp ${excluded_regions} /scratch/ || error "failed to download ${excluded_regions}"
excluded_regions_basename=`basename "${excluded_regions}"`

# Step 1: samtools sort by name with -n option; then fixmate; then sort by coordinates
#
echo " - samtools sort by name with -n option; then fixmate; then sort by coordinates"

samtools sort -@ 8 /scratch/${input_file_basename} --reference /scratch/$reference_file_basename -n -l 9 -m 1G | samtools fixmate -@ 8 --reference ${reference_file_basename} - - | samtools sort -@ 8 --reference ${reference_file_basename} -l 9 -m 1G -O CRAM -o ${sample}_fixmate.sorted.cram -

# exit if the fixmated and sorted cram file doesn't exist
#
if [ ! -f /scratch/${sample}_fixmate.sorted.cram ]; then
    echo "The File ${sample}_fixmate.sorted.cram not found!"
    exit 1
fi

# Step 2: generate index file for the sorted fixmate cram file
#
echo " - samtools index the _fixmate.sorted.cram"
samtools index -c /scratch/${sample}_fixmate.sorted.cram || error "samtools index for ${sample}_fixmate.sorted.cram failed!"

# Step 3: run CNVScanner on autosome and chrX and chrY individually
#
echo " - CNVScanner for autosome only"
grep -v -e X -e Y /scratch/${chr_list_basename} > /scratch/autosome.bed
mkdir -p /scratch/autosomes
echo ""

cnvscanner -i /scratch/${sample}_fixmate.sorted.cram -o /scratch/autosomes -R /scratch/$reference_file_basename -e /scratch/$excluded_regions_basename -V $reference_version -r /scratch/autosome.bed -m 1 -T 2 -w /scratch/$equal_window_basename -N $sample -S $equal_window_size || error "cnvscanner run for the autosomes on ${sample}_fixmate.sorted.cram failed!"

echo " - CNVScanner for chrX only"
grep X /scratch/$chr_list_basename > /scratch/chrX.bed
mkdir -p /scratch/chrX
echo ""

cnvscanner -i /scratch/${sample}_fixmate.sorted.cram -o /scratch/chrX -R /scratch/$reference_file_basename -e /scratch/$excluded_regions_basename -V $reference_version -r /scratch/chrX.bed -m 1 -T 2 -w /scratch/$equal_window_basename -N $sample -S $equal_window_size || error "cnvscanner run for the chrX on ${sample}_fixmate.sorted.cram failed!"

echo " - CNVScanner for chrY only"
grep Y /scratch/$chr_list_basename > /scratch/chrY.bed
mkdir -p /scratch/chrY
echo ""

cnvscanner -i /scratch/${sample}_fixmate.sorted.cram -o /scratch/chrY -R /scratch/$reference_file_basename -e /scratch/$excluded_regions_basename -V $reference_version -r /scratch/chrY.bed -m 1 -T 2 -w /scratch/$equal_window_basename -N $sample -S $equal_window_size || error "cnvscanner run for the chrY on ${sample}_fixmate.sorted.cram failed!"

# now combine autosomes/chrX/chrY results together
#
cat /scratch/autosomes/${sample}_fixmate.sorted.cram.cnv.txt <(grep -v start /scratch/chrX/${sample}_fixmate.sorted.cram.cnv.txt) <(grep -v start /scratch/chrY/${sample}_fixmate.sorted.cram.cnv.txt) > /scratch/${sample}_fixmate.sorted.cram.cnv.txt

bcftools concat /scratch/autosomes/${sample}_fixmate.sorted.cram.cnv.vcf /scratch/chrX/${sample}_fixmate.sorted.cram.cnv.vcf /scratch/chrY/${sample}_fixmate.sorted.cram.cnv.vcf -o /scratch/${sample}_fixmate.sorted.cram.cnv.vcf

bgzip /scratch/${sample}_fixmate.sorted.cram.cnv.vcf
tabix -p vcf /scratch/${sample}_fixmate.sorted.cram.cnv.vcf.gz

mv /scratch/autosomes/${sample}_fixmate.sorted.cram.WGS_Coverage_Summary_Report.txt /scratch/${sample}_autosome_fixmate.sorted.cram.WGS_Coverage_Summary_Report.txt

# clean-up before copying
#
rm -rf /scratch/*bed
rm -rf /scratch/*fa
rm -rf /scratch/*fai
rm -rf /scratch/autosomes
rm -rf /scratch/chrX
rm -rf /scratch/chrY

# remove the older cram/bam files
#
echo " - remove all temp CRAM files"
rm -rf /scratch/${sample}_sorted.cram /scratch/${sample}_fixmate.cram /scratch/$input_file_basename /scratch/*bam

###############################
# copy files to s3
#
echo " - copying output files to s3"

aws s3 cp /scratch/ "${s3_output_dir_url}" --recursive --exclude "*" --include "*cnv.txt" --include "*cnv.vcf.gz*" --include "*.WGS_Coverage_Summary_Report.txt"

# complete
#
date
