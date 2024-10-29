<img src="images/BCM-HGSC-Logo.png" width=250>

---

### Table of Contents

- [Description](#description)
- [Build CNVScanner](#build-cnvscanner)
- [How To Use](#how-to-use)
- [Contact Info](#Contact-Info)

---

## Description

CNVScanner is a stand-alone tool to detect copy number variations (CNVs) in both autosomal and sex chromosomes from short-read sequencing data using multiple strategies. The coverage depth is used to determine CNV type after mitigating potential sequencing artifacts using a binning tactic. CNVScanner generates putative CNVs by applying statistical significance tests to bins and combining neighboring bins with similar coverage. Paired-end read information adds an orthogonal layer of support, enhancing the confidence and accuracy of the detected CNVs. Additionally, split reads are crucial for precisely pinpointing the breakpoints of variations, contributing to more accurate mapping and characterization.

CNVScanner has an accuracy, recall and F1 score surpassing 90% when evaluated against the NIST benchmark for deletion and tandem duplications of the HG002 dataset in high confidence regions. Moreover, our CNVScanner's proficiency extends to accurately identifying approximately two-thirds of the total CNV boundaries within 10 base pairs resolution when compared to those of the NIST benchmark truthset, with one-third correctly locating the exact boundary positions. Additionally, CNVScanner demonstrated an 87% (54/62) sensitivity when compared to published results for sample BAB3596 from Baylor Genetics (Liu, P., et al., DOI: 10.1016/j.cell.2017.01.037)

CNVScanner performs  effectively on variants 1 kb and larger for a 30x genome, but this roughly scales with depth of coverage. 

[Back To The Top](#Table-of-Contents)

## Build CNVScanner

See [INSTALL](INSTALL) for complete details. Please download the [release tarballs](https://github.com/phuang7008/CNVScanner/releases) of your choice and build CNVScanner using the following steps:

    wget https://github.com/phuang7008/Scandium/releases/download/CNVScanner_vXXX/CNVScanner-vXXX.tar.gz
    tar zxvf CNVScanner-vXXX.tar.gz 

    build htslib (see htslib install instruction at the htslib website)
    
    autoreconf -i
    ./configure
    make
    make install or make prefix=<your dir choice> install'

[Back To The Top](#Table-of-Contents)

## How To Use

For CNVScanner v1 series, here is an example of the run command: 

    cnvscanner -i input_bam -o output_dir -R reference -e exclude_region -V reference_version -r chromosome_list_to_be_processed -m 3 -T 12 -N sample_name -w genome_equal_window_bedfile -S minimal_CNV_length -B searching_distance_for_breakpoints

For CNVScanner v2 series, here is an example of the run command: 

    cnvscanner -i input_bam -o output_dir -R reference -e exclude_region -V reference_version -r chromosome_list_to_be_processed -m 3 -T 12 -N sample_name -w genome_equal_window_bedfile -S minimal_CNV_length -B searching_distance_for_breakpoints -M lowmappability_bedfile -L GC_below_25%_bedfile -G GC_above_85%_bedfile

the exclude_region: Ns regions; segdup >= 10,000bps; tandom repeats >= 10,000bps

genome_equal_window_bedfile: can be generated using bedtools makewindows

lowmappability_bedfile: can be downloaded from GA4GH (needs to be sorted and merged)

GC_below_25%_bedfile: can be downloaded from GA4GH   (needs to be sorted and merged)

GC_above_85%_bedfile: can be downloaded from GA4GH   (needs to be sorted and merged)

[Back To The Top](#Table-of-Contents)

## Contact Info

If users encounter any issues, please contact pemhuang@gmail.com

[Back To The Top](#Table-of-Contents)
