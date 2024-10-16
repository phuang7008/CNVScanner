<img src="images/BCM-HGSC-Logo.png" width=250>

### Table of Contents

- [Description](#description)
- [Installation](#installation)
- [How To Use](#how-to-use)
- [Contact Info](#Contact-Info)

---

## Description

CNVScanner is a software program that was developed at Baylor College of Medicine Human Genome Sequencing Center. It is used to detect copy number variations (CNVs) including deletions (DELs) and duplications (DUPs) in both autosomal and sex chromosomes for individual samples using multiple strategies.

[Back To The Top](README.md)

---

## Installation

#### Requirements

Building CNVScanner requires a few programs and libraries to be present.

Check to ensure the system provides the followings

    GNU make
    Both C and C++ compiler (e.g. gcc/g++ or related modules)

The build requires autotools scripts like,

    autoheader
    autoconf
    autoreconf

The following external libraries are needed for the CNVScanner installation

    libbz2
    libcrypto
    libcurl
    libdl
    liblzma
    libm
    libpthread
    librt
    libssl
    libz

#### Building Configure

Before the installation of CNVScanner, users need to first install 
    [htslib](https://github.com/samtools/htslib)
package

CNVScanner has been tested against htslib from version 1.10 to 1.19.1, but Do Not use v1.12 as there is a bug in v1.12.

Important Note:
htslib must be installed inside cnvscanner directory. Please create a 'htslib' 
sub-directory like cnvscanner/htslib if the 'htslib' sub-directory doesn't exist.
The 'htslib' sub-directory name should not include any version numbers.

Please follow the htslib 'INSTALL' instructions all the way to the final 
step 'make install' such as

    make prefix=<path to cnvscanner/htslib directory> install

Once htslib is successfully installed, users can proceeds to build CNVScanner.
Please download any release package you prefer and decompress the package

    tar zxvf CNVScanner-release-version.tar.gz
    cd CNVScanner-release-version

Run the following to create a brand new configure for your system

    autoreconf -i

Once users have generated the configure file, then do the following

    ./configure
    make
    make install

The './configure' generates Makefiles for the compilation.

The 'make' compiles source code in both src/ and lib/
and generates executables.

The 'make install' command installs the compiled executables into the /usr/local
directory. Users can change the installation location by adding --prefix=DIR
option to the ./configure run or via 'make prefix=DIR install'.

#### NOTES

[Back To The Top](README.md)

## How To Use

Here is an example how to run cnvscanner

    cnvscanner -i input_bam -o output_dir -R reference -e exclude_region -V reference_version -r chromosome_list_to_be_processed -m 3 -T 2 -N sample_name -w genome_equal_window_bedfile -S minimal_CNV_length -B searching_distance_for_breakpoints -M lowmappability_bedfile -L GC_below_25%_bedfile -G GC_above_85%_bedfile

the exclude_region: Ns regions; segdup >= 10,000bps; tandom repeats >= 10,000bps

genome_equal_window_bedfile: can be generated using bedtools makewindows

lowmappability_bedfile: can be downloaded from GA4GH

GC_below_25%_bedfile: can be downloaded from GA4GH

GC_above_85%_bedfile: can be downloaded from GA4GH

[Back To The Top](README.md)

## Contact Info

If users encounter any issues, please contact pemhuang@gmail.com

[Back To The Top](README.md)
