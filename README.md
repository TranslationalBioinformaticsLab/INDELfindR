# INDELfindR

## About
INDELfindR is an R based command line tool for detecting simple and complex insertion deletion (INDEL) which outputs indel calls in a VCF v.4.3 file compatible with downstream analysis and annotation tools.

## Installation and Setup

#### Requirements

INDELfindR is compatible with Linux, MacOS, and Windows. 

*Note:* INDELfindR does not support parallel processing on Windows. When running on Windows, be sure to use the flag `-nc 1`  to instruct the script to use one core.

Depends: R 4.1.0

We have developed a Docker image to run INDELfindR. To use Docker to setup and run INDELfindR, see the bottom of this README document. 

To install and set up INDELfindR without Docker, follow the instructions below.

INDELfindR has the following R package dependencies:

The following packages can be installed from CRAN:

```
argparse/2.1.3
plyr/1.8.6
tidyverse/1.3.1
inline/0.3.19
bettermc/1.1.2
purrr/0.3.4
dplyr/1.0.8
```

To install the required dependency versions use the devtools R package:

```
if(!require(devtools)) install.packages("devtools")
require(devtools)
install_version("<each_package>", version = "<version_number")
```

Several packages need to be installed via the Bioconductor 3.14 release.

First, install Bioconductor v.3.14

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

Then install the following packages with Bioconductor 3.14:

BiocManager::install(c("GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38","bamsignals",GenomicAlignments))

The package versions which are installed via Bioconductor 3.14 in the command above should be as follows:
```
GenomicFeatures/1.46.1
BSgenome.Hsapiens.UCSC.hg38/1.4.4
BSgenome.Hsapiens.UCSC.hg19/1.4.3
bamsignals/1.26.0
GenomicAlignments/1.30.0
```

## Quickstart

INDELfindR accepts indexed BAM files as input. INDELfindR only accepts human DNA sequence data which has been aligned to the hg38 reference genome version.

After installing the INDELfindR R package, download the [indelfindr.R script](https://github.com/TranslationalBioinformaticsLab/INDELfindR/blob/main/R/indelfindr.R) and [functions.R script](https://github.com/TranslationalBioinformaticsLab/INDELfindR/blob/main/R/functions.R) and store them in the same directory. 
You can then run INDELfindR from the command line by navigating to the directory containing the running script and running the following command:

```
Rscript indelfindr.R -a <indexed_bam.bam> (options)
```

## Usage Instructions

```
usage: indelfindr.R [-h] [-v] [-q] -a filename [-f number] [-b number]
                    [-l number] [-nr number] [-mq number] [-t filename]
                    [-nc number] [-p] [-vaf number] [-dp number] [-o filepath]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print extra output [default]
  -q, --quietly         Print little output
  -a filename, --alignment_file filename
                        Filepath to alignment file in bam format
  -f number, --flanking_region_length number
                        Minimum number of `=` or `X` CIGAR operators required
                        to flank a string of nearby `I` `D` operators in order
                        for a call to be made (default = 10).
  -b number, --bam_bin_size number
                        Length of non-overlapping sliding windows to use to
                        extract overlapping reads from bam file and store in
                        memory at a time (default = 10000000). Note: Larger windows
                        require loading more reads into memory at one time,
                        while smaller windows require longer runtime due to
                        evaluating more reads which overlap neighboring
                        sliding windows and removing more duplicate read calls
  -l number, --min_indel_length number
                        Minimum number of ref or alt allele basepairs for an
                        indel to be called (default = 3)
  -nr number, --number_reads number
                        Minimum number of supporting reads for a single given
                        indel to be called (default = 4)
  -mq number, --mapq_filter number
                        Minimum read mapping quality required to evaluate a
                        read (default = 20).
  -t filename, --target_regions filename
                        File path to .bed file containing regions in which to
                        perform variant calling. Chromosome syntax of 
                        bed file must match the aligned reads (ex. `Chr1`, `chr1`, or `1`).
  -nc number, --number_cores number
                        Number of cores to use for fork/join parallel
                        computation during indel calling and filtering
                        (default = 2).
  -p, --primary_chromosomes
                        Call indels in primary chromosomes only ignoring ALT
                        contig alignments (chr1-22,X,Y,M only) (default =
                        True).
  -vaf number, --vaf_filter number
                        Minimum variant allele freqency required for an indel
                        to be reported (default = 0.01).
  -dp number, --read_depth_filter number
                        Minimum indel range read depth required for an indel
                        to be reported (default = 10).
  -o filepath, --outname filepath
                        Define the output directory path
  -gv hgversion, --genome_version hgversion
                        Specify the genome version used to align your bam file (hg19 or hg38 only - default hg38)
```

#### Run INDELfindR with Demo Data
We provide a test bam file for running with your INDELfindR setup. We also provide two output files containing the expected output for you to compare your test run results against.

To run INDELfindR with the demo data we provide, download both the bam file, named EGFR_mutations_indelfindr_demo.bam and the .bai index file, named EGFR_mutations_indelfindr_demo.bam.bai, from the /demo/ directory of our INDELfindR Github repository. Make sure the .bai index file is downloaded after downloading the bam file since the program performs an internal check to make sure that the index file creation date is after the creation data of the bam file. 

Once you've installed the indelfindr R package following the above instructions, you can download the INDELfindR executable Rscript stored in the /INDELfindR_script/ directory of the INDELfindR Github repository, and run it on the command line. To run the demo, run the following command within the directory containing the indelfindr.R Rscript:

```
Rscript indelfindr -a /<your_path_to_file>/EGFR_mutations_indelfindr_demo.bam -t /<your_path_to_file>/EGFR_regions_of_interest.txt -o EGFR_mutations_indelfindr_demo

```

The output files for this analysis that we provide for your comparison are stored inside the same /demo/ directory and are named EGFR_mutations_indelfindr_demo.vcf and
EGFR_mutations_indelfindr_demo.indel.cigars.csv.

If the results you recieve are the same as those in our provided demo results files, you are ready to run your analysis!

## Running INDELfindR with ANNOVAR Variant Annotation

We provide a shell script template that you may use to call indels with INDELfindR and annotate the resulting VCF output file with ANNOVAR databases on Unix systems. You may download this shell script template [here](https://github.com/TranslationalBioinformaticsLab/INDELfindR/blob/main/demo/indelfindr_example_shell_script.sh).

## ANNOVAR Setup for INDELfindR

For users' convenience we've outlined a quickstart tutorial to set up the ANNOVAR databases used in our provided INDELfindR with ANNOVAR shell script.

### Set up COSMIC, ClinVar, and dpsnp databases:

Install ANNOVAR and prepare ANNOVAR databases using this tutorial:
https://annovar.openbioinformatics.org/en/latest/user-guide/startup/

### Set Up Some Suggested Annovar Databases

Note: avsnp150 (below) and is a left aligned dbsnp database for use with left aligned vcf. 

avsnp150  was chosen from ANNOVAR, instead of downloading dbSNP build 155 and preparing a database from scratch, since it includes annovar indexes for faster annotation.

```bash

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20210501 humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

```

### Custom Database Build: Current COSMIC v94 database

Since COSMIC70 is the latest version available through Annovar directly, set up v94 database.

See this ink for complete instructions: https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#cosmic-annotations

I've included my commands below but I would use the link above for complete instructions.

#### Get Download Authorization Code
echo "youremail@institution.edu:your_email_password" | base64

(this produces a base64 code which you can copy and paste into the commands below)

#### Download coding variants

curl -H "Authorization: Basic Z2VvcmdlX3RvbGxlZnNvbkBicm93bi5lZHU6TmV3cG9ydDk0IQo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/CosmicCodingMuts.vcf.gz

curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v94/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1633980713&Signature=2PL9NJwDJpqbZflB2zeo1V%2BNDrU%3D" > CosmicCodingMuts.vcf

#### Download non-coding variants
curl -H "Authorization: Basic Z2VvcmdlX3RvbGxlZnNvbkBicm93bi5lZHU6TmV3cG9ydDk0IQo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/CosmicNonCodingVariants.vcf.gz

curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v94/VCF/CosmicNonCodingVariants.vcf.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1633981220&Signature=X8JpUPgshA3PtIesTXZejadaNEw%3D" > CosmicNonCodingMuts.vcf.gz

#### Download mutant export
curl -H "Authorization: Basic Z2VvcmdlX3RvbGxlZnNvbkBicm93bi5lZHU6TmV3cG9ydDk0IQo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/CosmicMutantExport.tsv.gz

curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v94/CosmicMutantExport.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1633981344&Signature=RtjmsnVfZreiHz7f5nPA%2BsW13%2Fo%3D" > CosmicMutantExport.tsv.gz

#### Download noncoding mutant export

curl -H "Authorization: Basic Z2VvcmdlX3RvbGxlZnNvbkBicm93bi5lZHU6TmV3cG9ydDk0IQo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/CosmicNCV.tsv.gz

curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v94/CosmicNCV.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1634050162&Signature=bfa1NUzyHH1uzprzRDiE1AUKthk%3D" > CosmicNCV.tsv.gz

#### Prepare COSMIC database for ANNOVAR

```bash

perl ./prepare_annovar_user.pl -dbtype cosmic /path_to_your_directory/cosmic_database/CosmicMutantExport.tsv -vcf /path_to_your_directory/cosmic_database/CosmicCodingMuts.vcf.gz > ./hg38_cosmic94.txt

```

###### Coding
```
perl ./prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg38_cosmic94_coding.txt 
```
##### Noncoding
```
perl ./prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingMuts.vcf > hg38_cosmic94_noncoding.txt
```


## Running INDELfindR with Docker

### Installing Docker
You can run INDELfindR using the Docker images without installing INDELfindR, R, and it's dependencies. 

To run INDELfindR from a Docker image, first [install Docker](https://docs.docker.com/install/).

Once Docker is installed, double-click the Docker.app in the Applications folder to start Docker Desktop. A whale icon will appear in the top status bar to indicate that Docker is running and accessible from the terminal. You can quit Docker once you are finished using INDELfindR by clicking the Docker whale icon in the top status bar and clicking "Quit Docker Desktop."

### Running INDELfindR with Docker

Once Docker is running and you may want to increase the amount of resources that are available to Docker if you are running on a local computer.
To do this, click on the Docker whale icon, navigate to Preferences > Resources (Advanced) and increase memory to at least 4GB and CPUs to 4.

You can run INDELfindR by running the Docker commands below in the Mac/Linux terminal or Windows PowerShell. All INDELfindR command line arguments may be specified after the  call to `indelfindr` in the following `docker run` commands.

Example run on Mac or Linux:

```
docker run --cpuset-cpus="0-1" -v "$PWD":/data indelfindr -a /data/test.sorted.bam -o test.sorted -b 100000000 (options)
```

Example run on Windows:
*Reminder*: Since parallel processing is not supported on Windows, --number_cores is set to 1 . (In `docker run`, `--cpuset-cpus="0"` specifies  the first core - core "0")
```
docker run --cpuset-cpus="0" -v ${pwd}:/data indelfindr -a /data/test.sorted.bam -o test.sorted -nc 1 -b 100000000 (options)
```

### Running INDELfindR with Singularity on Shared Servers

Build a Singlularity image using the INDELfindR image from DockerHub:

```
singularity build indelfindr_test.simg docker://uzunlab/indelfindr:1.0
```

You may need to bind the data directory paths to your Singularity path to be able to read in and write data to your filesystem. You need to replace the paths to the data directories to reflect your own file system paths.

```
export SINGULARITY_BINDPATH="/gpfs/scratch,/gpfs/data" singularity shell indelfindr_test
```

Run INDELfindR Singlularity Image:

```
vcf_prefix="sample_1"
dir="/<data_directory_path>/"

singularity run indelfindr_test.simg -a $dir/bams/${vcf_prefix}.sorted.bam --bam_bin_size 100000000 --outname $dir/variant_calling/${vcf_prefix}.indelfindr.singularity.test -nc 1
```



