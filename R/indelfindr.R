#!/usr/bin/env Rscript

#############################################################################################
#
# Script to call indels using refined cigar strings in non-overlapping sliding windows
# genome wide or across the user specified genomic regions and output indel calls in VCF 4.3.
#
#############################################################################################

logo <- r"{
             _   _    ____  U _____ u  _       _____              _   _    ____    ____
    ___     | \ |"|  |  _"\ \| ___"|/ |"|     |" ___|    ___     | \ |"|  |  _"\U |  _"\ u
   |_"_|   <|  \| |>/| | | | |  _|" U | | u  U| |_  u   |_"_|   <|  \| |>/| | | |\| |_) |/
    | |    U| |\  |uU| |_| |\| |___  \| |/__ \|  _|/     | |    U| |\  |uU| |_| |\|  _ <
  U/| |\u   |_| \_|  |____/ u|_____|  |_____| |_|      U/| |\u   |_| \_|  |____/ u|_| \_\
.-,_|___|_,-.||   \\,-.|||_   <<   >>  //  \\  )(\\,-.-,_|___|_,-.||   \\,-.|||_   //   \\_
 \_)-' '-(_/ (_")  (_/(__)_) (__) (__)(_")("_)(__)(_/ \_)-' '-(_/ (_")  (_/(__)_) (__)  (__)

}"

cat("\n \n \n#############################################################################################\n")

cat(logo)

cat("
#############################################################################################\n
\n
Thank you for using INDELfindR:\nA tool for detecting somatic insertion and deletion events\nwritten in R v4.1.0\n
\n
#############################################################################################\n \n")

#############################################################################################
#
# Define Args:
#
#############################################################################################

suppressMessages(library(argparse))

# create parser object
parser <- ArgumentParser()

# define required parameters and options
parser$add_argument("-v", "--verbose", action="store_true", default=F,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
parser$add_argument("-a", "--alignment_file", type="character",required=T,
                    help="Filepath to alignment file in bam format",
                    metavar="filename")
parser$add_argument("-f", "--flanking_region_length", type="integer", default=10,
                    help="Minimum number of `=` or `X` CIGAR operators required to flank a string of nearby `I` `D` operators in order for a call to be made (default = 10).",
                    metavar="number")
parser$add_argument("-b", "--bam_bin_size", type="integer", default=10000000,
                    help="Length of non-overlapping sliding windows to use to extract overlapping reads from bam file and store in memory at a time (default = 10000000). Larger windows require loading more reads into memory at one time, while smaller windows require longer runtime due to evaluating more reads which overlap neighboring sliding windows and removing more duplicate read calls",
                    metavar="number")
parser$add_argument("-l", "--min_indel_length", type="integer", default=3,
                    help="Minimum number of ref or alt allele basepairs for an indel to be called (default = 3)",
                    metavar="number")
parser$add_argument("-nr", "--number_reads", type="integer", default=4,
                    help="Minimum number of supporting reads for a single given indel to be called (default = 4)",
                    metavar="number")
parser$add_argument("-mq", "--mapq_filter", type="integer", default=20,
                    help="Minimum read mapping quality required to evaluate a read (default = 20).",
                    metavar="number")
parser$add_argument("-t", "--target_regions", type="character", default=F,
                    help="File path to .bed file containing regions in which to perform variant calling. Chromosome syntax must match bed file (ex. `Chr1`, `chr1`, or `1`).",
                    metavar="filename")
parser$add_argument("-nc", "--number_cores", type="integer", default=2,
                    help="Number of cores to use for fork/join parallel computation during indel calling and filtering (default = 2).",
                    metavar="number")
parser$add_argument("-p", "--primary_chromosomes", action="store_true", default=T,
                    help="Call indels in primary chromosomes only ignoring ALT contig alignments (chr1-22,X,Y,M only) (default = True)."
)
parser$add_argument("-vaf", "--vaf_filter", type="double", default=0.01,
                    help="Minimum variant allele freqency required for an indel to be reported (default = 0.01).",
                    metavar="number")
parser$add_argument("-dp", "--read_depth_filter", type="double", default=10,
                    help="Minimum indel range read depth required for an indel to be reported (default = 10).",
                    metavar="number")
parser$add_argument("-o", "--outname", default=getwd(),
                    help="Define the output directory path",metavar="filepath")

args <- parser$parse_args()

bamPath <- args$alignment_file
bam_region_bin_size <- args$bam_bin_size
verbose_arg <- args$verbose
flanking_region_length <- args$flanking_region_length
target_regions <- args$target_regions
number_cores <- args$number_cores
primary_chromosomes <- args$primary_chromosomes
min_indel_length <- args$min_indel_length
mapq_threshold <- args$mapq_filter
min_supporting_reads <- args$number_reads
min_vaf <- args$vaf_filter
min_read_depth <- args$read_depth_filter
outname <- args$outname

#############################################################################################
#
# #To run during dev for debugging:
#
#############################################################################################
#
# bamPath <- "/Users/George/indel_detection_tool_project/benchmarking/indelfindr_results/subset_exome_tumor_validate.sorted.bam"
# bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/EGFR_mutations_reference_dwgsim.sorted.bam"
# #bamPath <- "/Users/George/indel_detection_tool_project/benchmarking/indelfindr_results/miss_validate_igv/miss_5.sorted.bam"
# bam_region_bin_size <- 100000000 #dev-on
# verbose_arg <- FALSE
# flanking_region_length <- 10
# target_regions <- F
# target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt"
# #target_regions <- "/Users/George/indel_detection_tool_project/indelfindr/dev_scripts/complex_indel_bam.bed"
# number_cores <- 2 # Make default 2
# primary_chromosomes <- T
# min_indel_length <- 3
# mapq_threshold <- 20
# min_supporting_reads <- 4
# min_vaf <- 0.01
# min_read_depth <- 10
# #zero_based <- F
# outname <- "/Users/George/indel_detection_tool_project/test_output_dir/test_vcf_output"
# Example command line run
# Rscript indelfindr.R -a /Users/George/indel_detection_tool_project/bam_files_for_testing/EGFR_mutations_reference_dwgsim.sorted.bam -t /Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt

#############################################################################################
#
# Load Dependency Packages and Helper Script
#
#############################################################################################

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
functions.name <- file.path(script.basename, "functions.R")
#print(paste("Sourcing",functions.name,"from",script.name))
source(functions.name)

suppressMessages(library(GenomicAlignments))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(inline))
suppressMessages(library(bettermc))
suppressMessages(library(stringr))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

#############################################################################################
#
# Run INDELfindR Wrapper Function
#
#############################################################################################

run_indelfindr(bamPath,bam_region_bin_size,verbose_arg,flanking_region_length,target_regions,number_cores,primary_chromosomes,min_indel_length,mapq_threshold,min_supporting_reads,min_vaf,min_read_depth,outname)
