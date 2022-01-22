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
Thank you for using indelFindR:\nA tool for detecting complex insertion deletion mutations \nfrom WGS alignment data written in R v4.1.1\n
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

# specify our desired options 
# by default ArgumentParser will add an help option 
# turn off for dev:

parser$add_argument("-v", "--verbose", action="store_true", default=F,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-a", "--alignment_file", type="character",required=T, 
                    help="Filepath to alignment file in bam format",
                    metavar="filename")
parser$add_argument("-f", "--flanking_region_length", type="integer", default=10, 
                    help="Minimum number of `=` or `X` CIGAR operators required to flank a string of nearby `I` `D` operators in order for a call to be made",
                    metavar="number")
parser$add_argument("-b", "--bam_bin_size", type="integer", default=100000, 
                    help="Length of non-overlapping sliding windows to use to extract overlapping reads from bam file and store in memory at a time (default = 100000). Larger windows require loading more reads into memory at one time, while smaller windows require longer runtime due to evaluating more reads which overlap neighboring sliding windows and removing more duplicate read calls",
                    metavar="number")
parser$add_argument("-l", "--min_indel_length", type="integer", default=2, 
                    help="Minumum number of ref or alt allele basepairs for an indel to be called (default = 5 nucleotides)",
                    metavar="number")
parser$add_argument("-nr", "--number_reads", type="integer", default=2, 
                    help="Minumum number of supporting reads for a single given indel to be called",
                    metavar="number")
parser$add_argument("-mq", "--mapq_filter", type="integer", default=20, 
                    help="Minumum read mapping quality required to evaluate a read.",
                    metavar="number")
parser$add_argument("-t", "--target_regions", type="character", default=F, 
                    help="File path to .bed file containing regions in which to perform variant calling. Chromosome syntax must match bed file (ex. `Chr1`, `chr1`, or `1`).",
                    metavar="filename")
parser$add_argument("-nc", "--number_cores", type="integer", default=2, 
                    help="Number of cores to use for fork/join parallel computation during indel calling and filtering.",
                    metavar="number")
parser$add_argument("-p", "--primary_chromosomes", action="store_true", default=T, 
                    help="Call indels in primary chromosomes only ignoring ALT contig alignments (chr1-22,X,Y,M only)."
)
parser$add_argument("-vaf", "--vaf_filter", type="double", default=0.01, 
                    help="Minumum variant allele freqency required for an indel to be reported",
                    metavar="number")
parser$add_argument("-dp", "--read_depth_filter", type="double", default=20, 
                    help="Minumum indel range read depth required for an indel to be reported",
                    metavar="number")
# parser$add_argument("-z", "--zero_based", action="store_true",
#                     help="convert variant calls from zero based to one based")
parser$add_argument("-o", "--outname", default=getwd(),
                    help="Define the output directory path")

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
#zero_based <- args$zero_based
outname <- args$outname

# To run during dev:
# bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/EGFR_mutations_reference_dwgsim.sorted.bam"
# bam_region_bin_size <- 100000 #dev-on
# verbose_arg=FALSE
# flanking_region_length <- 10
# target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt"
# number_cores <- 8 # Make default 2
# primary_chromosomes <- T
# min_indel_length <- 3
# mapq_threshold <- 20
# min_supporting_reads <- 4
# min_vaf <- 0.01
# min_read_depth <- 10
# #zero_based <- F
# outname <- "/Users/George/indel_detection_tool_project/test_output_dir/test_vcf_output"

# bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/bam_for_testing_1656854-1657354.bam"
# bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/EGFR_mutations_reference_dwgsim.sorted.bam"
# #bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/test_read.bam"
# #max_window_length <- 5 #dev-on
# bam_region_bin_size <- 100000 #dev-on
# verbose_arg=FALSE
# #window_load_length <- 500
# flanking_region_length <- 10
# #sliding_window_size <- 20
# #Optional parameters for testing:
# target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt"
# target_regions <- F
# #target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/debug_10_27.txt"
# number_cores <- 10 # Make default 2
# primary_chromosomes <- T
# min_indel_length <- 3
# mapq_threshold <- 20
# min_supporting_reads <- 1
# min_vaf <- 0.01
# min_read_depth <- 20

#############################################################################################
#
# Load Packages
#
#############################################################################################

suppressMessages(library(GenomicAlignments))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38)) 
suppressMessages(library(tictoc)) 
suppressMessages(library(bamsignals)) 
suppressMessages(library(parallel)) 
suppressMessages(library(inline))
suppressMessages(library(bettermc))

#############################################################################################
#
# Load Functions:
#
#############################################################################################

source("/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/functions.R")
#source("/Users/George/indel_detection_tool_project/indelfindr/src/functions.R")

#############################################################################################
#
# Load Utility Data
#
#############################################################################################

options(scipen=20)

# get all hg38 chromosome lengths
chg38_chromosome_lengths <- getChromInfoFromUCSC("hg38") # all hg38 chr/altcontig lengths  (595 chromosomes/contigs in hg38 annotation)

# load hg38 reference genome sequence data (688MB) - add option to load this optionally, if user opts to use BBtools Reformat tool, or other aligner, for extended cigar string, it will cut down on runtime and not require loading reference genome sequence into memory
hg38_genome <- BSgenome.Hsapiens.UCSC.hg38

#############################################################################################
#
# Running Script:
#
#############################################################################################

# define data table for reporting all refined cigar string frames in sliding windows

# master_indel_record_table <- c("chr","sliding_reference_start_record",
#                                "sliding_reference_end_record",
#                                "exploded_refined_cigar_string",
#                                "collapsed_refined_cigar_string",
#                                "ref_allele",
#                                "alt_allele",
#                                "strand",
#                                "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))


master_indel_record_table <- c("chr","start_pos",
                               "end_pos",
                               "refined_cigar_string",
                               "collapsed_cigar_string",
                               "reference_allele",
                               "alt_allele",
                               "strand",
                               "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))

# load target regions

if (target_regions != FALSE) {
  
  # load target regions
  target_regions_table <- read.table(target_regions,sep="\t",col.names = c("chr","start","stop"))
  chr_in_bam <- unique(target_regions_table$chr)
  
} else {
  
  # get chr names in bam file to search
  p = ScanBamParam(what=c("rname", "pos"))
  
  chr_in_bam <- na.omit(unique(as.data.frame(scanBam(bamPath, param=p))$rname))
}

# Subset primary chromosomes only for extracting reads for indel calling (optional)
if (primary_chromosomes == T){
  chr_in_bam <- get_primary_chroms(chr_in_bam)
}

for (each_chromosome in chr_in_bam){
  
  message(paste("Analyzing chromosome:",each_chromosome))
  
  # Define chr subset reference sequence
  hg38_genome_chr_subset <- hg38_genome[[each_chromosome]] 

    if (target_regions != FALSE){
    
    for (regions_table_row in 1:nrow(subset(target_regions_table,chr=chr_in_bam))){
      
      min_pos <- subset(target_regions_table,chr=chr_in_bam)$start[regions_table_row]
      max_pos <- subset(target_regions_table,chr=chr_in_bam)$stop[regions_table_row]
      
      # get total number reads 
      # Check if this is needed
      #overlap.counts <- suppressMessages((bamCount(bamPath,GRanges(Rle(each_chromosome,1), IRanges(start=min_pos, end=max_pos)))))
      
      target_bin_size <- bam_region_bin_size
      
      if (max_pos-min_pos < bam_region_bin_size){
        target_bin_size <- max_pos-min_pos
      }
      
      intervals <- seq_with_uneven_last(from=min_pos,to=max_pos,by=target_bin_size)
      
      # define coordinates of sliding windows
      sliding_windows_per_bam_region=data.frame(chr=each_chromosome,
                                                start=intervals[-(length(intervals))],
                                                end=intervals[-1])
      
      # Initialize dataframe for collecting each bam region indel calls
      # per_bam_region_indel_records <- c("chr","sliding_reference_start_record",
      #                                   "sliding_reference_end_record",
      #                                   "exploded_refined_cigar_string",
      #                                   "collapsed_refined_cigar_string",
      #                                   "ref_allele",
      #                                   "alt_allele",
      #                                   "strand",
      #                                   "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))
      # 
      per_bam_region_indel_records <- c("chr",
                                        "start_pos",
                                        "end_pos",
                                        "refined_cigar_string",
                                        "collapsed_cigar_string",
                                        "reference_allele",
                                        "alt_allele",
                                        "strand",
                                        "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))
      
      # Find indels in each bam region using all available cores in parallel

        # per_bam_region_indel_records <-
        #   do.call(
        #     rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))
      
      # Test without calling handelrs  
      per_bam_region_indel_records <-
          do.call(
            rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F))
      
      # Add each chromosome results to the master results table
      master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
      
      # Find indels in each bam region using all available cores in parallel
      # if (min_indel_length >= 2){
      #   
      #   per_bam_region_indel_records <-
      #     do.call(
      #       rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))
      # } else {
      #   print("This is working")
      #   per_bam_region_indel_records <-
      #     do.call(
      #       rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))
      # }
      
      ## Add each chromosome results to the master results table
      #master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
      
      wait()
      
     } # end each targeted region in bed file iteration
    
  } else if (target_regions == FALSE){ #end targeted region mode conditional 
    
    # do all analysis per chromosome, per bam region chunk alignment region, per read in alignment region
    
    size_chr <- chg38_chromosome_lengths[which(chg38_chromosome_lengths[,1]==each_chromosome),2]
    
    # p_for_chr_extract <- ScanBamParam(which=GRanges(
    #   Rle(each_chromosome),
    #   IRanges(1, size_chr)),
    #   what=c("pos"))
    # 
    # max_pos<- max(na.omit(unlist(scanBam(bamPath, param=p_for_chr_extract)))) # add na omit for finding max read start when there are NA positions returned - find out why and decide if should remove them. 
    # min_pos<- min(na.omit(unlist(scanBam(bamPath, param=p_for_chr_extract)))) # add na omit for finding max read start when there are NA positions returned - find out why and decide if should remove them. 
    # 
    # get total number reads 
    #overlap.counts <- suppressMessages((bamCount(bamPath,GRanges(Rle(each_chromosome,1), IRanges(start=min_pos, end=max_pos)))))
    
    target_bin_size <- bam_region_bin_size
    
    # if (max_pos-min_pos < bam_region_bin_size){
    #   target_bin_size <- max_pos-min_pos
    # }
    
    intervals <- seq_with_uneven_last(from=1,to=size_chr,by=target_bin_size)
    
    if (length(intervals)==1){
      p_for_read_length_extract <- ScanBamParam(which=GRanges(
        Rle(each_chromosome), 
        IRanges(1, size_chr)),
        what=c("pos","qwidth"))
      
      read_length <- unlist(scanBam(bamPath, param=p_for_read_length_extract))[2]
      intervals <- c(intervals[1],intervals[1]+read_length)
      
      rm(p_for_read_length_extract)
    }
    
    # define coordinates of sliding windows
    sliding_windows_per_bam_region=data.frame(chr=each_chromosome,
                                              start=intervals[-(length(intervals))],
                                              end=intervals[-1])
    
    #rm(p_for_chr_extract)
    
    # Initialize dataframe for collecting each bam region indel calls
    per_bam_region_indel_records <- c("chr",
                                      "start_pos",
                                      "end_pos",
                                      "refined_cigar_string",
                                      "collapsed_cigar_string",
                                      "reference_allele",
                                      "alt_allele",
                                      "strand",
                                      "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))
    
    # Find indels in each bam region using all available cores in parallel

      # per_bam_region_indel_records <-
      #   do.call(
      #     rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))

    # Find indels in each bam region using all available cores in parallel
 
    # Test without calling handelrs  
    per_bam_region_indel_records <-
      do.call(
        rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F))

    master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)

    # # Find indels in each bam region using all available cores in parallel
    # if (min_indel_length >= 2){
    #   
    #   per_bam_region_indel_records <-
    #     do.call(
    #       rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))
    #   # Add each chromosome results to the master results table
    #   
    # } else {
    #   print("This is working")
    #   per_bam_region_indel_records <-
    #     do.call(
    #       rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F,mc.stdout=c("output"),mc.warnings=c("output")))
    # }
    #master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
    
    wait()
  } # end genome wide mode conditional
  
} # end chr iteration

# Process results

#write a wrapper function around these processing steps to handle cases with no indels found, generate warnings etc.

# remove reads which cover multiple bam regions, were extracted, and searched more than once
dup_indices <- duplicated(master_indel_record_table[,c("read_name","chr","start_pos","end_pos","alt_allele")])
master_indel_record_table_no_dup_reads <- master_indel_record_table[!dup_indices,]

# get supporting read counts for each unique indel candidate
collapsed_read_counts_with_strand <- rename(dplyr::count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele,strand), FREQ = n)

collapsed_read_counts <- rename(count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele), FREQ = n)

# collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_stand_bias_pval_vaf_and_dp(collapsed_read_counts,master_indel_record_table_no_dup_reads)
collapsed_read_counts_with_strand_correction_pvalue_for_writing <- suppressMessages(calculate_strand_bias_pval_vaf_and_dp_parallel(collapsed_read_counts,master_indel_record_table_no_dup_reads,number_cores))

# filter one read calls
collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered <- collapsed_read_counts_with_strand_correction_pvalue_for_writing[which(collapsed_read_counts_with_strand_correction_pvalue_for_writing$FREQ>min_supporting_reads & collapsed_read_counts_with_strand_correction_pvalue_for_writing$VAF>=min_vaf  & collapsed_read_counts_with_strand_correction_pvalue_for_writing$DP >= min_read_depth),]

# adjust start_pos to reflect padding base
collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos <- as.integer(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos)
collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos <- as.integer(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos) # no longer needed after ref and query fix # - 1

# # convert variant coordinates from zero-based to 1-based coordinate system (ex. if RefGene was used)
# if (zero_based == T){
#   collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos <- collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos + 1
#   collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos <- collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos + 1
# }

message(paste0("Saving ",nrow(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered), " indels to output filepath: ",outname,"."))

# # make directory
# suppressWarnings(
#   "suppressWarnings("In dir.create(file.path("/Users/George/indel_detection_tool_project/test_output_dir/")) :
#   '/Users/George/indel_detection_tool_project/test_output_dir' already exists")"
# )

#dir.create(file.path(outname)) # suppress warnings

# write out VCF file
write.vcf(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered,phased=FALSE,paste0(outname,".vcf"))

#write out data table with cigar string
write.table(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered,paste0(outname,".csv"),sep=",",row.names=F)

message("Run Complete")
