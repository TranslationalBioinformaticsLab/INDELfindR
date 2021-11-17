#############################################################################################
#
# Script to call indels using refined cigar strings in non-overlapping sliding windows
# genome wide or across the user specified genomic regions and output indel calls in VCF 4.3.
#
#############################################################################################

#############################################################################################
#
# Load Packages
#
#############################################################################################
suppressMessages(library(GenomicAlignments))
suppressMessages(library(argparse))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38)) 
suppressMessages(library(tictoc)) 
suppressMessages(library(bamsignals)) 
suppressMessages(library(parallel)) 
suppressMessages(library(inline))

#############################################################################################
#
# Define Args:
#
#############################################################################################

#bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/bam_for_testing_1656854-1657354.bam"
bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/EGFR_mutations_reference_dwgsim.sorted.bam"
#bamPath <- "/Users/George/indel_detection_tool_project/bam_files_for_testing/test_read.bam"
#max_window_length <- 5 #dev-on
bam_region_bin_size <- 5000 #dev-on
verbose_arg=FALSE
#window_load_length <- 500
flanking_region_length <- 10
#sliding_window_size <- 20
#Optional parameters for testing:
target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt"
#target_regions <- F
#target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/debug_10_27.txt"
number_cores <- 10
primary_chromosomes <- T


#############################################################################################
#
# Load Functions:
#
#############################################################################################

source("/Users/George/indel_detection_tool_project/indelfindr/src/functions.R")

#############################################################################################
#
# Load Utility Data
#
#############################################################################################

options(scipen=20)

# get all hg38 chromosome lengths
chg38_chromosome_lengths <- getChromInfoFromUCSC("hg38") # all hg38 chr lengths  (595 chrs in hg38 annotation including alts)

# load hg38 reference genome sequence data (688MB) - add option to load this optionally if user doesn't opt to use BBtools Reformat tool, or other aligner, for extended cigar string
hg38_genome <- BSgenome.Hsapiens.UCSC.hg38

#############################################################################################
#
# Running Script:
#
#############################################################################################

# define data table for reporting all refined cigar string frames in sliding windows

master_indel_record_table <- c("chr","sliding_reference_start_record",
                               "sliding_reference_end_record",
                               "exploded_refined_cigar_string",
                               "collapsed_refined_cigar_string",
                               "ref_allele",
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
  
  #head(as.data.frame(scanBam(bamPath,param=p)))
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
      overlap.counts <- suppressMessages((bamCount(bamPath,GRanges(Rle(each_chromosome,1), IRanges(start=min_pos, end=max_pos)))))
      
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
      per_bam_region_indel_records <- c("chr","sliding_reference_start_record",
                                        "sliding_reference_end_record",
                                        "exploded_refined_cigar_string",
                                        "collapsed_refined_cigar_string",
                                        "ref_allele",
                                        "alt_allele",
                                        "strand",
                                        "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))
     
      # Find indels in each bam region using all available cores in parallel

      per_bam_region_indel_records <-
        do.call(
          rbind, mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = T))
      
      # Add each chromosome results to the master results table
      master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
      
      wait()
      
     } # end each targeted region in bed file iteration
    
  } else if (target_regions == FALSE){ #end targeted region mode conditional 
    
    # do all analysis per chromosome, per bam region chunk alignment region, per read in alignment region
    
    size_chr <- chg38_chromosome_lengths[which(chg38_chromosome_lengths[,1]==each_chromosome),2]
    
    p_for_chr_extract <- ScanBamParam(which=GRanges(
      Rle(each_chromosome),
      IRanges(1, size_chr)),
      what=c("pos"))
    
    # unlist(scanBam(bamPath, param=p_for_chr_extract))
    
    max_pos<- max(na.omit(unlist(scanBam(bamPath, param=p_for_chr_extract)))) # add na omit for finding max read start when there are NA positions returned - find out why and decide if should remove them. 
    min_pos<- min(na.omit(unlist(scanBam(bamPath, param=p_for_chr_extract)))) # add na omit for finding max read start when there are NA positions returned - find out why and decide if should remove them. 
    
    # get total number reads 
    overlap.counts <- suppressMessages((bamCount(bamPath,GRanges(Rle(each_chromosome,1), IRanges(start=min_pos, end=max_pos)))))
    
    target_bin_size <- bam_region_bin_size
    
    if (max_pos-min_pos < bam_region_bin_size){
      target_bin_size <- max_pos-min_pos
    }
    
    intervals <- seq_with_uneven_last(from=min_pos,to=max_pos,by=target_bin_size)
    
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
    
    rm(p_for_chr_extract)
    
    # Initialize dataframe for collecting each bam region indel calls
    per_bam_region_indel_records <- c("chr","sliding_reference_start_record",
                                      "sliding_reference_end_record",
                                      "exploded_refined_cigar_string",
                                      "collapsed_refined_cigar_string",
                                      "ref_allele",
                                      "alt_allele",
                                      "strand",
                                      "read_name")  %>% purrr::map_dfc(setNames, object = list(as.character()))
    
    # Find indels in each bam region using all available cores in parallel
    
    per_bam_region_indel_records <-
      do.call(
        rbind, mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores))
    # Add each chromosome results to the master results table
    
    master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
    
    wait()
  } # end genome wide mode conditional
  
} # end chr iteration

# Process results

# remove reads which cover multiple bam regions, were extracted, and searched more than once
dup_indices <- duplicated(master_indel_record_table[,c("read_name","chr","start_pos","end_pos","alt_allele")])
master_indel_record_table_no_dup_reads <- master_indel_record_table[!dup_indices,]

# get supporting read counts for each unique indel candidate
collapsed_read_counts_with_strand <- rename(dplyr::count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele,strand), FREQ = n)

collapsed_read_counts <- rename(count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele), FREQ = n)

#collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_stand_bias_pval_vaf_and_dp(collapsed_read_counts,master_indel_record_table_no_dup_reads)
collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_strand_bias_pval_vaf_and_dp_parallel(collapsed_read_counts,master_indel_record_table_no_dup_reads,number_cores)

# write out VCF file
write.vcf(collapsed_read_counts_with_strand_correction_pvalue_for_writing,phased=FALSE,"/Users/George/indel_detection_tool_project/test_vcf_file_output/test_vcf_output.vcf")

#write out data table with cigar string
write.table(collapsed_read_counts_with_strand_correction_pvalue_for_writing,"/Users/George/indel_detection_tool_project/test_output/test_output_table.csv",sep=",",row.names=F)
