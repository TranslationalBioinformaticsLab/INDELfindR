#############################################################################################
#
# Script to retrieve alignment data in single base pair sliding windows across aligned regions
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
# Outline:
#
#############################################################################################

# Outline of Indel Detection Tool Script Code Logic

# ### Initialize the Run
# 1. Load dependency packages silently
# 2. Load and parse arguments from command line
# 3. Load hg38 human reference genome
# 
# ## Run in either whole genome mode or targeted region mode:
# 
# ### Targeted Region Mode:
# 1. Load BED file with targeted regions
# 2. Determine what chromosomes are present in the BED file 
# 3. Initialize the `master_indel_record_table` with placeholders for the required indel record fields including genomic position, collapsed refined cigar string, exploded cigar string, alt, ref, strand, read name
# 4. For each targetted region, split it into ranges small to load the BAM file
# 5. Initialize the `per_bam_region_indel_records` tibble with placeholders for the required indel record fields including genomic position, collapsed refined cigar string, exploded cigar string, alt, ref, strand, read name
# 6. Iterate through each read in the bam region
# 7. Refine the cigar string to replace M operators with match (=) or mismatch (X) by comparing the query sequence to the reference sequence at that position and explode the cigar string into long format `10D5M` = `DDDDDDDDDDMMMMM`
# 8. At start of per read iteration, initialize indel counters and flags: 
#   1. `match_operator_counter`=0
# 2. `indel_candidate_container`=c()
# 3. `consecutive_indel_operator_flag`=F
# 9. Iterating over each cigar operator of the read. If the operator is a match operator `=` add 1 to the `match_operator_counter` and continue. *Note: A caveat to this algorithm is that indels which begin in the first 10 bases of the read will not be recognized, since we require a 10bp matching operator flanking region on either end of the indel candidate*
#   10. If the `match_operator_counter` >= 10 and the current operator is equal to I,D, or X , reset the `match_operator_counter` to zero, add the current operator to the `indel_candidate_container`, and set `consecutive_indel_operator_flag` == T and move to the next operator 
# 11. If the current operator is equal to `I,D or X` and `consecutive_indel_operator_flag` == T, add it to the `indel_candidate_container`
# 12. If the current operator is equal to `=`, add it to the `indel_candidate_container` and add 1 to the `match_operator_counter` and continue.
# 13. If the `match_operator_counter` => 10 and "consecutive_indel_operator_flag" == T
# 1. Set `consecutive_indel_operator_flag`=F
# 2. Remove the last 10 "=" operators from the `indel_candidate_container`
# 3. Get query seq cooridnates using current operator number - 10 as the end, and start is equal to end minus the length of the indel_candidate_container
# 4.  If `I or D` are in the `indel_candidate_container`, retrieve the following information from the read object using the candidate indel coordinates and assemble into a `candidate_indel_record` tibble row.:
#   1. genomic position
# 2. collapsed refined cigar string
# 3. exploded cigar string
# 4. alt
# 5. ref
# 6. strand 
# 7. read name
# 
# 5.  Add the `candidate_indel_record` to the `per_bam_region_indel_records` tibble 
# 6.  Set  `indel_candidate_container`=c()
# 14. Continue until the end of the read, then move to the next read
# 15. Once all of the reads have been processed in a bam region, add the `per_bam_region_indel_records` to the `master_indel_record_table`
# 16. Move on to the next BAM region
# 17. Repeat the per read searching in the next region
# 18. Once all bam regions are searched, Remove any duplicate rows with matching read name, genomic position, and alt alleles (not cigar string or ref to save resources, since alt allele and position are enough)
# 19. For each unique set of genomic positions and alt alleles, get all read names and save as vector. Append the vector to the `master_indel_record_table`
# 20. Count number of matching rows to get per read frequency for + and - strand separately
# 21. Join the tables by genomic position and alt alleles retaining + strand freq and - strand freq columns
# 22. For each indel record in the joined table, retrieve the strand info and read name of all overlapping reads withing that genomic range. Filter out reads which exist in that indel's sorted read name vector to get the reads not supporting the indel. Count the number of non-supporting reads which are + or - stranded. 
# 23. Assemble 2x2 contingency table using the number of +/- strand reads (in columns) carrying either reference or variant alleles (in rows).
# 24. Run Fishers Exact test to determine pvalue of strand bias. 
# 25. Use Fishers test outcome to rescale indel quality scores (as is done in variant score recalibration VSQR - recommended by GATK) or to hard filter indel from the output table.
# 26. Write output table as VCF.

#############################################################################################
#
# Define Args:
#
#############################################################################################

bamPath <- "/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/test_data/fastq/SRR12520438_1.sorted.bam"
#max_window_length <- 5 #dev-on
bam_region_bin_size <- 100000 #dev-on
verbose_arg=FALSE
#window_load_length <- 500
flanking_region_length <- 10
#sliding_window_size <- 20
#Optional parameters for testing:
#target_regions <- "/Users/George/indel_detection_tool_project/data_for_testing/EGFR_regions_of_interest.txt"
target_regions <- F
#number_cores <- 2
# number_cores <- 4
#number_cores <- 10
number_cores <- 8
# number_cores <- detectCores()
primary_chromosomes <- T

#############################################################################################
#
# Load Functions:
#
#############################################################################################

source("/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/functions.R")
 
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

start_time <- Sys.time()

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

if (primary_chromosomes == T){
  chr_in_bam <- get_primary_chroms(chr_in_bam)
}

for (each_chromosome in chr_in_bam){
  
 # each_chromosome <- "chr21"
  
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
     
      avail_cores <- number_cores
      per_bam_region_indel_records <-
        do.call(
          rbind, mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F))
      
      # Add each chromosome results to the master results table
      master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
      
      # clean up zombie processes
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
        rbind, mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region,per_bam_region_indel_records,mc.cores=number_cores,mc.preschedule = F))
    # Add each chromosome results to the master results table
    
    master_indel_record_table <- rbind(master_indel_record_table,per_bam_region_indel_records)
    
    # clean up zombie processes
    wait()
    
    write_per_chromosome_indel_table(master_indel_record_table,each_chromosome)
    
  } # end genome wide mode conditional
  
} # end chr iteration

# Process results

end_time <- Sys.time() - start_time

message(paste0("Completed in",end_time,"seconds."))

# remove reads which cover multiple bam regions, were extracted, and searched more than once
dup_indices <- duplicated(master_indel_record_table[,c("read_name","chr","start_pos","end_pos","alt_allele")])
master_indel_record_table_no_dup_reads <- master_indel_record_table[!dup_indices,]

# get supporting read counts for each unique indel candidate
collapsed_read_counts_with_strand <- rename(dplyr::count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele,strand), FREQ = n)

collapsed_read_counts <- rename(count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele), FREQ = n)

collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_stand_bias_pval_vaf_and_dp(collapsed_read_counts,master_indel_record_table_no_dup_reads)

end_processing_time <- Sys.time() - end_time

message(paste0("completed processing call data in",end_time-end_processing_time,"seconds"))
# write out VCF file
write.vcf(collapsed_read_counts_with_strand_correction_pvalue_for_writing,phased=FALSE,"/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/test_vcf_output_chr21_24core_nopre.vcf")

#write out data table with cigar string
write.table(collapsed_read_counts_with_strand_correction_pvalue_for_writing,"/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/test_output_chr21_24core_nopre.csv",sep=",",row.names=F)
