# local source path: /Users/George/indel_detection_tool_project/scripts/indel_detection_tool_functions.R

#############################################################################################
# 
# Store Functions for Main Rscript
# 
#############################################################################################

# Function 1) Convert bam to long format

convert_cigar_to_verbose_string <- function(each_cigar_string){
  
  cigar_operators <- unlist(str_split(each_cigar_string, "(\\d+)"))[-1]
  number_op_reps <- head(unlist(str_split(each_cigar_string, "[a-zA-Z]+")),-1)
  
  extended_string <- c()
  
  for (number_operators in 1:length(cigar_operators)){
    
    operator_reps <- rep(cigar_operators[number_operators],number_op_reps[number_operators])
    extended_string <- append(extended_string,operator_reps)
    
  }
  
  return(extended_string)
}

# Function 2.1) Improved Function to Refine M operators as match or mismatch using reference genome hg38

refine_cigar_string <- function(exploded_cigar_string,query_sequence_string,read_pos,query_read_length,each_chromosome,reference_sequence){

  exploded_cigar_string_for_refinement <- exploded_cigar_string
  
  reference_base_counter <- 0
  query_base_counter <- 0
  
  reference_sequence_testing <- unlist(str_split(toString(reference_sequence),""))
  query_sequence_string_testing <- unlist(str_split(toString(query_sequence_string),""))
  
  for (each_operator in 1:length(exploded_cigar_string_for_refinement)){
    
    if (exploded_cigar_string_for_refinement[[each_operator]]=="M"){
      
      reference_base_counter <- reference_base_counter + 1
      query_base_counter <- query_base_counter + 1
      #print("tick1")
      
      if (reference_sequence_testing[reference_base_counter] == query_sequence_string_testing[query_base_counter]){
        # print("match")
        exploded_cigar_string_for_refinement[[each_operator]] = "="
        # print("tick2")
        
      } else {
        # print("mismatch")
        exploded_cigar_string_for_refinement[[each_operator]] = "X" 
        # print("tick3")
        
      }
      
    } else if (exploded_cigar_string_for_refinement[[each_operator]] == "I"){
      query_base_counter <- query_base_counter + 1
      # print("tick4")
      
    } else if (exploded_cigar_string_for_refinement[[each_operator]] == "D"){
      reference_base_counter <- reference_base_counter + 1
      #  print("QUERY not added")
      
    }
    
  }
  
  # Debug view:
  # print(paste("Cigar operator:",exploded_cigar_string_for_refinement[each_operator]))
  # print(paste("Ref counter:",reference_base_counter))
  # print(paste("Query counter:",query_base_counter))
  
  return(exploded_cigar_string_for_refinement)
}

# Function 3: Translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

translate_cigar_index_to_ref_and_query <- function(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string){
  
  cigar_coord_start <- min(cigar_coords)
  cigar_coord_end <- max(cigar_coords)
  
  cigar_coords_for_query_end <- max(cigar_coords_for_query)
  
  exploded_cigar_string_record <- refined_cigar_string[cigar_coord_start:cigar_coord_end]
  
  D_indices <- which(exploded_cigar_string_record=="D")
  I_indices <- which(exploded_cigar_string_record=="I")
  equals_indices <- which(exploded_cigar_string_record=="=")
  X_indices <- which(exploded_cigar_string_record=="X")
  
  #define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices <- sort(c(D_indices,X_indices,equals_indices))
  reference_allele_record <- reference_sequence_for_translation[reference_indices]
  
  if (length(reference_allele_record) == 0){
    reference_allele_record = reference_sequence_for_translation[1]
  }
  
  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =
  alternate_indices <- sort(c(I_indices,X_indices,equals_indices))
  alternate_allele_record <- query_sequence_string_for_translation[alternate_indices]
 
  #reference_allele_record <- reference_sequence[cigar_coord_start:cigar_coord_end]
  #alternate_allele_record <- query_sequence_string[cigar_coord_start:cigar_coord_end]
  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")
  
  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.

}

# Function 3.1: Whole_read translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

translate_cigar_index_to_ref_and_query_whole_read <- function(reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string){
  
  exploded_cigar_string_record <- refined_cigar_string
  
  D_indices <- which(exploded_cigar_string_record=="D")
  I_indices <- which(exploded_cigar_string_record=="I")
  equals_indices <- which(exploded_cigar_string_record=="=")
  X_indices <- which(exploded_cigar_string_record=="X")
  
  #define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices <- sort(c(D_indices,X_indices,equals_indices))
  reference_allele_record <- reference_sequence_for_translation[reference_indices]
  
  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =
  #alternate_indices <- sort(c(I_indices,X_indices,equals_indices))
  alternate_allele_record <- query_sequence_string_for_translation
  
  #reference_allele_record <- reference_sequence[cigar_coord_start:cigar_coord_end]
  #alternate_allele_record <- query_sequence_string[cigar_coord_start:cigar_coord_end]
  cigar_string_record <- paste(refined_cigar_string,collapse=",")
  
  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.
  
}

# Function 3.2: Sliding, algo dev: Translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

translate_cigar_index_to_ref_and_query_algo_dev_sliding <- function(sliding_window_per_read_start,sliding_window_per_read_end,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string){
  
  cigar_coord_start <- sliding_window_per_read_start
  cigar_coord_end <- sliding_window_per_read_end
  
  exploded_cigar_string_record <- refined_cigar_string[cigar_coord_start:cigar_coord_end]
  
  D_indices <- which(exploded_cigar_string_record=="D")
  I_indices <- which(exploded_cigar_string_record=="I")
  equals_indices <- which(exploded_cigar_string_record=="=")
  X_indices <- which(exploded_cigar_string_record=="X")
  
  #define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices <- sort(c(D_indices,X_indices,equals_indices))
  reference_allele_record <- reference_sequence_for_translation[reference_indices]
  
  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =
  alternate_indices <- sort(c(I_indices,X_indices,equals_indices))
  alternate_allele_record <- query_sequence_string_for_translation[alternate_indices]
  
  #reference_allele_record <- reference_sequence[cigar_coord_start:cigar_coord_end]
  #alternate_allele_record <- query_sequence_string[cigar_coord_start:cigar_coord_end]
  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")
  
  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.
  
}

# Function 3.3: Translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

translate_cigar_index_to_ref_and_query_v2 <- function(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string){
  
  cigar_coord_start <- min(cigar_coords)
  cigar_coord_end <- max(cigar_coords)
  
  cigar_coords_for_query_end <- cigar_coords_for_query[length(cigar_coords_for_query)]
  cigar_coords_for_query_start <- cigar_coords_for_query[1]
  
  exploded_cigar_string_record <- refined_cigar_string[cigar_coord_start:cigar_coord_end]
  
  D_indices <- which(exploded_cigar_string_record=="D")
  I_indices <- which(exploded_cigar_string_record=="I")
  equals_indices <- which(exploded_cigar_string_record=="=")
  X_indices <- which(exploded_cigar_string_record=="X")
  
  #define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices <- sort(c(D_indices,X_indices,equals_indices))
  reference_allele_record <- reference_sequence_for_translation[reference_indices]
  
  if (length(reference_allele_record) == 0){
    reference_allele_record = reference_sequence_for_translation[1]
  }
  
  
  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =
  #alternate_indices <- sort(c(I_indices,X_indices,equals_indices))
  #alternate_allele_record <- query_sequence_string_for_translation[alternate_indices]
  
  alternate_allele_record <- query_sequence_string_for_translation
  
  if (cigar_coords_for_query_end < cigar_coords_for_query_start){
    alternate_allele_record <- alternate_allele_record[1]
  }
  
  #reference_allele_record <- reference_sequence[cigar_coord_start:cigar_coord_end]
  #alternate_allele_record <- query_sequence_string[cigar_coord_start:cigar_coord_end]
  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")
  
  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.
  
}

# Function 4: Convert exploded CIGAR string to normal shortened format CIGAR string

unexplode_cigar_string <- function(exploded_cigar_string){
  
  #use run length encoding to summarize exploded cigar string
  operator_lengths <- rle(exploded_cigar_string)[[1]]
  operator_values <- rle(exploded_cigar_string)[[2]]
  
  cigar_record <- ""
  for (operator_num in 1:length(operator_values)){
    opertator_substring <- str_c(operator_lengths[operator_num],operator_values[operator_num])
    cigar_record <- str_c(cigar_record,opertator_substring)
  }
  
  return(cigar_record)
}

# Function 5: Get Indel Read Depth

getIndelCoverage <- function(indel_chr,indel_start,indel_end,bamPath){
  
  indel_read_depth <- countOverlaps(GRanges(paste0(indel_chr,":",indel_start,"-",indel_end)),bamfile <- readGAlignments(bamPath),minoverlap = as.numeric((indel_end-indel_start)+1))
  
  return(indel_read_depth)
}

# Function 6: Get Variant Allele Frequency

getVaf <- function(num_supporting_reads,indel_read_depth){
  
  vaf <- num_supporting_reads/indel_read_depth
  return(vaf)
  
}

# Function 7

calculate_strand_bias_pval_vaf_and_dp_parallel <- function(collapsed_read_counts,master_indel_record_table_no_dup_reads,number_cores) {
  
  cols_to_add <- c("VAF","DP","strand_bias_pval","ref_forward","ref_rev","alt_forward","alt_rev")  %>% purrr::map_dfc(setNames, object = list(as.character()))
  
  cols_to_add <-
    do.call(
      rbind, bettermc::mclapply(1:nrow(collapsed_read_counts),calculate_strand_bias_pval_vaf_and_dp_parallel_each_indel,collapsed_read_counts,master_indel_record_table_no_dup_reads,mc.cores=number_cores,mc.preschedule = F))
  
  colnames(cols_to_add) <- c("VAF","DP","strand_bias_pval","ref_forward","ref_rev","alt_forward","alt_rev")
  
  collapsed_read_counts_with_strand_correction_pvalue_for_writing <- cbind(collapsed_read_counts,cols_to_add)
  
  return(collapsed_read_counts_with_strand_correction_pvalue_for_writing)
  wait()
  
}

# Function 7.1: Get strand bias pval, vaf, and dp for each indel (used by the parallelized function calculate_strand_bias_pval_vaf_and_dp_parallel())

calculate_strand_bias_pval_vaf_and_dp_parallel_each_indel <- function(each_unique_indel_candidate,collapsed_read_counts,master_indel_record_table_no_dup_reads){
  
  each_unique_indel_record <- collapsed_read_counts[each_unique_indel_candidate,]
  all_supporting_read_matches <- match_df(master_indel_record_table_no_dup_reads,each_unique_indel_record)
  indel_reads_vector <- all_supporting_read_matches[,c("read_name","strand")]
  
  variant_strand_counts <- table(indel_reads_vector$strand)
  
  strand_check_chr <- each_unique_indel_record$chr
  strand_check_start_pos <- as.numeric(each_unique_indel_record$start_pos)
  strand_check_end_pos<- as.numeric(each_unique_indel_record$end_pos)
  
  # set params for retrieving reads to extract non-supporting read strandedness counts
  strand_check_params=ScanBamParam(simpleCigar=FALSE, 
                                   which=GRanges(seqnames=Rle(strand_check_chr), 
                                                 ranges=IRanges(strand_check_start_pos, strand_check_end_pos)), 
                                   mapqFilter=mapq_threshold,
                                   what=c("rname","strand"), # make mapqFilter an arg to parse
                                   scanBamFlag(isSecondaryAlignment=FALSE,
                                               isNotPassingQualityControls=FALSE,
                                               isDuplicate=FALSE, #user would need bam with marked duplicates for this to be useful. Can test with BAM with marked duplicates. This could cause issues with amplicon seq
                                               isUnmappedQuery=FALSE))
  
  gal_read_counts <- readGAlignments(bamPath, param=strand_check_params, use.names=TRUE)
  non_supporting_reads <- gal_read_counts[!(rownames(mcols(gal_read_counts)) %in% indel_reads_vector$read_name),]
  
  reference_strand_counts <- table(mcols(non_supporting_reads)$strand)
  
  # make contingency tables
  contingency_table <- rbind(reference_strand_counts[c("+","-")],variant_strand_counts[c("+","-")])
  
  ref_forward <- reference_strand_counts["+"]
  ref_rev <- reference_strand_counts["-"]
  alt_forward <- variant_strand_counts["+"]
  alt_rev <- variant_strand_counts["-"]
  
  contingency_table[is.na(contingency_table)] <- 0
  
  ref_forward[is.na(ref_forward)] <- 0
  ref_rev[is.na(ref_rev)] <- 0
  alt_forward[is.na(alt_forward)] <- 0
  alt_rev[is.na(alt_rev)] <- 0
  
  pvalue <- fisher.test(contingency_table)$p.value
  
  dp_value <- length(gal_read_counts)
  vaf_value <- (contingency_table[2,1]+contingency_table[2,2])/dp_value
  
  row_to_add_to_cols_to_add <- c(vaf_value,dp_value,pvalue,ref_forward,ref_rev,alt_forward,alt_rev)
  
  return(row_to_add_to_cols_to_add)
}


# Function 8: Write indel calls to VCF file output

write.vcf <- function(master_indel_call_table,phased,output_filename){
  
  results_table <- master_indel_call_table

  file <- output_filename
  
  # assign column name
  sample_name <- "EGFR_mutation_bam"
  assign("SAMPLE_ID", sample_name)
  
  if (phased == TRUE){
    phase_seperator <- "|"
  } else if (phased == FALSE){
    phase_seperator <- "/"
  }
  
  genotype_col<- rep(paste0(".",phase_seperator,"."),nrow(results_table))
  
  format_col_pre <- data.frame(cbind(results_table$DP,genotype_col))
  
  format_col <- unlist(lapply(1, function(x) paste(format_col_pre[,"genotype_col"],format_col_pre[,"V1"], sep=":")))
  
  af_col <- unlist(lapply(results_table$VAF, function(x) paste("AF=",x, sep="")))
  strand_bias_col <- unlist(lapply(results_table$strand_bias_pval, function(x) paste("SB=",x, sep="")))
  
  ref_fwd_count_col <- unlist(lapply(results_table$ref_forward, function(x) paste("RF=",x, sep="")))
  ref_rev_count_col <- unlist(lapply(results_table$ref_rev, function(x) paste("RR=",x, sep="")))
  alt_fwd_count_col <- unlist(lapply(results_table$alt_forward, function(x) paste("VF=",x, sep="")))
  alt_rev_count_col <- unlist(lapply(results_table$alt_rev, function(x) paste("VR=",x, sep="")))

  INFO_col <- unlist(lapply(1, function(x) paste(af_col,strand_bias_col,ref_fwd_count_col,ref_rev_count_col,alt_fwd_count_col,alt_rev_count_col, sep=";")))
  
  results_table$reference_allele[results_table$reference_allele == ""] <- NA
  results_table$alt_allele[results_table$alt_allele == ""] <- NA
  
  CHROM_record <- results_table$chr
  POS_record	 <- results_table$start_pos
  ID_record <- rep(".",nrow(results_table))
  REF_record <- results_table$reference_allele
  ALT_record <- results_table$alt_allele
  QUAL_record <- rep(".",nrow(results_table))
  FILTER_record <- rep(".",nrow(results_table))
  INFO_record <- INFO_col
  FORMAT_record <- rep("GT:DP")
  SAMPLE_DATA_record<- format_col
  
  vcf_data_lines <- tibble(CHROM_record,POS_record,ID_record,REF_record,ALT_record,QUAL_record,FILTER_record,INFO_record,FORMAT_record,SAMPLE_DATA_record)
  
  colnames(vcf_data_lines) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",SAMPLE_ID)
  
  vcf_data_lines_sorted <- vcf_data_lines[with(vcf_data_lines, order(CHROM,POS)),]
  
  cat(file=file, '##fileformat=VCFv4.3\n##filedate=',as.character(Sys.Date()),'\n##source="Produced with write.vcf() of indelfindR R package v0.0.1 for complex indel detection"\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n##INFO=<ID=SB,Number=A,Type=Float,Description="Strand Bias Pvalue">\n##INFO=<ID=RF,Number=A,Type=Float,Description="Reference Forward Count">\n##INFO=<ID=RR,Number=A,Type=Float,Description="Reference Reverse Count">\n##INFO=<ID=VF,Number=A,Type=Float,Description="Alternate Forward Count">\n##INFO=<ID=VR,Number=A,Type=Float,Description="Alternate Reverse Count">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#')
  cat(file=file, paste0(colnames(vcf_data_lines_sorted), collapse="\t"),"\n", append=TRUE)
  write.table(vcf_data_lines_sorted, file=file,
              quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t", na=paste("."))
  
}

# Function 9: Get intervals with last interval included when it is not within the by interval

seq_with_uneven_last <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# Function 10: Run algorithm on each read and collect indel calls to add to master table

run_algo_on_one_read_explicit_args <- function(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,each_chromosome,min_indel_length){
  
  # define algo counters
  match_operator_counter=0
  indel_candidate_container=c()
  consecutive_indel_operator_flag=F
  number_deletions_encountered <- 0
  
  # iterate through each cigar operator per read
  for (each_operator in 1:length(refined_cigar_string)){
    
    #print(indel_candidate_container)
    
    operator <- refined_cigar_string[each_operator]
    
    # add exception for soft clipping, consider S as match operators 
    if (operator == "S"){
      operator <- "="
    }
    
    # if the operator matches the reference and there is not a current indel candidate
    if (operator == "=" && consecutive_indel_operator_flag == F){
      
      match_operator_counter = match_operator_counter + 1
      #print("ONE")
      
      # if the first variant operator is encountered: 
    } else if (operator %in% c("I","D") && match_operator_counter >= flanking_region_length && consecutive_indel_operator_flag == F){
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO")
      
      # Catch indels which begin at very start of read and dont have 10bp flanking prior. (Ff the first variant operator is encountered before the 10th bp in the read)
    } else if (operator %in% c("I","D") && each_operator < flanking_region_length && consecutive_indel_operator_flag == F) {
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO.5")
      
      # If there is an immediately adjacent variant operator
    } else if (operator %in% c("I","D") && (match_operator_counter < flanking_region_length) && consecutive_indel_operator_flag == T){
      
      indel_candidate_container <- c(indel_candidate_container,operator)
      
      match_operator_counter <- 0
      
      #print("THREE")
      
      # count match operators after encountering variant operators
    } else if ((operator == "=" && consecutive_indel_operator_flag == T) | (operator == "X" && consecutive_indel_operator_flag == T)){
      indel_candidate_container <- c(indel_candidate_container,operator)
      match_operator_counter = match_operator_counter + 1
      
      #print("FOUR")
    } 
    
    # Check for complete indel match: if there are more than 10 match operators after a candidate indel is found. See below for exception for indel reaching end of read.
    if (match_operator_counter >= flanking_region_length && consecutive_indel_operator_flag == T){
      
      # reset the consecutive_indel_operator_flag to false, since string of closely spaced operators is broken
      consecutive_indel_operator_flag = F
      
      #print("FIVE")
      
      #print(indel_candidate_container)
      
      first_remove <- length(indel_candidate_container)-(flanking_region_length-1)
      last_remove <- length(indel_candidate_container)
      
      indel_candidate_container <- indel_candidate_container[-(first_remove:last_remove)]
      
      # calculate number_deletions_per_candidate and add to number_delections_encountered
      number_deletions_per_candidate <- sum(indel_candidate_container == "D")
      number_deletions_encountered <- number_deletions_encountered + number_deletions_per_candidate
      
      # begin conditional min_indel_length conditional      
      if (length(indel_candidate_container) >= min_indel_length){
        
        cigar_end <- each_operator - flanking_region_length
        cigar_start <- cigar_end-(length(indel_candidate_container)-1)
        cigar_coords <- cigar_start:cigar_end
        
        # Define indel records and add indel candidate record to per bam region table
        cigar_end_for_query <- each_operator-flanking_region_length-number_deletions_encountered
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)
        cigar_coords_for_query<- cigar_start:cigar_end_for_query
        
        reference_start_record <- read_pos+cigar_start-number_leading_softclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container)-1
        chr_record <- each_chromosome
        
        cigar_start_for_reference <- cigar_start-number_leading_softclips-1
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container)-1
        
        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] #can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)
        
        candidate_indel_record <-tibble(chr=chr_record,
                                        start_pos=reference_start_record,
                                        end_pos=reference_end_record,
                                        refined_cigar_string=toString(exploded_cigar_string_record),
                                        collapsed_cigar_string=cigar_string_record,
                                        reference_allele=reference_allele_record,
                                        alt_allele=alternate_allele_record,
                                        strand = read_strand,
                                        read_name = read_name_record)
        
        #add indel record to per region table
        per_bam_region_indel_records <- rbind(per_bam_region_indel_records,candidate_indel_record)
        
        # clear the indel_candidate_container
        indel_candidate_container=c()
        
      } else {
        # don't save indel and keep moving on nothing
        indel_candidate_container=c()
      } # end conditional min_indel_length conditional
      
      
      # add exception if the indel candidate runs all the way into the end of the read
    } else if (each_operator == length(refined_cigar_string) && consecutive_indel_operator_flag == T){
      
      # reset the consecutive_indel_operator_flag to false, since string of closely spaced operators is broken
      consecutive_indel_operator_flag = F
      
      #print("SIX")
      
      # remove flanking "=" or "X" operators
      
      candidate_rle <- rle(indel_candidate_container)
      
      pattern <- c("D","I")
      
      # if candidate ends in "i" or "D" remove nothing
      if (candidate_rle[2]$values[length(candidate_rle[2]$values)] %in% pattern){
        
        num_to_remove <- 0
        
      } else {
        
        # get max index of I or D in RLE values
        operator_to_cut_to <- max(which(candidate_rle[2]$values=="I" | candidate_rle[2]$values== "D"))
        # define which rle values to keep (all those between D or I)
        operators_before_operator_to_cut_to <- 1:operator_to_cut_to
        # count number X or = operators to remove
        num_to_remove <- sum(candidate_rle$lengths[-c(1:operator_to_cut_to)])
        
        last_remove <- length(indel_candidate_container)-(num_to_remove-1)
        first_remove <- length(indel_candidate_container)
        
        # trim indel candidate to remove non indel operators
        indel_candidate_container <- indel_candidate_container[-(first_remove:last_remove)]
        
      }
      
      number_deletions_per_candidate <- sum(indel_candidate_container == "D")
      number_deletions_encountered <- number_deletions_encountered + number_deletions_per_candidate
      
      # begin conditional min_indel_length conditional      
      if (length(indel_candidate_container) >= min_indel_length){
        
        #print(indel_candidate_container)
        
        cigar_end <- each_operator-num_to_remove
        cigar_start <- cigar_end-(length(indel_candidate_container)-1)
        cigar_coords <- cigar_start:cigar_end
        
        cigar_end_for_query <- each_operator-num_to_remove-number_deletions_encountered
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)
        cigar_coords_for_query<- cigar_start:cigar_end_for_query
        
        # Define indel records and add indel candidate record to per bam region table
        reference_start_record <- read_pos+cigar_start-number_leading_softclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container)-1
        
        chr_record <- each_chromosome

        cigar_start_for_reference <- cigar_start-number_leading_softclips-1
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container)-1
        
        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] # can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)
        
        candidate_indel_record <-tibble(chr=chr_record,
                                        start_pos=reference_start_record,
                                        end_pos=reference_end_record,
                                        refined_cigar_string=toString(exploded_cigar_string_record),
                                        collapsed_cigar_string=cigar_string_record,
                                        reference_allele=reference_allele_record,
                                        alt_allele=alternate_allele_record,
                                        strand = read_strand,
                                        read_name = read_name_record)
        
        #add indel record to per region table
        per_bam_region_indel_records <- rbind(per_bam_region_indel_records,candidate_indel_record)
        
        # clear the indel_candidate_container
        indel_candidate_container=c()
        
      } else {
        # end conditional min_indel_length conditional      
        indel_candidate_container=c()
        
      }
      
    }
    
  } # end each operator iteration
  
  return(per_bam_region_indel_records)
  
}

# Function 11 write intermediate chrom calls table:

write_per_chromosome_indel_table <- function(master_indel_record_table,each_chromosome){
  
  # remove reads which cover multiple bam regions, were extracted, and searched more than once
  dup_indices <- duplicated(master_indel_record_table[,c("read_name","chr","start_pos","end_pos","alt_allele")])
  master_indel_record_table_no_dup_reads <- master_indel_record_table[!dup_indices,]
  
  # get supporting read counts for each unique indel candidate
  collapsed_read_counts_with_strand <- rename(dplyr::count(master_indel_record_table_no_dup_reads[-1,], chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele,strand), FREQ = n)
  
  collapsed_read_counts <- rename(count(master_indel_record_table_no_dup_reads[-1,], chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele), FREQ = n)
  
  collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_stand_bias_pval_vaf_and_dp(collapsed_read_counts,master_indel_record_table_no_dup_reads)
  
  # write out VCF file
  write.vcf(collapsed_read_counts_with_strand_correction_pvalue_for_writing,phased=FALSE,paste0("/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/test_vcf_output_SRR12520438_1.sorted_",each_chromosome,".vcf"))
  
  #write out data table with cigar string
  write.table(collapsed_read_counts_with_strand_correction_pvalue_for_writing,paste0("/gpfs/data/dgamsiz/Uzun_Lab/gtollefs/indel_detection_project/running_directory/test_output_SRR12520438_1.sorted_,",each_chromosome,".csv"),sep=",",row.names=F)
  
}

# Function 12: Extract all reads per bam region and run algo

run_algo_all_reads_each_bam_region <- function(row_num,per_bam_region_indel_records){
  
  bam_region_number <- row_num
  bam_region_chr <- each_chromosome
  bam_region_start <- sliding_windows_per_bam_region$start[row_num]
  bam_region_end <-  sliding_windows_per_bam_region$end[row_num]
  
  # apply filters:
  
  # 1) isSecondaryAlignment=FALSE
  # 2) isNotPassingQualityControls=FALSE
  # 3) isDuplicate=FALSE
  # 4) isUnmappedQuery=FALSE
  # use mapq filer=20 arbitrarily until determine which to use
  
  parameters=ScanBamParam(simpleCigar=FALSE, 
                          which=GRanges(Rle(c(sliding_windows_per_bam_region[row_num,1])), 
                                        IRanges(c(sliding_windows_per_bam_region[row_num,2]), c(sliding_windows_per_bam_region[row_num,3]))), 
                          mapqFilter=mapq_threshold,
                          what=c("rname","pos","strand","cigar","seq"), # make mapqFilter an arg to parse
                          #what=c("rname","pos","strand","cigar","seq","qual"), # make mapqFilter an arg to parse
                          scanBamFlag(isSecondaryAlignment=FALSE,
                                      isNotPassingQualityControls=FALSE,
                                      isDuplicate=FALSE, #user would need bam with marked duplicates for this to be useful. Can test with BAM with marked duplicates. This could cause issues with amplicon seq
                                      isUnmappedQuery=FALSE))
  
  gal <- readGAlignments(bamPath, param=parameters, use.names=TRUE)
  
  if (length(gal) >= 1){ # minimum bam region coverage 1 read
    
    for (read_num in 1:length(gal)){
      
      #print(read_num)
      
      # get cigar string, check if indel operators are in the read, and decide to do work
      cigar_string <- mcols(gal)$cigar[[read_num]]
      
      indel_operators <- c("I","D")
      pattern = paste(indel_operators, collapse="|")
      
      if (grepl(pattern, cigar_string)==T){
        
        #message(paste("This read has indel operators:",cigar_string))
        
        # if indel is present in cigar string, define read other bam fields
        query_sequence_string <- mcols(gal)$seq[[read_num]]
        #read_qual_string <- mcols(gal)$qual[[read_num]]
        read_pos <- mcols(gal)$pos[[read_num]]
        query_read_length <- length(query_sequence_string)
        read_strand <- mcols(gal)$strand[[read_num]]
        read_name_record <- (names(gal)[[read_num]])
        
        # Remove
        # if (verbose_arg==T){
        #   
        #   # print debugging info
        #   message("Current read:")
        #   message(read_name_record)
        #   message(each_chromosome)
        #   message(read_pos)
        #   
        # }
        
        # Run function to convert cigar to long format
        exploded_cigar_string <- convert_cigar_to_verbose_string(cigar_string)
        
        # determine how much reference sequence to extract by counting the D's
        number_deletions <- sum(exploded_cigar_string == "D")
        number_insertions <- sum(exploded_cigar_string == "I")
        
        operator_values <- rle(exploded_cigar_string)[[2]]
        operator_lengths <- rle(exploded_cigar_string)[[1]]
        
        if (operator_values[1] == "S"){
          number_leading_softclips <- operator_lengths[1]
        } else {
          number_leading_softclips <- 0
        }
        
        # get reference sequence for alignment region
        reference_sequence <- Views(hg38_genome_chr_subset, start=read_pos, end=read_pos+query_read_length+number_deletions)[[1]]
        
        # refine cigar string 
        refined_cigar_string <- refine_cigar_string(exploded_cigar_string,query_sequence_string,read_pos,query_read_length,each_chromosome,reference_sequence)
        
        # Run algorithm on read and collect indel calls
        per_bam_region_indel_records <- run_algo_on_one_read_explicit_args(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,each_chromosome,min_indel_length)
        
      } # end conditional if read has any I or D operators
      
    } # end per read iteration
    
    rm(gal)
    
  } # end coverage > 2 conditional
  
  return(per_bam_region_indel_records)
} # end run algo all reads function

# Function 11.1: Wrap run_algo_all_reads_each_bam_region in a tryCatch function for writing errors to log.

safe_run_algo_all_reads_each_bam_region <- function(row_num,per_bam_region_indel_records,log.path){
  tryCatch({
    
    #Function here
    bam_region_results <- run_algo_all_reads_each_bam_region(row_num,per_bam_region_indel_records)
    
  }, error = function(err.msg){
    # Add error message to the error log file
    write(paste0(toString(err.msg),"encountered processing read number:",read_num," in bam region number",row_num," at ",bam_region_chr,":",read_pos," - read name:",read_name_record), log.path, append=TRUE)
  })
  
  return(bam_region_results)
  
}


# Function 12: Define C function to wait on child processes to get rid of zombie processes

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

# Function 13: Get primary chromosomes in bam file
get_primary_chroms <- function(chr_in_bam){
  
if (startsWith(as.character(chr_in_bam[1]),"chr")){
  primary_chromosome_list_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",'chr10',"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
} else if(startsWith(as.character(chr_in_bam[1]),"Chr")){
  primary_chromosome_list_Chr <- c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9",'Chr10',"Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21","Chr22","ChrX","ChrY","ChrM")
} else {
  primary_chromosome_list_nochr <- c("1","2","3","4","5","6","7","8","9",'10',"11","12","13","14","15","16","17","18","19","20","21","22","X","Y","M")
}

chr_in_bam <- intersect(chr_in_bam,primary_chromosome_list_chr)

return(chr_in_bam)
}


# Function 14: Get ref and alt alleles, with extra conditional for adding padding base for simple insertion or deletions (only added to single bp indels per vcf 4.3 specifications)

translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding <- function(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string,padding_base){
  
  cigar_coord_start <- min(cigar_coords)
  cigar_coord_end <- max(cigar_coords)
  
  cigar_coords_for_query_end <- cigar_coords_for_query[length(cigar_coords_for_query)]
  cigar_coords_for_query_start <- cigar_coords_for_query[1]
  
  exploded_cigar_string_record <- refined_cigar_string[cigar_coord_start:cigar_coord_end]
  
  D_indices <- which(exploded_cigar_string_record=="D")
  I_indices <- which(exploded_cigar_string_record=="I")
  equals_indices <- which(exploded_cigar_string_record=="=")
  X_indices <- which(exploded_cigar_string_record=="X")
  
  # define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices <- sort(c(D_indices,X_indices,equals_indices))
  reference_allele_record <- reference_sequence_for_translation[reference_indices]
  
  if (length(reference_allele_record) == 0){
    reference_allele_record = reference_sequence_for_translation[1]
  }
  
  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =

  alternate_allele_record <- query_sequence_string_for_translation
  
  if (cigar_coords_for_query_end < cigar_coords_for_query_start){
    alternate_allele_record <- alternate_allele_record[1]
  }
  
  #test to see if it is a simple insertion
  # if (length(D_indices)==1 & length(I_indices)==0){
  #   reference_allele_record <- paste0(alternate_allele_record,reference_allele_record)
  # } else if (length(I_indices)==1 & length(D_indices)==0){
  #   alternate_allele_record <- paste0(alternate_allele_record,reference_allele_record)
  # }
  
  # add preceding padding base
   reference_allele_record <- paste0(padding_base,reference_allele_record)
   alternate_allele_record <- paste0(padding_base,alternate_allele_record)
  # 
  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")
  
  return(list(reference_allele_record,alternate_allele_record,cigar_string_record)) # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.
  
}


# Function 15: Run algorithm on each read and collect indel calls to add to master table. Version for adding padding to simple indels.

run_algo_on_one_read_explicit_args_with_simple_indel_padding <- function(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,each_chromosome,min_indel_length){
  
  # define algo counters
  match_operator_counter=0
  indel_candidate_container=c()
  consecutive_indel_operator_flag=F
  number_deletions_encountered <- 0
  
  # iterate through each cigar operator per read
  for (each_operator in 1:length(refined_cigar_string)){
    
    #print(indel_candidate_container)
    
    operator <- refined_cigar_string[each_operator]
    
    # add exception for soft clipping, consider S as match operators 
    if (operator == "S"){
      operator <- "="
    }
    
    # if the operator matches the reference and there is not a current indel candidate
    if (operator == "=" && consecutive_indel_operator_flag == F){
      
      match_operator_counter = match_operator_counter + 1
      #print("ONE")
      
      # if the first variant operator is encountered: 
    } else if (operator %in% c("I","D") && match_operator_counter >= flanking_region_length && consecutive_indel_operator_flag == F){
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO")
      
      # Catch indels which begin at very stard of read and dont have 10bp flanking prior. (Ff the first variant operator is encountered before the 10th bp in the read)
    } else if (operator %in% c("I","D") && each_operator < flanking_region_length && consecutive_indel_operator_flag == F) {
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO.5")
      
      # If there is an immediately adjacent variant operator
    } else if (operator %in% c("I","D") && (match_operator_counter < flanking_region_length) && consecutive_indel_operator_flag == T){
      
      indel_candidate_container <- c(indel_candidate_container,operator)
      
      match_operator_counter <- 0
      
      #print("THREE")
      
      # count match operators after encountering variant operators
    } else if ((operator == "=" && consecutive_indel_operator_flag == T) | (operator == "X" && consecutive_indel_operator_flag == T)){
      indel_candidate_container <- c(indel_candidate_container,operator)
      match_operator_counter = match_operator_counter + 1
      
      #print("FOUR")
    } 
    
    # Check for complete indel match: if there are more than 10 match operators after a candidate indel is found. See below for exception for indel reaching end of read.
    if (match_operator_counter >= flanking_region_length && consecutive_indel_operator_flag == T){
      
      # reset the consecutive_indel_operator_flag to false, since string of closely spaced operators is broken
      consecutive_indel_operator_flag = F
      
      #print("FIVE")
      
      #print(indel_candidate_container)
      
      first_remove <- length(indel_candidate_container)-(flanking_region_length-1)
      last_remove <- length(indel_candidate_container)
      
      indel_candidate_container <- indel_candidate_container[-(first_remove:last_remove)]
      
      # calculate number_deletions_per_candidate and add to number_delections_encountered
      number_deletions_per_candidate <- sum(indel_candidate_container == "D")
      number_deletions_encountered <- number_deletions_encountered + number_deletions_per_candidate
      
      # begin conditional min_indel_length conditional      
      if (length(indel_candidate_container) >= min_indel_length){
        
        cigar_end <- each_operator - flanking_region_length
        cigar_start <- cigar_end-(length(indel_candidate_container)-1)
        cigar_coords <- cigar_start:cigar_end
        
        # Define indel records and add indel candidate record to per bam region table
        cigar_end_for_query <- each_operator-flanking_region_length-number_deletions_encountered
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)
        cigar_coords_for_query<- cigar_start:cigar_end_for_query
        
        reference_start_record <- read_pos+cigar_start-number_leading_softclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container)-1
        chr_record <- each_chromosome

        cigar_start_for_reference <- cigar_start-number_leading_softclips-1
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container)-1
        
        padding_base <- reference_sequence[cigar_start_for_reference-1]
          
        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string,padding_base)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] #can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)
        
        candidate_indel_record <-tibble(chr=chr_record,
                                        start_pos=reference_start_record,
                                        end_pos=reference_end_record,
                                        refined_cigar_string=toString(exploded_cigar_string_record),
                                        collapsed_cigar_string=cigar_string_record,
                                        reference_allele=reference_allele_record,
                                        alt_allele=alternate_allele_record,
                                        strand = read_strand,
                                        read_name = read_name_record)
        
        #add indel record to per region table
        per_bam_region_indel_records <- rbind(per_bam_region_indel_records,candidate_indel_record)
        
        # clear the indel_candidate_container
        indel_candidate_container=c()
        
      } else {
        # don't save indel and keep moving on nothing
        indel_candidate_container=c()
      } # end conditional min_indel_length conditional
      
      
      # add exception if the indel candidate runs all the way into the end of the read
    } else if (each_operator == length(refined_cigar_string) && consecutive_indel_operator_flag == T){
      
      # reset the consecutive_indel_operator_flag to false, since string of closely spaced operators is broken
      consecutive_indel_operator_flag = F
      
      #print("SIX")
      
      #print(indel_candidate_container)
      
      # remove flanking "=" or "X" operators
      
      candidate_rle <- rle(indel_candidate_container)
      
      pattern <- c("D","I")
      
      # if candidate ends in "i" or "D" remove nothing
      if (candidate_rle[2]$values[length(candidate_rle[2]$values)] %in% pattern){
        
        num_to_remove <- 0
        
      } else {
        
        # get max index of I or D in RLE values
        operator_to_cut_to <- max(which(candidate_rle[2]$values=="I" | candidate_rle[2]$values== "D"))
        # define which rle values to keep (all those between D or I)
        operators_before_operator_to_cut_to <- 1:operator_to_cut_to
        # count number X or = operators to remove
        num_to_remove <- sum(candidate_rle$lengths[-c(1:operator_to_cut_to)])
        
        last_remove <- length(indel_candidate_container)-(num_to_remove-1)
        first_remove <- length(indel_candidate_container)
        
        # trim indel candidate to remove non indel operators
        indel_candidate_container <- indel_candidate_container[-(first_remove:last_remove)]
        
      }
      
      number_deletions_per_candidate <- sum(indel_candidate_container == "D")
      number_deletions_encountered <- number_deletions_encountered + number_deletions_per_candidate
      
      # begin conditional min_indel_length conditional      
      if (length(indel_candidate_container) >= min_indel_length){
        
        #print(indel_candidate_container)
        
        cigar_end <- each_operator-num_to_remove
        cigar_start <- cigar_end-(length(indel_candidate_container)-1)
        cigar_coords <- cigar_start:cigar_end
        
        cigar_end_for_query <- each_operator-num_to_remove-number_deletions_encountered
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)
        cigar_coords_for_query<- cigar_start:cigar_end_for_query
        
        # Define indel records and add indel candidate record to per bam region table
        reference_start_record <- read_pos+cigar_start-number_leading_softclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container)-1
        
        chr_record <- each_chromosome
        cigar_start_for_reference <- cigar_start-number_leading_softclips-1
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container)-1
        
        padding_base <- reference_sequence[cigar_start_for_reference-1]
        
        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string,padding_base)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] # can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)
        
        candidate_indel_record <-tibble(chr=chr_record,
                                        start_pos=reference_start_record,
                                        end_pos=reference_end_record,
                                        refined_cigar_string=toString(exploded_cigar_string_record),
                                        collapsed_cigar_string=cigar_string_record,
                                        reference_allele=reference_allele_record,
                                        alt_allele=alternate_allele_record,
                                        strand = read_strand,
                                        read_name = read_name_record)
        
        #add indel record to per region table
        per_bam_region_indel_records <- rbind(per_bam_region_indel_records,candidate_indel_record)
        
        # clear the indel_candidate_container
        indel_candidate_container=c()
        
      } else {
        
        # end conditional min_indel_length conditional      
        indel_candidate_container=c()
        
      }
      
    }
    
  } # end each operator iteration
  
  return(per_bam_region_indel_records)
  
}

# Function 16: Extract all reads per bam region and run algo. Version to add simple indel padding.

run_algo_all_reads_each_bam_region_with_simple_indel_padding <- function(row_num,per_bam_region_indel_records){
  
  bam_region_number <- row_num
  bam_region_chr <- each_chromosome
  bam_region_start <- sliding_windows_per_bam_region$start[row_num]
  bam_region_end <-  sliding_windows_per_bam_region$end[row_num]
  
  # apply filters:
  
  # 1) isSecondaryAlignment=FALSE
  # 2) isNotPassingQualityControls=FALSE
  # 3) isDuplicate=FALSE
  # 4) isUnmappedQuery=FALSE
  # use mapq filer=20 arbitrarily until determine which to use
  
  parameters=ScanBamParam(simpleCigar=FALSE, 
                          which=GRanges(Rle(c(sliding_windows_per_bam_region[row_num,1])), 
                                        IRanges(c(sliding_windows_per_bam_region[row_num,2]), c(sliding_windows_per_bam_region[row_num,3]))), 
                          mapqFilter=mapq_threshold,
                          what=c("rname","pos","strand","cigar","seq"), # make mapqFilter an arg to parse
                          #what=c("rname","pos","strand","cigar","seq","qual"), # make mapqFilter an arg to parse
                          scanBamFlag(isSecondaryAlignment=FALSE,
                                      isNotPassingQualityControls=FALSE,
                                      isDuplicate=FALSE, #user would need bam with marked duplicates for this to be useful. Can test with BAM with marked duplicates. This could cause issues with amplicon seq
                                      isUnmappedQuery=FALSE))
  
  gal <- readGAlignments(bamPath, param=parameters, use.names=TRUE)
  
  if (length(gal) >= 1){ # minimum bam region coverage 1 read
    
    for (read_num in 1:length(gal)){
      
      #print(read_num)
      
      # get cigar string, check if indel operators are in the read, and decide to do work
      cigar_string <- mcols(gal)$cigar[[read_num]]
      
      indel_operators <- c("I","D")
      pattern = paste(indel_operators, collapse="|")
      
      if (grepl(pattern, cigar_string)==T){
        
        #message(paste("This read has indel operators:",cigar_string))
        
        # if indel is present in cigar string, define read other bam fields
        query_sequence_string <- mcols(gal)$seq[[read_num]]
        #read_qual_string <- mcols(gal)$qual[[read_num]]
        read_pos <- mcols(gal)$pos[[read_num]]
        query_read_length <- length(query_sequence_string)
        read_strand <- mcols(gal)$strand[[read_num]]
        read_name_record <- (names(gal)[[read_num]])
        
        if (verbose_arg==T){
          
          # print debugging info
          message("Current read:")
          message(read_name_record)
          message(each_chromosome)
          message(read_pos)
          
        }
        
        # Run function to convert cigar to long format
        exploded_cigar_string <- convert_cigar_to_verbose_string(cigar_string)
        
        # determine how much reference sequence to extract by counting the D's
        number_deletions <- sum(exploded_cigar_string == "D")
        number_insertions <- sum(exploded_cigar_string == "I")
        
        operator_values <- rle(exploded_cigar_string)[[2]]
        operator_lengths <- rle(exploded_cigar_string)[[1]]
        
        if (operator_values[1] == "S"){
          number_leading_softclips <- operator_lengths[1]
        } else {
          number_leading_softclips <- 0
        }
        
        # get reference sequence for alignment region
        reference_sequence <- Views(hg38_genome_chr_subset, start=read_pos, end=read_pos+query_read_length+number_deletions)[[1]]
        
        # refine cigar string 
        refined_cigar_string <- refine_cigar_string(exploded_cigar_string,query_sequence_string,read_pos,query_read_length,each_chromosome,reference_sequence)
        
        # Run algorithm on read and collect indel calls
        per_bam_region_indel_records <- run_algo_on_one_read_explicit_args_with_simple_indel_padding(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,each_chromosome,min_indel_length)
        
      } # end conditional if read has any I or D operators
      
    } # end per read iteration
    
    rm(gal)
    
  } # end coverage > 2 conditional
  
  return(per_bam_region_indel_records)
} # end run algo all reads function


