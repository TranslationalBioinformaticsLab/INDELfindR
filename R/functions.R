# local source path: /Users/George/indel_detection_tool_project/scripts/indel_detection_tool_functions.R

#############################################################################################
#
# Store Functions for Main Rscript
#
#############################################################################################

# Function 1) Convert bam to long format

#' convert_cigar_to_verbose_string
#' @importFrom stringr str_split
#' @param each_cigar_string each_cigar_string
#'
#' @return extended_string
#' @export
#'
convert_cigar_to_verbose_string <- function(each_cigar_string){

  cigar_operators <- unlist(str_split(each_cigar_string, "(\\d+)"))[-1]
  number_op_reps <- utils::head(unlist(str_split(each_cigar_string, "[a-zA-Z]+")),-1)

  extended_string <- c()

  for (number_operators in 1:length(cigar_operators)){

    operator_reps <- rep(cigar_operators[number_operators],number_op_reps[number_operators])
    extended_string <- append(extended_string,operator_reps)

  }

  return(extended_string)
}

# Function 2.1) Improved Function to Refine M operators as match or mismatch using reference genome hg38

#' refine_cigar_string
#'
#' @param exploded_cigar_string exploded_cigar_string
#' @param query_sequence_string query_sequence_string
#' @param read_pos read_pos
#' @param query_read_length query_read_length
#' @param each_chromosome each_chromosome
#' @param reference_sequence_test reference_sequence
#' @param number_leading_softclips  number_leading_softclips
#' @param number_leading_hardclips number_leading_hardclips
#'
#' @return exploded_cigar_string_for_refinement
#' @export
#'
refine_cigar_string <- function(exploded_cigar_string,query_sequence_string,read_pos,query_read_length,each_chromosome,reference_sequence_test,number_leading_softclips,number_leading_hardclips){

  reference_sequence <- reference_sequence_test[-1]
  exploded_cigar_string_for_refinement <- exploded_cigar_string

  reference_base_counter <- number_leading_softclips
  query_base_counter <- number_leading_softclips

  reference_sequence_testing <- unlist(str_split(toString(reference_sequence),""))
  query_sequence_string_testing <- unlist(str_split(toString(query_sequence_string),""))

  start_of_search <- number_leading_softclips+number_leading_hardclips+1

  for (each_operator in start_of_search:length(exploded_cigar_string_for_refinement)){

    if (exploded_cigar_string_for_refinement[[each_operator]]=="M"){

      reference_base_counter <- reference_base_counter + 1
      query_base_counter <- query_base_counter + 1
      #print("tick1")

      if (reference_sequence_testing[reference_base_counter] == query_sequence_string_testing[query_base_counter]){
        #print("match")
        exploded_cigar_string_for_refinement[[each_operator]] = "="
        #print("tick2")

      } else {
        #print("mismatch")
        exploded_cigar_string_for_refinement[[each_operator]] = "X"
        #print("tick3")
      }

    } else if (exploded_cigar_string_for_refinement[[each_operator]] == "I"){
      query_base_counter <- query_base_counter + 1

      #print("tick4")

    } else if (exploded_cigar_string_for_refinement[[each_operator]] == "D"){
      reference_base_counter <- reference_base_counter + 1
      #print("QUERY not added")

    }

    # Debug view:
    # print(paste("Cigar operator:",exploded_cigar_string_for_refinement[each_operator]))
    # print(paste("Ref counter:",reference_base_counter))
    # print(paste("Query counter:",query_base_counter))
    # print(paste("Ref base:",reference_sequence_testing[reference_base_counter]))
    # print(paste("Query base:",query_sequence_string_testing[query_base_counter]))

  }

  return(exploded_cigar_string_for_refinement)
}

# Function 3.1: Whole_read translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

#' translate_cigar_index_to_ref_and_query_whole_read
#'
#' @param reference_sequence_for_translation reference_sequence_for_translation
#' @param query_sequence_string_for_translation query_sequence_string_for_translation
#' @param refined_cigar_string refined_cigar_string
#'
#' @return list(reference_allele_record,alternate_allele_record,cigar_string_record)
#' @export
#'
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
  alternate_allele_record <- query_sequence_string_for_translation

  cigar_string_record <- paste(refined_cigar_string,collapse=",")

  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.

}

# Function 3.2: Sliding, algo dev: Translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

#' translate_cigar_index_to_ref_and_query_algo_dev_sliding
#'
#' @param sliding_window_per_read_start sliding_window_per_read_start
#' @param sliding_window_per_read_end sliding_window_per_read_end
#' @param reference_sequence_for_translation reference_sequence_for_translation
#' @param query_sequence_string_for_translation query_sequence_string_for_translation
#' @param refined_cigar_string refined_cigar_string
#'
#' @return reference_allele_record alternate_allele_record cigar_string_record
#' @export
#'
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

  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")

  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.

}

# Function 3.3: Translate coordinates of windows with significant concentration of variant cigar operators to significant cigar string and ref, alt alleles

#' translate_cigar_index_to_ref_and_query_v2
#'
#' @param cigar_coords cigar_coords
#' @param cigar_coords_for_query cigar_coords_for_query
#' @param reference_sequence_for_translation reference_sequence_for_translation
#' @param query_sequence_string_for_translation query_sequence_string_for_translation
#' @param refined_cigar_string refined_cigar_string
#'
#' @return list(reference_allele_record,alternate_allele_record,cigar_string_record)
#' @export
#'
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

  alternate_allele_record <- query_sequence_string_for_translation

  if (cigar_coords_for_query_end < cigar_coords_for_query_start){
    alternate_allele_record <- alternate_allele_record[1]
  }

 cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")

  return(list(reference_allele_record,alternate_allele_record,cigar_string_record))
  # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.

}

# Function 4: Convert exploded CIGAR string to normal shortened format CIGAR string

#' unexplode_cigar_string
#'
#' @param exploded_cigar_string exploded_cigar_string
#'
#' @return cigar_record
#' @export
#'
unexplode_cigar_string <- function(exploded_cigar_string){

  #use run length encoding to summarize exploded cigar string
  operator_lengths <- base::rle(exploded_cigar_string)[[1]]
  operator_values <- base::rle(exploded_cigar_string)[[2]]

  cigar_record <- ""
  for (operator_num in 1:length(operator_values)){
    opertator_substring <- stringr::str_c(operator_lengths[operator_num],operator_values[operator_num])
    cigar_record <- stringr::str_c(cigar_record,opertator_substring)
  }

  return(cigar_record)
}

# Function 5: Get Indel Read Depth

#' getIndelCoverage
#'
#' @param indel_chr indel_chr
#' @param indel_start indel_start
#' @param indel_end indel_end
#' @param bamPath bamPath
#'
#' @return indel_read_depth
#' @export
#'
getIndelCoverage <- function(indel_chr,indel_start,indel_end,bamPath){

  indel_read_depth <- GenomicRanges::countOverlaps(GenomicRanges::GRanges(paste0(indel_chr,":",indel_start,"-",indel_end)),bamfile <- GenomicAlignments::readGAlignments(bamPath),minoverlap = as.numeric((indel_end-indel_start)+1))

  return(indel_read_depth)
}

# Function 6: Get Variant Allele Frequency

#' getVaf
#'
#' @param num_supporting_reads num_supporting_reads
#' @param indel_read_depth indel_read_depth
#'
#' @return vaf
#' @export
#'
getVaf <- function(num_supporting_reads,indel_read_depth){

  vaf <- num_supporting_reads/indel_read_depth
  return(vaf)

}

# Function 7

#' calculate_strand_bias_pval_vaf_and_dp_parallel
#' @import dplyr
#' @param collapsed_read_counts collapsed_read_counts
#' @param master_indel_record_table_no_dup_reads master_indel_record_table_no_dup_reads
#' @param number_cores number_cores
#' @param mapq_threshold mapq_threshold
#' @param bamPath bamPath
#' @return collapsed_read_counts_with_strand_correction_pvalue_for_writing
#' @export
#'
calculate_strand_bias_pval_vaf_and_dp_parallel <- function(collapsed_read_counts,master_indel_record_table_no_dup_reads,number_cores,mapq_threshold,bamPath) {

  cols_to_add <- c("VAF","DP","strand_bias_pval","ref_forward","ref_rev","alt_forward","alt_rev")  %>% purrr::map_dfc(stats::setNames, object = list(as.character()))

  cols_to_add <-
    do.call(
      rbind, bettermc::mclapply(1:nrow(collapsed_read_counts),calculate_strand_bias_pval_vaf_and_dp_parallel_each_indel,collapsed_read_counts,master_indel_record_table_no_dup_reads,mapq_threshold,bamPath,mc.cores=number_cores,mc.preschedule = F))

  colnames(cols_to_add) <- c("VAF","DP","strand_bias_pval","ref_forward","ref_rev","alt_forward","alt_rev")

  collapsed_read_counts_with_strand_correction_pvalue_for_writing <- cbind(collapsed_read_counts,cols_to_add)

  return(collapsed_read_counts_with_strand_correction_pvalue_for_writing)
  wait()

}

# Function 7.1: Get strand bias pval, vaf, and dp for each indel (used by the parallelized function calculate_strand_bias_pval_vaf_and_dp_parallel())

#' calculate_strand_bias_pval_vaf_and_dp_parallel_each_indel
#'
#' @param each_unique_indel_candidate each_unique_indel_candidate
#' @param collapsed_read_counts collapsed_read_counts
#' @param master_indel_record_table_no_dup_reads master_indel_record_table_no_dup_reads
#' @param mapq_threshold mapq_threshold
#' @param bamPath bamPath
#'
#' @return row_to_add_to_cols_to_add
#' @export
#'
calculate_strand_bias_pval_vaf_and_dp_parallel_each_indel <- function(each_unique_indel_candidate,collapsed_read_counts,master_indel_record_table_no_dup_reads,mapq_threshold,bamPath){

  each_unique_indel_record <- collapsed_read_counts[each_unique_indel_candidate,]
  all_supporting_read_matches <- plyr::match_df(master_indel_record_table_no_dup_reads,each_unique_indel_record)
  indel_reads_vector <- all_supporting_read_matches[,c("read_name","strand")]

  variant_strand_counts <- table(indel_reads_vector$strand)

  strand_check_chr <- each_unique_indel_record$chr
  strand_check_start_pos <- as.numeric(each_unique_indel_record$start_pos)
  strand_check_end_pos<- as.numeric(each_unique_indel_record$end_pos)

  # set params for retrieving reads to extract non-supporting read strandedness counts
  strand_check_params=Rsamtools::ScanBamParam(simpleCigar=FALSE,
                                   which=GenomicRanges::GRanges(seqnames=S4Vectors::Rle(strand_check_chr),
                                                 ranges=IRanges::IRanges(strand_check_start_pos, strand_check_end_pos)),
                                   mapqFilter=mapq_threshold,
                                   what=c("rname","strand"), # make mapqFilter an arg to parse
                                   Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE,
                                               isNotPassingQualityControls=FALSE,
                                               isDuplicate=FALSE, #user would need bam with marked duplicates for this to be useful. Can test with BAM with marked duplicates. This could cause issues with amplicon seq
                                               isUnmappedQuery=FALSE))

  gal_read_counts <- GenomicAlignments::readGAlignments(bamPath, param=strand_check_params, use.names=TRUE)
  non_supporting_reads <- gal_read_counts[!(rownames(S4Vectors::mcols(gal_read_counts)) %in% indel_reads_vector$read_name),]

  reference_strand_counts <- table(S4Vectors::mcols(non_supporting_reads)$strand)

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

  pvalue <- stats::fisher.test(contingency_table)$p.value

  dp_value <- length(gal_read_counts)
  vaf_value <- (contingency_table[2,1]+contingency_table[2,2])/dp_value

  row_to_add_to_cols_to_add <- c(vaf_value,dp_value,pvalue,ref_forward,ref_rev,alt_forward,alt_rev)

  return(row_to_add_to_cols_to_add)
}

# Function 8: Write indel calls to VCF file output

#' write.vcf
#'
#' @param master_indel_call_table master_indel_call_table
#' @param phased phased
#' @param output_filename output_filename
#'
#' @export
#'
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
  utils::write.table(vcf_data_lines_sorted, file=file,
              quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t", na=paste("."))

}

# Function 9: Get intervals with last interval included when it is not within the by interval

#' Title
#'
#' @param from from
#' @param to to
#' @param by by
#'
#' @return vec to
#' @return to
#' @export
#'
seq_with_uneven_last <- function (from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( utils::tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# Function 12: Define C function to wait on child processes to get rid of zombie processes

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- inline::cfunction(body=code, includes=includes, convention='.C')

# Function 13: Get primary chromosomes in bam file
#' get_primary_chroms
#'
#' @param chr_in_bam chr_in_bam
#'
#' @return chr_in_bam
#' @export
#'
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

#' translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding
#'
#' @param cigar_coords cigar_coords
#' @param cigar_coords_for_query cigar_coords_for_query
#' @param reference_sequence_for_translation reference_sequence_for_translation
#' @param query_sequence_string_for_translation query_sequence_string_for_translation
#' @param refined_cigar_string refined_cigar_string
#'
#' @return reference_allele_record alternate_allele_record cigar_string_record
#' @export
#'
translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding <- function(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string){

  cigar_coord_start <- min(cigar_coords)
  cigar_coord_end <- max(cigar_coords)

  cigar_coords_for_query_end <- cigar_coords_for_query[length(cigar_coords_for_query)]
  cigar_coords_for_query_start <- cigar_coords_for_query[1]

  exploded_cigar_string_record <- refined_cigar_string[cigar_coord_start:cigar_coord_end]

  D_indices <- which(exploded_cigar_string_record=="D")+1
  I_indices <- which(exploded_cigar_string_record=="I")+1
  equals_indices <- which(exploded_cigar_string_record=="=")+1
  X_indices <- which(exploded_cigar_string_record=="X")+1

  # define Reference seq:
  # build reference sequence retrieval coords with D, X, and =
  reference_indices_pre <- sort(c(1,D_indices,X_indices,equals_indices))
  number_ref_bases_to_get <- length(reference_indices_pre)
  reference_indices <- c(1:number_ref_bases_to_get)
  reference_allele_record <- reference_sequence_for_translation[reference_indices]

  if (length(reference_allele_record) == 0){
    reference_allele_record = reference_sequence_for_translation[1]
  }

  #define Alternate seq:
  # Build query read sequence retrieval coords to with I, X, and =
  query_sequence_string_for_translation[1] <- reference_sequence_for_translation[1]
  alternate_allele_record <- query_sequence_string_for_translation
  cigar_string_record <- paste(refined_cigar_string[cigar_coord_start:cigar_coord_end],collapse=",")

  return(list(reference_allele_record,alternate_allele_record,cigar_string_record)) # to access returned objects, use "reference_allele_record <- unlist(results)[[1]]", etc.

}

# Function 15: Run algorithm on each read and collect indel calls to add to master table. Version for adding padding to simple indels.

#' run_algo_on_one_read_explicit_args_with_simple_indel_padding
#'
#' @param per_bam_region_indel_records per_bam_region_indel_records
#' @param refined_cigar_string refined_cigar_string
#' @param flanking_region_length flanking_region_length
#' @param query_sequence_string query_sequence_string
#' @param reference_sequence reference_sequence
#' @param read_pos read_pos
#' @param query_read_length query_read_length
#' @param read_strand read_strand
#' @param read_name_record read_name_record
#' @param number_leading_softclips number_leading_softclips
#' @param number_leading_hardclips number_leading_hardclips
#' @param each_chromosome each_chromosome
#' @param min_indel_length min_indel_length
#'
#' @return per_bam_region_indel_records
#' @export
#'
run_algo_on_one_read_explicit_args_with_simple_indel_padding <- function(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,number_leading_hardclips,each_chromosome,min_indel_length){

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
    if (operator == "S" | operator == "H"){
      operator <- "="
    }

    # if the operator matches the reference and there is not a current indel candidate
    if ((operator == "=" | operator == "X" ) && consecutive_indel_operator_flag == F){

      match_operator_counter = match_operator_counter + 1
      #print("ONE")

      # if the first variant operator is encountered:
    } else if (operator %in% c("I","D") && match_operator_counter >= flanking_region_length && consecutive_indel_operator_flag == F){
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO")

      # Catch indels which begin at very start of read and dont have 10bp flanking prior. (If the first variant operator is encountered before the 10th bp in the read)
    } else if (operator %in% c("I","D") && match_operator_counter < flanking_region_length && consecutive_indel_operator_flag == F) {
      match_operator_counter <- 0
      indel_candidate_container <- c(indel_candidate_container,operator)
      consecutive_indel_operator_flag <- T
      #print("TWO.5")

      # If there is an immediately adjacent variant operator
    } else if (operator %in% c("I","D","X") && (match_operator_counter < flanking_region_length) && consecutive_indel_operator_flag == T){

      indel_candidate_container <- c(indel_candidate_container,operator)

      match_operator_counter <- 0

      #print("THREE")

      # count match operators after encountering variant operators
    } else if (operator == "=" && consecutive_indel_operator_flag == T){

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
        cigar_end_for_query <- each_operator-flanking_region_length-number_deletions_encountered-number_leading_hardclips
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)-1 # dev change add -1
        cigar_coords_for_query <- cigar_start:cigar_end_for_query

        reference_start_record <- read_pos+cigar_start-number_leading_softclips-number_leading_hardclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container) # dev change (remove)
        chr_record <- each_chromosome

        cigar_start_for_reference <- cigar_start-number_leading_softclips-number_leading_hardclips
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container) # dev change remove

        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] #can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)

        candidate_indel_record <-c(chr=chr_record,
                                   start_pos=reference_start_record,
                                   end_pos=reference_end_record,
                                   refined_cigar_string=toString(exploded_cigar_string_record),
                                   collapsed_cigar_string=cigar_string_record,
                                   reference_allele=reference_allele_record,
                                   alt_allele=alternate_allele_record,
                                   strand = as.character(read_strand),
                                   read_name = read_name_record)

        #add indel record to per region table
        per_bam_region_indel_records <- bind_rows(per_bam_region_indel_records,candidate_indel_record)

        # clear the indel_candidate_container
        indel_candidate_container=c() # dev off temp

      } else {
        # don't save indel and keep moving on nothing
        indel_candidate_container=c() # dev off temp
      } # end conditional min_indel_length conditional


      # add exception if the indel candidate runs all the way into the end of the read
    } else if (each_operator == length(refined_cigar_string) && consecutive_indel_operator_flag == T){

      # reset the consecutive_indel_operator_flag to false, since string of closely spaced operators is broken
      consecutive_indel_operator_flag = F

     # print("SIX")
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

        cigar_end_for_query <- each_operator-num_to_remove-number_deletions_encountered-number_leading_hardclips
        cigar_start_for_query <- cigar_end_for_query-(length(indel_candidate_container)-1-number_deletions_per_candidate)-1 # dev change
        cigar_coords_for_query<- cigar_start:cigar_end_for_query

        # Define indel records and add indel candidate record to per bam region table
        reference_start_record <- read_pos+cigar_start-number_leading_softclips-number_leading_hardclips-2
        reference_end_record <- reference_start_record+length(indel_candidate_container) #-1 # dev change

        chr_record <- each_chromosome
        cigar_start_for_reference <- cigar_start-number_leading_softclips-number_leading_hardclips
        cigar_end_for_reference  <- cigar_start_for_reference + length(indel_candidate_container) #-1 # dev change

        reference_sequence_for_translation <- reference_sequence[cigar_start_for_reference:cigar_end_for_reference]
        query_sequence_string_for_translation <- query_sequence_string[cigar_start_for_query:cigar_end_for_query]
        indel_record_results_list <- translate_cigar_index_to_ref_and_query_v2_with_simple_indel_padding(cigar_coords,cigar_coords_for_query,reference_sequence_for_translation,query_sequence_string_for_translation,refined_cigar_string)
        reference_allele_record <- toString(unlist(indel_record_results_list)[[1]])
        alternate_allele_record <- toString(unlist(indel_record_results_list)[[2]])
        exploded_cigar_string_record <- str_split(unlist(indel_record_results_list)[[3]],",")[[1]] # can convert this to regular condensed format cigar string
        cigar_string_record <- unexplode_cigar_string(exploded_cigar_string_record)

        candidate_indel_record <-c(chr=chr_record,
                                   start_pos=reference_start_record,
                                   end_pos=reference_end_record,
                                   refined_cigar_string=toString(exploded_cigar_string_record),
                                   collapsed_cigar_string=cigar_string_record,
                                   reference_allele=reference_allele_record,
                                   alt_allele=alternate_allele_record,
                                   strand = as.character(read_strand),
                                   read_name = read_name_record)


        #add indel record to per region table
        per_bam_region_indel_records <- bind_rows(per_bam_region_indel_records,candidate_indel_record)

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

#' run_algo_all_reads_each_bam_region_with_simple_indel_padding
#'
#' @param row_num row_num
#' @param per_bam_region_indel_records per_bam_region_indel_records
#' @param mapq_threshold mapq_threshold
#' @param bamPath bamPath
#' @param each_chromosome each_chromosome
#' @param sliding_windows_per_bam_region sliding_windows_per_bam_region
#' @param verbose_arg verbose_arg
#' @param hg38_genome_chr_subset hg38_genome_chr_subset
#' @param flanking_region_length flanking_region_length
#' @param min_indel_length min_indel_length
#'
#' @return per_bam_region_indel_records
#' @export
#'
run_algo_all_reads_each_bam_region_with_simple_indel_padding <- function(row_num,per_bam_region_indel_records,mapq_threshold,bamPath,each_chromosome,sliding_windows_per_bam_region,verbose_arg,hg38_genome_chr_subset,flanking_region_length,min_indel_length){

  #print(row_num)
  bam_region_number <- row_num
  bam_region_chr <- each_chromosome
  bam_region_start <- sliding_windows_per_bam_region$start[row_num]
  bam_region_end <-  sliding_windows_per_bam_region$end[row_num]

  # apply filters:
  # 1) isSecondaryAlignment=FALSE
  # 2) isNotPassingQualityControls=FALSE
  # 3) isDuplicate=FALSE
  # 4) isUnmappedQuery=FALSE

  parameters=Rsamtools::ScanBamParam(simpleCigar=FALSE,
                          which=GenomicRanges::GRanges(S4Vectors::Rle(c(sliding_windows_per_bam_region[row_num,1])),
                                        IRanges::IRanges(c(sliding_windows_per_bam_region[row_num,2]), c(sliding_windows_per_bam_region[row_num,3]))),
                          mapqFilter=mapq_threshold,
                          what=c("rname","pos","strand","cigar","seq"), # make mapqFilter an arg to parse
                          #what=c("rname","pos","strand","cigar","seq","qual"), # make mapqFilter an arg to parse
                          Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE,
                                      isNotPassingQualityControls=FALSE,
                                      isDuplicate=FALSE, #user would need bam with marked duplicates for this to be useful. Can test with BAM with marked duplicates. This could cause issues with amplicon seq
                                      isUnmappedQuery=FALSE))

  gal <- GenomicAlignments::readGAlignments(bamPath, param=parameters, use.names=TRUE)

    indel_operators <- c("I","D")
    pattern = paste(indel_operators, collapse="|")

    gal_with_indels <- gal[which(grepl(pattern,S4Vectors::mcols(gal)$cigar)==T),]

    if (length(gal_with_indels) >= 1){

    for (read_num in 1:length(gal_with_indels)){

     #print(read_num) #debugging line

       cigar_string <- gal_with_indels@cigar[[read_num]] # improved replacement
       query_sequence_string <- gal_with_indels@elementMetadata@listData$seq[[read_num]]# improved replacement
       read_pos <- gal_with_indels@elementMetadata@listData$pos[[read_num]]#  improved replacement
       query_read_length <- length(query_sequence_string)

        read_strand <- gal_with_indels@strand[read_num]@values # improved replacement
        read_name_record <- gal_with_indels@NAMES[read_num]

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

        operator_values <- base::rle(exploded_cigar_string)[[2]]
        operator_lengths <- base::rle(exploded_cigar_string)[[1]]

        if (operator_values[1] == "S"){
          number_leading_softclips <- operator_lengths[1]
          number_leading_hardclips <- 0
        } else if (operator_values[1] == "H"){
          number_leading_hardclips <- operator_lengths[1]
          number_leading_softclips <- 0
        } else {
          number_leading_softclips <- 0
          number_leading_hardclips <- 0
        }

        # if (operator_values[length(operator_values)] == "S"){
        #   number_tailing_softclips <- operator_lengths[length(operator_values)]
        #   number_tailing_hardclips <- 0
        # } else if (operator_values[length(operator_values)] == "H"){
        #   number_tailing_hardclips <- operator_lengths[length(operator_values)]
        #   number_tailing_softclips <- 0
        # } else {
        #   number_tailing_softclips <- 0
        #   number_tailing_hardclips <- 0
        # }

        # get reference sequence for alignment region
        reference_sequence <- IRanges::Views(hg38_genome_chr_subset, start=read_pos-1, end=read_pos+query_read_length+number_deletions)[[1]]

        # refine cigar string
        refined_cigar_string <- refine_cigar_string(exploded_cigar_string,query_sequence_string,read_pos,query_read_length,each_chromosome,reference_sequence,number_leading_softclips,number_leading_hardclips)

        # Run algorithm on read and collect indel calls
        per_bam_region_indel_records <- run_algo_on_one_read_explicit_args_with_simple_indel_padding(per_bam_region_indel_records,refined_cigar_string,flanking_region_length,query_sequence_string,reference_sequence,read_pos,query_read_length,read_strand,read_name_record,number_leading_softclips,number_leading_hardclips,each_chromosome,min_indel_length)

    } # end per read iteration

    rm(gal_with_indels)
    rm(gal)

  } # end coverage > 2 conditional

  return(per_bam_region_indel_records)
} # end run algo all reads function

# Function 17: Filter read records for all indels calls per chromosome (rather than all at once at the end of the script, causing a bottleneck on large bam files)

#' filter_and_annotate_calls
#'
#' @param master_indel_record_table master_indel_record_table
#' @param mapq_threshold mapq_threshold
#' @param bamPath bamPath
#' @param number_cores number_cores
#' @param min_supporting_reads min_supporting_reads
#' @param min_vaf min_vaf
#' @param min_read_depth min_read_depth
#' @param each_chromosome each_chromosome
#'
#' @return collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered
#' @export
#'
filter_and_annotate_calls <- function(master_indel_record_table,mapq_threshold,bamPath,number_cores,min_supporting_reads,min_vaf,min_read_depth,each_chromosome){

  if (nrow(master_indel_record_table) != 0){

    dup_indices <- duplicated(master_indel_record_table[,c("read_name","chr","start_pos","end_pos","alt_allele")])
    master_indel_record_table_no_dup_reads <- master_indel_record_table[!dup_indices,]

    # get supporting read counts for each unique indel candidate
    collapsed_read_counts_with_strand <- rename(dplyr::count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele,strand), FREQ = n)

    collapsed_read_counts <- rename(count(master_indel_record_table_no_dup_reads, chr,start_pos,end_pos,refined_cigar_string,collapsed_cigar_string,reference_allele,alt_allele), FREQ = n)

    # collapsed_read_counts_with_strand_correction_pvalue_for_writing <- calculate_stand_bias_pval_vaf_and_dp(collapsed_read_counts,master_indel_record_table_no_dup_reads)
    collapsed_read_counts_with_strand_correction_pvalue_for_writing <- suppressMessages(calculate_strand_bias_pval_vaf_and_dp_parallel(collapsed_read_counts,master_indel_record_table_no_dup_reads,number_cores,mapq_threshold,bamPath))

    # filter by minimum number of supporting reads
    collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered <- collapsed_read_counts_with_strand_correction_pvalue_for_writing[which(collapsed_read_counts_with_strand_correction_pvalue_for_writing$FREQ>min_supporting_reads & collapsed_read_counts_with_strand_correction_pvalue_for_writing$VAF>=min_vaf  & collapsed_read_counts_with_strand_correction_pvalue_for_writing$DP >= min_read_depth),]

    # adjust start_pos to reflect padding base
    collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos <- as.integer(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$end_pos)
    collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos <- as.integer(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered$start_pos) # no longer needed after ref and query fix # - 1

    message(paste0("Found ",nrow(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered), " indels in chromosome: ",each_chromosome,"."))

    return(collapsed_read_counts_with_strand_correction_pvalue_for_writing_filtered)
  }

}


#' run_indelfindr
#'
#' @param bamPath bamPath
#' @param bam_region_bin_size bam_region_bin_size
#' @param verbose_arg verbose_arg
#' @param flanking_region_length flanking_region_length
#' @param target_regions target_regions
#' @param number_cores number_cores
#' @param primary_chromosomes primary_chromosomes
#' @param min_indel_length min_indel_length
#' @param mapq_threshold mapq_threshold
#' @param min_supporting_reads min_supporting_reads
#' @param min_vaf min_vaf
#' @param min_read_depth min_read_depth
#' @param outname outname
#'
#' @export
#'
run_indelfindr <- function(bamPath,bam_region_bin_size,verbose_arg,flanking_region_length,target_regions,number_cores,primary_chromosomes,min_indel_length,mapq_threshold,min_supporting_reads,min_vaf,min_read_depth,outname){

  options(scipen=20)

  # get all hg38 chromosome lengths
  chg38_chromosome_lengths <- GenomeInfoDb::getChromInfoFromUCSC("hg38") # all hg38 chr/altcontig lengths  (595 chromosomes/contigs in hg38 annotation)

  # load hg38 reference genome sequence data (688MB) - Future dev idea: add option to load this optionally, if user opts to use BBtools Reformat tool, or other aligner, for extended cigar string, it will cut down on runtime and not require loading reference genome sequence into memory but this would require we dev a handler for refined cigar bam input.
  hg38_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  # load hg19 reference genome sequence data (677MB) -
  #hg19_genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  #############################################################################################
  #
  # Running Script:
  #
  #############################################################################################

  # define data table for reporting all refined cigar string frames in sliding windows
  master_indel_record_table <- c("chr","start_pos",
                                 "end_pos",
                                 "refined_cigar_string",
                                 "collapsed_cigar_string",
                                 "reference_allele",
                                 "alt_allele",
                                 "strand",
                                 "read_name")  %>% purrr::map_dfc(stats::setNames, object = list(as.character()))

  # Subset primary chromosomes only for extracting reads for indel calling (optional)
  if (primary_chromosomes == T){
    #chr_in_bam <- get_primary_chroms(chr_in_bam)
    chr_in_bam <-  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",'chr10',"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
  }

  if (target_regions != FALSE) {
    # load target regions
    target_regions_table <- read.table(target_regions,sep="\t",col.names = c("chr","start","stop"))
    chr_in_bam <- unique(target_regions_table$chr)
  }

  for (each_chromosome in chr_in_bam){

    message(paste("Analyzing chromosome:",each_chromosome))

    # Define chr subset reference sequence
    hg38_genome_chr_subset <- hg38_genome[[each_chromosome]]

    if (target_regions != FALSE){

      for (regions_table_row in 1:nrow(subset(target_regions_table,chr=chr_in_bam))){

        min_pos <- subset(target_regions_table,chr=chr_in_bam)$start[regions_table_row]
        max_pos <- subset(target_regions_table,chr=chr_in_bam)$stop[regions_table_row]

        target_bin_size <- bam_region_bin_size

        if (max_pos-min_pos < bam_region_bin_size){
          target_bin_size <- max_pos-min_pos
        }

        intervals <- seq_with_uneven_last(from=min_pos,to=max_pos,by=target_bin_size)

        # define coordinates of sliding windows
        sliding_windows_per_bam_region=data.frame(chr=each_chromosome,
                                                  start=intervals[-(length(intervals))],
                                                  end=intervals[-1])

        per_bam_region_indel_records <- c("chr",
                                          "start_pos",
                                          "end_pos",
                                          "refined_cigar_string",
                                          "collapsed_cigar_string",
                                          "reference_allele",
                                          "alt_allele",
                                          "strand",
                                          "read_name")  %>% purrr::map_dfc(stats::setNames, object = list(as.character()))

        # Find indels in each bam region using all available cores in parallel
        per_bam_region_indel_records <-
          do.call(
            rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mapq_threshold,bamPath,each_chromosome,sliding_windows_per_bam_region,verbose_arg,hg38_genome_chr_subset,flanking_region_length,min_indel_length,mc.cores=number_cores,mc.preschedule = F))

        per_chrom_filtered_calls <- filter_and_annotate_calls(per_bam_region_indel_records,mapq_threshold,bamPath,number_cores,min_supporting_reads,min_vaf,min_read_depth,each_chromosome)

        ## Add each chromosome results to the master results table
        master_indel_record_table <- rbind(master_indel_record_table,per_chrom_filtered_calls)

        wait()

      } # end each targeted region in bed file iteration

    } else if (target_regions == FALSE){ #end targeted region mode conditional

      # do all analysis per chromosome, per bam region chunk alignment region, per read in alignment region
      size_chr <- chg38_chromosome_lengths[which(chg38_chromosome_lengths[,1]==each_chromosome),2]
      target_bin_size <- bam_region_bin_size
      intervals <- seq_with_uneven_last(from=1,to=size_chr,by=target_bin_size)

      # define coordinates of sliding windows
      sliding_windows_per_bam_region=data.frame(chr=each_chromosome,
                                                start=intervals[-(length(intervals))],
                                                end=intervals[-1])

      # Initialize dataframe for collecting each bam region indel calls
      per_bam_region_indel_records <- c("chr",
                                        "start_pos",
                                        "end_pos",
                                        "refined_cigar_string",
                                        "collapsed_cigar_string",
                                        "reference_allele",
                                        "alt_allele",
                                        "strand",
                                        "read_name")  %>% purrr::map_dfc(stats::setNames, object = list(as.character()))

      # Find indels in each bam region using all available cores in parallel
      per_bam_region_indel_records <-
        do.call(
          rbind, bettermc::mclapply(1:nrow(sliding_windows_per_bam_region),run_algo_all_reads_each_bam_region_with_simple_indel_padding,per_bam_region_indel_records,mapq_threshold,bamPath,each_chromosome,sliding_windows_per_bam_region,verbose_arg,hg38_genome_chr_subset,flanking_region_length,min_indel_length,mc.cores=number_cores,mc.preschedule = F))

      # Apply filters to indel call list per chromosome results table
      per_chrom_filtered_calls <- filter_and_annotate_calls(per_bam_region_indel_records,mapq_threshold,bamPath,number_cores,min_supporting_reads,min_vaf,min_read_depth,each_chromosome)

      master_indel_record_table <- rbind(master_indel_record_table,per_chrom_filtered_calls)

      wait()

    } # end genome wide mode conditional

  } # end chr iteration

  message(paste0("Saving ",nrow(master_indel_record_table), " indels to output filepath: ",outname,"."))

  # # make directory
  # suppressWarnings(
  #   "suppressWarnings("In dir.create(file.path("/Users/George/indel_detection_tool_project/test_output_dir/")) :
  #   '/Users/George/indel_detection_tool_project/test_output_dir' already exists")"
  # )

  #dir.create(file.path(outname)) # suppress warnings

  # write out VCF file
  write.vcf(master_indel_record_table,phased=FALSE,paste0(outname,".indelfindr.vcf"))

  #write out data table with cigar string
  write.table(master_indel_record_table,paste0(outname,".indelfindr.cigars.csv"),sep=",",row.names=F)

  message("Run Complete")

}
