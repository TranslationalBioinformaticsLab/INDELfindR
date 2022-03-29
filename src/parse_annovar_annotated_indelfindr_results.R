#!/usr/bin/env Rscript

######################################################################
#
# Script to Parse VCF FORMAT info from ANNOVAR output.txt file with 
# INDELfindR results
#
######################################################################

suppressMessages(library(argparse))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option 
parser$add_argument("-d", "--results_directory", type="character",required=T, 
                    help="Filepath to directory containing the annovar output .txt file for parsing. This is only for the .txt output produced by our provided shell script for combined INDELfindR and ANNOVAR analysis.",
                    metavar="filepath")
parser$add_argument("-p", "--prefix", type="character",required=T, 
                    help="Specify sample prefix of annovar output .txt file from combined INDELfindR and ANNOVAR shell script.",
                    metavar="filename")

args <- parser$parse_args()

read_in_results_directory <- args$results_directory
file_prefix <- args$prefix

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(Rsamtools))

############################################################
#
# Define Parsing Function
#
############################################################

parse_and_clean_indelfindr_annovar_text_file_from_vcfinput_run <- function(annovar_inputvcf_textfile_output_filename){
  
  annovar_table_pre <- read.table(annovar_inputvcf_textfile_output_filename,sep="\t",fill=T)
  colnames(annovar_table_pre) <- annovar_table_pre[1,]
  annovar_table <- annovar_table_pre[-1,-c(18:20)]
  
  base::colnames(annovar_table)[18:25]<- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  annovar_table_parsed_pre <- within(annovar_table, INFO<-data.frame(do.call('rbind', strsplit(as.character(INFO), ';', fixed=TRUE))))
  
  annovar_table_parsed_pre$INFO$X1 <- unlist(lapply(strsplit(as.character(annovar_table_parsed_pre$INFO$X1), "="), '[', 2))
  annovar_table_parsed_pre$INFO$X2 <- unlist(lapply(strsplit(as.character(annovar_table_parsed_pre$INFO$X2), "="), '[', 2))
  
  annovar_table_parsed_info_only <- annovar_table_parsed_pre[,25]
  annovar_table_no_parsed_info <- annovar_table_parsed_pre[,c(1:24)]
  annovar_table_no_parsed_info_end <- annovar_table_parsed_pre[,c(26:27)]
  
  base::colnames(annovar_table)[18:27]<- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample1")
  
  annovar_table_parsed_pre_sample <- within(annovar_table, Sample1<-data.frame(do.call('rbind', strsplit(as.character(Sample1), ':', fixed=TRUE))))
  annovar_table_parsed_info_only_sample <- annovar_table_parsed_pre_sample[,"Sample1"]
  #annovar_table_no_parsed_info_sample <- annovar_table_parsed_pre_sample[,-which(colnames(annovar_table_parsed_pre_sample)=="Sample1")]
  parsed_col_headers_sample <- unlist(strsplit(as.character(annovar_table[1,"FORMAT"]),":",fixed=TRUE))
  colnames(annovar_table_parsed_info_only_sample) <- parsed_col_headers_sample
  
  parsed_col_headers <- c("AF","SB")
  colnames(annovar_table_parsed_info_only) <- parsed_col_headers
  annovar_table_parsed <- cbind(annovar_table_no_parsed_info,annovar_table_parsed_info_only,annovar_table_parsed_info_only_sample)
  
  colnames(annovar_table_parsed)[27:30] <- c("Ref_Fwd_Reads","Ref_Rev_Reads","Alt_Fwd_Reads","Alt_Rev_Reads")
  return(annovar_table_parsed)
}

############################################################
#
# INDELfindR ANNOVAR .txt output parsing.
#
############################################################

# Define Filepaths
results_dir <- read_in_results_directory
indelfindr_annovar_file <- file_prefix

# Run parsing function: 
indelfindr_parsed_annovar_table <- parse_and_clean_indelfindr_annovar_text_file_from_vcfinput_run(paste0(results_dir,indelfindr_annovar_file,".txt"))

# Save parsed output
write.table(indelfindr_parsed_annovar_table,paste0(results_dir,file_prefix,".parsed.txt"),sep="\t",row.names = F,col.names = T)
