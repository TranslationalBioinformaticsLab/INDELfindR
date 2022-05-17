#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=28G
#SBATCH -t 12:00:00
#SBATCH -N 1

# replace the module names with your modules if using a shared Linux server
module load R/4.1.0
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
module load annovar

vcf_list="EGFR_mutations_indelfindr_demo"
dir="/your-path/indelfindr/demo/"

for vcf_prefix in $vcf_list

do

Rscript indelfindr.R --alignment_file $dir/${vcf_prefix}.bam --target_regions EGFR_regions_of_interest.txt --bam_bin_size 10000000 --flanking_region_length 10 --number_cores 4 --primary_chromosomes --min_indel_length 3 --mapq_filter 20 --number_reads 4 --vaf_filter 0.01 --read_depth_filter 10 --outname $dir/${vcf_prefix}

# To run, un-comment the below ANNOVAR command using databases which are described in the INDELfindR README Annotation with ANNOVAR section.
#table_annovar.pl $dir/${vcf_prefix}.indelfindr.vcf /gpfs/data/dgamsiz/Uzun_Lab/gtollefs/annovar/humandb -buildver hg38 -out $dir/${vcf_prefix}.indelfindr -remove -protocol refGene,avsnp150,clinvar_20210501,cosmic94_coding -operation g,f,f,f -nastring . -polish -vcfinput

done