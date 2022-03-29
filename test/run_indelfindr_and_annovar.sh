#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=28G
#SBATCH -t 12:00:00

module load R/4.1.0
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
module load annovar

bam_list="any_bam_sample_prefix"
dir="/path/to/your/indexed/bam/file"

for bam_prefix in $bam_list

do

Rscript indelfindr.R --alignment_file $dir/${bam_prefix}.bam --bam_bin_size 10000000 --flanking_region_length 10 --number_cores 4 --primary_chromosomes --min_indel_length 3 --mapq_filter 20 --number_reads 4 --vaf_filter 0.01 --read_depth_filter 10 --outname $dir/variant_calling/${bam_prefix}.indelfindr

table_annovar.pl $dir/variant_calling/${bam_prefix}.indelfindr.vcf /gpfs/data/dgamsiz/Uzun_Lab/gtollefs/annovar/humandb -buildver hg38 -out $dir/variant_calling/${bam_prefix}.indelfindr -remove -protocol refGene,avsnp150,clinvar_20210501,cosmic94_coding -operation g,f,f,f -nastring . -polish -vcfinput

Rscript --results_directory $dir/variant_calling/ --prefix ${bam_prefix}.indelfindr.hg38_multianno

done
