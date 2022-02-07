# INDELfindR

## About
INDELfindR is an R based command line tool for detecting simple and complex insertion deletion (INDEL) variants which outputs variant calls in a VCF v4.3 file compatible with downstream analysis and annotation tools.

## Setup

#### Requirements

INDELfindR is compatible with Linux, Unix, and MacOS.

INDELfindR requires R version >= 4.1.0

#### Installation

Install the INDELfindR R package from CRAN to install INDELfindR.

```
install.packages("INDELfindR")

# test installation
library(INDELfindR)
```

## Quickstart

After installing the INDELfindR R package, run INDELfindR from the command line by calling:

```
Rscript indelfindr.R -b <indexed_bam.bam> (...)
```

Please see our documentation(hyperlink here) for full use instructions and parameter option definitions.

## Usage Instructions

```
usage: indelfindr.R [-h] [-v] [-q] -a filename [-f number] [-b number]
                    [-l number] [-nr number] [-mq number] [-t filename]
                    [-nc number] [-p] [-vaf number] [-dp number] [-o OUTNAME]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print extra output [default]
  -q, --quietly         Print little output
  -a filename, --alignment_file filename
                        Filepath to alignment file in bam format
  -f number, --flanking_region_length number
                        Minimum number of `=` or `X` CIGAR operators required
                        to flank a string of nearby `I` `D` operators in order
                        for a call to be made
  -b number, --bam_bin_size number
                        Length of non-overlapping sliding windows to use to
                        extract overlapping reads from bam file and store in
                        memory at a time (default = 100000). Larger windows
                        require loading more reads into memory at one time,
                        while smaller windows require longer runtime due to
                        evaluating more reads which overlap neighboring
                        sliding windows and removing more duplicate read calls
  -l number, --min_indel_length number
                        Minumum number of ref or alt allele basepairs for an
                        indel to be called (default = 5 nucleotides)
  -nr number, --number_reads number
                        Minumum number of supporting reads for a single given
                        indel to be called
  -mq number, --mapq_filter number
                        Minumum read mapping quality required to evaluate a
                        read.
  -t filename, --target_regions filename
                        File path to .bed file containing regions in which to
                        perform variant calling. Chromosome syntax must match
                        bed file (ex. `Chr1`, `chr1`, or `1`).
  -nc number, --number_cores number
                        Number of cores to use for fork/join parallel
                        computation during indel calling and filtering.
  -p, --primary_chromosomes
                        Call indels in primary chromosomes only ignoring ALT
                        contig alignments (chr1-22,X,Y,M only).
  -vaf number, --vaf_filter number
                        Minumum variant allele freqency required for an indel
                        to be reported
  -dp number, --read_depth_filter number
                        Minumum indel range read depth required for an indel
                        to be reported
  -o OUTNAME, --outname OUTNAME
                        Define the output directory path
```

## Setting Up the Annotation Databases

Describe Annovar Db prep here


## Installing with Docker
Describe Docker installation here