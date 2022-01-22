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