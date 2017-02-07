#!/bin/bash
#
#  CAP - CLIP Analysis Pipeline
# 
#  Simon H. Rasmussen and Mireya Plass
#  Bioinformatics Centre
#  University of Copenhagen
#
######################################################################################

BASE=$(pwd)
file7=$BASE"/resources/SRR248532_chr4.fastq"
adapt="TCGTATGCCGTCTTCTGCTTG"

# Run the pipeline
scripts/CLAP.sh $file7 $adapt "" 0 1 2 1 1 2 0 test_CLAP 1 8 1
