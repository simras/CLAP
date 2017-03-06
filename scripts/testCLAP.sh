#!/bin/bash
#
#  CLAP - CLIP Analysis Pipeline
# 
#  Simon H. Rasmussen and Mireya Plass
#  Bioinformatics Centre
#  University of Copenhagen
#
######################################################################################
#  Dataset: Reads that map to chr4 from the PAR-CLIP dataset SRR248532
#
######################

BASE=$(pwd)
file7=$BASE"/resources/SRR248532_chr4.fastq"
adapt="TCGTATGCCGTCTTCTGCTTG"

#  Run the pipeline
#  Should be run from the CLAP/ base Directory

scripts/CLAP.sh $file7 $adapt "" 0 1 2 1 1 2 0 test_CLAP 1 8 1
