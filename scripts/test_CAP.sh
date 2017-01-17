#!/bin/bash
#
#  CAP - CLIP Analysis Pipeline
# 
#  Simon H. Rasmussen and Mireya Plass
#  Bioinformatics Centre
#  University of Copenhagen
#
######################################################################################


#Rclip 1 , NNNGGTTNN , iCLIP 1
#Rclip 12, NNNCCACNN , iCLIP 2
#Rclip 13, NNNCGGANN , RNAseq 1
#Rclip 16, NNNTTAANN , RNAseq 2
BASE=$(pwd)
file1=$BASE"/data/test_IMP_iCLIP_1.fastq"
file2=$BASE"/data/test_IMP_iCLIP_2.fastq"
file3=$BASE"/data/test_IMP_iCLIP_3.fastq"
file4=$BASE"/data/test_IMP_iCLIP_4.fastq"
adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
scripts/CAP.sh $file1 $adapt GGTT 5 1 1 1 1 1 1 test_iCLIP 0 1 1


