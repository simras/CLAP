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
file5=$BASE"/data/test_IMP_iCLIP_5.fastq"
adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
scripts/CAP.sh $file5 $adapt GGTT 5 1 1 1 1 2 1 test_iCLIP 1 8 1


