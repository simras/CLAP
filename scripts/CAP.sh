#!/bin/bash
#
#  CAP - CLIP Analysis Pipeline
# 
#  Simon H. Rasmussen and Mireya Plass
#  Bioinformatics Centre
#  University of Copenhagen
#
####################################################################################################
# COMMENTS:
# Phread 33 option not tested. UPDATE: I wonder if it is or not...
# Adaptor alphabet {ACGTacgt} no other IUPAC like N in the sequence.
# How do we share all of our custom data?
# How do we integrate peak calling?
# What should we do with genome annotation, index and sequences.
# Does our genome have Y-chromosome? it should have to be general.
# I wonder about the fundamental setup in this script. Is it for analyzing one or multiple datasets at a time???
####################################################################################################

trap "" 1
if [ $# -lt 14 ]
  then
    echo "ARGUMENTS:" 
    echo "\$1: Filename"
    echo "\$2: Remove adapters?" 
    echo "    0: No"
    echo "    1: use file"
    echo "    If all datasets have the same 3'adaptor just input the adaptor sequence, ex: ACCTGCA..."
    echo "\$3: Sequence fixed barcode" 
    echo "\$4: length of random barcode"
    echo "\$5: Remove duplicates?"
    echo "    0: No" 
    echo "    1: Yes"  
    echo "\$6: Type of analysis" 
    echo "    1: fixed barcode, random barcodes"
    echo "    2: no fixed barcode, no random barcodes"
    echo "    3: only fixed barcode, no random barcodes"
    echo "    4: no fixed barcode, only random barcodes"
    echo "\$7: UCSC Custom Tracks (bed tracks)"
    echo "    0: No UCSC custom tracks" 
    echo "    1: UCSC custom tracks" 
    echo "\$8: Stranded protocol?"
    echo "    0: Strandless" 
    echo "    1: Stranded" 
    echo "\$9: Index"
    echo "    1: Genome index" 
    echo "    2: Genome index + exon junction index" 
    echo "\$10: Model"
    echo "    0: Model T>C conversions (PAR-CLIP), conversion prob 0.125" 
    echo "    1: No model (RNA-Seq, iCLIP or HITS-CLIP)" 
    echo "\$11: Output name"
    echo "\$12: Quality scores"
    echo "    0: Phread 64" 
    echo "    1: Phread 33" 
    echo "\$13: Number of threads?"
    echo "    Input number of threads"
    echo "\$14: Peak calling?"
    echo "    0: No"
    echo "    1: Yes"
    
    exit 0
fi

set -x
set -e

if [ ${12} -eq 0 ]
then 
    # Trimming option
    Qtype=64  # illumina Phred+64
    #Qtype=65  # illumina Phred+64 BGI

    # BWA option
    phreadOpt=" -I "
else
    # Trimming option
    Qtype=32  # illumina Phred+33

    # BWA option
    phreadOpt=""
fi

if [ ${13} -lt 1 ]
then
    # Number of threads
    threads=1
else
    threads=${13}
fi
   
barcode=$3
trim=$4

# Q = 1 - P, Q is the fraction you expect not to map due to contamination, see Kerpedjiev et al 2014.
P=0.5

PP=0.99
PP_exon=0.99
pid=$$

# Base dir
BASE=$(pwd)
scripts=$BASE"/scripts"
# Configurations
# PSSM location
bwa=$BASE/../bwa-pssm/bwa
pyicos=$BASE/../pyicos/pyicoclip

if [ ${#12} -gt 0 ]
then
    outFolder=$(pwd)/${11}"_"$pid
else
    outFolder=$(pwd)/${f1%.*}"_"$pid
fi

mkdir $outFolder
trap "echo $outFolder" 0
#exec 1> $outFolder/mapping_pipe.log
#exec 2>&1

# Mapping index location
idx1=$BASE"/resources/hg19.fa"
#idx2=$BASE"/indexes/bwa-pssm-exon/exons_ENS_70_collapsed_newIDs_2.fa"
idx3=$BASE"/resources/ensembl70_ej.fa"	    

# Min read length - species dependent
length=20

#path="/home/mplass/miRNA_project"
#seqs="/seqdata/krogh/mplass/miRNA_project/sequences/"

# Treads
T=$threads
# heap size
H=400
# seed length, seeding disabled
L=1024
# Error prob binomial model see original BWA paper
E=0.04

# Filenames 
file_list=$1
FILE=$(basename $file_list)
FILE_NAME=$(dirname $file_list)
#for FILE in $file_list
#do
    
    if [ $6 -eq 2 -o $6 -eq 4 ]
    then
	if [ "$barcode" == "" ]
	then
	    barcode="nobc"
	fi
    fi
    ftmp=$barcode"_assigned_"$FILE
    f2=$barcode"_noADPT_"$FILE
    f3=$barcode"_noDUPS_"$FILE
    final_bed=$barcode"_merged_"${FILE%.*}.bed
    clusters=$barcode"_"${FILE%.*}"_clusters".bed
    
    if [ $6 -eq 1 -o $6 -eq 3 ]
    then
	# Remove barcode and select data of max edit dist 1, offset is at what position the fixed barcode starts XXXTTGTXX_READ_ > pirmer XXACAAXXX
	echo "Remove barcodes and pick out data"
	cd $outFolder
	cat $FILE_NAME/$FILE|$scripts/get_data_from_barcode.py -p $barcode -o 3 -m 1 1> $outFolder/$ftmp
	cd ..
    else
	echo "Data has no fixed barcodes"
	cp $FILE_NAME/$FILE $outFolder/$ftmp
    fi
    echo "adapter length" ${#2}
    if [ ${#2} -gt 0 ]
    then
	# Remove Adaptors
	echo "Remove Adaptors"
	adapterS=$2
	cd $outFolder
	$scripts/adapterHMM_multiPY.py -q $Qtype -s $adapterS -f $outFolder/$ftmp -o $outFolder/$f2 -l $(($length+$trim)) -p $threads
	cd ..
    else
	echo "Don't remove Adaptors"
	mv $outFolder/$ftmp $outFolder/$f2
    fi
    
    # Remove duplicates
    if [ $6 -eq 1 -o $6 -eq 4 ]
    then
	echo "Remove duplicates with random barcodes"
    else
	echo "Remove duplicates without random barcodes"
	trim=0
    fi
    
    if [ $5 -eq 1 ]
    then
	cd $outFolder
	$scripts/duplication_barcode.py -f $outFolder/$f2 -t $trim > $outFolder/$f3
	cd ..
    else
	mv $outFolder/$f2 $outFolder/$f3
    fi
    
    # cut suffix .fastq
   # FILE=${FILE%.*}
    map_reads=${FILE%.*}
    echo "Map"
    if [ $9 -ne 2 ]
    then
	if [ ${10} -eq 0 ]
	then
	    # Map with error model
	    $bwa pssm -n $E -l $L -m $H -t $T -P $P $phreadOpt -G $BASE/cluster/error_model_125.txt $idx1 $FILE_NAME/$map_reads.fastq  > $outFolder/$map_reads.sai
	else
	    $bwa pssm -n $E -l $L -m $H -t $T -P $P $phreadOpt $idx1 $FILE_NAME/$map_reads.fastq  > $outFolder/$map_reads.sai
	fi
	$bwa samse -f $outFolder/$map_reads.sam $idx1 $outFolder/$map_reads.sai $FILE_NAME/$map_reads.fastq
	$scripts/mappingStats.sh $outFolder/$map_reads.sam $PP y > $outFolder/$map_reads.mappingstats.txt
    fi
    if [ $9 -eq 2 ]
    then
	# Exon-junction mapping
	map_exon=$map_reads.unmapped
	$BASE/scripts/select_unmapped_reads.pl -f $FILE_NAME/$map_reads.fastq -s $outFolder/$map_reads.sam  > $outFolder/$map_exon.fastq
    fi
    if [ ${10} -eq 0 ]
    then
	nice $p
	$bwa pssm -n $E -l $L -m $H -t $T -P $P $phreadOpt -G $BASE/cluster/error_model_125.txt $idx3 $outFolder/$map_exon.fastq > $outFolder/$map_exon.sai
    else
	nice $bwa pssm -n $E -l $L -m $H -t $T -P $P $phreadOpt $idx3 $outFolder/$map_exon.fastq > $outFolder/$map_exon.sai
    fi
    #    nice $pssm pssm -n $E -l $L -m $H -t $T $dindex $outFolder/$FILE.unmapped.fastq > $outFolder/$FILE.unmapped.sai
    $bwa samse -f $outFolder/$map_exon.sam $idx3 $outFolder/$map_exon.sai $outFolder/$map_exon.fastq
    $scripts/mappingStats.sh $outFolder/$map_exon.sam $PP y > $outFolder/$map_exon.mappingstats.txt
done

# Make bed-files 
# Genome
$BASE/scripts/parse_sam_files_final.pl -p $idx1 -t pssm -c $PP -f $outFolder/$map_reads.sam -o $outFolder/$map_reads.bed -a

# Exon junctions
$BASE/scripts/parse_sam_files_ej_final.pl -p $idx3 -t pssm -c $PP -f $outFolder/$map_exon.sam -o $outFolder/$map_exon.bed -a

# peak Calling
exon_annot=$BASE"/resources/ensembl70.long.exons"
exon_annot="/seqdata/krogh/mplass/RBPs/FDR/ensembl70.long.exons"
peak_in=.bed
peak_out=.pk
if [ ${14} -eq 0 ]
then
    
    cp $peak_in $peak_out
else
    if [ $8 -eq 0 ]
    then
	$pyicos $peak_in $peak_out -f bed -F bed_pk --region $exon_annot --stranded
    else
	# not stranded
	$pyicos $peak_in $peak_out -f bed -F bed_pk --region $exon_annot
    fi
#    # not stranded
#    $pyicos RBP.all.bed RBP.pk -f bed -F bed_pk --region $exon_annot

fi
if [ $9 -eq 1 ]
then
    # rename genome reads file
    cat $outFolder/$map_reads.bed  > $outFolder/$final_bed
fi

if [ $9 -eq 2 ]
then
    # Merge reads from genome and exon junctions
    cat $outFolder/$map_reads.bed $outFolder/$map_exon.bed > $outFolder/$final_bed
fi

if [ $7 -eq 1 ]
then
    #    echo "Printing UCSC custom tracks..."
    cat $clusters |sort -s -k6,6 -k1,1 -k2,2g -k3,3g | $scripts/cluster_reads.py -2 -u -w - 1> /dev/null 2> "UCSC_-_"${f1%.*}".track"
    cat $clusters |sort -s -k6,6 -k1,1 -k2,2g -k3,3g | $scripts/cluster_reads.py -2 -u -w + 1> /dev/null 2> "UCSC_+_"${f1%.*}".track"
    if [ $8 -eq 1 ]
    then
	echo "Printing UCSC custom tracks for stranded protocol..."
	cat $outFolder/$final_bed|sort -s -k6,6 -k1,1 -k2,2g -k3,3g|$scripts/cluster_reads.py -2 1> $outFolder/$clusters
	cat $outFolder/$clusters | /home/mplass/programs/pk2bedGraph_info.pl -c 1,2,7,6 -s - -n "-_strand_"$name -d "Clusters on minus-strand" > $outFolder/"UCSC_m_"${f1%.*}".track"
	cat $outFolder/$clusters | /home/mplass/programs/pk2bedGraph_info.pl -c 1,2,7,6 -s + -n "+_strand_"$name -d  "Clusters on plus-strand" > $outFolder/"UCSC_p_"${f1%.*}".track"
	
	gzip $outFolder/"UCSC_m_"${f1%.*}".track"
	gzip $outFolder/"UCSC_p_"${f1%.*}".track"
    else
	echo "Printing UCSC custom track for strandless protocol..."
	cat $outFolder/$final_bed |sort -s -k1,1 -k2,2g -k3,3g | $scripts/cluster_reads.py -s "+" -2 1> $outFolder/$clusters
	cat $outFolder/$clusters | /home/mplass/programs/pk2bedGraph_info.pl -c 1,2,7,6 -s "+" -n $name -d "Clusters on both strands" > $outFolder/"UCSC_mp_"${f1%.*}".track"
	
	gzip $outFolder/"UCSC_-+_"${f1%.*}".track"
    fi
else
    echo "Do not print UCSC custom tracks..."
    cat $outFolder/$final_bed |sort -s -k6,6 -k1,1 -k2,2g -k3,3g | $scripts/cluster_reads.py -2 1> $outFolder/$clusters
fi

#done
echo "Files are in:"
echo $outFolder
