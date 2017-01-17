#!/usr/bin/perl -w

####################################################################################################
## ATTENTION: reads that are correctly mapped but not split mapped will be discarded ###############
####################################################################################################

use strict;
use lib '/home/mplass/modules/';
use Sequence;
use Alignment;
use General;
use Getopt::Long;

my $file   ="";
my $remove ="";
my $type;
my $fasta;
my $cutoff = 0.5;
my $maxcutoff = 1;
my $bedout = "out.bed";
my $help;
my $no_seq;

&GetOptions('f:s'   => \$file,      ## query file with the sam file
	    't=s'   => \$type,    
	    'p=s'   => \$fasta,     ## path where the chr files are
	    'c:f'   => \$cutoff,    ## cut off to define the uniquely mapped reads
	    'm:f'   => \$maxcutoff, ## max cut off to select reads
	    'o:s'   => \$bedout,    ## bed output file
	    'r:s'   => \$remove,    ## list of IDs to be disregarded from the analysis
	    'a'     => \$no_seq,    ## flag saying if the sequence and the CIGAR alignment should be printed
	    'help'  => \$help,      ## help
	   );

## define options
if (defined($help) || !($fasta||$type)){
    print STDERR "
    USAGE ./parse_sam_files_ej_final.pl -p path 

      -p:     file where to find the sequence of the exon junctions (compulsory)  

      -t:     type of file to anaylize. select pssm or bwa
      
    OPTIONS:
      
      -f:     file containing the alignment of reads in sam format
              Default: reads from STDIN
       
      -c:     posterior probability cut-off to define reads as uniquely mapped
              Default: 0.5

      -m:     maximum posterior probability cut-off to define reads as uniquely mapped
              Default: 1
       
      -r:     file containing the read IDs to be disregarded from the analysis
              Default: none
      
      -a:     print sequences and extended cigar in the stats file
              Default: no printing

      -o:     name of the output bed file
              Default: out.bed

      -help:  shows this information\n\n";
    exit (1);
}

## define reading filehandle
my $fh;    ## filehandle
if ($file eq ""){
    $fh =*STDIN;
}
else{
    open ($fh, "<$file") or die ("cannot open $file for reading");
}

## define hash with IDs to be removed
my %remove;
unless ($remove  eq ""){
    %remove = %{General::read_table($remove,0)};
}


## save the sequence of the multifasta file in a hash
my %fasta;
my $id;

unless (defined($no_seq)){
    open (ARX, "<$fasta") || die ("cannot open $fasta for reading");
    while (<ARX>){
	chomp;
	if ($_=~m/>/){
	    $_=~s/>//;
	    $id = $_;
	    $fasta{$id}="";
	}
	else{
	    $fasta{$id}.=$_;
	}
    }
    close (ARX);
}


## open OUTPUT files
open (BED, ">$bedout") || die ("cannot open $bedout for writing");

### read and process the SAM file
while (<$fh>){
    chomp;
    my @line = split;
    my $mapping_quality= $line[4];
    my $pp;
    next if $_=~m/^@/;             ## discard header lines
    next if ($remove{$line[0]});   ## discard reads in the file of "reads to discard"

    next unless ($line[1] == 0 || $line[1]== 16);            ## we only keep reads mapped 

    if ($type eq "pssm"){
	my @v;
	if ($line[-1] =~m/PP/){
	    @v = split (/:/, $line[-1]);
	}
	elsif ($line[-2]=~m/PP/){
	    @v = split (/:/, $line[-2]);
	}
	else{
	    @v = split (/:/, $line[-3]);
	}
	$pp= $v[2];
	next if ($pp <= $cutoff || $pp > $maxcutoff);

    }
    else{
	next if $mapping_quality <= $cutoff;
    }


    my $readname= $line[0];
    my $chr = $line[2];
	
    ### define target coordinates in the genome (close notation) ###
    my $tstart =  $line[3];      ## target start	
    next if ($line[3] == 0);

    my $cigar_line = Alignment::parse_cigar ($line[5]);    ## cigar_line (full)
    my $target_len = scalar (my @match=$cigar_line=~m/(S|D|M)/g);	
    if ($cigar_line=~m/^(S+)/){  ## recalculate "real" target start after soft clipping
	my $start= length ($1);
	$tstart -=$start;
    }
    
    my $tend = $tstart+$target_len-1;  ## target end
    ## define the strand (of the mapping)
    my $strand = "+";
    if ($line[1] == 16){
	$strand = "-";
    } 
    ### genome coordinates
    my @genomic = split (":", $line[2]);
    my $gstrand= $genomic[1];

    ### final strand
    my $final_strand;
    if ($gstrand eq $strand){
	$final_strand = "+";
    }
    else{
	$final_strand = "-";
    }
    
    ## we get the genomic coordinates and the half read, qualities and cigar line
    my $gstart1;
    my $gend1;
    my $gstart2="";
    my $gend2="";
    
    ## if the read is not split mapped, we can discard it.
    ## how do we know is not split mapped?
    # 1.- the mapping start is after the length of the first exon
    # 2.- the mapping end is before the length of the first exon

    my $exon1len = $genomic[3]+1-$genomic[2];

    if ($tstart > $exon1len || $tend <= $exon1len){
	next;
    }
    ###  chr8:-:74742616:74742707:74722765:74722864
    ### chr16:-:70292883:70292982:70292021:70292120
    if ($gstrand eq "-"){
	$gstart1= $genomic[3]+1 - $tend;
	$gend1 =  $genomic[3]+1 - $tstart;
	my $offset= $genomic[2]-$gstart1;

	if ($offset > 0){
	    $gstart1 = $genomic[2];
	    $gend2= $genomic[5];
	    $gstart2 = $genomic[5]-$offset+1;
	}
    }
    else{
	$gstart1 = $genomic[2] -1 + $tstart;
	$gend1 =   $genomic[2] -1 + $tend;
	my $offset= $gend1-$genomic[3];

	if ($offset > 0){
	    $gend1   = $genomic[3];
	    $gstart2 = $genomic[4];
	    $gend2   = $genomic[4]+$offset -1;
	}
    }
    
    
    ### get the read sequence. The sequence of the read is always reported to the reference
    if (defined($no_seq)){
	print BED join ("\t", $genomic[0],$gstart1,$gend1,$readname."-1",$pp,$final_strand),"\n";
	print BED join ("\t", $genomic[0],$gstart2,$gend2,$readname."-2",$pp,$final_strand),"\n";
    }
    else{
	## get the read sequence
	my $readseq   = uc($line[9]);
	my $quality_scores= $line[10];
	
	### get the exon-junction sequence ###
	my $chr_file=$line[2];
	my $tseq = uc(substr ($fasta{$chr_file},$tstart-1, ($tend+1-$tstart)));
	
	if (length ($tseq) <  $tend+1 -$tstart){
	    print STDERR "SEQUENCE SHORTER THAN INPUT: $line[0]\t$chr_file\t$tstart\t$tend\t+\n";
	    print STDERR "TSEQ: $tseq\n";
	    my $diff = ($tend+1 -$tstart)-length ($tseq);
	    $tseq .= "N"x$diff;
	}
	
	## if the read is mapped antisense to the mRNA, we get the complement
	## otherwise the mutation patterns will be wrong
   	
	if ($strand eq "-"){              
	    $readseq =~tr/ATGC/TACG/;           ## complement read
	    $tseq =~tr/ATGC/TACG/;              ## complement target sequence
	}
	
        ### recover the alignment between the read and the target read 
	my @aln = Alignment::get_alignment_from_cigar($readseq,$tseq,$cigar_line);
    
	### make the extended cigar_line
	my $newcigar =Alignment::get_extended_cigar (@aln, $cigar_line);
    
	### my get quality scores
	my @stats = get_aligned_stats ($aln[0],$aln[1],$newcigar,$quality_scores,$final_strand);
	
	
	### divide the stats array
	my $len1 = $gend1+1-$gstart1;  ## genomic length of first chunk
	my $len2 = $gend2+1-$gstart2;  ## genomic length of first chunk
	
	my $len_count=$len1;
	
	if ($len1 +$len2 > length ($stats[1])){ ## we have indels in the alignment
	    my $len_count=0;
	    my @cigar_aln= split (//, $stats[1]);
	    
	    while ($len_count < $len1){
		if ($cigar_aln[$len_count] eq "I"){
		    $len1++;
		    $len_count++;
		}
		else{
		    $len_count++;
		}
	    }
	    $len_count++;
	}
	
	my @stats1;
	my @stats2;
	
	for (my $i = 0; $i < @stats -1 ; $i++){          ## we don't want to modify the insertion count
	    my $slice1= substr($stats[$i],0,$len_count);
	    my $slice2= substr($stats[$i],$len_count);
	    if ($strand eq "-"){
		$slice1 = reverse($slice1);
		$slice2 = reverse($slice2);
	    }
	    if ($final_strand eq "-"){
		$slice1 = reverse($slice1);
		$slice2 = reverse($slice2);
		if ($i == 0){
		    $slice1 =~tr/ATGC/TACG/;
		    $slice2 =~tr/ATGC/TACG/;
		}
	    }
	    push (@stats1, $slice1);
	    push (@stats2, $slice2);
	}

	## we need to divide the insertions in 2 groups and recalculate their position 
	## according to which half they are in

	my @ins1;
	my @ins2;
	my $ins1 = "-";
	my $ins2 = "-";
	
	my $l1= length ($stats1[0]);
	my $l2= length ($stats2[0]);
	
	if ($stats[-1] ne "-"){
	    my @ins = split (/;/, $stats[-1]);
	    for (my $i =0 ; $i < @ins; $i++){
		if ($ins[$i] <= $len_count){
		    push (@ins1, $ins[$i]);
		}
		else{
		    push (@ins2, $ins[$i]-$l1);
		}
	    }
	    
	    for (my $i = 0; $i < @ins1; $i++){
		if ($gstrand eq "-"){
		    $ins1[$i] = $l1 - $ins1[$i];
		}
	    }
	    
	    for (my $i = 0; $i < @ins2; $i++){
		if ($gstrand eq "-"){
		    $ins2[$i] = $l2 - $ins2[$i];
		}
	    }
	    
	    if ($gstrand eq "-"){
		@ins1 = reverse(@ins1);
		@ins2 = reverse(@ins1);
	    }
	    
	    if (@ins1 > 0){
		$ins1= join (";", @ins1);
	    }
	    if (@ins2 > 0){
		$ins2= join (";", @ins2);
	    }
	}
	### print BED FILE
	### format chr start end name score strand
	### score: number of tc mutations
	
	print BED join ("\t", $genomic[0],$gstart1,$gend1,$readname."-1",$pp,$final_strand,@stats1,$ins1),"\n";
	print BED join ("\t", $genomic[0],$gstart2,$gend2,$readname."-2",$pp,$final_strand,@stats2,$ins2),"\n";

    }
}

close $fh unless $file eq "";
close (BED);
    
    



##################################################
sub get_aligned_stats{
    my @query  = @{$_[0]};
    my @target = @{$_[1]};
    my @cig  = split (//, $_[2]);
    my @qual = split (//, $_[3]);
    my $strand = $_[4];
    my @alcig;
    my @alqual;
    my @fquery; ## final query removing "insertions"
     my @insertions;
    
    for (my $i = 0; $i < @query; $i++){
	
	if ($query[$i] eq "-"){     ## deletion in the read
	    push (@alcig, shift (@cig));
	    push (@alqual,"-");
	    push (@fquery, $query[$i]);
	}
	elsif($target[$i] eq "-"){  ## insertion in the read
	    shift (@cig);           # remove one position from both the read and the quality
	    shift (@qual);
	}	
	else{
	    push (@alcig, shift (@cig));
	    push (@alqual,shift (@qual));
	    push (@fquery, $query[$i]);
	}
    }
    my $ins="-";
    if (@insertions > 0){
	$ins = join (";", @insertions);
    }
  
    my $query = join ("", @fquery);
    my $cig = join ("", @alcig);
    my $qual = join ("", @alqual);

    if (length ($query) != length ($cig) || length ($query) !=  length ($qual)){
	print STDERR "wrong lengths\n";
	exit (1);
    }
    
    return ($query,$cig,$qual,$ins);
}

