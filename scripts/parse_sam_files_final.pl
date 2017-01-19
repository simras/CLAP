#!/usr/bin/perl -w

use strict;
use lib './resources/';
use Pipeline_functions;
use Getopt::Long;

my $file   ="";
my $remove ="";
my $type;
my $path;
my $cutoff = 0.5;
my $maxcutoff = 1;
my $bedout = $$."out.bed";
my $batchsize=10000;
my $help;
my $no_seq;

&GetOptions('f:s'   => \$file,      ## query file with the sam file
	    't=s'   => \$type,
	    'p=s'   => \$path,      ## path where the chr files are
	    'c:f'   => \$cutoff,    ## cut off to define the uniquely mapped reads
	    'm:f'   => \$maxcutoff, ## max cut off to select reads
	    'o:s'   => \$bedout,    ## bed output file
	    'r:s'   => \$remove,    ## list of IDs to be disregarded from the analysis
	    'b:i'   => \$batchsize, ## number of lines that will be processed at the same time
	    'a'     => \$no_seq,    ## flag saying if the sequence, the qualities and the CIGAR alignment should be printed
	    'help'  => \$help,      ## help
	   );

## define options
if (defined($help) || !($path || $type)){
    print STDERR "
    USAGE ./parse_sam_files.pl -p path -t pssm

      -p:     index file

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
      
      -a:     only produce the bed file and do not recover the sequence
              Default: recover sequence
      
      -b:     number of records to process at a time
              Default: 10000

      -o:     name of the output bed file
              Default: out.bed

      -help:  shows this information\n\n";
    exit (1);
}

#print STDERR $batchsize, "\n";
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
    %remove = %{Pipeline_functions::read_table($remove,0)};
}

## open OUTPUT file
open (BED, ">$bedout") || die ("cannot open $bedout for writing");


my @batcharray;
my $batchcounter=0;
my %chr_sizes;

### read and process the SAM file
while (<$fh>){
    chomp;
    my @line = split;
    my $mapping_quality= $line[4];
    my $pp;
    if ($_=~m/^\@SQ/){
	my $chr_name= $line[1];
	$chr_name=~s/SN://;
	my $chr_length= $line[2];
	$chr_length=~s/LN://;
	$chr_sizes{$chr_name}=$chr_length;
    }

    next if $_=~m/^@/;             ## discard header lines
    next if ($remove{$line[0]});   ## discard reads in the file of "reads to discard"

    next unless ($line[1] == 0 || $line[1]== 16);  ## discard reads that are not mapped

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
    my $cigar_line = Pipeline_functions::parse_cigar ($line[5]);    ## cigar_line (full)
    my $target_len = scalar (my @match=$cigar_line=~m/(S|D|M)/g);	

    if ($cigar_line=~m/^(S+)/){  ## recalculate "real" target start after soft clipping
	my $start= length ($1);
	$tstart -=$start;
    }
    my $tend = $tstart+$target_len-1;  ## target end
    
    ## discard reads that are outside the chr
    next if ($tstart < 1 || $tend >$chr_sizes{$chr});
    
    ## define the strand
    my $strand = "+";
    if ($line[1]==16){
	$strand = "-";
    }

    ## get the read sequence
    my $readseq   = $line[9];
    my $quality_scores= $line[10];
    if ($strand eq "-"){
	$readseq =~tr/ATGC/TACG/;
	$readseq =reverse ($readseq);
	$cigar_line = reverse ($cigar_line);
	$quality_scores= reverse ($quality_scores);
    } 
    
    ### we push the info we need in @tosave
    my @tosave = ($readname,$chr,$tstart,$tend,$strand,$target_len,$cigar_line,$mapping_quality,$readseq,$pp,$quality_scores);
    push (@batcharray, \@tosave);
    $batchcounter++;
    
    if ($batchcounter== $batchsize){
	&process_batch (\@batcharray,\*BED, $path,$no_seq);
	@batcharray=();
	$batchcounter=0;
    }
}

close $fh unless $file eq "";

if (scalar (@batcharray) > 0){
    &process_batch (\@batcharray,*BED,$path,$no_seq);
}

close (BED);
    

### FUNCTIONS ################################################################

## pas as parameters the pointers of BED, TC and ST file pointers
sub process_batch{
    my ($batch, $bed_fh,  $index,$noseq) = @_;
    my @batch = @{$batch};
    my @seqs;
    my $seqfile;

    unless  (defined($noseq)){
	$seqfile = recover_seqs($batch,$index);
    
	# ($readname,$chr,$tstart,$tend,$strand,$target_len,$cigar_line,$mapping_quality,$readseq,$pp);
	open (SEQ, "<$seqfile") || die ("cannot open $seqfile for reading");
	while (<SEQ>){
	    chomp;
	    my @line = split;
	    push (@seqs, uc($line[-1]));
	}
	close (SEQ);
    }
    
    for (my $i = 0; $i < @batch; $i++){
	my ($readname,$chr,$tstart,$tend,$strand,$target_len,$cigar_line,$mapping_quality,$readseq,$pp,$quality_scores) = @{$batch[$i]};
	
	unless (defined($noseq)){
	    my $tseq = $seqs[$i];
	
	    ### recover the alignment between the read and the target read 
	    my @aln = Pipeline_functions::get_alignment_from_cigar($readseq,$tseq,$cigar_line);
	
	    ### make the extended cigar_line
	    my $newcigar =Pipeline_functions::get_extended_cigar (@aln, $cigar_line);
	
	    ### my get quality scores
	    my @stats = get_aligned_stats ($aln[0],$aln[1],$newcigar,$quality_scores,$strand);
	    
	    ### print BED FILE
	    ### format chr start end name score strand sequence new_cigar qualitie
	    print $bed_fh join ("\t", $chr,$tstart,$tend,$readname,$pp,$strand,@stats),"\n";
	}
	else{
	    print $bed_fh join ("\t", $chr,$tstart,$tend,$readname,$pp,$strand),"\n";
	}
    }
    unless (defined($noseq)){
	unlink ($seqfile);
    }
}

##################################################
sub recover_seqs{
    my @batch = @{$_[0]};;
    my $index = $_[1];
    my $bedfile = "file.".$$.".bed";
    
    ## print bed file
    open (ARX, ">$bedfile") || die ("cannot open $bedfile for writing");
    for (my $i = 0; $i < @batch; $i++){
	print ARX join ("\t", $batch[$i][1], $batch[$i][2]-1, $batch[$i][3],".",".",$batch[$i][4]), "\n";
    }
    close (ARX);

    my $outfile = Pipeline_functions::getseq_bedtools($bedfile,$index);
    unlink ($bedfile);
    return $outfile;
}

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
	    shift (@cig);           ## remove one position from both the read and the quality
	    shift (@qual);
	    push (@insertions,scalar (@alcig));
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
    
    if ($strand eq "-"){
	$qual = reverse ($qual);
	$cig =reverse($cig);
	$query=~tr/ATGC/TACG/;
	$query =reverse($query);
	if (@insertions > 0){
	    for (my $i = 0; $i < @insertions; $i++){
		$insertions[$i]= length ($qual)-$insertions[$i];
	    }
	    @insertions= reverse (@insertions);
	    $ins = join (";", @insertions);	
	}
    }
    return ($query,$cig,$qual,$ins);
}
