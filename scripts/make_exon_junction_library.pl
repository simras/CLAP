#!/usr/bin/perl -w
use strict;

if (@ARGV != 2){
    print "introduce the annotation file and the genome fasta file\n";
    print "The exon junction  library will be printed to standard out\n";

    print "USAGE: ./make_exon_junction_library.pl annotation.all.txt hg19.fa\n\n";
    exit (1);
}

## save exons
my %tr;

open (ARX, "<$ARGV[0]") || die ("cannot open $ARGV[0] for reading");
while (<ARX>){
    chomp;
    my @line = split;
    my @c = @line[1..2];
    ## {chr}{strand}{gene}{transcript} (exon1, exon2, ...)
    push (@{$tr{$line[0]}{$line[3]}{$line[6]}{$line[7]}}, \@c);
}
close (ARX);

## define pairs of exons

my %junctions;

foreach my $chr (keys %tr){
    foreach my $strand (keys %{$tr{$chr}}){
	foreach my $gene (keys %{$tr{$chr}{$strand}}){
	    foreach my $tr (keys %{$tr{$chr}{$strand}{$gene}}){
		my @exons = sort {$a->[0] <=> $b->[0]} @{$tr{$chr}{$strand}{$gene}{$tr}};
		next if (@exons < 2);
		for (my $i = 0; $i < @exons-1; $i++){
		    my @e1 = ($exons[$i][0],$exons[$i][1]);
		    my @e2 = ($exons[$i+1][0],$exons[$i+1][1]);
		    if ($exons[$i][1]-$exons[$i][0] >= 99){
			@e1 = ($exons[$i][1]-99,$exons[$i][1]);
		    }
		    if ($exons[$i+1][1]-$exons[$i+1][0] >= 99){
			@e2 = ($exons[$i+1][0],$exons[$i+1][0]+99);
		    }
		    my @pair = (@e1,@e2);
		    my $jun = "$chr:$strand:$e1[1]:$e2[0]";
		    push (@{$junctions{$jun}},\@pair);
		}
	    }
	}
    }
}

# get the sequence of the exon pairs and print it to a multifasta file

my $genome_file = $ARGV[1];
my %final_set;
my @final_ids;
my %segments;

my $bedfile="annotation".$$.".bed";
open (T, ">$bedfile") || die ("cannot open $bedfile for writing");

foreach my $id (keys %junctions){
    my ($chr,$strand,$s,$e) = split (":", $id);
    my @pairs = @{$junctions{$id}};
    my @segments=@{$pairs[0]};
    if (@pairs > 1){
	for (my $i = 1; $i < @pairs; $i++){
	    if ($pairs[$i][0] < $segments[0]){
		$segments[0] = $pairs[$i][0];
	    }
	    if ($pairs[$i][3] > $segments[3]){
		$segments[3] = $pairs[$i][3];
	    }
	}
    }
    # Convert to closed coordinates
    $segments[0] = $segments[0] - 1;
    $segments[2] = $segments[2] - 1;

    # Convert mitochodrial chromosome annotation
    $chr=~s/MT/M/;

    print T join ("\t", $chr,$segments[0],$segments[1],".",".",$strand), "\n";
    print T join ("\t", $chr,$segments[2],$segments[3],".",".",$strand), "\n";
    
    my $id1=join ("_", $chr,$segments[0],$segments[1],".",".",$strand);
    my $id2=join ("_", $chr,$segments[2],$segments[3],".",".",$strand);

    push (@final_ids, $id1,$id2);

    my $junction_id =join (":",$chr,$strand,@segments);
    if ($strand eq "-"){
	@{$final_set{$id}} = ($id2,$id1);
    }
    else{
	@{$final_set{$id}} = ($id1,$id2);
    }
}
close (T);


my $seqfile= "out.".$$.".tab";
my @seqs;
my $command = "../bedtools2/bin/bedtools getfasta -fi $genome_file -bed $bedfile -tab -fo $seqfile -s";

system ("$command");
unlink ($bedfile);

open (SEQ, "<$seqfile") || die ("cannot open $seqfile for reading");

while (<SEQ>){
    chomp;
    my @line = split;
    push (@seqs, uc($line[-1]));
}
close (SEQ);

### join name and sequence
for (my $i = 0; $i < @seqs; $i++){
    $segments{$final_ids[$i]} = $seqs[$i];
}

## print out fasta file 
foreach my $k (keys %final_set){
    my $seq1 = $segments{$final_set{$k}[0]};
    my $seq2 = $segments{$final_set{$k}[1]};

    print ">$k\n$seq1$seq2\n";
}

		     


  
