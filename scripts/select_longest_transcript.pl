#!/usr/bin/perl -w

use strict;

if (@ARGV < 2){
    print "introduce the all file and all cds, 3utr and 5utr files\n";
    print "USAGE: ./select_longest_transcript.pl ensembl70.all.txt ensembl70.cds.txt ensembl70.3utr.txt ensembl70.5utr.txt\n";
    exit (1);
}

## gets the longest protein_coding

my $longest= get_longest ($ARGV[0]);


foreach my $file (@ARGV){
    my $out = $file;
    $out=~s/.txt/.long.txt/;
    &select_longest_tr ($longest,$file,$out);
}

### FUNCTIONS ######################################################################################
sub select_longest_tr{
    my %long = %{$_[0]};
    my $file = $_[1];
    my $outfile= $_[2];
    
    open (ARX, "<$file") || die ("cannot open $file for reading");
    open (OUT, ">$outfile") || die ("cannot open $file for writing");
     while (<ARX>){
	chomp;
	my @line = split;
	if ($long{$line[7]}){
	    print OUT "$_\n";
	}
     }
    close (ARX);
    close (OUT);
}

####################################################################################################
sub get_longest{
    my $file = shift;
    my %hash;
    my %longest;
    my %cds;
    
    open (ARX, "<$file") || die ("cannot open $file for reading");
    while (<ARX>){
	chomp;
	my @line = split;
	next unless ($line[9] eq "protein_coding");
	unless($hash{$line[6]}{$line[7]}){
	    $hash{$line[6]}{$line[7]} = $line[5];
	}
	if ($line[5] > $hash{$line[6]}{$line[7]}){
	    $hash{$line[6]}{$line[7]}  = $line[5];
	}
    }
    close (ARX);
    
    foreach my $gene (keys %hash){
	my $id;
	my $len = 0;
	
	foreach my $tr (keys %{$hash{$gene}}){
	    if ($hash{$gene}{$tr} > $len){
		$len = $hash{$gene}{$tr};
		$id= $tr;
	    }
	}
	$longest{$id}= $len;
    }
    return \%longest;
}
	
	    
