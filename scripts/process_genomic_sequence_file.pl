#!/usr/bin/perl -w
#
#   Simon H. Rasmussen
#   Bioinformatics Centre
#   Unversity of Copenhagen
############################   

use strict;
use warnings;

if (@ARGV != 2){
    print "Removed non-wanted chromosomes. Introduce the fasta file with the genomic sequence for all chromosomes and chromosome names that match, which you wish to extract in a komma-separated list e.g. 1,2, ... ,X,Y,MT \n";
    print "\nUSAGE. ./process_genomic_sequence.pl chromosome_sequences.fa <chromosome-names> \n\n";
    exit (1);
}

my @m_list = split(",",$ARGV[1]);
my %genes;
my $p = 0;
my @l;

open (ARX, "<$ARGV[0]") || die ("cannot open $ARGV[0] for reading");
while (<ARX>){
    #if ($_ =~ m/^>[0-9MTXY]+ /){
    chomp;
    @l = split;
    #print "XX $l[0], \n";
    if($_ =! m/^>/){
	if($p){
	    # Sequence lines we wanna print
	    print "$l[0]\n"; 
	}
    }elsif(match($l[0], @m_list)){
	# records  we wanna print
	$p = 1;
	print "$l[0]\n";
    }elsif ($_ =~ m/^>/){
	# records we don't wanna print
	$p=0;
    }   
}
close (ARX);


sub match {
    # Will match the list with chromosome names 
    # and the record names in the fasta file
    
    my $name = $_[0];    
    my @l = splice(@_,0,$#_ + 1);
    foreach my $ll (@l){
	if("$name" eq ">$ll"){
	    return 1;
        }
    }
    return 0
}
