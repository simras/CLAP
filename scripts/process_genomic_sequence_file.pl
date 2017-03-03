#!/usr/bin/perl -w
#
#   Simon H. Rasmussen
#   Bioinformatics Centre
#   Unversity of Copenhagen
############################
#   Developed for Ensembl ver. 87
#   Tested on species: hsa, mmu, dme and cel   

use strict;
use warnings;

if (@ARGV != 3){
    print "Removed non-wanted chromosomes. Introduce the fasta file with the genomic sequence for all chromosomes and chromosome names that match, which you wish to extract in a komma-separated list e.g. 1,2, ... ,X,Y,MT \n";
    print "\nUSAGE. ./process_genomic_sequence.pl chromosome_sequences.fa <chromosome-names> \n\n";
    exit (1);
}
my $mito = $ARGV[2];
my @m_list = split(",",$ARGV[1]);
my %genes;
my $p = 0;
my @l;

open (ARX, "<$ARGV[0]") || die ("cannot open $ARGV[0] for reading");
while (<ARX>){
    #if ($_ =~ m/^>[0-9MTXY]+ /){
    chomp;
    @l = split;
    if($_ !~ m/^>/){
	#print "XX1 $l[0]\n";
	if($p == 1){
	#print "XX2 $l[0]\n";
	    
	    # Sequence lines we wanna print
	    print "$l[0]\n"; 
	}#else{
	 #   print "XX3 $l[0]\n";
	#}	
    }elsif(match($l[0], @m_list)){
	# records  we wanna print
	$p = 1;
        # rename mitochondrial chromosomes for Human, Drosophila and C. Elegans
	$l[0] =~ s/$mito/M/g;
	#$l[0] =~ s/dmel_mitochondrion_genome/M/g;
	#$l[0] =~ s/MtDNA/M/g;

	# rename chromosome records to >chr1, chr2, ... , from >1, >2,...
	(my $chr = $l[0]) =~ s/>/>chr/g;
	print "$chr\n";
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
	    #print "$name == >$ll\n";
	    return 1;
        }
    }
    return 0
}
