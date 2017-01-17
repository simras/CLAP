#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $help;
my $fastq;
my $sam;

&GetOptions('h'   => \$help,    ## help
	    'f=s' => \$fastq,   ## fastq file
	    's=s' => \$sam,     ## sam file
    );

if (defined($help) || !($fastq ||$sam)){
    print "USAGE: ./select_unmapped_reads.pl -s file.sam -f file.fastq\n";
    exit (1);
}

### sam file
my $sam_command;
if (substr ($sam, -3) eq "bz2"){
    $sam_command = "bzcat $sam";
}
else{
    $sam_command = "cat $sam";
}

my %reads;
open (OUT, "$sam_command|") || die ("cannot open $sam_command");
while (<OUT>){
    next if $_=~m/^@/;
    chomp;
    my @line = split;
    if ($line[1] == 4 || $line[1] == 20){
	$reads{$line[0]}= 1;
    }
}
close (OUT);

#print "UNMAPPED READS: ", scalar (keys %reads), "\n";
### fastq file
my $fastq_command;
if (substr ($fastq, -3) eq "bz2"){
    $fastq_command = "bzcat $fastq";
}
else{
    $fastq_command = "cat $fastq";
}

my $counter = 0;
my $id;
open (OUT, "$fastq_command|") || die ("cannot open $fastq_command");
while (<OUT>){
    chomp;
    if ($counter % 4 == 0){
	my @line = split;
	$id = $line[0];
	$id =~s/@//;
	$id =~s/\/1$//;
	$id =~s/\/2$//;
    }
    if ($reads{$id}){
	print "$_\n";
    }
    $counter++;
}
close (OUT);
	

