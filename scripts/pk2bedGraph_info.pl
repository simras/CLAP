#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $cols;
my $strand;
my $track_name;
my $track_description;
my $help;

&GetOptions('c=s'   => \$cols,      
	    's=s'   => \$strand,
	    'n=s'   => \$track_name,
	    'd=s'   => \$track_description,
	    'help'  => \$help,      ## help
	   );

if (defined($help) || !($cols||$strand||$track_name||$track_description)){
    print STDERR "
    USAGE cat file.pk | ./pk2bedGraph_info.pl  -c 1,2,3,4 -n NAME -s + -d clusters of protein X

      -c:     columns where to find chr, start, pk_line, strand

      -s:     strand: + or -
      
      -n:     track_name (no spaces)

      -d:     track description

      -help\n\n";
    exit(1);
}


my $start;
my $end;
my ($chrp,$startp,$pkp,$strandp) = split (/,/, $cols);
$chrp--;
$startp--;
$pkp--;
$strandp--;

my $color;
if ($strand eq "+"){
    $color="38,139,3";
}
else{
    $color="238,0,0";
}


print "track type=bedGraph name=$track_name description=\"$track_description\" color=$color\n";

while (<STDIN>){
    chomp;
    my @line = split;
    next unless $line[$strandp] eq $strand;
    my @positions = split/\|/, $line[$pkp];
    $start = $line[$startp]-1;
    $end= $start;
    
    foreach my $p (@positions){
	my @p = split/:/, $p;
	$start = $end;
	$end= $start+$p[0];
	print "$line[$chrp]\t$start\t$end\t$p[1]\n"; ## chr,start,end,number of reads
    }
}


