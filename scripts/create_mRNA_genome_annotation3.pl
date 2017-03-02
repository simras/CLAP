#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 2){
    print "introduce the gtf file with the ensembl annotation and the output prefix\n";
    print "\nUSAGE. ./create_mRNA_genome_annotation.pl ensembl.gff output\n\n";
    exit (1);
}

my $outname= $ARGV[1];
my %genes;
open (ARX, "<$ARGV[0]") || die ("cannot open $ARGV[0] for reading");
while (<ARX>){
    next if ($_=~m/^#/);
    chomp;
    my @l = split;
    next unless ($l[2] eq "exon" || $l[2] eq "CDS");
   
    my $exon = read_ensembl87_gtf($_);

    # test if coding or non-coding exon
    if ($exon -> {"type"} eq "exon"){
	push (@{$genes{$exon->{"chr"}}{$exon->{"transcript_id"}}{"exons"}}, $exon);
    }
    elsif ($exon -> {"type"} eq "CDS"){
	push (@{$genes{$exon->{"chr"}}{$exon->{"transcript_id"}}{"CDS"}}, $exon);
    }
}
close (ARX);

my $all  = $outname.".all.txt";
my $all_introns=  $outname.".introns.all.txt";
my $cds  = $outname.".cds.txt";
my $cds_introns=  $outname.".introns.cds.txt";
my $utr5 = $outname.".5utr.txt";
my $utr3 = $outname.".3utr.txt";

open (ALL, ">$all") || die ("cannot $all for writing");
open (CDS, ">$cds") || die ("cannot $cds for writing");
open (UTR5, ">$utr5") || die ("cannot $utr5 for writing");
open (UTR3, ">$utr3") || die ("cannot $utr3 for writing");
open (ALLINT, ">$all_introns") || die ("cannot $all_introns for writing");
open (CDSINT, ">$cds_introns") || die ("cannot $cds_introns for writing");


## we reconstruct the transcripts.
foreach my $chr (keys %genes){
    foreach my $transcript (keys %{$genes{$chr}}){
	## sort the exons
	my @exons = sort{$a ->{"start"} <=> $b ->{"start"}}@{$genes{$chr}{$transcript}{"exons"}};
	
	my $gene_start = $exons[0]  -> {"start"};
	my $gene_end   = $exons[-1] -> {"end"};
	my $strand = $exons[0] ->{"strand"};
	if ($strand eq "-"){
	    @exons = sort{$b ->{"end"} <=> $a ->{"end"}}@exons;
	    $gene_start = $exons[0]->{"end"};
	    $gene_end   = $exons[-1]->{"start"};
	}
	### PRINT ALL ######################################################################################
	&print_all (\@exons, *ALL);
	&print_introns (\@exons, *ALLINT);
	
	####################################################################################################
	
	next unless $genes{$chr}{$transcript}{"CDS"};
	
	my $utr5_start = "-";
	my $utr5_end   = "-";
	my $utr3_start = "-";
	my $utr3_end   = "-";
	
	## we get the coding region of a gene
	my @cds =  sort{$a ->{"start"} <=> $b ->{"start"}}@{$genes{$chr}{$transcript}{"CDS"}};

	## we define the coding region start and the coding region end
	my $cds_start  = $cds[0]    -> {"start"};
	my $cds_end    = $cds[-1]   -> {"end"};
	
	if ($strand eq "-"){
	    @cds =  sort{$b ->{"end"} <=> $a ->{"end"}}@cds;
	    $cds_start = $cds[0]->{"end"}; 
	    $cds_end   = $cds[-1]->{"start"};
	}
	
	
	my $i = 0;
	my $done =0;

#	print $exons[$i] -> {"start"}, "\t", $exons[$i] -> {"end"}, "\t", $exons[$i]->{"strand"}, "\t",$start, "\t", $end, "\t", $exons[$i] -> {"gene_id"}, "\t",$exons[$i] -> {"transcript_id"},"\t", $exons[$i]->{"gene_name"}, "\t", $exons[$i]->{"biotype"}, "\n";

	if ($cds_start != $gene_start){   ## the transcript has a 5'UTR annotated 
	    if ($strand eq "+"){
		$utr5_start = $gene_start;
		$utr5_end   = $cds_start-1;
	    }
	    else{
		$utr5_start = $gene_start;
		$utr5_end   = $cds_start+1;
	    }
	    
#	    ## check if 5'utr end is well defined
	    while ($i < @exons && $done == 0){
		if ($cds[0]->{"start"} >=$exons[$i]->{"start"}  && $cds[0]->{"start"} <= $exons[$i]->{"end"}){
		    if ($cds[0]->{"start"} == $exons[$i]->{"start"} && $strand eq "+"){
			$utr5_end = $exons[$i-1] ->{"end"};
		    }
		    elsif($cds[0]->{"end"} == $exons[$i]->{"end"} && $strand eq "-"){
			$utr5_end = $exons[$i-1] ->{"start"};
		    }
		    $done = 1;
		}
		$i++;
	    }
	}
	if ($cds_end != $gene_end){     ## the transcript has a 3'UTR annotated 
	    if ($strand eq "+"){
		$utr3_start = $cds_end+1;
		$utr3_end   = $gene_end;
	    }
	    else{
		$utr3_start = $cds_end-1;
		$utr3_end   = $gene_end;
	    }
	    
#	    ## check if 3'utr end is well defined
	    $done = 0;
	    $i = 0;
	    while ($i < @exons && $done == 0){
		if ($cds[-1]->{"end"} >=$exons[$i]->{"start"}  && $cds[-1]->{"end"} <= $exons[$i]->{"end"}){
		    if ($cds[-1]->{"end"} == $exons[$i]->{"end"} && $strand eq "+"){
			$utr3_start = $exons[$i+1] ->{"start"};
		    }
		    elsif($cds[-1]->{"start"} == $exons[$i]->{"start"} && $strand eq "-"){
			$utr3_start = $exons[$i+1] ->{"end"}; ## this exon does not exit
		    }
		    $done = 1;
		}
		$i++;
	    }
	    
	}
	
        ### PRINT STUFF ######################################################################################
	if ($utr5_end ne "-"){
	    &print_func (\@exons, $utr5_start,$utr5_end, *UTR5);
	}
	if ($cds_start ne "-"){
	    &print_all (\@cds,*CDS);
	}
	if ($cds_start ne "-"){
	    &print_introns (\@cds,*CDSINT);
	}
	if ($utr3_end ne "-"){
	    &print_func (\@exons, $utr3_start, $utr3_end,*UTR3);
	}
    }
}
	
close (ALL);
close (CDS);
close (UTR5);
close (UTR3);
close (ALLINT);
close (CDSINT);


sub print_all{
    my @exons = @{$_[0]};
    my $fh = $_[1];
    my ($start,$end);
    for (my $i = 0;$i < @exons; $i++){
	if ($i == 0){
	    $start = 1;
	}
	else{
	    $start= $end+1;
	}
	$end = $start + abs($exons[$i] -> {"end"}-$exons[$i] -> {"start"});

	my $printchr;
	if ($exons[0]->{"chr"} eq "MT"){
	    $printchr="chrM";
	}
	elsif ($exons[0] ->{"chr"} =~/chr/){
	    $printchr=$exons[0]->{"chr"};
	}
	else{
	    $printchr="chr".$exons[0]->{"chr"};
	}
	
	print $fh $printchr, "\t", $exons[$i] -> {"start"}, "\t", $exons[$i] -> {"end"}, "\t", $exons[$i]->{"strand"}, "\t",$start, "\t", $end, "\t", $exons[$i] -> {"gene_id"}, "\t",$exons[$i] -> {"transcript_id"},"\t", $exons[$i]->{"gene_name"}, "\t", $exons[$i]->{"biotype"}, "\n";

    }
}

sub print_introns{
    my @exons = @{$_[0]};
    my $fh = $_[1];
    my ($start,$end);

    @exons = sort{$a ->{"start"} <=> $b ->{"start"}}@exons;
    
    my $printchr;
    if ($exons[0]->{"chr"} eq "MT"){
	$printchr="chrM";
    }
    elsif ($exons[0] ->{"chr"} =~/chr/){
	    $printchr=$exons[0]->{"chr"};
	}
    else{
	$printchr="chr".$exons[0]->{"chr"};
    }
   
    ## we define the introns;
    for (my $i = 0; $i < @exons -1; $i++){
	my ($start,$end);
	$start = 1;
	$end = ($exons[$i+1] -> {"start"}-1) +1 - ($exons[$i] -> {"end"}+1);
   
	
	print $fh $printchr, "\t", $exons[$i] -> {"end"}+1, "\t", $exons[$i+1] -> {"start"}-1, "\t", $exons[$i]->{"strand"}, "\t",$start, "\t", $end, "\t",$exons[$i] -> {"gene_id"}, "\t",$exons[$i] -> {"transcript_id"},"\t", $exons[$i]->{"gene_name"}, "\t", $exons[$i]->{"biotype"}, "\n";
    }
}
	
sub print_func{
    my ($e, $rstart,$rend,$fh) = @_;
    my @exons = @{$e};
    my $start="-";
    my $end ="-";
    my $started= 0;
    my $ended=0;
    my $done = 0;

    my $printchr;
    if ($exons[0]->{"chr"} eq "MT"){
	$printchr="chrM";
    }
    elsif ($exons[0] ->{"chr"} =~/chr/){
	    $printchr=$exons[0]->{"chr"};
	}
    else{
	$printchr="chr".$exons[0]->{"chr"};
    }
    
    for (my $i = 0;$i < @exons && ($started== 0 || $ended ==0); $i++){
	my $gstart = $exons[$i]->{"start"};
	my $gend =   $exons[$i] ->{"end"};
	if ($gend eq "-" || $gstart eq "-"){
	    print $exons[$i]->{"transcript_id"}, "\n";
	    exit (1);
	}
	if ($rstart >= $exons[$i] -> {"start"} && $rstart <= $exons[$i] ->{"end"}){ ## the exon contains the start
	    $start = 1;
	    $started =1;
	    if  ($exons[$i]->{"strand"} eq "+"){
		$gstart = $rstart;
	    }
	    else{
		$gend = $rstart;
	    }
	}
	elsif ($started == 1){
	    $start = $end+1;
	}
	if ($rend >= $exons[$i] -> {"start"} && $rend <= $exons[$i] ->{"end"}){ ## the exon contains the end
	    $ended = 1;
	    if  ($exons[$i]->{"strand"} eq "+"){
		$gend = $rend;
		if ($start eq "-"){
		    print "the exon contains the end: +strand \n", join(" ",keys $exons[$i],"\n",values $exons[$i]), "\n";
		    exit (1);
		}
		$end = $start + abs ($gend-$gstart);

	    }
	    else{
		$gstart = $rend;
		if ($start eq "-"){
		    print "the exon contains the end -strand: \n", join(" ",keys $exons[$i],"\n",values $exons[$i]), "\n";
		    exit (1);
		}
		$end = $start + abs ($gend-$gstart);
	    }
	}
	elsif($started == 1){  ## the exon does not contains the end
	    $end = $start + abs($gend-$gstart);
	}

	## we only need to print results when $start and $end are defined and 
	if ($started == 1 && $done== 0){
	    print $fh $printchr, "\t", $gstart, "\t", $gend, "\t", $exons[$i] -> {"strand"}, "\t",$start, "\t", $end, "\t", $exons[$i] -> {"gene_id"}, "\t",$exons[$i] -> {"transcript_id"},"\t", $exons[$i]->{"gene_name"}, "\t", $exons[$i]->{"biotype"}, "\n";
	    if ($ended == 1){
		$done= 1;
	    }
	}
    }
}

sub read_ensembl87_gtf{
    
    my $string = shift;
    my %hash;
    my @line = split("\t",$string);
    
    $hash{"chr"} = $line[0];
    $hash{"type"} = $line[2];
    $hash{"start"} = $line[3];
    $hash{"end"} = $line[4];
    $hash{"score"} = $line[5];
    $hash{"strand"} = $line[6];
    $hash{"frame"} = $line[7];
 #   print $string, "\n";
    my $parent_line = $line[8];
    $hash{"gene_id"} = $string; 
    $hash{"gene_id"} =~ /(gene_id.")([A-Z0-9]+)/;
    $hash{"gene_id"} = $2;
    $hash{"gene_name"} = $string; 
    $hash{"gene_name"} =~ /(gene_name.")([_A-Za-z0-9\.\/-]+)/;
    $hash{"gene_name"} = $2;
# maybe it should be transcript_biotype instead
    $hash{"biotype"} = $string;
    $hash{"biotype"} =~ /(gene_biotype.")([A-Za-z0-9_]+)/;
    $hash{"biotype"} = $2;
    $hash{"transcript_id"} = $string;
    $hash{"transcript_id"} =~ /(transcript_id.")([A-Z0-9]+)/;
    $hash{"transcript_id"} = $2;
#    print Dumper(\%hash);
    return \%hash;
}
sub read_ensembl86_gtf{
    #use Data::Dumper;
   #my $string = shift;
    # my @line = split;
    my @line = split;
    my %hash;
 #   my @line = split($string);
    $hash{"chr"} = $line[0];
    $hash{"type"} = $line[2];
    $hash{"start"} = $line[3];
    $hash{"end"} = $line[4];
    $hash{"score"} = $line[5];
    $hash{"strand"} = $line[6];
    $hash{"frame"} = $line[7];
   
    for (my $i = 8; $i < @line-1; $i++){
    if ($line[$i] eq "gene_id"){
        $hash{"gene_id"} = $line[$i+1];
        $hash{"gene_id"} =~s/\"|;//g;
        $i++;
    }
    elsif ($line[$i] eq "transcript_id"){
        $hash{"transcript_id"} = $line[$i+1];
        $hash{"transcript_id"} =~s/\"|;//g;
        $i++;
    }
    elsif ($line[$i] eq "gene_name"){
        $hash{"gene_name"} = $line[19];
        $hash{"gene_name"} =~s/\"|;//g;
        $i++;
    }
    elsif ($line[$i] eq "gene_biotype"){
        $hash{"biotype"} = $line[$i+1];
        $hash{"biotype"} =~s/\"|;//g;
        $i++;
    }
    }
#    print Dumper(\%hash);
    return \%hash;
}
