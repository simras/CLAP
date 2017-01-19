package Pipeline_functions;

use strict;

sub read_table{
    my $file = shift;
    my $pos = shift;
    my %hash;

    open (ARX, "<$file") || die ("cannot open $file for reading"); 
    while (<ARX>){
	next if ($_=~m/^#/);
	chomp;
	my @line = split;
	my $id = splice (@line,$pos,1);
	if (@line == 0){
	    $hash{$id}= 1;
	}
	elsif (@line== 1){
	    $hash{$id}= $line[0];
	}
	else{
	    $hash{$id}= \@line;
	}
    }
    close (ARX);
    return \%hash;
}


sub getseq_bedtools{
    my ($bedfile,$genomefile)= @_;
    my $outfile= "out.".$$.".tab";
    my $command = "bedtools getfasta -fi $genomefile -bed $bedfile -tab -fo $outfile -s";
    system ("$command");
    return $outfile;
}


sub parse_cigar{
    my $line = $_[0];
    my @matches =$line=~m/(\d+[SMID])/g;
    my $alignment="";
    foreach my $m (@matches){
	if ($m=~m/(\d+)([SMID])/){
	    my $num = $1;
	    my $let = $2;
	    $alignment .= $let x $num;
	}
    }
    return $alignment;
}


sub get_alignment_from_cigar{
    my ($qseq,$tseq,$cigar) = @_;
    my @qseq = split(//, $qseq);
    my @tseq = split(//, $tseq);
    my @cigar = split(//, $cigar);
    
    my @alnqseq;
    my @alntseq;
    
    foreach my $pos (@cigar){
	if ($pos eq "M" || $pos eq "S"){
	    push (@alnqseq, shift (@qseq));
	    push (@alntseq, shift (@tseq));
	}
	elsif ($pos eq "I"){
	    push (@alnqseq, shift (@qseq));
	    push (@alntseq, "-");
	}
	elsif ($pos eq "D"){
	    push (@alnqseq, "-");
	    push (@alntseq, shift (@tseq));
	}
    }
    return (\@alnqseq,\@alntseq);
}


sub get_extended_cigar {
    my @qseq = @{$_[0]};
    my @tseq = @{$_[1]};
    my @cigar = split(//, $_[2]);
    
    my @newcigar;
    my %conversion = ("AC" => "a",
		      "AG" => "b",
		      "AT" => "c",
		      "CA" => "d",
		      "CG" => "e",
		      "CT" => "f",
		      "GA" => "g",
		      "GC" => "h",
		      "GT" => "i",
  		      "TA" => "j",
		      "TC" => "k",
 		      "TG" => "l",
		      "AN" => "x",
    		      "CN" => "x",
		      "GN" => "x",
  		      "TN" => "x",
  		      "NN" => "x",
		      "NA" => "x",
    		      "NC" => "x",
		      "NG" => "x",
  		      "NT" => "x");
    
    for (my $i = 0; $i < @cigar; $i++){
	if ($cigar[$i] eq "M"){
	    if ($qseq[$i] eq $tseq[$i]){
		push (@newcigar,$cigar[$i]);
	    }
	    else{
		my $key = $qseq[$i].$tseq[$i];
		push (@newcigar, $conversion{$key});
	    }
	}
	else{
	    push (@newcigar, $cigar[$i]);
	}
    }

    my $n = join ("", @newcigar);

    return $n;
}

1;
