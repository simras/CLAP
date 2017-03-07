#!/usr/bin/perl -w

use strict;

if (@ARGV < 1){
    print "introduce the file with the annotation\n";
    exit (1);
}

my %trans;

open (ARX, "<$ARGV[0]") || die ("cannot open $ARGV[0] for reading");
while (<ARX>){
    chomp;
    my @line = split;
    
    ##chr strand gene_name
    @{$trans{$line[0]}{$line[3]}{$line[6]}{"attributes"}} =($line[7],$line[8],$line[9]); ## tr_id gene_name gene_type
    ## exons
    my @exon = ($line[1],$line[2],$line[4],$line[5]);
    push (@{$trans{$line[0]}{$line[3]}{$line[6]}{"exons"}}, \@exon);


    ## define start
    unless ($trans{$line[0]}{$line[3]}{$line[6]}{"start"}){
	$trans{$line[0]}{$line[3]}{$line[6]}{"start"} = $line[1];
    }
    elsif ($line[1] < $trans{$line[0]}{$line[3]}{$line[6]}{"start"}){	
	$trans{$line[0]}{$line[3]}{$line[6]}{"start"} = $line[1];
    }

    ## define end
    unless ($trans{$line[0]}{$line[3]}{$line[6]}{"end"}){
	$trans{$line[0]}{$line[3]}{$line[6]}{"end"} = $line[2];
    }
    elsif ($line[2] > $trans{$line[0]}{$line[3]}{$line[6]}{"end"}){	
	$trans{$line[0]}{$line[3]}{$line[6]}{"end"} = $line[2];
    }

    ## define length
    $trans{$line[0]}{$line[3]}{$line[6]}{"length"} = $line[5];
    
}
close (ARX);


foreach my $chr (keys %trans){
    foreach my $strand (keys %{$trans{$chr}}){
	my @genes;

	## save genes
	foreach my $id (keys %{$trans{$chr}{$strand}}){
	    my @g = ($id, $trans{$chr}{$strand}{$id}{"start"},$trans{$chr}{$strand}{$id}{"end"});
	    push (@genes, \@g);
	}

	## sort genes
	@genes = sort{$a->[1]<=>$b->[1] || $a->[2] <=>$b->[2]}@genes;

	## we make a 
	my @clusters;
	my @c = @{$genes[0]};
	push (@c, 1);         ## count of number of genes
	push (@clusters,\@c);

	for (my $i = 1; $i < @genes; $i++){
	    if ($genes[$i][1] <= $clusters[-1][2]){
		$clusters[-1][3]++;   ## we have a gene more
		
		if ($genes[$i][2] > $clusters[-1][2]){   ## we make the cluster longer
		    $clusters[-1][2] = $genes[$i][2];
		}
		
		$clusters[-1][0]= join (":", $clusters[-1][0],$genes[$i][0]); ## add the Ids
	    }
	    else{
		my @c = @{$genes[$i]};
		push (@c, 1);
		push (@clusters,\@c);
	    }
	}

	for (my $i = 0; $i < @clusters; $i++){
	    my $id;
	    if ($clusters[$i][-1] ==1){
		$id = $clusters[$i][0];
	    }
	    else{
		my @ids = split (/:/, $clusters[$i][0]);
#		print join ("\t", @ids), "\n";
		my $maxlen=0;
		
		foreach my $p (@ids){
		    if ($trans{$chr}{$strand}{$p}{"length"} > $maxlen){
			$maxlen = $trans{$chr}{$strand}{$p}{"length"};
			$id = $p;
		    }
		}
	    }
#	    print "ID: $id\n";
	    my @exons = @{$trans{$chr}{$strand}{$id}{"exons"}};
	    for (my $z = 0; $z < @exons; $z++){
		my @toprint = ($chr,$exons[$z][0],$exons[$z][1],$strand,$exons[$z][2],$exons[$z][3],$id, @{$trans{$chr}{$strand}{$id}{"attributes"}});
		print join ("\t", @toprint), "\n";
	    }
	}
    }   
}
	
