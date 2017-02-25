CLAP relies on many dependencies, setting it up may require customization on your platform. Here The known errors will be discussed.

  scripts/make_exon_junction_library.pl

Error:

Can't exec "bedtools": No such file or directory at scripts/make_exon_junction_library.pl line 102.
cannot open out.7673.tab for reading at scripts/make_exon_junction_library.pl line 105.

Solution:
Either bedtoll is not installed properly or this line requires an explicit path to the binary of bedtool.
  
  my $command = "bedtools getfasta -fi $genome_file -bed $bedfile -tab -fo $seqfile -s";
  
It is located in 
bedtools2/bin
