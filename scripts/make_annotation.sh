#!/bin/bash
#
#  Simon H. Rasmussen
#  Bioinformatics Centre
#  University of Copenhagen
#
#
set -x
set -e

# Ensembl release version
ver=87

base_URL="ftp://ftp.ensembl.org/pub/release-"$ver

echo "unfortunately the names in the databases are not standardized in ver 87"
echo "So, if wget does not retrieve the file, find name of species in" $base_URL/fasta/$species/dna/
echo "and" $base_URL/gtf/species/

species="homo_sapiens"
name_gtf=$base_URL/gtf/$species/Homo_sapiens.GRCh38.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "MT"
rm resources/ensembl.$species.$ver.gtf
# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,Y,X,MT" "MT" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai
# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*

species="mus_musculus"
name_gtf=$base_URL/gtf/$species/Mus_musculus.GRCm38.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "MT"

# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,Y,X,MT" "MT" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai
# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*

species="caenorhabditis_elegans"
name_gtf=$base_URL/gtf/$species/Caenorhabditis_elegans.WBcel235.87.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "MtDNA"

# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "I,II,III,IV,V,X,MtDNA" "MtDNA" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai
# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*

# retrieve data
species="drosophila_melanogaster"
name_gtf=$base_URL/gtf/$species/Drosophila_melanogaster.BDGP6.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "dmel_mitochondrion_genome"

# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "2L,2R,3L,3R,4,X,Y,dmel_mitochondrion_genome" "dmel_mitochondrion_genome" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai
# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*

# Homo Sapiens
#1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,Y,X,MT
# Mus Musculus
#1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,Y,X,MT
# C. Elegans
# I,II,III,IV,V,X,MtDNA
# Drosophila Melanogastor
# 2L,2R,3L,3R,4,X,Y,dmel_mitochondrion_genome
# More species
# rattus_
# sacaramyces_cerevisie
# ftp://ftp.ensembl.org/pub/release-87/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.87.gtf.gz

# retrieve data
species="saccharomyces_cerevisiae"
name_gtf=$base_URL/gtf/$species/Saccharomyces_cerevisiae.R64-1-1.87.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "Mito"

# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito" "Mito" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai

# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*

species="rattus_norvegicus"
name_gtf=$base_URL/gtf/$species/Rattus_norvegicus.Rnor_6.0.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz

wget -O resources/ensembl.$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/ensembl.$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

# Create annotation files
scripts/create_mRNA_genome_annotation3.pl resources/ensembl.$species.$ver.gtf resources/ensembl.$species.$ver "MT"

# Process genomic sequence file
scripts/process_genomic_sequence_file.pl resources/$species.$ver.fa "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,X,Y,MT" "MT" > resources/ensembl.$species.$ver.fa
rm resources/$species.$ver.fa

# Make exon junction sequence library
scripts/make_exon_junction_library.pl resources/ensembl.$species.$ver.all.txt resources/ensembl.$species.$ver.fa > resources/ensembl.$species.$ver.ej.fa
rm resources/ensembl.$species.$ver.fa.fai
# make compressed archive
tar -czvf resources.$species.$ver.tar.gz resources/ensembl.$species.$ver.*
rm resources/ensembl.$species.$ver.*
