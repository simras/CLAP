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
name_dna=$base_URL/fasta/$species/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget -O resources/$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

species="mus_musculus"
name_gtf=$base_URL/gtf/$species/Mus_musculus.GRCh38.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Mus_musculus.GRCh38.dna.primary_assembly.fa.gz

wget -O resources/$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

species="caenorhabditis_elegans"
name_gtf=$base_URL/gtf/$species/Caenorhabditis_elegans.WBcel235.87.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

wget -O resources/$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz

species="drosophila_melanogaster"
name_gtf=$base_URL/gtf/$species/Drosophila_melanogaster.BDGP6.87.chr.gtf.gz
name_dna=$base_URL/fasta/$species/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz

wget -O resources/$species.$ver.gtf.gz $name_gtf
wget -O resources/$species.$ver.fa.gz $name_dna
gunzip resources/$species.$ver.gtf.gz
gunzip resources/$species.$ver.fa.gz



#species="Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
#species="Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz"

