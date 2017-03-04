CLIP Analysis Pipeline - CLAP, is used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data. The family of protocols where cross-linking immunoprecipitation coupled with high-throughput sequencing are referred to as CLIP. They are used identify binding sites of RNA-binding proteins and exist in the different flavours mentioned above. The data produced by these protocols have different characteristics which are account for in this pipeline such that it can be used to process the above three versions and combination of them like iCLIP using 4SU nucleosides (PAR-iCLIP).

## 1. SYSTEM REQUIREMENTS
A relatively powerful PC, minimun 4 GB memory and with Linux installed. The pipeline is tested on Ubuntu 16.04, in principle it should run on a standard Linux setup with BASH Shell and some version of awk installed. We will not guarantee that the
pipeline works on Apple computers, it could as they are based on FreeBSD, but we provide no support to the end of making it run on a Mac.

## 2. INSTALLATION AND CONFIGURATION
1. Clone repository, find a suitable directory and use git to dowload the repository

        git clone https://github.com/simras/CLAP/CLAP.git

1. Install Python (link https://www.python.org/)

2. Install bwa-pssm (https://github.com/pkerpedjiev/bwa-pssm)

3. Install bedTools (http://bedtools.readthedocs.io/en/latest/)

4. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo)

5. Download mapping indexes and other files from (https://sid.erda.dk/share_redirect/F1j2xb0jdB), copy to the folder CLAP, unpack and and merge with folder CLAP/resources (it should happen automatically)

6. Set paths in scripts/CLAP.sh 
Mapper: bwa-pssm (set path of executable) <BR>

## 3. FURTHER CONFIGURATION
Currently the pipeline is set up with an hg19 assemly and a processed ENSEMBL annotation. If one wished to analyze data from a different species or use a different annotation it has to be integrated following a number of steps. The scripts we provide assumes an ENSEMBL annotation GTF file, it will most likely not work with other types of anotation.

### Updating annotation and mapping indexes:
The script "scripts/make_annotation.sh" contains commands to download and process annotation and sequence files from Ensembl version 87. As Ensembl alters their data formats slightly across versions these script may need to be updated, they have been test on selected versions back to version 70. To create annotation for other species or other versions, configure the script by changing lines

        #Ensembl Version
        ver=87
        ...
        # Species name
        species="homo_sapiens"

You run the script like this
        
        scripts/make_annotation.sh

After creating new annotation it is neccessay to configure the pipeline and create new BWT-indexes for BWA-PSSM. This will be described in the following.

### Configure the pipeline such that your annotation and sequence file will be used (open scripts/CLAP.sh)

Change lines:

        idx1=$BASE"/resources/hg19.fa"
        exon_annot=$BASE"/resources/ensembl70.all.long_nooverlap.txt"
        idx3=$BASE"/resources/ensembl70_ej.fa"

Assuming the annotation and sequence files have been moved to resource folder, swap lines by

        idx1=$BASE"/resources/ensembl.mus_musculus.87.fa"
        exon_annot=$BASE"/resources/ensembl.mus_musculus.87.all.txt"
        idx3=$BASE"/resources/ensembl.mus_musculus.87_ej.fa"

Processed and tested annotation and sequence files can be found here

#### Share-links
<b><i>Homo sapians</i></b><BR>
https://sid.erda.dk/share_redirect/e0a4KihL8C
        
<b><i>Mus musculus</i></b><BR>
https://sid.erda.dk/share_redirect/BESvrerfdQ

<b><i>Drosophila melanogaster</i></b><BR>
https://sid.erda.dk/share_redirect/EIZDWhOQkU

<b><i>Caenorhabditis elegans</i></b><BR>
https://sid.erda.dk/share_redirect/c8tlU073R9

<b><i>Rattus norvegicus</i></b><BR>
https://sid.erda.dk/share_redirect/BB5p5bQ1lR

<b><i>Saccharomyces cerevisiae</i></b><BR>
https://sid.erda.dk/share_redirect/aHhDj6y8pW

They are retrieved to the CLAP-dir by

        wget <share-link>
        # unpack the compressed file
        tar -xvzf <tar.gz-file>
        
## 4. TEST-EXAMPLE
To test that everything works, run:

        scripts/testCLAP.sh

It maps, does peak calling and produces a UCSC custom track from reads that map to chr4 in the PAR-CLIP dataset SRR248532.

## 5. USAGE
All scripts are provided as are and will not be maintained or supported.<BR>
To get a help menu run:
        
        scripts/CLAP.sh

14 options have to be specified in the sequence presented in the help menu.<BR>
ARGUMENTS:<BR>

    $1: Filename
    $2: Remove adapters?
        0: No
        Input the adaptor sequence, ex: ACCTGCA...
    $3: Sequence fixed barcode
    $4: length of random barcode
    $5: Remove duplicates?
        0: No
        1: Yes
    $6: Type of analysis
        1: fixed barcode, random barcodes
        2: no fixed barcode, no random barcodes
        3: only fixed barcode, no random barcodes
        4: no fixed barcode, only random barcodes
    $7: UCSC Custom Tracks (bed tracks)
        0: No UCSC custom tracks
        1: UCSC custom tracks
    $8: Stranded protocol?
        0: Strandless
        1: Stranded
    $9: Index
        1: Genome index
        2: Genome index + exon junction index
    $10: Model
        0: Model T>C conversions (PAR-CLIP), conversion prob 0.125
        1: No model (RNA-Seq, iCLIP or HITS-CLIP)
    $11: Output name
    $12: Quality scores
        0: Phread 64
        1: Phread 33
    $13: Number of threads?
        Input number of threads
    $14: Peak calling?
        0: No
        1: Yes

Example analyses:

    PAR-CLIP (substitution model and no barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 0 Analysis_name 1 8 1

    HITS-CLIP (no substitution model and no barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 1 Analysis_name 1 8 1

    iCLIP (with multiplexing and duplication barcodes)
    scripts/CLAP.sh fastq-file TCGTATGCCGTCTTCTGCTTG GGTT 5 1 1 1 1 2 1 Analysis_name 1 8 1

The default substitution model has a T to C conversion rate at 12,5 %. A substitution model with different conversion probability can be created with the script scripts/mk_errorModel.py or the more general script where conversions from and to any nucleotide can be specified (See the repository of BWA-PSSM). <BR>

## 6. HOW TO CITE<BR>
M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. PLOS Computational Biology (in review)<BR>

## 7. LICENSE<BR>
Copyright (c) 2017, Simon H. Rasmussen. The software is open source and released under the MIT license.
