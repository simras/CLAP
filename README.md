CLIP Analysis Pipeline - CLAP, is used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data. The family of protocols where cross-linking immunoprecipitation coupled with high-throughput sequencing are referred to as CLIP. They are used identify binding sites of RNA-binding proteins and exist in the different flavours mentioned above. The data produced by these protocols have different characteristics which are account for in this pipeline such that it can be used to process the above three versions and combination of them like iCLIP using 4SU nucleosides (PAR-iCLIP).

## 1. SYSTEM EQUIREMENTS
A relatively powerful PC, minimun 4 GB memory and with Linux installed. The pipeline is tested on Ubuntu 16.04, in principle it should run on a standard Linux setup with BASH Shell and some version of awk installed. We will not guarantee that the
pipeline works on Apple computers, it could as they are based on FreeBSD, but we provide no support to the end of making it run on a Mac.

## 2. INSTALLATION AND CONFIGURATION
1. Install Python (link https://www.python.org/)

2. Install bwa-pssm (https://github.com/pkerpedjiev/bwa-pssm)

3. Install bedTools (http://bedtools.readthedocs.io/en/latest/)

4. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo)

5. Download mapping indexes and other files from (https://sid.erda.dk/share_redirect/F1j2xb0jdB), copy to the folder CLAP, unpack and and merge with folder CLAP/resources (it should happen automatically)

6. Set paths in scripts/CLAP.sh 
Mapper: bwa-pssm (set path of executable) <BR>

## 3. TEST-EXAMPLE
To test that everything works, run:

        CLAP/testCLAP.sh

It maps, does peak calling and produces a UCSC custom track from reads that map to chr4 in the PAR-CLIP dataset SRR248532.

## 4. USAGE
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

## 5. HOW TO CITE<BR>
M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. PLOS Computational Biology (in review)<BR>

## 6. LICENSE<BR>
Copyright (c) 2017, Simon H. Rasmussen. The software is open source and released under the MIT license.
