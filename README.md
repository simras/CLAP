CLAP - Pipeline used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data.  

## 1. SYSTEM EQUIREMENTS
Minimun 4 GB memory.
Linux shell with BASH and some version of awk installed

## 2. INSTALL AND CONFIGURE
1. Install Python (link https://www.python.org/)

2. Install bwa-pssm (https://github.com/pkerpedjiev/bwa-pssm)

3. Install bedTools (http://bedtools.readthedocs.io/en/latest/)

4. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo)

5. Download mapping indexes and other files from (https://sid.erda.dk/share_redirect/FOATbg5v14), copy to the folder CLAP, unpack and and merge with folder CLAP/resources

6. Set paths in scripts/CLAP.sh 
Mapper: bwa-pssm (set path of executable)
Mapping index: Set minimum read length for the species in question

## 3. USAGE
All scripts are provided as are and will not be maintained or supported.
To get a help menu run:
scripts/CLAP.sh

14 options have to be specified in the sequence presented in the help menu.
ARGUMENTS:
$1: Filename
$2: Remove adapters?
    0: No
    If all datasets have the same 3'adaptor just input the adaptor sequence, ex: ACCTGCA...
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

Example runs:
PAR-CLIP (substitution model and no barcodes)
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 0 <Analysis_name> 1 8 1

HITS-CLIP (no substitution model and no barcodes)
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 1 <Analysis_name> 1 8 1

iCLIP (with multiplexing and duplication barcodes)
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG GGTT 5 1 1 1 1 2 1 <Analysis_name> 1 8 1



## 4. HOW TO CITE
M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. PLOS Computational Biology (in review)

## LICENSE
Copyright (c) 2017, Simon H. Rasmussen. The software is open source and released under the MIT license.
