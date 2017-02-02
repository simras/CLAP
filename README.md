CLAP - Pipeline used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data. The family of protocols where cross-linking immunoprecipitation coupled with high-throughput sequencing are called CLIP. They are used identify binding sites of RNA-binding proteins and exist in the different flavours mentioned above. The data produced by these protocols have different characteristics which are account for in this pipeline such that it can be used to process the above three versions and combination of them like iCLIP using 4SU nucleosides (PAR-iCLIP).

## 1. SYSTEM EQUIREMENTS
A relatively powerful PC, minimun 4 GB memory and with Linux installed. The pipeline is tested on Ubuntu 16.04, in principle it should run on a standard Linux setup with BASH Shell and some version of awk installed.

## 2. STEPS OF INSTALLATION AND CONFIGURATION
1. Install Python (link https://www.python.org/)

2. Install bwa-pssm (https://github.com/pkerpedjiev/bwa-pssm)

3. Install bedTools (http://bedtools.readthedocs.io/en/latest/)

4. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo)

5. Download mapping indexes and other files from (https://sid.erda.dk/share_redirect/FOATbg5v14), copy to the folder CLAP, unpack and and merge with folder CLAP/resources

6. Set paths in scripts/CLAP.sh 
Mapper: bwa-pssm (set path of executable) <BR>
Mapping index: Set minimum read length for the species in question

## 3. USAGE
All scripts are provided as are and will not be maintained or supported.<BR>
To get a help menu run:<BR>
scripts/CLAP.sh<BR>

14 options have to be specified in the sequence presented in the help menu.<BR>
ARGUMENTS:<BR>
$1: Filename<BR>
$2: Remove adapters?<BR>
&nbsp;&nbsp;&nbsp;&nbsp;0: No<BR>
&nbsp;&nbsp;&nbsp;&nbsp;If all datasets have the same 3'adaptor just input the adaptor sequence, ex: ACCTGCA...<BR>
$3: Sequence fixed barcode<BR>
$4: length of random barcode<BR>
$5: Remove duplicates?<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: No<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: Yes<BR>
$6: Type of analysis<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: fixed barcode, random barcodes<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    2: no fixed barcode, no random barcodes<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    3: only fixed barcode, no random barcodes<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    4: no fixed barcode, only random barcodes<BR>
$7: UCSC Custom Tracks (bed tracks)<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: No UCSC custom tracks<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: UCSC custom tracks<BR>
$8: Stranded protocol?<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: Strandless<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: Stranded<BR>
$9: Index<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: Genome index<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    2: Genome index + exon junction index<BR>
$10: Model<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: Model T>C conversions (PAR-CLIP), conversion prob 0.125<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: No model (RNA-Seq, iCLIP or HITS-CLIP)<BR>
$11: Output name<BR>
$12: Quality scores<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: Phread 64<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: Phread 33<BR>
$13: Number of threads?<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    Input number of threads<BR>
$14: Peak calling?<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    0: No<BR>
&nbsp;&nbsp;&nbsp;&nbsp;    1: Yes<BR>

Example runs:<BR>
PAR-CLIP (substitution model and no barcodes)<BR>
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 0 <Analysis_name> 1 8 1<BR>

HITS-CLIP (no substitution model and no barcodes)<BR>
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG "" 0 1 2 1 1 2 1 <Analysis_name> 1 8 1<BR>

iCLIP (with multiplexing and duplication barcodes)<BR>
scripts/CLAP.sh <fastq-file> TCGTATGCCGTCTTCTGCTTG GGTT 5 1 1 1 1 2 1 <Analysis_name> 1 8 1<BR>

## 4. HOW TO CITE<BR>
M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. PLOS Computational Biology (in review)<BR>

## LICENSE<BR>
Copyright (c) 2017, Simon H. Rasmussen. The software is open source and released under the MIT license.
