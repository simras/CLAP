CAP - Pipeline used to analyse CLIP-seq (specifically PAR-CLIP, HITS-CLIP and iCLIP) data.  

## 1. SYSTEM EQUIREMENTS
Minimun 4 GB memory.
Linux shell with BASH and some version of awk installed

## 2. INSTALL AND CONFIGURE
1. Install Python (link https://www.python.org/)

2. Install bwa-pssm (https://github.com/pkerpedjiev/bwa-pssm)

3. Install bedTools (http://bedtools.readthedocs.io/en/latest/)

4. Install pyicos (https://bitbucket.org/regulatorygenomicsupf/pyicoteo)

5. Download mapping indexes and other files from (https://sid.erda.dk/share_redirect/FOATbg5v14), copy to the folder CAP, unpack and and merge with folder CAP/resources

6. Set paths in scripts/CAP.sh 
   Mapper: bwa-pssm (find path of executable)
   Mapping index: Set minimum read length for the species in question


## 3. USAGE
All scripts are provided as are and will not be maintained or supported.

## 4. HOW TO CITE

M Plass, SH Rasmussen and A Krogh. Highly accessible AU-rich regions in 3′ untranslated regions are hotspots for binding of proteins and miRNAs. PLOS Computational Biology (in review)

## LICENSE
Copyright (c) 2017, Simon H. Rasmussen. The software is open source and released under the MIT license.
