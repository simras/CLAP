#!/usr/bin/python
#   cutIT.py
#   Example cat seqs.fastq | cutIT.py 
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#

from types import *


def trimQs(fastq,ntrim,trim5,trim3,minl,fivePrime,cl,qtype):
    '''
    Trims low quality positions (only 1-5 ) from the 3' end the 5'
    
    '''
    import os
    import sys
    import gzip
    quals = ""
    quals = chr(qtype + 1) + chr(qtype + 2) + chr(qtype + 3) + chr(qtype + 4) + chr(qtype + 5)
    print >> sys.stderr, "trimQs.py: Phread", qtype
    print >> sys.stderr, "trimQs.py: remove bases with qualities,", quals
    line = ""
    file1 = sys.stdin
    i = 0
    line1 = '#'
    q = 0
    qq = 0
    while True:
        line1 = file1.readline()
        if line1 == '':
            break
        full_ID = line1
        
        #print >> sys.stderr, full_ID
        try:
            ID = (full_ID.split()[0])[1:]
        except:
            print >> sys.stderr, "trimQs.py: Empty line"
            next
        #        ID_rest = (full_ID.split())[1:]
        line1 = file1.readline()
        seq = line1
        line1 = file1.readline()
        l3 = line1
        line1 = file1.readline()
        Qs = line1.strip()
        lseq = len(Qs)
        start = lseq - trim3
        if fivePrime:
            j = trim5
            start = j
            if (lseq - start - 1) >= minl and not ntrim:
                while (lseq - j) > minl and Qs[j] in quals:
                    j = j + 1
                    start = j
                if start != trim5:
                    qq = qq + 1
            #print "start",start
            if (lseq - start - trim3) >= minl and not full_ID == "":
                if cl and "length" in full_ID:
                    firstLine = ("".join(full_ID.split("length")[0:-1])[1:]) + "length=" + str(lseq - start - trim3)
                else:
                    if cl:
                        firstLine = full_ID.strip()[1:] + " length=" + str(lseq - start - trim3)
                    else:
                        firstLine = full_ID.strip()[1:] + " length=" + str(lseq)
                print ">" + firstLine
                print seq[start:(lseq-trim3)].strip()
                print "+" + firstLine
                print Qs[start:(lseq-trim3)].strip()
                
            else:
                q = q + 1
        else:
            # Trim from 3'UTR
            if start >= minl and not ntrim:
                j = lseq - trim3 - 1
                
                while j >= 0 and Qs[j] in quals:
                    
                    j = j - 1
                start = j + 1
                if start != lseq - trim3:
                    qq = qq + 1
            if (start - trim5) >= minl and full_ID != "":
                if "length" in full_ID:
                    firstLine = ("".join(full_ID.split("length")[0:-1])[1:]) + "length=" + str(start-trim5)
                else:
                    if cl:
                        firstLine = full_ID.strip()[1:] + " length=" + str(start-trim5)
                    else:
                        firstLine = full_ID.strip()[1:] + " length=" + str(start)
                print ">" + firstLine
                print seq[trim5:start]
                #print >> sys.stderr,str(start), str(lseq-trim3)
                #print >> sys.stderr, seq[trim5:start]
                #print >> sys.stderr, seq,
                #print >> sys.stderr, Qs
                #print >> sys.stderr, ""
                print "+" + firstLine
                print Qs[trim5:start]
            else:
                q = q + 1
        i = i + 1
    print >> sys.stderr, "trimQs.py:",q,"reads filtered cause they were too short"
    print >> sys.stderr, "trimQs.py:",qq,"reads truncated due to low qualities: "

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", action="store", type="string", dest="fastqF",default="", help="input fastq file")
#    parser.add_option("-o", action="store", type="string", dest="outfastqF",default="", help="input fastq file")
    parser.add_option("-q", action="store", type="int",    dest="qtype",   default=35, help="type quality scores default illumina Phred+33 <35>, other common Phred+64 <66>, Sanger <33>, Solexa+64 <59>")
    parser.add_option("-t", action="store", type="int",    dest="trim5",   default=0, help="trim n chars from 5' end")
    parser.add_option("-s", action="store", type="int",    dest="trim3",   default=0, help="trim n chars from 3' end")
    parser.add_option("-l", action="store", type="int",    dest="minl",   default=20, help="length cutoff")
    parser.add_option("-n", action="store_true",    dest="notrim",   default=False, help="Don't trim") 
    parser.add_option("-5", action="store_true",    dest="fiveprime",   default=False, help="trim qualities from 5' end")
    parser.add_option("-c", action="store_false",    dest="changelength",   default=True, help="change length property") 

    (options, args) = parser.parse_args()

    trimQs(options.fastqF,options.notrim,options.trim5,options.trim3,options.minl,options.fiveprime,options.changelength,options.qtype)
