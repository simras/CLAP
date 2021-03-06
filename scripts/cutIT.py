#!/usr/bin/python
#   cutIT.py
#   Example cat seqs.fastq | cutIT.py 
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen

from types import *

def clearFile(fn):
    import os
    import sys 
    
    cwd = os.getcwd()
    if fn != "":
        file1 = open(cwd + "/" + fn,'r')
        sin = False
    else:    
        sin = True
    i = 0
    j = 0
    second = False
    l = '#'
    while l != "":
        if sin:
            sys.stdin.readline()
        else:
            l = file1.readline()

        if i % 2 == 0 and len(l.split()) > 1 or second and len(l.split()) > 1:
            second = True
            j = j + 1
            continue
        second = False
        if len(l.split()) == 1:
            print l,
        else:
            print l,
        i = i + 1
    print >> sys.stderr, "Number of : ", j

        
def cutIT_fasta(fastq,predfile,cutoff):
    import os
    import sys
    line = ""
    if fastq != "":
        file1 = open(fastq,'r')    
        lines1 = file1.readlines()
    else:
        lines1 = sys.stdin.readlines()   
    file2 = open(predfile,'r')
    lines2 = file2.readlines()
    numLines = len(lines2)
    for i in range(numLines):
        full_ID = lines1[2*i]
        ID = int(full_ID.split(".")[1])
        seq = lines1[2*i+1]
        l = lines2[i].split()
        ID2   = int(l[0].split(".")[1])
        start = int(l[1])
        if ID == ID2:
            if not full_ID == "" and not seq[0:start-1] == "" and start > cutoff:
                print full_ID,
                print seq[0:start-1]

def cutIT_fastq(fastq,predfile,cutoff,prime):
    import os
    import sys
    import gzip
    line = ""
    if predfile != "":
        file2 = open(predfile,'r')
        sin = False
    else:    
        sin = True

    if fastq[-2:] == "gz":
        file1 = gzip.open(fastq,'r')
    else:
        file1 = open(fastq,'r')
    
    i = 0
    k = 0
    line1 = '#'
    line2 = '#'
    next = True
    while True:
        if sin:
            if next:
                #Read cutfile
                cutline = sys.stdin.readline()
        else:
            if next:
                #Read cutfile
                cutline = file2.readline()
        line1 = file1.readline()
        if cutline.strip() == '':
            break
        full_ID = line1
        
        ID = (full_ID.split()[0])[1:]
        line1 = file1.readline()
        seq = line1
        line1 = file1.readline()
        l3 = line1
        line1 = file1.readline()
        l4 = line1
        
        l = cutline.split()
        ID2 = l[0]
        start = int(l[1])
        if prime:
            ll = len(seq)
        if ID == ID2:
            next = True 
            if not prime and not full_ID == "" and seq[0:start-1] != "" and start > cutoff:
                if "length" in full_ID:
                    firstLine = "".join(full_ID.split("length")[0:-1]) + "length=" + str(start-1)
                else:
                    firstLine = full_ID.strip() + " length=" + str(start-1)
               
                print firstLine
                print seq[0:(start-1)].strip()
                print "+" + firstLine[1:]
                print l4[0:(start-1)].strip()
                i = i + 1
            elif prime and not full_ID == "" and seq[ll - start - 2:] != "" and start + 1 >= cutoff and start  < ll:
                if "length" in full_ID:
                    oldll = int(full_ID.split("length=")[-1].split()[0])
                    firstLine = "".join(full_ID.split("length")[0:-1]) + "length=" + str(start + 1)
                else:
                    firstLine = full_ID.strip() + " length=" + str(start+1)
               
                print firstLine
                print seq[ll - start - 2:].strip()
                print "+" + firstLine[1:]
                print l4[ll - start - 2:].strip()
                i = i + 1
            else:
                k = k + 1
        else:
            k = k + 1
            next = False
    
    print >> sys.stderr,"cutIT.py: Sequences filtered out",k
    print >> sys.stderr,"cutIT.py: Sequences printed" ,i 
         
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", action="store", type="string", dest="fastqF",default="", help="input fastq file") 
    parser.add_option("-p", action="store", type="string", dest="predF",default="", help="input pred file") 
    parser.add_option("-c", action="store", type="int",    dest="cut",   default=0,  help="length cutoff") 
    parser.add_option("-5", action="store_true", dest="prime",   default=False,  help="5prime adapter?") 

    (options, args) = parser.parse_args()

    if options.fastqF[-5:] == "fastq" or options.fastqF[-2:] == "gz":
        cutIT_fastq(options.fastqF,options.predF,options.cut,options.prime)
    else:
        cutIT_fasta(options.fastqF,options.predF,options.cut)
