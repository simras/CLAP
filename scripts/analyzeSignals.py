#!/usr/bin/python
#   agrep.py
#   Example cat seqs.fa | agrep.py -p gcttcatagccccttcccaat -m 3
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   gcttcatagccccttcccaat

from types import *

def analyzeSignals(fn):
    import os
    import sys 
    from operator import itemgetter

    cwd = os.getcwd()
    if fn != "":
        file1 = open(cwd + "/" + fn,'r')
        lines = file1.readlines()
    else:    
        lines = sys.stdin.readlines()
    i = 0
    j = 0
    start = -1
    lst = []
    first = True
    for l in lines:
        if l[0] == ">":
            i = i + 1
            if not first:
                print ID,start
                if pred == 0:
                    j = j + 1
            ID = l[1:].rstrip()
            pred = 0
            first = False
            firstSig = True
        elif l[0:4] == "%GFF":
            if firstSig:
                tmp_pred = float(l.split()[6]) + .15
            else:
                tmp_pred = float(l.split()[6])
            if tmp_pred > pred:
                pred = tmp_pred
                start = int(l.split()[4])
                firstSig = False
    print >> sys.stderr, "analyzeSignals.py: Number of predictions with no Signal: ", j, i#, lst[0]
#    print sorted(lst, key=itemgetter(1))[0:10]
#    print sorted(lst, reverse=True,key=itemgetter(1))[0:10]

def analyzeSignals_viterbi(fn,prime):
    import os
    import sys
    from operator import itemgetter

    cwd = os.getcwd()
    if fn != "":
        file1 = open(cwd + "/" + fn,'r')
        sin = False
    #    lines = file1.readlines()
    else:    
        sin = True
    i = 0
    j = 0
    start = -1
    #lst = []
    first = True
    l = '#'
    ID = ""
    printID = True
    ldist = {}
    while True:
#        print "#1" + l + "#2"
        if sin:
            l = sys.stdin.readline()            
        else:
            l = file1.readline()
        if l == '':
            break
        if first and l.strip() == "":
            printID = False
        if l[0] == ">":
            printID = True
            i = i + 1
            if not first:
                print ID,start
                if ldist.has_key(lseq - start + 1):
                    ldist[lseq - start + 1] = ldist[lseq - start + 1] + 1
                else:
                    ldist[lseq - start + 1] = 1
            ID = l[1:].split()[0].rstrip()
            first = False
        elif l[0:4] == "%len":
            lseq = int(l.split("len ")[1])
            #print lseq
        elif l[0:7] == "%pred V":
            try:
                if not prime:
                    start = int(l.split("s")[1].split()[0])
                else:
                    # start is actually the end of the adapter
                    start = lseq - int(l.split("b")[-1].split()[0])
            except:
                
                if not prime:
                    start = int(l.split("b")[-1].split()[-1]) + 1
                else:
                    start = 0
                j = j + 1
            printID = True
        
    if start != -1 and printID:
        print ID,start
    for k,v in sorted(ldist.items(),key=itemgetter(0)):
        if k == 0:
            print >> sys.stderr, "analyzeSignals.py: Sequences with no adaptor:",v
        else:
            print >> sys.stderr, "analyzeSignals.py: Adaptor length",k,"number",v
                        
    print >> sys.stderr, "analyzeSignals.py: Sequences with adapter", i-j, "all sequences",i
#    print sorted(lst, key=itemgetter(1))[0:10]
#    print sorted(lst, reverse=True,key=itemgetter(1))[0:10]
        

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()

    parser.add_option("-p", action="store", type="string", dest="predF",default="", help="input pred file") 
    parser.add_option("-5", action="store_true", dest="prime",   default=False,  help="5prime adapter?") 


    (options, args) = parser.parse_args()

    analyzeSignals_viterbi(options.predF,options.prime)
