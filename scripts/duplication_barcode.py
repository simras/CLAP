#!/usr/bin/python
#   duplication_barcode.py
#   Example cat seqs.fa | duplication_barcode.py 
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   gcttcatagccccttcccaat

from types import *

def expected_num_errors_in_barcode(quals):
    ssum=0
    qs = "!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
    for q in quals:
        #print q, 10**(qs.index(q)/-10)
        ssum = ssum + (10**(qs.index(q)/-10))
    return ssum
#expected_num_errors_in_barcode("#2&7J")

def eval(fn,trim):
    import os
    import sys 
    from operator import itemgetter

    if fn != "":
        sin = False
        file1 = open(fn,'r')
    else:
        sin = True
    rec = []
    fHash = {}
    t = 0
    i = 0
    tt = 0
    ttt = 0
    First = True
    while True:
        if sin:
            l = sys.stdin.readline()
        else:
            l = file1.readline()
        if (i + 1) % 1000000 == 0:
            if (i + 1) %   10000000 == 0:
                print  >> sys.stderr, str((i/4000000.0)) + " M reads"
            else:
                print  >> sys.stderr, ".",
            
        if i % 4 == 0:
            if not First:                 
                if fHash.has_key(recs[1]):
                    tmp = fHash[recs[1]].split("\t")
                    c = int(tmp[0])
                    c = c + 1
                    oldrecs = tmp[1:]
                    if ord((oldrecs[1])[-2]) > ord((recs[3])[-2]):
                        fHash[recs[1]] = str(c) + "\t" + oldrecs[0] + "\t" + oldrecs[1]
                    else:
                        fHash[recs[1]] = str(c) + "\t" + recs[0] + "\t" + recs[3]         
                    tt = tt + 1
                else:
                    fHash[recs[1]] = "0" + "\t" + recs[0] + "\t" + recs[3]
                    t = t + 1
            recs = rec
            rec = []

        if i > 3:
            First = False
        if not l:
            break
        rec.append(l)
        i = i + 1
    if fHash.has_key(recs[1]):
        tmp = fHash[recs[1]].split("\t")
        c = int(tmp[0])
        c = c + 1
        oldrecs = tmp[1:]
        if ord((oldrecs[1])[-2]) > ord((recs[3])[-2]):
            fHash[recs[1]] = str(c) + "\t" + oldrecs[0] + "\t" + oldrecs[1]
        else:
            fHash[recs[1]] = str(c) + "\t" + recs[0] + "\t" + recs[3]         
        tt = tt + 1
    else:
        fHash[recs[1]] = "0" + "\t" + recs[0] + "\t" + recs[3]
        t = t + 1
    jj = 0
    jjj = 0
    rr = 0
    for k,v in fHash.iteritems():
        tmp = v.split("\t")
       
        qualities =  (tmp[2])[(trim):]
        c = int(tmp[0])
        if trim > 0:
            barcode = k[0:trim]
        sequence = k[trim:]
        ID = "".join(tmp[1])
        if trim > 0:
            print str(ID.strip()) + "_BC=" + barcode +"_duplications=" + tmp[0] + "\n" + sequence,"+\n" + qualities,
        else:
            print str(ID.strip()) + "duplications=" + tmp[0] + "\n" + sequence,"+\n" + qualities,
        if c > 0:
            jjj = jjj + 1
        rr = rr + c
        jj = jj + 1
        
    print >> sys.stderr, "duplication_barcode.py: Total reads", i/4, " unique reads ",jj," affect by duplication ",jjj
    
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-f", action="store", type="string", dest="predF",default="", help="input pred file")
    parser.add_option("-t", action="store", type="int", dest="trim",default=0, help="Trim") 
    
    
    (options, args) = parser.parse_args()
    
    eval(options.predF,options.trim)

    
