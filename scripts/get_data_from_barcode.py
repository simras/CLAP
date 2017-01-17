#!/usr/bin/python
#   agrep.py
#   Example cat seqs.fa | agrep.py -p gcttcatagccccttcccaat -m 3
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   gcttcatagccccttcccaat


def eval(filename,barcode,offset,allowed_mismatch):
    import os
    import sys
    first = False
    if filename != "":
        file1 = open(filename,'r')
        first = True
    rec = []
    recs = []
    t = 0
    i = 0
    First = True
    len_barcode = len(barcode)
    while True:
        if first:
            l = file1.readline()
        else:
            l = sys.stdin.readline()
        if i % 4 == 0 and not First:
            seq = recs[1]
            j = offset
            subseq = seq[j:(j+len_barcode)]
            mis = 0
            k = 0
            for k in range(len_barcode):
                if subseq[k] != barcode[k]:
                    mis = mis + 1
            if mis <= allowed_mismatch:
                sequence = (recs[1])[0:offset] + (recs[1])[(offset+len_barcode):]
                try: 
                    ID = strip(recs[0].split("length=")[0])
                    leng = str(len(sequence))
                except:
                    ID = recs[0]
                    leng = ""

                print ID + leng ,sequence,ID.replace("@","+"),(recs[3])[0:offset] + (recs[3])[(offset+len_barcode):],
            recs = []
        if (i + 1) % 1000000 == 0:
            print  >> sys.stderr, ".",
            if (i + 1) % 10000000 == 0:
                print  >> sys.stderr, str(((i+1)/4000000.0)) + " M reads"
        i = i + 1 
        if i > 3:
            First = False
        if not l:
            break
        recs.append(l)
    #seq = recs[1]
    #j = offset
    #subseq = seq[j:(j+npatt)]
    #mis = 0
    #k = 0
    #for k in range(npatt):
    #    if subseq[k] != patt[k]:
    #        mis = mis + 1
    #if mis <= amm:
    #    print recs[0],recs[1],recs[2],recs[3],
            
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-f", action="store", type="string", dest="file",default="", help="input pred file")
    parser.add_option("-p", action="store", type="string", dest="barcode",default="", help="Barcode")
    parser.add_option("-o", action="store", type="int", dest="offset",default=0, help="offset")
    parser.add_option("-m", action="store", type="int", dest="allow_mismatch",default=1, help="allowed mismatches")
    
    (options, args) = parser.parse_args()
    
    eval(options.file,options.barcode,options.offset,options.allow_mismatch)

    
