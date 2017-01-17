#!/usr/bin/python
#   mk_model.py
#   Example mk_model.py cgtcgtggacacactt
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   

def mk_model(fiveprime,seq,long_reads):
    import os
    seq = seq.upper()
    if long_reads:
        after =  "after:0.01"
        before = "before"
    else:
        after =  "bg:0.01"
        before = "bg"

    length = str(len(seq))
    if fiveprime:
        print "# Adapter HMM - A profile HMM of a 5' adaptor sequence in short read sequenzing"
    else:
        print "# Adapter HMM - A profile HMM of a 3' adaptor sequence in short read sequenzing"
    print "header { alphabet ACGT;wildcards N; } # The alphabet"
    print 
    print "#begin state"
    print "begin {trans " + before + " s1 d1;}"
    print
    if seq[0].lower() != "n":
        print "s1   {trans s2   i1:0.001      d2:0.001;  letter " + seq[0] + ":0.97:10;label s; signal X:0:" + length + ";}"
    else:
        print "s1   {trans s2   i1:0.001      d2:0.001;label s; signal X:0:" + length + ";}"
    for i in range(2,len(seq)):
        base = seq[i-1]
        n = i
        if True:
            if i < 3:
                p="0.97"
            else:
                p = "0.99"
            if fiveprime:
                after = "bg"
            else:
                after = ""
        else:
            p = "0.99"
            if long_reads:
                after =  "after:0.01"
                #before = "before"
            else:
                after = "bg:0.01"
                #before = "bg"

        if base.lower() != "n":
            print "s" + str(i) + "  {trans s" + str(i + 1) + "  i" + str(i) + ":0.001 " + after + "  d" + str(i + 1) + ":0.01; letter " + base + ":" + p + ":10;label s;}"
        else:
            print "s" + str(i) + "  {trans s" + str(i + 1) + "  i" + str(i) + ":0.001 " + after + "  d" + str(i + 1) + ":0.01;label s;}"
    n = n + 1
    if seq[-1].lower() != "n":
        if long_reads:
            print "s" + str(n) + "  {trans s1 d1 " + " after" + "             ; letter " + seq[-1] + ":0.99:10;label s;}"
        else:
            print "s" + str(n) + "  {trans s1 d1 " + " bg" + "             ; letter " + seq[-1] + ":0.99:10;label s;}"
    else:
        if long_reads:
            print "s" + str(n) + "  {trans s1 d1 " + " after" + "             ;label s;}"
        else:
            print "s" + str(n) + "  {trans s1 d1 " + " bg" + "             ;label s;}"
    print 
    print "# insertions"
    for i in range(1,len(seq)):
        print "i" + str(i) + " {trans   i" + str(i) + ":0.001  s" + str(i + 1) + "; label s;}"
    print
    print "# deletions"
    for i in range(1,len(seq)):
        print "d" + str(i) + "   {trans d" + str(i + 1) + ":0.001   i" + str(i) + ":0.001  s" + str(i + 1) +";   letter NULL;}"
    n = i + 1
    if long_reads:
        print "d" + str(n) + "  {trans s1 after          ;   letter NULL;}"
    else:
        print "d" + str(n) + "  {trans s1 bg             ;   letter NULL;}"
    
    print
    print "# Background"
    
    if fiveprime:
        print " bg{trans bg;order 2;label b;}"
    else:
        if long_reads:
            print "before {trans s1:0.01 before   ; order 2; label b;}"
            print "after  {trans s1:0.00001 after ; order 2; label b;}"
        else:
            print " bg{trans s1 bg;order 2;label b;}"
    #print "#e {letter NULL; end 1; }"
    print
    print "# End of model"

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-s", action="store", type="string", dest="seqs",default="", help="input model file")
    parser.add_option("-l", action="store_true", dest="long_reads",default=False, help="Set flag if the adaptor in many cases would have sequence after that could be modeled")
    parser.add_option("-5", action="store_true", dest="fiveprime",default=False, help="remove 5' adapter")
#   parser.add_option("-n", action="store", type="int",    dest="numl",   default=0,  help="number of lines") 


    (options, args) = parser.parse_args()

    mk_model(options.fiveprime,options.seqs,options.long_reads)
