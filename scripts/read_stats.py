#!/usr/bin/python
# 
#  Generate mapping statistics given a bwa-PSSM sam-file
#  Implements parsing of Posterior probability and filtering
#



def filter_reads(fname,t,t2,print_IDs,junction_reads):
    import sys
    u = 0
    m = 0
    cp = 0
    cm = 0
    mp = 0
    mm = 0
    i = 0
    a = 0
    for l in sys.stdin:
        
        if l[0] == "@":
            continue
        PP = 0
        a = a + 1
        try:
            ID = l.split()
            read_ID = ID[0]
            bit = int(ID[1])
            r_len = int(len(ID[9]))

            #  Only check certain bit-flags
                    
            # Exon location
            Exon_loc = ID[2]
            rel_start = int(ID[3])

            try:
                rest = l.split("PP:f:")[1]
                PP = float(rest.split()[0])
                m = m + 1
            except:
                PP=0
            revbit = bin(bit)[::-1]
            
            if int(revbit[2]) == 1:
                # Unmapped reads
                u = u + 1
            elif bit == 0:
                # Mapped to the plus strand
                if junction_reads:
                    chr_ID = ID[2]
                    strand = chr_ID.split(":")[1]
                    if strand == "-":
                        if PP > t:
                            # confidently mapped
                            cm = cm + 1
                            
                        else:
                            # multimapped
                            mm = mm + 1
                    else:
                        if PP > t:
                            # confidently mapped
                            cp = cp + 1
                            
                        else:
                            # multimapped
                            mp = mp + 1
                else:
                    if PP > t:
                        # confidently mapped
                        cp = cp + 1
                        #
                    else:
                        # multimapped
                        mp = mp + 1
            elif int(revbit[4]) == 1:
                # Mapped to the minus strand
                if PP > t:
                    # confidently mapped
                    cm = cm + 1
                    
                else:
                    # multimapped
                    mm = mm + 1
            else:
                print >>sys.stderr, "read_stats.py: Warning: Unexpected input"
        except:
            None
    try:
        print >>sys.stderr, "read_stats.py: Confidently mapped ",round(float(cm + cp)/a,2),"Confidently mapped plus ",round(float(cp)/a,2),"Confidently mapped minus ",round(float(cm)/a,2), " Multiple mapped ",round(float(mm + mp)/a,2)," unmapped ",round(float(u)/a,2)
        print "Confidently_Mapped_Minus","Confidently_Mapped_Plus","Multiple_Mapped_Minus","Multiple_Mapped_Plus", "Unmapped_reads","All_reads"
        print cm, cp, mp, mm, u, a
        if (cm + cp + mm + mp + u) != a:
            print >>sys.stderr, "Warning: ",a - (cm + cp + mm + mp + u)," reads were not counted."
    except ZeroDivisionError:
        
        print >>sys.stderr, "No reads mapped"
        print "No reads mapped"
        print cm,cp,mp,mm,u,a
        if (cm + cp + mm + mp + u) != a:
            print >>sys.stderr, "Warning: ",a - (cm + cp + mm + mp + u)," reads were not counted."

if __name__ == "__main__":                                                                                                                                         
    from optparse import OptionParser  
    parser = OptionParser()
    
    parser.add_option("-f", action="store", type="string", dest="f",default="", help="Input file.")
    parser.add_option("-t", action="store", type="float", dest="t",default="0.99", help="PP Threshold conf map vs. mult map")
    parser.add_option("-l", action="store", type="float", dest="l",default="0.01", help="PP low bound")
    parser.add_option("-m", action="store", type="int", dest="m",default="1", help="Mode 1: filter reads, 2: select reads and 3: count correctly mapped.")
    parser.add_option("-i", action="store_true", dest="i",default=False, help="Print IDs of selected reads.")
    parser.add_option("-j", action="store_true", dest="j",default=False, help="Reads mapped to exon junction index.")
    
    (options, args) = parser.parse_args()
        
#    countR(options.f,options.t,options.e,options.l,options.s)
    filter_reads(options.f,options.t,options.l,options.i,options.j)
