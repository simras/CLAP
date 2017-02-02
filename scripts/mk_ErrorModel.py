#!/usr/bin/python
#   mk_ErrorModel.py -m 0.125
#   Example 
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   


def wLines(mutP):
    # range of qualities 0...41
    for qual in range(42):
        for base in ["A","C","G","T"]:
            # quality base P(a|base) P(c|base) P(g|base) P(t|base)
            print qual, base, item(qual,"A",base,mutP),item(qual,"C",base,mutP),item(qual,"G",base,mutP),item(qual,"T",base,mutP)

def item(qual,base,obsBase,mutP):
    import math
    bg = [0.25,0.25,0.25,0.25]
    logC = 1/math.log(2)
    return logC * math.log(sumProb(qual,base,obsBase,mutP)/bg[base2I(base)])


def base2I(base):
    if   base == "A":
        return 0
    elif base == "C":
        return 1
    elif base == "G":
        return 2
    elif base == "T":
        return 3

def error(qual,f,to):
    qual = int(qual)
    quals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"

    if f == to:
        # no error
        return 1 - 10**(-qual/10.0)
    else:
        # an error
        return 10**(-qual/10.0)/3
    
def sumProb(qual,f,to,mutP):
    # sums the probabilities
    # regularizing probability
    pr = 0.0000001
    for a in ["A","C","G","T"]:
        pr = pr + tTOc(f,a,mutP) * error(qual,a,to)
    return pr

def tTOc(f,to,mutP):
    # Calc probability that the base to is a mutation of f. P(to|f)
    p = mutP
    if    f == "T" and to == "C":
        return p
    elif  f == "C" and to == "C":
        return 1 - p
    elif f == to:
        # No mutation
        return 1
    else:
        # No mutation and f and to are different that's impossible
        return 0


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-m", action="store", type="float",    dest="mutP",   default=0.125,  help="P(\"t\"|\"c\"): t to c mutation probability")


    (options, args) = parser.parse_args()


    wLines(options.mutP)
