#!/usr/bin/python
#   adapterHMM2_para.py
#   Example adapterHMM2_para.py -f infile.fastq -o out.fastq -l 20 -s TCGTATGCCGTCTTCTGCTTG -p 35 -n
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen

def run(orig_file,cut_file,cutoff,mp,model,seq,ntrim,trim,BS,prime,qtype):
    import os
    import subprocess 
    import sys

    cwd =  os.getcwd()
    adr = cwd + "/../scripts/"
    import random
    if cut_file == "outF":
        file_parts= orig_file.split(".")
        #print file_parts, file_parts[0:-1], file_parts[-1]
        cut_file = ".".join(file_parts[0:-1]) + "_noAdapt." + file_parts[-1]
        #print cut_file
    if prime:
        primeopt = "-5"
    else:
        primeopt = ""
    rannum = random.randint(0, 1000000)
    if seq != "":
        model = str(rannum) + "_adaptor_hmm.imod"
        os.system(adr + "mk_model.py " + primeopt + " -s " + seq + " > " + model)
    # Which platform?
    p1 = subprocess.Popen(["uname", "-s"],stdout=subprocess.PIPE)
    kernel = p1.communicate()[0].strip()
    print "adapterHMM_multiPY.py: Kernel",kernel
    if kernel  == "Linux":
        decodeanhmm = "decodeanhmm"
    elif kernel == "Darwin":
        decodeanhmm = "decodeanhmm_mac"
    else:
        print >>sys.stderr, "adapterHMM_multiPY.py: C binaeies are only compiled for Linux and Mac Unix platforms. This exception can be manually overwritten"
    # overwrite by uncommenting one of the below lines
    # decodeanhmm = "decodeanhmm"
    # decodeanhmm = "decodeanhmm_mac"

    blockSize = str(BS)
    outpipe = " > "
    if model == "":
        model = adr + "AdapterHMM.imod"
    if orig_file[-2:] == "gz":
        # Construct pipeline command: unzip, convert to fasta, decode, cut, return fastq file and gzip
        cmd = "zcat " + orig_file + " | "+ adr +"trimQs.py " + primeopt + " -l " + str(cutoff) + " -q " + str(qtype) + " " + ntrim + trim +"| awk \'1 == NR % 4,2 == NR % 4\' |" + adr + "multiPY.py -e -p " + str(mp) + " -b " + str(blockSize) + " -l 2 -c \" " + adr + decodeanhmm  + " -v -PrintNumbers -modelfile " + model + "\"" + " 2> /dev/null | " + adr + "analyzeSignals.py | " + adr + "cutIT.py -f " + orig_file + " -c " + str(cutoff) + " | gzip" + outpipe + cut_file
    else:
        # Construct pipeline command: convert to fasta, decode, cut and return fastq file
        cmd = "cat " + orig_file + " | "+ adr +"trimQs.py " + primeopt + " -l " + str(cutoff) + " -q " + str(qtype) + " " + ntrim + trim + "| awk \'1 == NR % 4,2 == NR % 4\' |" + adr + "multiPY.py -e -p " + str(mp) + " -b " + str(blockSize) + " -l 2 -c \" " + adr + decodeanhmm + " -v -PrintNumbers -modelfile " +  model + "\"" + " 2> /dev/null | " + adr + "analyzeSignals.py " + primeopt + " | " + adr + "cutIT.py -f " + orig_file + " " + primeopt + " -c " + str(cutoff) + outpipe + cut_file
    #print cmd
    os.system(cmd)
        
if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", action="store", type="string", dest="fastqF",   default="", help="input fastq file", metavar="FILE") 
    parser.add_option("-m", action="store", type="string", dest="mod",      default="", help="input model file", metavar="FILE")
    parser.add_option("-s", action="store", type="string", dest="mk_model", default="", help="Make Model File based on given sequence. If you give a sequence the -m option will be ignored.")
    parser.add_option("-o", action="store", type="string", dest="outF",     default="outF", help="output file", metavar="FILE") 
#    parser.add_option("-a", action="store", type="string", dest="adapter", default="", help="adapter sequence") 
    parser.add_option("-l", action="store", type="int",    dest="cut",      default=20,  help="length cutoff") 
    parser.add_option("-q", action="store", type="int",    dest="qtype",   default=35, help="type quality scores default illumina Phred+33 <35>, other common Phred+64 <66>, Sanger <33>, Solexa+64 <59>")
    parser.add_option("-b", action="store", type="int",    dest="BS",       default=100000,  help="block size in lines") 
    parser.add_option("-p", action="store", type="int",    dest="mult_p",   default=1,  help="number of different processes") 
    parser.add_option("-n", action="store_true",           dest="notrim",   default=False, help="Don't trim base calls with low quality scores") 
    parser.add_option("-5", action="store_true",          dest="fiveprime",   default=False, help="5 prime adapter") 
    parser.add_option("-t", action="store", type="int",    dest="trim",     default=0, help="trim n bases from 5' end")
#    parser.add_option("-n", action="store", type="int",    dest="numl",   default=0,  help="number of lines") 


    (options, args) = parser.parse_args()
    if options.notrim:
        ntrim = " -n "
    else:
        ntrim = ""
    if options.trim > 0:
        trim = " -t " + str(options.trim)
    else:
        trim = ""
    #if options.mk_model != "":
    run(options.fastqF,options.outF,options.cut,options.mult_p,options.mod,options.mk_model,ntrim,trim,options.BS,options.fiveprime,options.qtype)
