#!/usr/bin/python
#   multiPY.py
#   Example cat seqs.fa | multiPY.py -c wc -p 3
#   By Simon H. Rasmussen
#   Bioinformatics Centre
#   University of Copenhagen
#   Implementation of multiprocess batch processing

import os
import sys
import threading

def multiPY(in_file,command, lines_in_record, procs,batchSize=100000,err=False):
    import subprocess
    import random
    ran = str(random.randint(0,10000))
    stopi = procs
    if in_file == "":
        fil = False
    else:
        fil = True
        inF = open(in_file,'r')
    stop = False
    stopp=[False]
    k = 0
    while not stop:
        t = []
        for j in range(procs):
            lines = ""
            ## Popen write out batchSize lines to file
            for l in range(lines_in_record * batchSize):
                if fil:
                    line = inF.readline()
                else:
                    line = sys.stdin.readline()
                if not line:
                    stop = True
                    stopi = j
                    break
                lines = lines + line
            if lines != "":
                t.append(threading.Thread(target=thread_it,args=[j,stopi,lines,command,lines_in_record,err,ran]))
                t[j].start()
            #print "Indices ",k, j
            if stop:
                break
        k = k + 1
        for j in range(len(t)):
            t[j].join()
            out_file_name = "out_file_" + ran + "_" + str(j) + ".txt" 
            out_file = open(out_file_name,'r')
            lines = out_file.readlines()
            
            for l in lines:
                print l,
            out_file.close()
            
#    print "Clean up"
    for j in range(procs):
        tmp_file_name = "tmp_file_" +  ran + "_" + str(j) + ".txt" 
        os.system("rm " + tmp_file_name + " 2> /dev/null")
        out_file_name = "out_file_" +  ran + "_" + str(j) + ".txt" 
        os.system("rm " + out_file_name + " 2> /dev/null")

def thread_it(j,stopi,lines,command,lines_in_record,err,ran):
    import subprocess
    import random
    if j > stopi:
        return
    out_file_name = "out_file_" + ran + "_" + str(j) + ".txt" 
    tmp_file_name = "tmp_file_" + ran + "_" + str(j) + ".txt"
    tmp_file = open(tmp_file_name,'w')
    tmp_file.writelines(lines)
    tmp_file.close()
    p = subprocess.Popen(["cat " + tmp_file_name + " |" + command + " > " + out_file_name], shell=True)
    p.wait()

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", action="store", type="string", dest="fname",default="", help="input file")
    parser.add_option("-l", action="store", type="int",    dest="lines",   default=0,  help="Number of lines in one recored") 
    parser.add_option("-b", action="store", type="int",    dest="batchSize",   default=100000,  help="Batch Size") 
    parser.add_option("-p", action="store", type="int",    dest="procs",   default=0,  help="Number of processes")
    #    parser.add_option("-p", action="store", type="string", dest="patt", default="", help="pattern") 
    #   parser.add_option("-l", action="store", type="string", dest="label", default="x", help="pattern")
#    parser.add_option("-o", action="store", type="string", dest="out_file", default="", help="output File")
    parser.add_option("-c", action="store", type="string", dest="command", default="", help="Command")
    parser.add_option("-e", action="store_true", dest="err", default=False, help="Suppress std. error")


    (options, args) = parser.parse_args()

    multiPY(options.fname,options.command,options.lines,options.procs,options.batchSize,options.err)



