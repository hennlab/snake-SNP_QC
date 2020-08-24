#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
import gzip
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-R","--report",type="string",dest="report",action="store",help="genome studio report to call from")
parser.add_option("-O","--output",type="string",dest="outputROOT",action="store",help="output root")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

if options.outputROOT == None:
    print "specify output root with -O"
    sys.exit()


outTFAM = open(options.outputROOT + ".tfam", 'w')
outTPED = open(options.outputROOT + ".tped", 'w')
        

if options.report.find(".gz") == -1: # report is not compressed
    # Iterate through Illumina GenomeStudio Report
    n = 0
    for line in open(options.report, 'r'):
        n += 1
        line = line.replace("\n", "")
        line = line.replace("\r", "")

        fields = line.split("\t")

        if n == 1: # header row -- write tfam info
            for i in range(3,len(fields), 3):
                k = fields[i].split(".")
                id = ".".join(k[:len(k)-1])
                out = [id, id, "0", "0", "-9", "-9"]
                outTFAM.write(" ".join(out) + "\n")

        else:
            snp = fields[0]
            chr = fields[1]
            pos = fields[2]

            out = [chr, snp, "0", pos]
            for i in range(3, len(fields), 3):
                gt = fields[i]
		print(i)
		print(n)
		print("gt0", gt[0])
		print("gt1", gt[1]) 
                if gt == "NC":
                    gt = "00"
                out.append(gt[0])
                out.append(gt[1].rstrip())
        
            outTPED.write(" ".join(out) + "\n")

else: # report is gzipped
    # Iterate through Illumina GenomeStudio Report
    n = 0
    for line in gzip.open(options.report):
        n += 1
        line = line.replace("\n", "")
        line = line.replace("\r", "")

        fields = line.split("\t")

        if n == 1: # header row -- write tfam info
            for i in range(3,len(fields), 3):
                k = fields[i].split(".")
                id = ".".join(k[:len(k)-1])
                out = [id, id, "0", "0", "-9", "-9"]
                outTFAM.write(" ".join(out) + "\n")

        else:
            snp = fields[0]
            chr = fields[1]
            pos = fields[2]

            out = [chr, snp, "0", pos]

            for i in range(3, len(fields), 3):
                gt = fields[i]
                if gt == "NC":
                    gt = "00"
                out.append(gt[0])
                out.append(gt[1])
        
            outTPED.write(" ".join(out) + "\n")

                    
    
