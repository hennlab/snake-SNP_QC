#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-R","--report",type="string",dest="report",action="store",help="genome studio report to call from")
parser.add_option("-T","--thresholds",type="string",dest="thresholds",action="store",help="thresholds file")
parser.add_option("-O","--output",type="string",dest="outputROOT",action="store",help="output directory")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

if options.outputROOT == None:
    print "specify output root with -O"
    sys.exit()

if options.thresholds == None:
    print "specify thresholds file with -T"
    sys.exit()


# Initialize output files
outTFAM = open(options.outputROOT + ".tfam", 'w')
outTPED = open(options.outputROOT + ".tped", 'w')

# Get thresholds for each site
thresholdsX = {}
thresholdsY = {}
n = 0
for line in open(options.thresholds, 'r'):
    n += 1
    line = line.replace("\n", "")
    if n == 1:
        continue
    else:
        fields = line.split("\t")
        snp = fields[0]

        if fields[1] == "NA" or fields[2] == "NA":
            Tx = fields[1]
            Ty = fields[2]
        else:
            Tx = float(fields[1])
            Ty = float(fields[2])
            
        thresholdsX[snp] = Tx # dicitonary with tx for each snp
        thresholdsY[snp] = Ty # dictionary with ty for each snp
                

# Iterate through Illumina GenomeStudio Report
if options.report.find(".gz") != -1:
    file = gzip.open(options.report)
else:
    file = open(options.report,'r')

n = 0
for line in file:
    n += 1
    line = line.replace("\n", "")
    line = line.replace("\r", "")

    fields = line.split("\t")

    if n == 1: # header row -- write tfam info
        for i in range(3,len(fields), 3):
            k = fields[i].split(".")
            id = ".".join(k[:len(k)-1])
            out = [id, id, "0", "0", "-9", "-9"] # don't know PID, MID, sex, or case/control status
            outTFAM.write(" ".join(out) + "\n")

    else:
        snp = fields[0]
        chr = fields[1]
        pos = fields[2]
        Tx = thresholdsX[snp] # get Tx and Ty from dictionary for that snp
        Ty = thresholdsY[snp]

        out = [chr, snp, "0", pos] # start output array with map info

        for i in range(3, len(fields), 3):
            # extract one sample's gt, normalize X and Y intensities
            gt = fields[i]
            x = fields[i + 1]
            y = fields[i + 2]

            if x != "NaN":
                x = float(x)
            if y != "NaN":
                y = float(y)

            if gt != "NC" or Tx == "NA" or Ty == "NA": # If not no call or no thresholds defined, output original call
                if gt == "NC":
                    out.append("0")
                    out.append("0")
                else:
                    out.append(gt[0])
                    out.append(gt[1])
        
            elif gt == "NC":
                ### AA is always quadrant 4 by definition, BB is always quadrant 1 by definition because X always tags A and Y always tags B

                if x == "NaN" or y == "NaN": # No intensity data -- assign as no call
                    out.append("0")
                    out.append("0")
                elif x >= Tx and y <= Ty: # lower right quadrant
                    out.append("A")
                    out.append("A")
                elif x <= Tx and y >= Ty: # upper left quadrant
                    out.append("B")
                    out.append("B")
                elif x > Tx and y > Ty: # upper right quadrant
                    out.append("A")
                    out.append("B")
                elif x < Tx and y < Ty: # lower left quadrant
                    out.append("0")
                    out.append("0")
                    
        outTPED.write(" ".join(out) + "\n") # write output to file

        
                    
