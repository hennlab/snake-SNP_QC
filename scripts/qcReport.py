#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# January 10, 2014

import sys
import gzip
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-R","--report",type="string",dest="report",action="store",help="genome studio report to call from")
parser.add_option("-C","--callrate",type="string",dest="callrate",action="store",help="sample call rate filter")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

if options.callrate == None:
    cr = 0.99
else:
    cr = float(options.callrate)

if options.report.find(".gz") != -1:
    file = gzip.open(options.report)
else:
    file = open(options.report,'r')

CR = {}
nsnps = 0
n = 0
for line in file:
    n += 1
    line = line.replace("\n", "")
    line = line.replace("\r", "")
        
    fields = line.split("\t")
    if n == 1:
        for i in range(3,len(fields),3):
            CR[i] = 0
    else:
        nsnps += 1
        for i in range(3,len(fields),3):
            gt = fields[i]
            if gt != "NC":
                CR[i] += 1


KEEP = {}
for i in CR:
    if CR[i]/float(nsnps) >= cr:
        KEEP[i] = ""
        KEEP[i + 1] = ""
        KEEP[i + 2] = ""


if options.report.find(".gz") != -1:
    file = gzip.open(options.report)
else:
    file = open(options.report,'r')

for line in file:    
    line = line.replace("\n", "")
    line = line.replace("\r", "")
        
    fields = line.split("\t")

    out = []
    for i in range(len(fields)):
        if i < 3 or i in KEEP:
            out.append(fields[i])
    print "\t".join(out)            
    

