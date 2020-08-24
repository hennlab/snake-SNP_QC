#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
from optparse import OptionParser
from calcMeanSD import *

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--betas",type="string",dest="betas",action="store",help="betas.txt file path")
parser.add_option("-R","--report",type="string",dest="report",action="store",help="GenomeStudio report file path")
parser.add_option("-I","--minint",type="string",dest="minIntensity",action="store",help="minimum mean int signal for comm hom cluster")
parser.add_option("-Z","--z",type="string",dest="z",action="store",help="z-score threshold")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

if options.z == None:
    options.z = 7

if options.betas == None:
    print "specify betas.txt file with -B"
    sys.exit()

z = int(options.z)

if options.minIntensity != None:
    options.minIntensity = float(options.minIntensity)
else:
    options.minIntensity = 0.2

### Print header line to std out
head = ["SNP", "Tx", "Ty"] 
print "\t".join(head)


### Parse betas.txt file
### Order is as follows:
### meanY ~ meanX
### meanX ~ meanY
### sdY ~ sdX
### sdX ~ sdY
beta0 = [] # list container for beta intercept
beta1 = [] # list container for beta of slope
for line in open(options.betas, 'r'):
    line = line.replace("\n", "")
    if line.find("Beta0") == -1:
        fields = line.split("\t")
        beta0.append(float(fields[1]))
        beta1.append(float(fields[2]))


### Find thresholds for each site
### AA is always quadrant 4 by definition, BB is always quadrant 1 by definition
### because X always tags A and Y always tags B

if options.report.find(".gz") != -1:
    file = gzip.open(options.report)
else:
    file = open(options.report,'r')

n = 0
for line in file:
    n += 1
    line = line.replace("\n", "")
    line = line.replace("\r", "")

    if n == 1: # skip header line
        continue
    else:
        fields = line.split("\t")

        # get snp name
        snp = fields[0]

        # Extract the mean and sd for each common allele homozygote clusters in the noise dimension

        X_AA = []
        Y_AA = []
        X_BB = []
        Y_BB = []
        nAA = 0
        nBB = 0
        genotypes = []

        for i in range(3, len(fields), 3):
            gt = fields[i]
            x = float(fields[i + 1])
            y = float(fields[i + 2])

            genotypes.append(gt)
            
            if gt == "AA": # make arrays of X and Y intensities for AA genotype
                X_AA.append(x)
                Y_AA.append(y)
                nAA += 1
                
            if gt == "BB": # make arrays of X and Y intensities for BB genotype
                X_BB.append(x)
                Y_BB.append(y)
                nBB += 1

        if nAA <= 2 and nBB <= 2: # Too few points in common allele homozygote cluster
            Tx = "NA"
            Ty = "NA"

        else:
            try:
                meanXAA, devXAA = calcMeanSD(X_AA)
                meanYAA, devYAA = calcMeanSD(Y_AA)
            except: # only 0 or 1 point in cluster -- mean and SD don't apply, so only need mean and SD of other cluster
                meanXBB, devXBB = calcMeanSD(X_BB)
                meanYBB, devYBB = calcMeanSD(Y_BB)

            try:
                meanXBB, devXBB = calcMeanSD(X_BB)
                meanYBB, devYBB = calcMeanSD(Y_BB)
            except: # only 0 or 1 point in cluster -- mean and SD don't apply, so only need mean and SD of other cluster
                meanXAA, devXAA = calcMeanSD(X_AA)
                meanYAA, devYAA = calcMeanSD(Y_AA)

            if nAA <= 2 and nBB <= 2: # Not enough points in common allele homozygote cluster
                Tx = "NA"
                Ty = "NA"

            else:
                # Calculate Thresholds depending on which homozygote cluster is tagging the common allele for that SNP
                if nAA >= nBB:
                    if meanXAA < options.minIntensity: # site has less than min. intensity to recall
                        Tx = "NA" # mark with "NA" so skip new genotype calls in zCall.py
                        Ty = "NA"
                    else:
                        Ty = meanYAA + z * devYAA
                        meanXBB = beta1[1]*meanYAA + beta0[1] # Solve for the mean of the minor allele hom. cluster based on betas and mean of common allele hom. cluster
                        devXBB = beta1[3]*devYAA + beta0[3] # Solve for the sd of the minor allele hom. cluster based on betas and sd of common allele hom. cluster
                        Tx = meanXBB + z * devXBB # Use inferred mean and sd to find Tx

                if nAA < nBB:
                    if meanYBB < options.minIntensity: # site has less than min. intensity to recall
                        Tx = "NA" # mark with "NA" so skip new genotype calls in zCall.py
                        Ty = "NA"
                    else:
                        Tx = meanXBB + z * devXBB
                        meanYAA = beta1[0] * meanXBB + beta0[0] # Solve for the mean of the minor allele hom. cluster based on betas and mean of common allele hom. cluster
                        devYAA = beta1[2] * devXBB + beta0[2] # Solve for the sd of the minor allele hom. cluster based on betas and sd of common allele hom. cluster
                        Ty = meanYAA + z * devYAA # Use inferred mean and sd to find Ty


        # Write thresholds to std out
        out = [snp, Tx, Ty]
        out = [str(o) for o in out]

        print "\t".join(out)


