#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
import gzip
from optparse import OptionParser
from calcMeanSD import *

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-R","--report",type="string",dest="report",action="store",help="GenomeStudio report file path")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

### Write header line
head = ["SNP", "meanX", "meanY", "sdX", "sdY", "nMinorHom", "nCommonHom"] 
print "\t".join(head) # Write header line

if options.report.find(".gz") != -1:
    file = gzip.open(options.report)
else:
    file = open(options.report,'r')

## Iterate over each SNP in genome studio report
n = 0
for line in file:
    n += 1
    line = line.replace("\n", "")
    line = line.replace("\r", "")

    if n == 1:
        continue
    else:
        fields = line.split("\t")

        # get snp name
        snp = fields[0]

        # extract genotypes into python list
        genotypes = [fields[i] for i in range(3, len(fields), 3)]

        # Get number of points in each genotype cluster
        nAA = genotypes.count("AA")
        nAB = genotypes.count("AB")
        nBB = genotypes.count("BB")
        nNC = genotypes.count("NC")
        nGenotypes = nAA + nAB + nBB
        nTotal = len(genotypes)
        

        # Calculate Missing Rate (ignore SNPs that have less than 99% call rate)
        if float(nGenotypes) / float(nTotal) < 0.99:
            continue

        # Make sure there are at least 10 points in each homozygote cluster
        if nAA < 10 or nBB < 10:
            continue

        # Calculate MAF
        if nAA > nBB:
            maf = (nAB + 2 * nBB) / float(2 * nTotal)
        elif nAA <= nBB:
            maf = (nAB + 2 * nAA) / float(2 * nTotal)

        # MAF check ( >5% MAF)
        if maf < 0.05:
            continue

        # Hardy-Weinberg Equilibrium Check (don't use site if p_hwe < 0.00001)
        chiCritical = 19.5 # p = 0.00001 for 1 DOF

        if nAA > nBB:
            p = 1.0 - maf
            q = maf        
            expAA = p**2 * nTotal
            expAB = 2 * p * q * nTotal
            expBB = q**2 * nTotal

        if nBB >= nAA:
            p = 1.0 - maf
            q = maf        
            expAA = q**2 * nTotal
            expAB = 2 * p * q * nTotal
            expBB = p**2 * nTotal

        chiSquare = ((nAA - expAA)**2 / float(expAA)) + ((nAB - expAB)**2 / float(expAB)) + ((nBB - expBB)**2 / float(expBB))
        if chiSquare > chiCritical:
            continue

        # Extract the mean and sd for each common allele homozygote clusters in the noise dimension

        X_AA = []
        Y_AA = []
        X_BB = []
        Y_BB = []

        for i in range(3, len(fields), 3):
            gt = fields[i]
            x = float(fields[i + 1])
            y = float(fields[i + 2])

            if gt == "AA": # make arrays of X and Y intensities for AA genotype
                X_AA.append(x)
                Y_AA.append(y)
            if gt == "BB": # make arrays of X and Y intensities for BB genotype
                X_BB.append(x)
                Y_BB.append(y)

        
        meanXAA, devXAA = calcMeanSD(X_AA)
        meanYAA, devYAA = calcMeanSD(Y_AA)
        
        meanXBB, devXBB = calcMeanSD(X_BB)
        meanYBB, devYBB = calcMeanSD(Y_BB)
        

        if meanXAA >= meanYAA: ## AA is in the lower right quadrant
            meanY = meanYAA
            devY = devYAA
            meanX = meanXBB
            devX = devXBB

        elif meanXAA < meanYAA: ## AA is in the upper left quadrant; However this should never be used because by definition AA is always in the lower right quadrant
            meanY = meanYBB
            devY = devYBB
            meanX = meanXAA
            devX = devXAA

        if nAA >= nBB:
            out = [snp, meanX, meanY, devX, devY, nBB, nAA] # output array
        elif nBB > nAA:
            out = [snp, meanX, meanY, devX, devY, nAA, nBB] # output array
        out = [str(o) for o in out]

        print "\t".join(out)

