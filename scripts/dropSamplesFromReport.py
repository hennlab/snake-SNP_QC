#! /usr/bin/python

import sys
import gzip

report = sys.argv[1]
dropSamples = sys.argv[2]

samples = {}
for line in open(dropSamples, 'r'):
    line = line.replace("\n", "")
    samples[line] = ""


if report.find(".gz") != -1:
    file = gzip.open(report)
else:
    file = open(report,'r')

n = 0
dropColumns = {}
for line in file:
    n += 1
    line = line.replace("\n", "")
    fields = line.split("\t")
    out = [fields[0],fields[1],fields[2]]
    if n == 1:
        for i in range(3, len(fields)):
            k = fields[i].split(".")
            id = ".".join(k[:len(k)-1])
            if id in samples:
                dropColumns[i] = ""
            else:
                out.append(fields[i])

    else:
        for i in range(3, len(fields)):
            if i not in dropColumns:
                out.append(fields[i])

    print "\t".join(out)

