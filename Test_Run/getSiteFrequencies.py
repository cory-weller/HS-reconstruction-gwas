#!/usr/bin/env python

import gzip

import sys

vcf_gzfile = sys.argv[1]

with gzip.open(vcf_gzfile, 'r') as infile:
    for line in infile:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            L = len(line.rstrip().split())
        else:
            splitLine = line.rstrip().split()
            contig = splitLine[0]
            position = splitLine[1]
            refCount = splitLine[9:].count("0")
            altCount = splitLine[9:].count("1")
            print('\t'.join([str(x) for x in [contig, position, refCount, altCount]]))
