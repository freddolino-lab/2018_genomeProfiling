#!/usr/bin/env python2

"""
This script is to extract sequence namess and sequences from FASTQ file and convert into a .csv file, output to screen

USAGE:  python2 fastqSeq2csv.py [filename|stdin]

EXAMPLE1: zcat ~/data/ss/0_seq/CviAII_ec_tn5_S2_R1.barcode.fastq.gz | python2 ssPytools/fastqSeq2csv.py

EXAMPLE2: python2 ssPytools/fastqSeq2csv.py ~/data/ss/0_seq/CviAII_ec_tn5_S2_R1.barcode.fastq 

ARGS: 
filename: input FASTQ file name
stdin: left empty if input FASTQ file is input via standard input (screen input)

INPUT: an unzipped .fastq file, or screen input of unzipped .fastq file

OUTPUT: a .csv file with the first column being the sequences names, the second column being the sequences; no header; printed to screen
"""

import sys

if len(sys.argv) > 1:
    filenameGiven = True
    filename = sys.argv[1]
else:
    filenameGiven = False

if filenameGiven:
    inputFile = open(filename)
else:
    inputFile = sys.stdin

with inputFile:
	for counter, line in enumerate(inputFile):
		if (counter % 4 == 0): # 'barcode' qname
			title = line.rstrip()
			truncated = title[1:].split(' ')[0]
			sys.stdout.write(truncated)
		if (counter % 4 == 1): # 'barcode' sequence
			sys.stdout.write(',' + line.rstrip() + '\n')

