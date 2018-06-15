#!/usr/bin/env python2

"""
This script is to: 
1) count barcode and output into a tab separated table
2) convert input FASTQ file to FASTA file

USAGE: python2 ss-03.barcodeCounting.py <inputFastq> <outputTsv>

DEPENDENCIES: Biopython

INPUT: a FASTQ file, not zipped

OUTPUT:
1) printed into given output file name: tab-delimited text file with first column being barcode sequences, and second column being counts in given FASTQ file. No header. 
2) printed into standard output (screen): converted FASTA file of all instances of barcodes, with read names and sequences.
"""

import sys

fq_input_name = sys.argv[1]
output_name = sys.argv[2]

import Bio
from Bio.SeqIO.QualityIO import FastqGeneralIterator

count_dict = dict()

with open(fq_input_name) as in_handle:
	for title, seq, qual in FastqGeneralIterator(in_handle):
		if (not('N' in seq)):
			print (">" + title)
			print (seq)
			if (seq in count_dict):
				count_dict[seq] += 1
			else:
				count_dict[seq] = 1

# redirect printing to file
orig_stdout = sys.stdout
f = open(output_name, "w")
sys.stdout = f

print("Sequence" + "\t" + "Count")
for seq, count in count_dict.items():
    print('{}\t{}'.format(seq, count))
