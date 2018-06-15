#!/use/bin/env python2

"""
This script is to merge barcode, UMI, and insertion location tables into one by sequence identifier of each paired-end read pair containing barcode and UMI, i.e. the `QNAME` in SAM file. The merged table is  written into csv file with all colums as: 
`barcode sequence,UMI sequence, RNAME, POS, STRAND`
with no header. RNAME, POS, STRAND are defied the same as in SAM format conventions. 
Stats for merging process are also printed to standard output (screen). 
USAGE:
ss-05-mergeTables.py <barcodeFilename> <umiFilename> <locFilename> <outputFilename>

EXMAPLE:
python2 ss-05-mergeTables.py Cvi.barcode.csv Cvi.umi.csv Cvi.genomic.csv Cvi.barcodeUmiGenomic.csv

INPUTS:
1) barcodeFilename: a csv file with two columns. First column: sequence identifiers; second column: barcode sequences.
2) umiFilename: a csv file with two columns. First column: sequence identifiers; second column: UMI sequences. 
3) locFilename: a csv file with four columns. First column: sequence identifiers; second column: reference names of aligned genomic sequences, e.g. genome or plasmids; third: position of alignment; forth: strand of alignment. 
4) outputFilename: a string as file name of output csv file.

OUTPUTS:
1) A csv file with information of `barcode sequence,UMI sequence, RNAME, POS, STRAND`, with no header. 
2) Standard output (screen print) of stats for merging process: 
"""

import pandas as pd
import sys

qnameBarcodeFilename = sys.argv[1]
qnameUmiFilename = sys.argv[2]
qnameLocFilename = sys.argv[3]
outputFilename = sys.argv[4]

# read in files
print('read in qname -> barcode')
qnameBarcode = pd.read_csv(qnameBarcodeFilename, sep = ',', header = None)
qnameBarcode.columns = ['qname', 'barcode']
print('barcode file lines (non-unique barcode counts): %d' % qnameBarcode.shape[0])
sys.stdout.flush() 

print('read in qname -> UMI')
qnameUmi = pd.read_csv(qnameUmiFilename, sep = ',', header = None)
qnameUmi.columns = ['qname', 'umi']
print('umi file lines (non-unique umi counts): %d' % qnameUmi.shape[0])
sys.stdout.flush() 

print('read in qname -> insertion location')
qnameLoc = pd.read_csv(qnameLocFilename, sep = ',', header = None)
qnameLoc.columns = ['qname', 'rname', 'pos', 'strand']
print('genomic file lines (non-unique genomic location counts): %d' % qnameLoc.shape[0])
sys.stdout.flush() 

# merge barcode and UMI
print('merge barcode and UMI')
qnameBarcodeUmi = pd.merge(qnameBarcode, qnameUmi, on = 'qname', how = 'inner')
print('Total barcode-umi applicable read counts: %d' % (qnameBarcodeUmi.shape[0]))
print(' -> reads with barcode but not UMI: %d' % (qnameBarcode.shape[0] - qnameBarcodeUmi.shape[0]))
print(' -> reads with UMI but not barcode: %d' % (qnameUmi.shape[0] - qnameBarcodeUmi.shape[0]))
sys.stdout.flush() 

# merge barcode-UMI and genomic locations
print('merge barcode-UMI and genomic locations')
qnameBarcodeUmiLoc = pd.merge(qnameBarcodeUmi, qnameLoc, on = 'qname', how = 'inner')
print('Total barcode-umi-loc applicable read counts: %d' % (qnameBarcodeUmiLoc.shape[0]))
print(' -> reads with locations but not barcode-UMI: %d' % (qnameLoc.shape[0] - qnameBarcodeUmiLoc.shape[0]))
print(' -> reads with barcode-UMI but not locations: %d' % (qnameBarcodeUmi.shape[0] -qnameBarcodeUmiLoc.shape[0]))
sys.stdout.flush() 

# remove qname column
barcodeUmiLocLite = qnameBarcodeUmiLoc.drop('qname', axis = 1)

# write to file
print('write to file')
sys.stdout.flush() 
barcodeUmiLocLite.to_csv(outputFilename, index = False, header = False)
