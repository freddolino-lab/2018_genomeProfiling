#!/usr/bin/env python2

import pandas as pd
import sys

"""
DO NOT USE THIS SCRIPT.

This script is to remove any barcode with any base of quality score below 30.

USAGE: python2 ss-054-filterBarcodeQuality.py <input_text_qname_to_quality.ssv> <input_text_qname_to_barcode.csv> <output_csv_filename.csv>

EXAMPLE: python2 ss-054-filterBarcodeQuality.py qnameQual.ssv qnameBarcode.csv barcode

ARGS: None

INPUTS:
1) input_text_qname_to_quality.ssv: plain text file with two columns separated by two spaces. First column: first line from barcode FASTQ file, i.e. the sequence identifier; second column: forth line from barcode FASTQ file, i.e. the barcode quality. Obtained from `awk` processed FASTQ file.
2)input_text_qname_to_barcode.csv: plain text file with two columns separated by a comma. First column: first line from barcode FASTQ file, i.e. the sequence identifier; second column: second line from barcode FASTQ file, i.e. the barcode sequences. Obtained from `awk` processed FASTQ file. 

OUTPUT:
output_csv_filename.csv: a text file of comma-separated-values with two columns: first column for sequence identifiers and second column for high quality barcode sequences. 
"""

def checkHighQual(qualString):
    highQualScores = '?@ABCDEFGHI'
    index = 0
    highQualFlag = True
    qualList = list(qualString)
    while(index < len(qualString) and highQualFlag):
        if (qualList[index] not in highQualScores):
            highQualFlag = False
        index += 1
    return(highQualFlag)

if __name__=='__main__':
    qnameQualFilename = sys.argv[1]
    qnameHighQual = []
    with open(qnameQualFilename) as qnameQualFile:
        for qnameQual in qnameQualFile:
            qnameQualInfo = qnameQual.rstrip().rsplit(' ')
            qname = qnameQualInfo[0]
            qual = qnameQualInfo[1]
            highQualFlag = checkHighQual(qual)
            if (highQualFlag):
                qnameHighQual.append(qname)
    
    # print(qnameHighQual)
    qnameHighQualDf = pd.DataFrame({'qname': qnameHighQual})
    # print(qnameHighQualDf)
    
    qnameBcFilename = sys.argv[2]
    qnameBc = pd.read_csv(qnameBcFilename, header = None)
    qnameBc.columns = ['qname', 'barcode']
    qnameBcHighQual = pd.merge(qnameBc, qnameHighQualDf, on = 'qname', how = 'right')
    qnameBcHighQual.to_csv(sys.argv[3], header = False, index = False, sep = ',')

