#!/usr/bin/env python2

import sys
import pandas as pd

"""
This script is to merge barcode count table with barcode insertion location information, and add GC content as an additional colum.n
USAGE:
python2 ss-07-mapCountsLoc.py <countTableFilename> <locationTableFilename> <outputTableFilename>
EXAMPLE: 
python2 ss-07-mapCountsLoc.py barcodeCountTable.csv Cvi.barcodeUmiGenomic.sorted.pyFiltered.csv countLocation.out.csv
INPUTS:
countTableFilename: a string of file name of a csv file of five columns, containing barcode sequences and corresponding counts in all samples (DNA_1, RNA_1, DNA_2,RNA_2). The file HAS a header, in the format of `barcode,sampleName1, sampleName2, sampleName3, sampleName4`. 
locationTableFilename: a string of file name of a csv file of five columns, containing: barcode sequences, UMI counts, insertion references, insertion positions, insertion strandness. NO HEADER.  
outputTableFilename: a string of file name for output file
OUTPUT:
A csv file containing 10 fields, containing: barcode sequences ,barcode counts in sample1,counts in sample2,counts in sample3, counts in sample4, UMI counts, insertion reference names,insertion positions, insertion strandness, GC contents of barcode sequences in form of percentage`. WITH HEADER. 
"""

def calcGCpercentage(barcode):
    """
    calculating GC content for single input barcode
    """
    barcodeList = list(barcode)
    totalLength = len(barcodeList)
    # count GC
    gcCount = 0
    for base in barcodeList:
        if (base == 'G' or base == 'C'):
            gcCount += 1
    gcPercentage = float(gcCount) / float(totalLength)
    return(gcPercentage)


if __name__ == '__main__': 
       
    # parse filenames
    countTableFilename = sys.argv[1]
    locTableFilename = sys.argv[2]
    outputTableFilename = sys.argv[3]
    # read files
    countTable = pd.read_csv(countTableFilename, sep = ',', header = 0)
    # print(countTable.head())
    locTable = pd.read_csv(locTableFilename, sep = ',', header = None)
    locTable.columns = ['barcode', 'umiCount', 'rname', 'pos', 'strand']
    # print(locTable.head())
    
    # merge barcode count table and loc table
    countLocTable = pd.merge(countTable, locTable, on = 'barcode', how = 'outer')
    countLocTable.sort_values(by = 'barcode', axis = 0, ascending = True, inplace = True)
    
    # calc GC content and add to data frame
    gc = countLocTable['barcode'].map(calcGCpercentage)
    countLocTable['gc'] = gc
    # print(countLocTable.head(50))
    
    # write to file
    countLocTable.to_csv(outputTableFilename, sep = ',', na_rep = '', header = True, index = False)
    
