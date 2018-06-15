#!/usr/bin/env python2
"""
This script calculate the Knockout effect of genome footprinting insertion
USAGE:
python2 ss-09-koEffectCircular.NCBI.py <ratioGR> <inputGFF> 
EXAMPLE: 
python2 ss-09-koEffectCircular.NCBI.py tmarchNormAverageIntCoord.gr Ecoli.NCBI.U00096.3.gene.20180112.gff > koEffect.filterSparseWindow.NCBI.csv
INPUT:
ratioGR: a .gr file with first column being coordinates and second column being corresponding transcriptional propensity measure, in this case RNA/DNA ratios quantile normalized and averaged between replicates, separated by a tab, with no header.  
inputGFF: a GFF3 file with gene annotation. Notice it should be compatible with genome version. 
OUTPUT:
Standard output (screen output): a csv file with the following fields: `start of gene annotation, end of gene annotation, strandness of gene, attribute of gene annotation, median ratio in upstream 500bp window, median ratio in downstream 500bp window, median ratio in first 500bp window, median ratio in last 500bp window`. File has header, and is comma-separated.  
"""

import sys
import numpy as np

GENOMESIZE = 4641652 # 1-based

def parseAttribute(attributeString):
    """
    parse the Attributes filed of GFF file (input as string), keeping only Gene or Synonym
    """
    # cleanedString = attributeString.replace()
    fields = attributeString.split(';')
    newStringList = []
    for field in fields:
        if (field.startswith('gene=') or field.startswith('locus_tag=')):
            newStringList += [field]
    newString = ';'.join(newStringList)
    return(newString)

def convertToDigitGrList(lineList):
    """
    convert GR file into a vector of length of GENOMESIZE storing coverage (GR value) for each position
    """
    digitGr = [np.nan] * (GENOMESIZE + 1) # 1-based
    for line in lineList:
        fields = line.strip().split('\t')
        digitGr[int(fields[0])] = [float(fields[1])]
    return (digitGr)

def convertToGrList(lineList):
    """
    conver GR file input (list of lines) into a list of list containing the position and coverage as inner-list
    """
    grDoubleList = [[-1, -1] for x in xrange(GENOMESIZE)] 
    # print(grDoubleList)
    for i,line in enumerate(lineList):
        fields = line.strip().split('\t')
        grDoubleList[i][0] = int(fields[0])
        grDoubleList[i][1] = float(fields[1])
    return (grDoubleList)

def parseGffInputLine(gffLine):
    """
    parse each line of GFF file from NCBI, return a list of interesting attributes
    """
    fields = gffLine.strip().split('\t')
    start = int(fields[3])
    end = int(fields[4])
    strand = str(fields[6])
    annotation = str(fields[8])
    return([start, end, strand, annotation])

def getWindowValues(start, end, gr):
    """
    retrieve all values within a window of start and end (half open) from list of list converted from GR file, return all values within the window as a list
    """
    # print('start: ' + str(start) + ' end: ' + str(end))
    i = 0
    windowValues = []
    while (i < len(gr)):
        if (gr[i][0] >= start and gr[i][0] < end):
            windowValues += [gr[i][1]]
        if(gr[i][0] >= end):
            return (windowValues)
        i += 1
    return(windowValues)
    

def getUpstreamValues(gr, start, end, strand, windowLength = 500):
    """
    get values upstream of feature (gene with start, end, and strand) of windowLength, be aware of the orientation (upstream of plus strand would be smaller coordinates and upstream of minus stand would be larger coordinates)
    """
    # input end is inclusive from gff file
    allWindowValues = []
    if (strand == '+'):
        upstreamUpBound = start - windowLength
    elif (strand == '-'):
        upstreamUpBound = end + windowLength
    else:
        sys.exit('unknown strandness in gr file')
    # catch if upstream is out of bound
    # inclusive starts and exclusive ends
    if (upstreamUpBound <= 0):
        allWindowValues = getWindowValues(start = 1, end = start, gr = gr) + getWindowValues(start = (GENOMESIZE + upstreamUpBound + 1), end = (GENOMESIZE + 1), gr = gr)
    elif (upstreamUpBound > GENOMESIZE):
        allWindowValues = getWindowValues(start = (end + 1), end = (GENOMESIZE + 1), gr = gr) + getWindowValues(start = 1, end = (upstreamUpBound - GENOMESIZE), gr = gr)
    else:
        if (strand == '+'):
            allWindowValues = getWindowValues(start = upstreamUpBound, end = start, gr = gr)
        elif (strand == '-'):
            allWindowValues = getWindowValues(start = (end + 1), end = (upstreamUpBound + 1), gr = gr)
    # print(allWindowValue)
    return(allWindowValues)

def getDownstreamValues(gr, start, end, strand, windowLength = 500):
    """
    get values downstream of feature (gene with start, end, and strand) of windowLength, be aware of the orientation (downstream of plus strand would be larger coordinates and downstream of minus stand would be smaller coordinates
    """
    # input end is inclusive from gff file
    allWindowValues = []
    if (strand == '+'):
        downstreamDownBound = end + windowLength # on-spot, inclusive
    elif (strand == '-'):
        downstreamDownBound = start - windowLength # on-spot, inclusive
    else:
        sys.exit('unknown strandness in gr file')
    # catch if downstream out of bound
    if (downstreamDownBound > GENOMESIZE):
        allWindowValues = getWindowValues(gr = gr, start = end + 1, end = GENOMESIZE + 1) + getWindowValues(gr = gr, start = 1, end = downstreamDownBound - GENOMESIZE + 1)
    elif (downstreamDownBound <= 0):
        allWindowValues = getWindowValues(gr = gr, start = GENOMESIZE + downstreamDownBound, end = GENOMESIZE + 1) + getWindowValues(gr = gr, start = 1, end = start + 1)
    else:
        if (strand == '+'):
            allWindowValues = getWindowValues(gr = gr, start = end + 1, end = downstreamDownBound + 1)
        else:
            allWindowValues = getWindowValues(gr = gr, start = downstreamDownBound, end = start)
    # print(allWindowValues)
    return(allWindowValues)

def getFirstWindowValues(gr, start, end, strand, windowLength = 500):
    """
    get values in the first windowLength of feature (gene with start, end, and strand), be aware of the orientation (first of plus strand would be smaller coordinates and first of minus stand would be larger coordinates
    """
    allWindowValues = []
    if (end - start + 1 >= windowLength): # input end is inclusive from gff file, if length between one and two time the windown size, take the overlapping window of the given length
        if (strand == '+'):
            allWindowValues = getWindowValues(gr = gr, start = start, end = start + windowLength) # end parameter exclusive
        elif (strand == '-'):
            allWindowValues = getWindowValues(gr = gr, start = end - windowLength + 1, end = end + 1)
        else:
            sys.exit('unknown strandness in gr file')
    else: 
        # halfPoint = int((end - start) / 2)
        # if (strand == '+'):
        #     allWindowValues = getWindowValues(gr = gr, start = start, end = start + halfPoint)
        # elif (strand == '-'):
        #     allWindowValues = getWindowValues(gr = gr, start = start + halfPoint + 1, end = end + 1)
        # else:
        #     sys.exit('unknown strandness in gr file')
        allWindowValues = []
    return(allWindowValues)    

def getLastWindowValues(gr, start, end, strand, windowLength = 500):
    """
    get values in the last windowLength of feature (gene with start, end, and strand), be aware of the orientation (last of plus strand would be larger coordinates and last of minus stand would be smaller coordinates
    """
    allWindowValues = []
    if (end - start + 1 >= windowLength): # input end is inclusive from gff file, if length between one and two time the windown size, take the overlapping window of the given length
        if (strand == '+'):
            allWindowValues = getWindowValues(gr = gr, start = end - windowLength + 1, end = end + 1) # end parameter exclusive 
        elif (strand == '-'):
            allWindowValues = getWindowValues(gr = gr, start = start, end = start + windowLength)
        else:
            sys.exit('unknown strandness in gr file')
    else:       
        # halfPoint = int((end - start) / 2)
        # if (strand == '+'):
        #     allWindowValues = getWindowValues(gr = gr, start = start + halfPoint + 1, end = end + 1)
        # elif (strand == '-'):
        #     allWindowValues = getWindowValues(gr = gr, start = start, end = start + halfPoint + 1)
        # else: 
        #     sys.exit('unknown strandness in gr file')
        allWindowValues = []
    return(allWindowValues)

def getWindowValuesMedian(allWindowValue, listLengthAtLeast = 10):
    """
    get the median value given a list of values, return nan for list with length smaller than listLengthAtList including empty list. i.e. return median of list when the list has at least listLengthAtLeast elements, otherwise return nan
    """
    if (len(allWindowValue) >= 10):
        windowMedian = np.nanmedian(allWindowValue)
    else:
        windowMedian = np.nan
    return (windowMedian)

if __name__ == '__main__':
    # read input
    ## GR file
    grInputFilename = sys.argv[1]
    grInputFile = open(grInputFilename)
    grInput = grInputFile.readlines()
    ## GFF file
    gffInputFilename = sys.argv[2]
    gffInputFile = open(gffInputFilename)
    gffInput = gffInputFile.readlines()
    
    # convert partial gr file to all locus full gr in a list
    grList = convertToGrList(grInput)
    # print(grList[100:200])
    # print(getWindowValues(start = 8000, end = 8100, gr = grList))
    # loop through annotation
    ## test
    gff = gffInput
    # print header
    print('start,end,strand,attribute,upstreamWindow,downstreamWindow,firstWindow,lastWindow')
    for line in gff:
        parsedGffLineFields = parseGffInputLine(line)
        # print(parsedGffLineFields)
        upstreamValues = getUpstreamValues(grList, parsedGffLineFields[0], parsedGffLineFields[1], parsedGffLineFields[2], windowLength = 500)
        # print(upstreamValues)
        upstreamMedian = getWindowValuesMedian(upstreamValues)
        downstreamValues = getDownstreamValues(grList, parsedGffLineFields[0], parsedGffLineFields[1], parsedGffLineFields[2], windowLength = 500)
        downstreamMedian = getWindowValuesMedian(downstreamValues)
        firstWindowValues = getFirstWindowValues(grList, parsedGffLineFields[0], parsedGffLineFields[1], parsedGffLineFields[2], windowLength = 500)
        # print(firstWindowValues)
        firstWindowMedian = getWindowValuesMedian(firstWindowValues)
        lastWindowValues = getLastWindowValues(grList, parsedGffLineFields[0], parsedGffLineFields[1], parsedGffLineFields[2], windowLength = 500)
        # print(lastWindowValues)
        lastWindowMedian = getWindowValuesMedian(lastWindowValues)
        # print the csv file line
        # fieldsStringList =  map(str, parsedGffLineFields)
        print(str(parsedGffLineFields[0]) + ',' + str(parsedGffLineFields[1]) + ',' + parsedGffLineFields[2] + ',\'' + parseAttribute(parsedGffLineFields[3]) + '\',' + str(upstreamMedian) + ',' + str(downstreamMedian) + ',' + str(firstWindowMedian) + ',' + str(lastWindowMedian))

