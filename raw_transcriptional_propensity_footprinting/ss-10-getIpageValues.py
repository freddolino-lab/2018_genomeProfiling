#!/usr/bin/env python2
"""
For each gene annotation, get log2 median value of RNA/DNA ratio within the gene and the flanking region (input length at EACH SIDE)
positioning is 1-based through out this script

USAGE:
python ss-10-getIpageValues.py <RNADNAratio.gr> <annotationNCBI.gff> <flankingLengthOneSide>

EXAMPLE:
python2 ss-10-getIpageValues.py normedRatios.gr U00096.3.gene.gff 2500 >10-iPage/ss-ipageInputLogMedianFlank.exp

INPUT
RNADNAratio.gr: a string of file name of the GR file of RNA/DNA ratios at different coordinates; GR files are tab separated two columns with first being coordinate values and second being corresponding ratios
annotationNCBI.gff: a string of file name of the annoation GFF3 file for genes on E. coli genome
flankingLengthOneSide: a integer of length of flanking region on each side of the annotated gene region
OUTPUT
"""

import sys
from numpy import median
from numpy import log2

GENOMESIZE = 4641652

def parseAttribute(attributeString):
    """
    get the b-number from one Attribute field of one line from GFF file in NCBI format, keeping only Gene or Synonym
    First written in `scott-5-koEffectCircular.py`.
    """
    # cleanedString = attributeString.replace()
    fields = attributeString.split(';')
    newStringList = []
    for field in fields:
        if (field.startswith('locus_tag=')):
            newStringList += [field[10:]]
    newString = ';'.join(newStringList)
    return(newString)

def parseGffInputLine(gffLine):
    """
    parse each line of GFF file from NCBI, return a list of interesting attributes
    First written in `scott-5-koEffectCircular.py`.
    """
    fields = gffLine.strip().split('\t')
    start = int(fields[3])
    end = int(fields[4])
    strand = str(fields[6])
    annotation = str(fields[8])
    return([start, end, strand, annotation])


def parseRatioListToDigitList(grList):
    """
    expand GR file containing position and ratio information, into a list with indices as positions on genome and elements as ratios
    """
    maxPos = GENOMESIZE
    grDigitList = [None] * (maxPos + 1)
    for line in grList:
        fields = line.strip().split('\t')
        pos = int(fields[0])
        value = float(fields[1])
        grDigitList[pos] = value
    return(grDigitList)

def getValueByInterval(start, end, flankingLen, ratio):
    """
    retrieve all ratio values given a start position, and end position, and a flanking length. Tolerate circular genome. 
    """
    newStart = start - flankingLen
    newEnd = end + flankingLen
    if (newStart < 0):
        leftRegion = ratio[(GENOMESIZE + newStart):(GENOMESIZE + 1)]
        newStart = 0
    else:
        leftRegion = []
    if (newEnd > GENOMESIZE):
        rightRegion = ratio[0:(newEnd - GENOMESIZE)]
        newEnd = GENOMESIZE + 1
    else: 
        rightRegion = []
    
    region = ratio[newStart:newEnd]
    values = leftRegion + region + rightRegion
    return(values)

def main(ratioFilename, gffFilename, flankingLen):
    # get gene information from GFF file
    gffFile = open(gffFilename)
    gff = gffFile.readlines()
    gffFile.close()
    gffList = []
    # list of list for GFF info, each element is a list as [strand, end, b-number]
    for i in range(0, len(gff)):
        fields = parseGffInputLine(gff[i])
        start = fields[0] # inclusive start
        end = fields[1] + 1 # exlusive end 
        strand = fields[2]
        annotation = parseAttribute(fields[-1])
        gffList.append([start, end, annotation])
    
    
    # get RNA/DNA ratio from GR file
    ratioFile = open(ratioFilename)
    ratioList = ratioFile.readlines()
    ratioFile.close()
    ratioPositionList = parseRatioListToDigitList(ratioList)
    # get corresponding expression value for each gene
    
    for i in range(0, len(gffList)):
        start = gffList[i][0]
        end = gffList[i][1]
        annotation = gffList[i][2]
        # get values in the flanked interval
        values = getValueByInterval(start = start, end = end, flankingLen = flankingLen, ratio = ratioPositionList)
        # remove None-s
        validValues = [x for x in values if x is not None]
        if (len(validValues) > 0):
             # get median of values, using function in Numpy
            metric = median(validValues)
            logMetric = log2(metric)
            print(annotation + '\t' + str(logMetric))
        

if __name__ == '__main__':
    ratioFilename = sys.argv[1]
    gffFilename = sys.argv[2]
    flankingLen = int(sys.argv[3])
    
    main(ratioFilename, gffFilename, flankingLen)
