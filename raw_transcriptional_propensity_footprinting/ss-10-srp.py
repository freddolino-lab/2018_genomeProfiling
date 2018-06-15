#!/usr/bin/env python2
"""
This script is to label gene-level transcriptional propensity by the genes' SRP status (whether or not the coding products were recognized by SRP).  
USAGE:
python2 ss-10-srp.py

INPUTs:
Input files are coded as inputs to main function.
srpFilename: a comma-separated list of SRP gene names and corresponding b-numbers
ipageValueFilename: a gr file (tab-separated two-column text file) for gene-level transcriptional propensity, with first column as b-numbers of genes, and second column as corresponding transcriptional propensity
gffFilename: a NCBI style GFF3 file annotating genes on E. coli genome

OUTPUT:
Screen output (standard output): a text file of comma-separated values with each column being: gene b-number, gene-level transcriptional propensity, gene start position, gene end position, gene strandness, gene SRP status
Standard error: statistics when searching SRP gene names and corresponding b-numbers across the value list
"""

from __future__ import print_function
from __future__ import division
import sys

def parseAttribute(attributeString):
    """
    get the b-number from one Attribute field of one line from GFF file in NCBI format, keeping only Gene or Synonym
    First written in `scott-5-koEffectCircular.py`.
    """
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

def parseSrpFile(srpLineList):
    """
    parse input list of SRP gene names, expected input to be 
    """
    genenameList = []
    bnumList = []
    for srp in srpLineList:
        fields = srp.strip().split(',')
        genename = fields[0]
        bnum = fields[1]
        genenameList.extend([genename])
        if (bnum is not ''):
            bnumList.extend([bnum])
    return((genenameList, bnumList))

def eprint(*args, **kwargs):
    """
    print to stderr
    credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)

def main(srpFilename, ipageValueFilename, gffFilename):
    srpFile = open(srpFilename)
    srpLineList = srpFile.readlines()
    srpFile.close()
    
    srpParsedTuple = parseSrpFile(srpLineList)
    srpGenenameList = srpParsedTuple[0]
    srpBnumList = srpParsedTuple[1]

    ipageValueFile = open(ipageValueFilename)
    ipageValueLineList = ipageValueFile.readlines()
    ipageValueFile.close()
    
    gffCoordDict = dict()
    with open(gffFilename) as gffFile:
        for line in gffFile:
            gffFields = parseGffInputLine(line)
            bnum = parseAttribute(gffFields[3])
            gffCoordDict[bnum] = gffFields[0:3]
    
    
    outputLineList = []
    srpFlag = dict.fromkeys(srpBnumList, 0)
    
    for ipageValueLine in ipageValueLineList:
        ipageFields = ipageValueLine.strip().split('\t')
        ipageBnum = ipageFields[0]
        ipageValue = float(ipageFields[1])
        posList = gffCoordDict[ipageBnum]
        posString = ','.join(map(str, posList))
        if (ipageBnum in srpBnumList): #any(i in a for i in b)
            tag = 'SRP'
            srpFlag[ipageBnum] = srpBnumList.count(ipageBnum)
        else:
            tag = 'nonSRP'
        outputLineList.extend([','.join([ipageBnum, str(ipageValue), posString, tag])])
    
    srpFound = srpFlag.values().count(1)
    srpNotFound = srpFlag.values().count(0)
    srpNotFoundNames = {k: v for k, v in srpFlag.iteritems() if v == 0}.keys()
    srpRepeatedNames = {k: v for k, v in srpFlag.iteritems() if v > 1}.keys()
    
    eprint('[Total genes with values]')
    eprint(str(len(ipageValueLineList)))
    eprint('[Total requested SRP genes with valid b-number]')
    eprint(str(len(srpBnumList)))
    eprint('[More than 1 SRP gene mapped to same b-number]')
    eprint(','.join(srpRepeatedNames))
    eprint('[Count of genes coding SRP found in value table]')
    eprint(str(srpFound))
    eprint('[Count of genes coding SRP NOT found in value table]')
    eprint(str(srpNotFound))
    eprint('[Names of genes coding SRP NOT found in value table]')
    eprint(','.join(srpNotFoundNames))
    
    for outputLine in outputLineList:
        print(outputLine)
    

if __name__ == '__main__':
    srpFilename = 'ss-list_SRP.genename.mapped.csv'
    ipageValueFilename = 'ss-ipageInputLogMedianFlank.exp'
    gffFilename = 'Escherichia_coli/MG1655/NCBI/U00096.3.gene.gff'
    
    main(srpFilename, ipageValueFilename, gffFilename)
