#!/usr/bin/env python2
"""
This script is to map gene names to b-number, using given GFF file

INPUT: NCBI-style GFF file (e.g. U00096.3.gene.gff), list of lines of gene names/synonyms

GFF sample line:
U00096.3	Genbank	gene	4386047	4389370	.	-	.	ID=gene4249;Dbxref=EcoGene:EG12478;Name=mscM;gbkey=Gene;gene=mscM;gene_biotype=protein_coding;gene_synonym=ECK4155,JW4120,yjeP;locus_tag=b4159

OUTPUT: b-number corresponds to each line of gene names/synonyms

USAGE:
python2 mapGenenameToBnum.py <GFF_filename> <genename_list_filename> 

EXAMPLE:
python2 ~/projects/genomeProfiling/ssPytools/mapGenenameToBnum.py /srv/diaorch/Escherichia_coli/MG1655/NCBI/U00096.3.gene.gff ss-list_SRP.genename.csv >ss-list_SRP.genename.mapped.csv 2>ss-list_SRP.genename.mapped.err
"""

from __future__ import print_function
import sys

def parseGffLine(gffLine, mode):
    """
    parse each line of GFF file, to get a tuple of (b_num, [names, synonyms])
    mode = [primaryName, synonym]
    """
    fields = gffLine.strip().split('\t')
    annotation = str(fields[8])
    attributeFields = annotation.split(';')
    bnum = None
    for attribute in attributeFields:
        if (attribute.startswith('locus_tag=')):
            bnum = attribute[10:]
        if (mode == 'primaryName' and attribute.startswith('Name=')):
            newName = attribute[5:]
        if (mode == 'synonym' and attribute.startswith('gene_synonym=')):
            newNameList=attribute[13:].split(',')
    # newName = ';'.join(newNameList)
    if (mode == 'primaryName'):
        return(bnum, newName)
    elif(mode == 'synonym'):
        return(bnum, newNameList)

def parsePrimaryNames(gffLineList):
    """
    parse GFF list to create dictionary for matching with primary gene names as keys and b-numbers as values
    """
    nameBnumDict = dict()
    for line in gffLineList:
        bnum, name = parseGffLine(gffLine = line, mode = 'primaryName')
        nameBnumDict[name] = bnum
    return(nameBnumDict)

def parseSecondarySynonyms(gffLinesList):
    """
    parse GFF list to create dictionary for all possible synonyms as keys and corresponding b-numbers as values
    """
    nameBnumDict = dict()
    for line in gffLinesList:
        bnum, synonymList = parseGffLine(gffLine = line, mode = 'synonym')
        geneDict = dict.fromkeys(synonymList, bnum)
        for k, v in geneDict.iteritems():
            if (k in nameBnumDict.keys()):
                nameBnumDict[k].extend([v])
            else:
                nameBnumDict[k] = [v]
    return(nameBnumDict)

def eprint(*args, **kwargs):
    """
    print to stderr
    credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)

def main(gffFilename, genenameFilename):
    
    gffFile = open(gffFilename)
    gffLineList = gffFile.readlines()
    gffFile.close()
    
    genenameFile = open(genenameFilename)
    genenameLineList = genenameFile.readlines()
    genenameFile.close()
    genenameList = map(str.strip, genenameLineList)
    
    convertedDict = dict.fromkeys(genenameList, None)
    
    primaryNameDict = parsePrimaryNames(gffLineList)
    secondaryNameDict = parseSecondarySynonyms(gffLineList)
    
    geneBnumList = []
    geneMappedType = []
    for genename in genenameList:
        # print(genename)
        if (genename in primaryNameDict.keys()):
            # print(primaryNameDict[genename])
            geneBnumList.extend([primaryNameDict[genename]])
            geneMappedType.extend(['primary'])
        elif (genename in secondaryNameDict.keys()):
            # print(secondaryNameDict[genename])
            if (len(secondaryNameDict[genename]) == 1):
                geneBnumList.extend(secondaryNameDict[genename])
                geneMappedType.extend(['synonym'])
            else:
                eprint('[WARNING] No primary gene name but multiple GFF synonyms record found for requested gene: ' + genename + ' - marked as conflict type and kept the smallest record')
                foundBnumList = secondaryNameDict[genename]
                bnumInts = map(lambda x: int(x[1:]), foundBnumList)
                selectedBnum = min(bnumInts)
                geneBnumList.extend(['b' + str(selectedBnum)])
                geneMappedType.extend('conflict')
        else:
            eprint('[WARNING] No primary gene name or synonym found for requested gene: ' + genename + ' - marked as notFound type and use empty string as b-number')
            geneBnumList.extend([''])
            geneMappedType.extend(['notFound'])
            
    
    for i in range(0, len(genenameList)):
        print(genenameList[i] + ',' + geneBnumList[i] + ',' + geneMappedType[i])
    

if __name__ == '__main__':
    gffFilename = sys.argv[1]
    genenameFilename = sys.argv[2]
    main(gffFilename, genenameFilename)
