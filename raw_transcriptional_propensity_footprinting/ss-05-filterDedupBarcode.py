#!/use/bin/env python2

"""
This script is to 1) filter barcodes with N in it, 2) filter barcodes with more than one insertion locations that are more than 2 nt apart, and 3) count UMIs for each surviving barcode. INPUT MUST BE SORTED BY LINES, ALPHABETICALLY. 

USAGE:
python2 ss-05-filterDedupBarcode.py <inputSortedFilename> <outputFilename>

EXAMPLE:
python2 ~/projects/genomeProfiling/ss-05-filterDedupBarcode.py Cvi.barcodeUmiGenomic.sorted.csv Cvi.barcodeUmiGenomic.filtered.csv

INPUT:
inputSortedFilename: a string of file name of a csv file containing the following six columns, separated by a single comma: barcode sequences, UMI sequences, insertion references, insertion positions, insertion strandness. No header for input. THE INPUT MUST BE SORTED BY LINES ALPHABETICALLY.
outputFilename: a string for file name of output csv file

OUTPUT:
1) A csv file containing the following four columns, in order: barcode sequences, unique UMI counts, insertion references, insertion locations, insertion strandness.  
2) Standard output (screen output) of stats of filtering process. 
"""

import sys
import pandas as pd

def getLocationInt(loc_string):
    loc_list = loc_string.split(',')
    pos = int(loc_list[1])
    return(pos)

def getMaxMinLocDictKey(loc_dict):
    locList = list(loc_dict.keys())
    locIntList = map(getLocationInt, locList)
    diff = max(locIntList) - min(locIntList)
    return(diff)

def getMostAbundantLoc(loc_dict):
    maxValue = -1
    maxKey = ''
    maxKeyInt = 4641653
    for key, value in loc_dict.iteritems():
        loc = getLocationInt(key)
        if (value >= maxValue and loc < maxKeyInt):
            maxValue = value
            maxKey = key
            maxKeyInt = loc
    return(maxKey)

if __name__ == '__main__': 
    inputFilename = sys.argv[1]
    outputFilename = sys.argv[2]
    
    outfile = open(outputFilename, 'w')
    neighbor_threshold = 2
    with open(inputFilename) as f: # read barcode-umi-loc table line-by-line
        # initialize previous info
        first_line = f.readline()
        previous_line = first_line
        previous_info = first_line.rstrip().rsplit(',')
        previous_barcode = previous_info[0]
        barcode_counter = 0
        multi_umi_counter = 0
        drop_multiple_loc = 0 # counter for barcodes with multiple insertion locations that are not neighbors
        drop_multiple_loc_ambiguous = 0
        neighbor_counter = 0 # counter for barcodes with multiple insertion nearby locations that are neighbor
        umi_set = set([previous_info[1]])
        loc_dict = {(str(previous_info[2]) + ',' + str(previous_info[3]) + ',' + str(previous_info[4])):1}
        # iterate through the rest of lines (second +)
        for line in f:
            current_info = line.rstrip().rsplit(',')
            current_barcode = current_info[0]
            current_umi = current_info[1]
            current_loc = str(current_info[2]) + ',' + str(current_info[3]) + ',' + str(current_info[4])
            if (line !=  previous_line):
                if (current_barcode != previous_barcode):
                    # new barcode in next line
                    barcode_counter += 1
                    if (len(loc_dict) == 1): # no multiple insertion locations
                        if (previous_barcode.find('N') == -1):
                            outfile.write(previous_barcode + ',' + str(len(umi_set)) + ',' + getMostAbundantLoc(loc_dict) + '\n')
                            if (len(umi_set) > 1): 
                                multi_umi_counter += 1
                    elif (getMaxMinLocDictKey(loc_dict) <= neighbor_threshold):
                        if (previous_barcode.find('N') == -1):
                            outfile.write(previous_barcode + ',' + str(len(umi_set)) + ',' + getMostAbundantLoc(loc_dict) + '\n')
                            if (len(umi_set) > 1):
                                multi_umi_counter += 1
                        neighbor_counter += 1
                    else: 
                        drop_multiple_loc += 1
                        if (previous_barcode.find('N') != -1):
                            drop_multiple_loc_ambiguous += 1
                    # reset 'previous' infos
                    previous_barcode = current_barcode
                    previous_info = current_info
                    previous_line = line
                    umi_set = set()
                    loc_dict = {}
                # update current sets
            umi_set.add(current_umi) # check for new UMI
            if (current_loc in loc_dict):
                loc_dict[current_loc] += 1
            else:
                loc_dict[current_loc] = 1
    
    # umi_set.add(current_umi) # check for new UMI
    # current_loc = str(current_info[2]) + ',' + str(current_info[3]) + ',' + str(current_info[4])
    # if (current_loc in loc_dict):
    #     loc_dict[current_loc] += 1
    # else:
    #     loc_dict[current_loc] = 1
    
        barcode_counter += 1
        if (len(loc_dict) == 1):
            if (len(umi_set) > 1):
                multi_umi_counter += 1
            outfile.write(current_barcode + ',' + str(len(umi_set)) + ',' + getMostAbundantLoc(loc_dict) + '\n')
        elif (getMaxMinLocDictKey(loc_dict) <= neighbor_threshold):
            if (len(umi_set) > 1):
                multi_umi_counter += 1
            neighbor_counter += 1
            outfile.write(current_barcode + ',' + str(len(umi_set)) + ',' + getMostAbundantLoc(loc_dict) + '\n')
        else:
            drop_multiple_loc += 1
            if (previous_barcode.find('N') != -1):
                drop_multiple_loc_ambiguous += 1
        

    print('Total barcode count before filtering: ' + str(barcode_counter))
    print('Barcode counts with multiple UMIs: ' + str(multi_umi_counter))
    print('Barcode inserted in nearby locations: ' + str(neighbor_counter))
    print('Dropped for multiple insertion locations: ' + str(drop_multiple_loc))
    print('Dropped for multiple insertion locations,ambiguous: ' + str(drop_multiple_loc_ambiguous))

    outfile.close()
