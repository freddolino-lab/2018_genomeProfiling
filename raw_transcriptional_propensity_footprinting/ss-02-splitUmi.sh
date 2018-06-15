#!/usr/bin/env bash

# this script is to retrieve the UMI on read 2
# UMI is at 
# ARGS: length of second part of Y-linker (the part that is at the 3' end of reads)

FILEPREFIXSET=('CviAII_Feb2018_S1' 'CviAII_Oct2017_S2' 'CviQI_Feb2018_S2' 'CviQI_Oct2017_S1')
YLINKER2LENSET=(6 4 6 4)

UMILEN="4"

YLINKER1='ACTACGCACGCGACGA'
YLINKER1LEN='16'

tempSuffix='_R2.umiYlinker.fastq.gz'

source ss-00-envSetUp.sh

for ((i=0; i<${#FILEPREFIXSET[@]};i++));
do
    filePrefix=${FILEPREFIXSET[i]}
    yLinker2Len=${YLINKER2LENSET[i]}
    filename="$rawDataSubPath""$filePrefix""$r2Suffix"
    printf '==Start of analysis==\n'
    printf '%s\n' "$filename"
    printf '=== Read Structure ===\n'
    printf '16 b Y-linker + 4 b UMI + %s b Y-linker\n' "$yLinker2Len"
    printf '======================\n'
    tmpName=${filename%$r2Suffix}
    filenameStem=${tmpName#$rawDataSubPath}
    # cut out the y-linker part 1 sequence as anchored 5' adapters (no indels) to get the last (4 + Y-linker_length) bases as UMI-Ylinker (see UMI struct info)
    printf ' + getting (UMI + Y-linker_part2)\n'
    "$CUTADAPT" --discard-untrimmed --no-indels -g ^"$YLINKER1" -o $splitSubPath$filenameStem$tempSuffix $filename
    # remove the last 6 (data from Feb2018) or 4 (data from Oct2017) part of remaining sequences to get UMI
    printf ' + getting UMI\n'
    "$CUTADAPT" -u -"$yLinker2Len" -o "$splitSubPath""$filenameStem""$umiSuffix" "$splitSubPath""$filenameStem""$tempSuffix"
    # delete temp file of construct+genomic sequence fastq.gz
    printf ' + cleaning up\n'
    rm $splitSubPath$filenameStem$tempSuffix
done

