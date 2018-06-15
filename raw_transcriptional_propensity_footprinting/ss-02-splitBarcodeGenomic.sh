#!/usr/bin/env bash

# This script is to use Cutadapt to modify Reads 1 to cut out and save the barcode as fastq.gz file while filtering ambiguous construct, and remove construct sequences from the remaining part of Read 1 to retrieve the genomic region corresponding to the barcode 
# DEPENDENCIES: need to run set-up script for env var
# INPUT: none in command line. Input files are looped through as in defined folder in set-up script
# OUTPUT: status in command line. Output files are saved to defined folder in set-up script
# USAGE: bash ss-02-splitBarcodeGenomic.sh

# dependency
source ss-00-envSetUp.sh
tempBarcodeSuffix='_R1.barcodeConstruct.fastq.gz'
tempGenomicSuffix='_R1.constructGenomic.fastq.gz'

# create split folder
mkdir $splitSubPath

CONSTRUCT_SEQ="GGCGCGCCTACCACCGTCAAAAAAAACGGCGCTTTTTAGCGCCGTTTTTATTTTTCAACCTTAGATGTGTATAAGAGACAG"
BARCODE_LEN=15
# CONSTRUCT_GENOMIC_LEN=111 # read 1 length 126, fist 15 b are barcode
GENOMIC_LEN=30

for filename in "$rawDataSubPath$dataPrefix"*"$r1Suffix";
do
    printf 'processing file: %s\n' "$filename"
    tmpName=${filename%$r1Suffix}
    filenameStem=${tmpName#$rawDataSubPath}
    # cut out to get the first 15 bases as barcode
    printf ' + splitting barcode from reads 1\n'
    printf '  + getting barcode+construct\n'
    ## remove genomic sequence from read 1
    "$CUTADAPT" -u -"$GENOMIC_LEN" -o "$splitSubPath""$filenameStem""$tempBarcodeSuffix" "$filename"
    printf '  + getting barcode\n'
    ## requiring trimmed counstruct as anchored 3' adapter
    "$CUTADAPT" --discard-untrimmed --no-indels -a "$CONSTRUCT_SEQ"$ -o "$splitSubPath""$filenameStem""$barcodeSuffix" "$splitSubPath""$filenameStem""$tempBarcodeSuffix" 
    rm "$splitSubPath""$filenameStem""$tempBarcodeSuffix" 
    ## remove the first 15 bases of barcode to get the construct+genomic part of sequences
    printf ' + splitting genomic sequences from reads 1\n'
    printf '  + getting construct+genomic\n'
    "$CUTADAPT" -u "$BARCODE_LEN" -o "$splitSubPath""$filenameStem""$tempGenomicSuffix" "$filename"
    ## remove construct and filter reads 
    printf '  + getting genomic\n'
    "$CUTADAPT" --discard-untrimmed --no-indels -g ^"$CONSTRUCT_SEQ" -o "$splitSubPath""$filenameStem""$genomicSuffix" "$splitSubPath""$filenameStem""$tempGenomicSuffix"
    # delete temp file of construct+genomic sequence fastq.gz
    rm "$splitSubPath""$filenameStem""$tempGenomicSuffix"
done
