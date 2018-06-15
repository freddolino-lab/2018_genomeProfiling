#!/usr/bin/env bash

# A brutal and ugly brute force way to deduplicate barcode-position pairs, using only unix utilites, and awk. This was intended to be a double check for Python deduplicating-and-filtering results.  

# Input starts from `Cvi.barcodeUmiGenomic.csv`, a text file of comma-separated-values, with each line containing five fields: barcode sequence, UMI sequence, insertion reference, insertion positions, insertion strandness. This file doesn't have to be sorted, as sorting is the first step below. 

# Notice there is no default file paths for reading and writing files. 

# sorting
sort Cvi.barcodeUmiGenomic.csv > Cvi.barcodeUmiGenomic.sorted.csv

# filtering, using Python script as described in documentation
# python2 /path/to/script/ss-05-filterDedupBarcode.py Cvi.barcodeUmiGenomic.sorted.csv Cvi.barcodeUmiGenomic.sorted.pyFiltered.csv >05.filter.log 2>05.filter.err

# filtering, using Unix utilities and awk
echo 'counting appearence of each exact same lines'
uniq -c Cvi.barcodeUmiGenomic.sorted.csv > Cvi.barcodeUmiGenomic.sorted.uniqC.ssv
echo 'extracting unique lines'
awk '{print $2}' Cvi.barcodeUmiGenomic.sorted.uniqC.ssv > Cvi.barcodeUmiGenomic.sorted.uniq.csv

echo 'calculating uniq barcode-loc pairs regardless of umi counts'
awk -F',' '{print $1"," $3"," $4"," $5}' Cvi.barcodeUmiGenomic.sorted.uniq.csv | sort | uniq -c >bcWithLoc.ssv
awk -F',' '{print $1"," $3"," $4"," $5}' Cvi.barcodeUmiGenomic.sorted.uniq.csv | sort | uniq >bcWithLoc.csv

echo 'counting uniq barcode-loc pairs with only one insertion site regardless of umi counts'
awk -F',' '{print $1}' bcWithLoc.csv | sort | uniq -c > bcWithLocCount.ssv
wc -l bcWithLocCount.ssv
echo 'counting uniq barcode with only 1 insertion sites'
awk '{if($1==1) print}' bcWithLocCount.ssv | wc -l
echo 'counting uniq barcode with multiple insertion sites'
awk '{if($1>1) print}' bcWithLocCount.ssv | wc -l

echo 'checking how many of the non-uniquely mapped barcode are on plasmid'
awk -F',' '{print $2"," $4"," $5"," $6}' py_multiLocInfo.csv | sort | uniq > multiLocInfoBc.csv
grep -v 'NC_000913.3' multiLocInfoBc.csv | wc -l
awk -F',' '{print $2"," $3}' py_multiLocInfo.csv | sort | uniq | awk -F',' '{print $1}' | sort | uniq -c | awk '{if($1>1) print}' | wc -l


echo 'filter out barcodes with multiple insertion sites (keeping only uniquely mapped barcodes)'
awk -F',' '{print $1"," $3"," $4"," $5}' Cvi.barcodeUmiGenomic.sorted.uniq.csv | sort | uniq | awk -F',' '{print $1}' | uniq -c | awk '{if($1==1) print $2}' > bcUmiWithUniqLoc.ssv

echo 'map barcodes to UMIs and locations'
join -1 1 -2 1 -t' ' <(sort bcUmiWithUniqLoc.ssv) <(sort Cvi.barcodeUmiGenomic.sorted.uniq.csv | sed 's/,/ /1') | sed 's/ /,/' > bcUmiLocWithUniqLoc.csv

echo 'count UMI for each unquely inserted barcode'
awk -F',' '{print $1"," $3"," $4"," $5}' bcUmiLocWithUniqLoc.csv | sort | uniq -c > umiCountWithUniqBcLoc.ssv

echo 'format into bacode-UMIcounts-location'
paste -d',' <(cat umiCountWithUniqBcLoc.ssv | awk '{print $2}' | awk -F',' '{print $1}') <(cat umiCountWithUniqBcLoc.ssv | awk '{print $1}') <(cat umiCountWithUniqBcLoc.ssv | awk '{print $2}' | awk -F',' '{print $2"," $3",", $4}') >umiCountWithUniqBcLoc.formated.csv

echo 'remove ambiguous barcodes (barcodes with N)'
awk -F',' '!($1 ~ /N/) {print}' umiCountWithUniqBcLoc.formated.csv > umiCountWithUniqBcLoc.formated.unambiguous.csv

# To this point, the same goal as the Python script has been achieved. File names are different.  

echo 'counting unambiguous barcodes on genome'
awk -F',' '{if($3=="NC_000913.3") print}' umiCountWithUniqBcLoc.formated.unambiguous.csv > umiCountWithUniqBcLoc.formated.unambiguous.genome.csv

echo 'counting unambiguous barcodes on genome with UMI>1'
awk -F',' '{if($2>1) print}' umiCountWithUniqBcLoc.formated.unambiguous.genome.csv | wc -l

