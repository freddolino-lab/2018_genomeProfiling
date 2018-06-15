# Analysis of transcriptional propensity footprinting

Author

Publication.

## Data analysis

### Environment variables set-up

Environment variables were set up in script `ss-00-envSetUp.sh` for shorthands for paths to data, programs, and scripts, and naming conventions of in/output files at each stage. Details of each variable were described in the script.  

### Analysis pipeline, description and source code

#### Raw data quality control

Raw sequencing data quality control was performed using `FastQC, v0.10.1` and `MultiQC, version 1.0 dev0`. Script contained one command for each QC step, not provided.

#### Footprinting of barcode inseriton positions

##### Read processing: splitting barcodes and corresponding genomic sequences, with UMIs

Barcodes and corresponding genomic sequences were split from Reads 1 from footprinting sequencing data, with requirement of a recognizable construct part in each read, using `cutadapt, version 1.8.1`. See bash script `ss-02-splitBarcodeGenomic.sh`.  

UMIs for corresponding barcode-genomic combinations were split from Reads 2 from footprinting sequencing data, with a recognizable construct part and no indels in first part of Y-linker sequence in each read, using `cutadapt, version 1.8.1`. See script `ss-02-splitUmi.sh`.  

Both scripts were intended to automatically loop through all files qualifying the reqired suffices under given directory. Suffix and directory requirements were set up in environmental variable set-up process, see script `ss-00-envSetUp.sh`.  

##### Read processing: quality trimming of genomic sequences

The genomic sequences matched to barcodes, as processing of Read 1 of footprinting data in previous step, were trimmed by quality, using `Trimmomatic, version 0.33`. with parameters `-phred33` and `TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20`. Notice there is no trimming on the leading end. Script contained one command for each sequence file, not provided.  

Quality control for examination of trimmed sequences was performed using `FastQC, v0.10.1` and `MultiQC, version 1.0 dev0`. Script contained one command for each QC step, not provided.  

The extracted barcodes from the two sequencing runs for both CviAII and CviQI (four samples in total) were pooled into a single barcode pool. Genomic sequences and UMIs were similarly treated. Script contained one command for each type of sequences (barcodes, genomic sequences, and UMIs), not provided.  

##### Identifying the insertion locations

Insertion positions of barcodes were identified by alignment of genomic sequences corresponding to the barcodes. Alignment reference contained sequences from MG1655 genome ( U00096.3), pBAD-Flp, pSAS31 (both with and without a 15b wildcard barcode), and pBT1-proD-mCherry sequences, using `bowtie2-build` indexer, version 2.1.0. Script contained one command, not provided.  

Alignment of genomic sequences was performed using Bowtie2 (version 2.1.0), using `preset` of `--very-sensitive`, for each trimmed genomic sequences files. Script contained one command for each trimmed sequence file, not provided. The script used was numbered as step 04, hence the missing number 04 in script numbering.  

###### Filtering low quality barcodes

Pooled barcodes were filtered to remove any barcodes with any base of quality score below 30. One implementation is shown in `ss-05-filterBarcodeQuality.py`.  

However, this way of implementation is a horrible idea - overcomplicated, too specific, and used inconsistent file formats. This step can be easily implemented using one FASTQ file of pooled barcodes. This script is for documentation purpose, as it is the implementation used in published analysis. **Do not** use this script in any new implementation.  

##### Merging barcode sequences, UMI sequences, and insertion genomic locations

Filtered barcode sequences, UMI sequences, and insertion positions (including reference name, position on reference, and strandness), are merged into one table, while keeping only entries with all three forms of informations. See script `ss-05-mergeTables.py`.  

##### Deduplicating barcode-position pairs

First sort the table of barcodes, UMI counts, and genomic positions, by lines containing all fields, alphabetically. This can be done using a `sort` in unix utilities.  

The deduplication process is to remove completely duplicated entries that have exactly the same barcode sequences, UMI sequences, and insertion positions. The resulted unique entries are also filtered to remove barcodes with:  

 1. any ambiguous base (N);  
 2. more than one insertion positions, unless the insertion locations are close to each other (position difference smaller or equal to 2). In this special case, the position with highest abundances (identified by counts of unique UMIs) are kept, while positions with lower abundances are filtered out.  

Then for each barcode-to-position combination, count the number of unique UMIs. (At this stage, all remaining UMIs are unique.) See script `ss-05-filterDedupBarcode.py`.  

Note that, for deduplicating process, the input csv file **must** be sorted by lines (i.e. entries containing all fields) beforehand.   

A checking Bash script that serves the same purpose are also provided, as `ss-05-filterDedupBarcodeBash.sh`.  

#### Quantitation of transcriptional propensity

##### Read processing: splitting barcodes from RNA and DNA count data for transcriptional propensity

Barcodes from RNA and DNA counting sequencing runs were retrieved using `cutadapt`, with construct sequence 'GGCGCGCCTACCACCGTCAAATC' as anchored 3' adapters, and parameters of `--discard-untrimmed --no-indels`. Script not provided.  

##### Barcode counting

Counts of barcodes of barcodes, as measures of barcode abundances in each RNA or DNA sample, were counted using a custom Python script. In this process, ambiguous barcodes were removed. See `ss-03-barcodeCounting.py`.  

#### Mapping of transcriptional propensities to insertion locations

Counts of barcodes in all samples are merged into one table, keeping only barcodes with valid count in at least one sample. Script not provided, numbered as step 06.  

The table with all counts of barcodes is then merged with the barcode-location table, keeping only barcodes with any counts and insertion position. See script `ss-07-mapCountsLoc.py`.  

### Profiling HU binding landscape

HU binding data was retrieved from NCBI SRA database: data for HupA as SRR353962 and HupB as SRR353967. Both data sets are single-ended ChIP-seq data. In summary, the reads were cut for sequencing adapters by Cutadapt, trimmed for quality by Trimmomatic, and than aligned using Bowtie2, as described in detail in script `ss-06-HU.sh`. The coverage of each position, defined as the counts of templates of reads that spanned the position, was calculated using an in-house Python script, see `ss-06-HUdataSpline.py`. The source code of Python module needed for spline normalization was provided in this repository, under directory `spline_normalization_module_source_code`. To run this as a module, environmental variables would need to be changed accordingly.  


### Quantification of knock-out effects of insertions

Knock-out effects of genes refer to the cases where barcode inserted into or near genes could interfere with genes' function. To quantify the potential knock-out effects, for each gene, the median RNA/DNA ratio in the following four windows are calculated:  

 1. upstream 500 bp of each gene;  
 2. downstream 500 bp of each gene;  
 3. first 500 bp of each gene; and  
 4. last 500 bp of each gene.  

Each median calculation required at least 10 insertions within the investigated window. Fist and last 500 bp windows were allowed to overlap, and genes with length less than 500 bp were left out. See script `ss-09-koEffectCircular.NCBI.py`.  

As input, the raw RNA/DNA ratios of two replicates were quantile normalized, and for each position, the normalized ratios, or the average of normalized ratios when more than one ratio identified for the position, was the representation of transcriptional propensity for a position. This allowed the input to have a unique representation of transcriptional propensity for each position footprinted on the genome. Output took form of a csv file with gene information (start, end, strand, attributes) and the statistics for the four windows described above.  

### Estimating gene-level transcriptional propensity for functional analysis

Gene-level transcriptional propensity was described as the log2 media of quantile normalized and averaged RNA/DNA ratios within a the region of a gene and a flanking region 2500 on each side of the gene. The input, quantile normalized and average RNA/DNA ratios were obtained in same way as in 'Quantification of knock-out effects of insertions'. The list of gene-level transcriptional propensity was calculated using script `ss-10-getIpageValues.py`.  

iPAGE analysis was performed in command-line, with parameters: 1) 9 bins, and 2) allowing dependency of GO terms. The script is not provided as it was an one-line run, see [iPAGE tool page](https://tavazoielab.c2b2.columbia.edu/iPAGE/) for detailed information and manuals for this tool.  

Genes were also categorized by whether or not the coding products were recognized by SRP, base on the gene names provided in previous study (Moffitt, Jeffrey R., et al., 2016) (see Methods in publication for details). Since the provided list was in gene names, the mapping of gene names to *E. coli* b-numbers were based on annotation from NCBI for genome version U00096.3. For gene names with multiple matches as "names" or "synonyms" of b-numbered genes, the following priority for the b-numbers was used (higher priority on the top):  
 + the searching gene names were annotated as primary name;  
 + the searching gene names were annotated as one of the synonyms;  
 + the smallest b-number when there were multiple matched synonyms, but no primary gene names.  

SRP gene names were mapped to b-numbers using script `ssPytools/mapGenenameToBnum.py`, and SRP gene information and transcriptional propensity table was summarized by script `ss-10-srp.py`.  

## License

Copyright 2018 Rucheng Diao, University of Michigan.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  

### Lience information as required by programs used
GNU bash, version 4.3.11 (1), licensed under GPLv3+.  
Python 2.7, [PSF license agreement](https://docs.python.org/2.7/license.html#psf-license-agreement-for-python-release).  
Biopython, version 1.68, licensed under BSD 3-Clause License. [Full license Agreement](https://github.com/biopython/biopython/blob/master/LICENSE.rst).  
