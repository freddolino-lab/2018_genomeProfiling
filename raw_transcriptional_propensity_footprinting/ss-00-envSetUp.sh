#!/usr/bin/env bash

# This script is an example to set up environment variables for analyzing genome profiling footprinting data.
# Environmental variables include paths to data and softwares, file name conventions, etc. 

# overall env vars
## Path for raw sequencing data
runDataPathOct2017='' # raw sequencing data for footprinting, 201710
runDataPathFeb2018='' # raw sequencing data for footprinting, 201802
transcriptDataPath2017='' # raw sequencing data for transcriptional propensity (RNA, DNA, and KanR retaining DNA samples)
# dataDirPathRun1='/home/diaorch/data/ss-Feb2018/'
dataDirPath='' # project folder for all outputs and intermediate results
dataPrefix='Cvi' # raw data prefix for sequencing data, requires renaming after symbolically linking raw sequencing data files from `runDataPath*` to `rawDataSubPath` directory below

## software paths
### Env var under this section are names to softwares used in analysis. When needed, i.e. the executable software files are not added to global environmental variable and cannot be called at all positions, the paths of the software should be added to the variable.
# For example, when cutadapt is not added to global environmental variable, the first line below should be `CUTADAPT='<path/to/>cutadapt'`
CUTADAPT='cutadapt' # path to executable cutadapt 
TRIMMOMATICSPATH='trimmomatic-0.33.jar' # path to executable Trimmomatic jar file
BOWTIE2BUILDPATH='bowtie2-build' # path to executable bowtie2-build
BOWTIE2PATH='bowtie2' # path to executable bowtie2 aligner
SAMTOOLSPATH='samtools' # path to executable

## in-house scripts
### Env var under section are paths and names to the custom generalized tools used in analysis. These scripts are uploaded under the same GitHub repository: `2018_genomeProfiling/raw_transcriptional_propensity_footprinting/ssPytools/`
projectPath='ssPytools/' # path to ssPytools, change when needed
SAM2FASTAPATH="$projectPath"sam2fasta.py # path to SAM-to-FASTA conversion script; see script docstring
FASTQ2CSVPATH="$projectPath"fastqSeq2csv.py # path to FASTQ-TO-CSV conversion script, csv file contains sequence names, and sequences; see script docstring

# analysis
## quality control
rawDataSubPath=$dataDirPath'00-seq/' # sub-directory for symbolic links to raw data
qcSubPath=$dataDirPath'01-qc/' # sub-directory for quality control outputs

## barcode-genomic splitting
splitSubPath=$dataDirPath'02-split/' # sub-directory for fastq files of split genomic sequences, barcodes, and UMIs
r1Suffix='_R1_001.fastq.gz' # raw data read 1 suffix
r2Suffix='_R2_001.fastq.gz' # raw data read 2 suffix
barcodeSuffix='_R1.barcode.fastq.gz' # split barcode file suffix
genomicSuffix='_R1.genomic.fastq.gz' # split footprint genomic sequence file suffix
### Sequences of identified constuct are defined in the script `ss-02-splitBarcodeGenomic.sh`.

## UMI retrieval
## UMI length, length of preceeding and follwing sequences, are defined in the script `ss-02-splitMiniBarcode.sh`.
umiSuffix='_R2.umi.fastq.gz' # split UMI sequence file suffix

## pooling sequences
poolSubPath="$dataDirPath"03-pool/ # sub-directory for pooled results of barcode sequences, footprint genomic sequences, and UMIs, from CviAII and CviQI digested, two batches (201710 and 201802) samples
poolSuffix='.pool.fastq.gz' # pooled data suffix

## quality trimming
trimSubPath="$dataDirPath"03-pool/ # sub-directory for quality trimming of footprint genomic sequences
trimSuffix='.trim.fastq.gz' # quality-trimmed genomic sequence file suffix

## alignment
refSubPath="$dataDirPath"04-ref/ # sub-direcotry for custom-built alignment references, including genome and plasmid sequences, see publication Method section
refFasta="$refSubPath"ss_ref.fasta # reference FASTA file suffix
refBt2Base="$refSubPath"ss_ref # Bowtie2 reference name base
alignSubPath="$dataDirPath"04-align/ # sub-directory for alignment results
alignedSuffix='.sam' # alignment result suffix

## summary of barcode, UMI, and genomic location table
sumSubPath="$dataDirPath"05-sum/ # sub-directory for pairing barcode sequences to genomic positions, i.e. the footprinting process

## transcriptional propensity
transcriptSubPath="$dataDirPath"06-transcriptional/ # path to master table of RNA and DNA counts of each barcode in replicate samples
transcriptDataSuffix='.fastq.gz' # suffix for split barcode file for RNA and DNA counts
transcriptCutSuffix='.barcode.fastq.gz' # suffix for split barcode file for RNA and DNA counts
transcriptCountTable='.count.tsv' # suffix of master table file of RNA and DNA counts in replicates

## megeing barcode counts and location and calculate GC content
mapSubPath="$dataDirPath"07-map/ # sub-directory for merging RNA and DNA counts to footprint of barcodes 
