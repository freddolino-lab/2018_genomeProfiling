# download data
~/packages/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch.2.8.2 SRR353967
~/packages/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch.2.8.2 SRR353962
mv ~/data/ncbi/sra/SRR35396* .

# Convert .sra to .fastq files
~/packages/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump.2.8.2 -O . SRR353962.sra
~/packages/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump.2.8.2 -O . SRR353967.sra

# Organizing files
mkdir 60_seq
mv S* 60_seq/

# QC
mkdir 61_qc
fastqc -o 61_qc/ 60_seq/SRR353962.fastq 60_seq/SRR353967.fastq

# Soft link
# ln -s 60_seq/SRR353962.fastq 60_seq/HupA.fastq
# ln -s 60_seq/SRR353967.fastq 60_seq/HupB.fastq
# Cutadapt cannot take soft link

# Cut and trim
mkdir 62_cutAndTrim

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o 62_cutAndTrim/HupA.cutadapt.fastq.gz 60_seq/SRR353962.fastq > 62-HupA.cutadapt.log 2> 62-HupA.cutadapt.err
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o 62_cutAndTrim/HupB.cutadapt.fastq.gz 60_seq/SRR353967.fastq > 62-HupB.cutadapt.log 2> 62-HupB.cutadapt.err

TrimmomaticSE -threads 8 -phred33 62_cutAndTrim/HupA.cutadapt.fastq.gz 62_cutAndTrim/HupA.cutadapt.trimmo.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >62-HupA.trimmo.log 2>62-HupA.trimmo.err
TrimmomaticSE -threads 8 -phred33 62_cutAndTrim/HupB.cutadapt.fastq.gz 62_cutAndTrim/HupB.cutadapt.trimmo.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >62-HupB.trimmo.log 2>62-HupB.trimmo.err

# QC
mkdir 63_qc
fastqc -o 63_qc/ 62_cutAndTrim/HupA.cutadapt.trimmo.fastq.gz 62_cutAndTrim/HupB.cutadapt.trimmo.fastq.gz

# HupA HupB paired-end adapter cutting
cutadapt -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -o 62_cutAndTrim/HupA.cutadapt.kmer.fastq.gz 62_cutAndTrim/HupA.cutadapt.fastq.gz > 62-HupA.cutadapt.kmer.log 2> 62-HupA.cutadapt.kmer.err
cutadapt -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -o 62_cutAndTrim/HupB.cutadapt.kmer.fastq.gz 62_cutAndTrim/HupB.cutadapt.fastq.gz > 62-HupB.cutadapt.kmer.log 2> 62-HupB.cutadapt.kmer.err

TrimmomaticSE -threads 8 -phred33 62_cutAndTrim/HupA.cutadapt.kmer.fastq.gz 62_cutAndTrim/HupA.cutadapt.kmer.trimmo.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >62-HupA.trimmo.kmer.log 2>62-HupA.trimmo.kmer.err
TrimmomaticSE -threads 8 -phred33 62_cutAndTrim/HupB.cutadapt.kmer.fastq.gz 62_cutAndTrim/HupB.cutadapt.kmer.trimmo.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >62-HupB.trimmo.kmer.log 2>62-HupB.trimmo.kmer.err

# QC
fastqc -o 63_qc/ 62_cutAndTrim/HupA.cutadapt.kmer.trimmo.fastq.gz 62_cutAndTrim/HupB.cutadapt.kmer.trimmo.fastq.gz

# Alignment
mkdir 64_align
INDEX_NAME="/srv/diaorch/Escherichia_coli/MG1655/Regulondb_20170908/u00096.3"

bowtie2 --very-sensitive -p 4 -x "$INDEX_NAME" -U 62_cutAndTrim/HupA.cutadapt.kmer.trimmo.fastq.gz -S 64_align/HupA.sam > 4-HupA_align.log 2> 4-HupA_align.err

bowtie2 --very-sensitive -p 4 -x "$INDEX_NAME" -U 62_cutAndTrim/HupB.cutadapt.kmer.trimmo.fastq.gz -S 64_align/HupB.sam > 4-HupB_align.log 2> 4-HupB_align.err

# SAM to BAM, sort and index BAM
samtools view -Sb 64_align/HupA.sam >64_align/HupA.bam
samtools view -Sb 64_align/HupB.sam >64_align/HupB.bam

rm 64_align/HupA.sam
rm 64_align/HupB.sam

samtools sort -@ 8 64_align/HupA.bam -o 64_align/HupA.sorted.bam
samtools sort -@ 8 64_align/HupB.bam -o 64_align/HupB.sorted.bam

samtools index 64_align/HupA.sorted.bam
samtools index 64_align/HupB.sorted.bam

# Coverage
mkdir 65_coverage
vi 65_coverage/Hup.coverage.cfg
echo $PYTHONPATH

mkdir 65_coverage/HupA
python /home/mbwolfe/src/RNA_IPOD_scripts/map_coverage.py 65_coverage/Hup.coverage.cfg 64_align/HupA.sorted.bam 65_coverage/HupA/

mkdir 65_coverage/HupB
python /home/mbwolfe/src/RNA_IPOD_scripts/map_coverage.py 65_coverage/Hup.coverage.cfg 64_align/HupB.sorted.bam 65_coverage/HupB/

# Peak calling
mkdir 66_peakCalling
mkdir 66_peakCalling/HupA
mkdir 66_peakCalling/HupB

# Coverage normalized
mkdir 67_coverageNormalized
mkdir 67_coverageNormalized/HupA
mkdir 67_coverageNormalized/HupB

python scott-67-splineCoverage.py

