######################################
# chIP-seq histone H3 acetylation WT #
######################################

# Downloading data from ncbi
prefetch -O ~/bioCourse/chip_seq/raw/ -v SRR8374225

# Converting downloaded file to fastq (not split as single ended)
fastq-dump --outdir ~/bioCourse/chip_seq/fastq/ ~/bioCourse/chip_seq/raw/SRR8374225.sra

# Quality checking the data using fastqc
fastqc -o ~/bioCourse/chip_seq/QA/ ~/bioCourse/chip_seq/fastq/SRR8374225.fastq

# trimming file
bbduk.sh in=./fastq/SRR8374225.fastq ref=./QA/adapters.fa out=./fastq/SRR8374225_trimmed.fastq

# Quality checking file after trimming
fastqc -o ~/bioCourse/chip_seq/QA/ ~/bioCourse/chip_seq/QA/SRR8374225_trimmed.fastq

# Download reference genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

# Creating bowtie2 indexed genome
bowtie2-build TAIR10_chr_all.fas tair10

# Aligning trimmed fastq to prepared genome
bowtie2 -x ~/bioCourse/ref/Arabidopsis_genome/tair10 -U ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed.fastq \
-S ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed.sam

# Filtering sam file
samtools view -Sh ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed.sam \
| grep -e "^@" -e "XM:i:[012][^0-9]" | grep -v "XS:i:" > ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered.sam

# convert sam to bam
samtools view -S -b ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered.sam \
> ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered.bam

# Sortng bam file
samtools sort ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered.bam \
-o ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted.bam

# Remove duplicates from PCR
samtools rmdup -s ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted.bam \
~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted_rmdup.bam

# Index the bam file to allow for visualisation in IGV
samtools index ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted_rmdup.bam

# Convert bam files to bed files
bedtools genomecov -ibam ~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted_rmdup.bam -bg > \
~/bioCourse/chip_seq/fastq/SRR8374225_trimmed_filtered_sorted_rmdup.bedgraph

############
# ATAC-seq #
############

# Downloading data from ncbi
prefetch -O ~/bioCourse/atac_seq/raw -v SRR4000468

# Generating fastq outputs from downloaded data --split-files as the data is paired
fastq-dump --outdir ~/bioCourse/atac_seq/fastq/ --split-files \
~/bioCourse/atac_seq/raw/SRR4000468.sra

# Quality checking data
fastqc -o ~/bioCourse/atac_seq/QA ~/bioCourse/atac_seq/fastq/*

# Trimming data
bbduk.sh in1=~/bioCourse/atac_seq/fastq/SRR4000468_1.fastq in2=~/bioCourse/atac_seq/fastq/SRR4000468_2.fastq \
ref=~/bioCourse/atac_seq/QA/adapters.fa out1=~/bioCourse/atac_seq/fastq/SRR4000468_1_trimmed.fastq \
out2=~/bioCourse/atac_seq/fastq/SRR4000468_2_trimmed.fastq

# Quality checking trimmed data
fastqc -o ~/bioCourse/atac_seq/QA ~/bioCourse/atac_seq/fastq/*trimmed*

# aligning trimmed reads to the reference genome
bowtie2 -x ~/bioCourse/ref/Arabidopsis_genome/tair10 -1 ~/bioCourse/atac_seq/fastq/SRR4000468_1_trimmed.fastq \
-2 ~/bioCourse/atac_seq/fastq/SRR4000468_1_trimmed.fastq -S ~/bioCourse/atac_seq/fastq/SRR4000468.sam

# convert sam to bam
samtools view -S -b ~/bioCourse/atac_seq/fastq/SRR4000468.sam \
> ~/bioCourse/atac_seq/fastq/SRR4000468.bam

# Sortng bam file
samtools sort ~/bioCourse/atac_seq/fastq/SRR4000468.bam \
-o ~/bioCourse/atac_seq/fastq/SRR4000468_sorted.bam

# Remove duplicates from PCR
samtools rmdup -s ~/bioCourse/atac_seq/fastq/SRR4000468_sorted.bam \
~/bioCourse/atac_seq/fastq/SRR4000468_sorted_rmdup.bam

# Indexing bam file to allow for visualisation in igv
samtools index ~/bioCourse/atac_seq/fastq/SRR4000468_sorted_rmdup.bam

# calling peaks
macs2 callpeak -t ~/bioCourse/ATAC_seq/SRR4000468_sorted_rmdup.bam -n SRR4000468_ -g 130e6 -q 0.05 --nomodel --extsize 200 -B --shift -100 --nolambda  -f BAM --call-summits

# Converting bam to bedgraph
bedtools genomecov -ibam ~/bioCourse/atac_seq/fastq/SRR4000468_sorted_rmdup.bam -bg > \
~/bioCourse/atac_seq/fastq/SRR4000468_sorted_rmdup.bed

###########
# RNA-seq # 
###########

# Downloading required rna seq data
echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/SRR8380412/SRR8380412_1.fastq.gz > ~/biocourse/rna_seq/fastq/files.txt
echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR838/002/SRR8380412/SRR8380412_2.fastq.gz >> ~/biocourse/rna_seq/fastq/files.txt
wget -i ~/biocourse/rna_seq/fastq/files.txt

# Unziping files
gunzip ~/biocourse/rna_seq/fastq/*

# Quality checking downloaded data
fastqc -o ~/bioCourse/rna_seq/QA ~/bioCourse/rna_seq/fastq/*

# Preparing the genome for hisat2
hisat2-build ~/bioCourse/rna_seq/hisat/tair10 ~/bioCourse/ref/raw/TAIR10_chr_all.fas

# Aligning reads to the reference genome
hisat2 -p 10 -x ~/bioCourse/ref/hisat/tair10 -1 ~/bioCourse/rna_seq/fastq/SRR8380412_1.fastq \
-2 ~/bioCourse/rna_seq/fastq/SRR8380412_2.fastq -S ~/bioCourse/rna_seq/fastq/SRR8380412.sam