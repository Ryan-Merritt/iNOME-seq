# creating required directories
cd ~/bioCourse/
mkdir ATAC_Seq
mkdir ATAC_Seq/QA
mkdir ATAC_Seq/data/1_raw
mkdir ATAC_Seq/data/2_trimmed
mkdir ATAC_Seq/data/3_mapped
mkdir ATAC_Seq/QA/1_raw
mkdir ATAC_Seq/QA/2_trimmed
mkdir ATAC_Seq/adapters

# Downloading required data
cd ATAC_Seq/data/1_raw
echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR400/009/SRR4000469/SRR4000469_1.fastq.gz > files.txt
echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR400/009/SRR4000469/SRR4000469_2.fastq.gz >> files.txt
wget -i files.txt

# Quality checking raw files
fastqc -o ../../QA/1_raw *fastq.gz

# Aligning trimmed reads to the reference genome
bowtie2 -p 10 -x ../../../ref/Arabidopsis_genome/tair10 -1 SRR4000469_1.fastq.gz \
-2 SRR4000469_2.fastq.gz -S ../3_mapped/SRR4000469.sam

# Converting sam to bam
cd ../3_mapped
samtools view -S -b SRR4000469.sam > SRR4000469.bam

# Sorting bam file
samtools sort SRR4000469.bam -o SRR4000469_sorted.bam

# Remove duplicates from pcr
samtools rmdup -s SRR4000469_sorted.bam SRR4000469_sorted_rmdup.bam

# Converting bam to bedgraph
bedtools genomecov -ibam SRR4000469_sorted_rmdup.bam -bg > SRR4000469_sorted_rmdup.bedgraph

# Calling peaks
macs2 callpeak -t SRR4000469_sorted_rmdup.bam -n SRR4000469_q005 -g 130e6 -q 0.05 --nomodel --extsize 200 -B --shift -100 --nolambda -f BAM --call-summits