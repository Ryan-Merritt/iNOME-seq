# creating required directories
cd ~/bioCourse/
mkdir ATAC_Seq_root
mkdir ATAC_Seq_root/QA
mkdir ATAC_Seq_root/data
mkdir ATAC_Seq_root/data/1_raw
mkdir ATAC_Seq_root/data/2_trimmed
mkdir ATAC_Seq_root/data/3_mapped
mkdir ATAC_Seq_root/QA/1_raw
mkdir ATAC_Seq_root/QA/2_trimmed
mkdir ATAC_Seq_root/adapters

# Downloading required data (ATAC-seq data root)
cd ATAC_Seq_root/data/1_raw
cat SRR4000472.txt | cut -f 10 | sed -n 2,2p > files.txt | sed -i 's/;/\n/g' files.txt | rm SRR4000472.txt
wget -i files.txt

# quality checking raw data
fastqc -o ../../QA/1_raw *fastq.gz

# Aligning data to the reference genome
bowtie2 -p 10 -x ../../../ref/Arabidopsis_genome/tair10 -1 SRR4000472_1.fastq.gz \
-2 SRR4000472_2.fastq.gz -S ../3_mapped/SRR4000472.sam

# Converting sam to bam
cd ../3_mapped
samtools view -S -b SRR4000472.sam > SRR4000472.bam

# Sorting bam file
samtools sort -n SRR4000472.bam -o SRR4000472_sorted.bam

# Removing pcr duplicates
samtools fixmate -m SRR4000472_sorted.bam SRR4000472_fixmate.bam
samtools sort SRR4000472_fixmate.bam -o SRR4000472_fixmate_sorted.bam
samtools markdup SRR4000472_fixmate_sorted.bam SRR4000472_rmdup.bam
rm SRR4000472_fixmate.bam SRR4000472_fixmate_sorted.bam SRR4000472_sorted.bam 
# Converting to bedgraph
bedtools genomecov -ibam SRR4000472_rmdup.bam -bg > SRR4000472_rmdup.bedgraph

# Calling peaks
macs2 callpeak -t SRR4000472_rmdup.bam -n SRR4000472_q005 -g 130e6 -q 0.05 --nomodel --extsize 200 -B --shift -100 --nolambda -f BAM --call-summits