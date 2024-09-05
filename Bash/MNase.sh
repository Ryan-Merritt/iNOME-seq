prefetch -v SRR957780
fastq-dump --outdir . /home/u1896702/ncbi/public/sra/SRR957780.sra
gzip SRR957780.fastq
mv SRR957780.fastq.gz mNASE.fastq.gz

trimmomatic SE -threads 40 -phred33 mNASE.fastq.gz \
mNASE_col_T.fastq.gz \
ILLUMINACLIP:/home/u1896702/adapters/se.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35

hisat2 -p 40 -x ~/ecotypes/annotations/Col/hisat2/Col \
-U mNASE_col_T.fastq.gz -S mNASE_col_T.sam

samtools view -S -b mNASE_col_T.sam > mNASE_col.bam
samtools sort -n mNASE_col.bam -o mNASE_col_sorted.bam
samtools fixmate -m mNASE_col_sorted.bam mNASE_col_fixmate.bam
samtools sort mNASE_col_fixmate.bam -o mNASE_col_fixmate_sorted.bam
samtools markdup mNASE_col_fixmate_sorted.bam mNASE_col_rmdup.bam

samtools view -bq 30 mNASE_col_rmdup.bam > mNASE_col_30.bam

samtools index mNASE_col_30.bam
samtools view -b mNASE_col_30.bam Chr1 Chr2 Chr3 Chr4 Chr5 > mNASE_col_30_edited.bam

samtools index mNASE_col_30_edited.bam
bamCoverage --normalizeUsing CPM --binSize 10 -b mNASE_col_30_edited.bam -o ./bw/mNASE_col_10.bw

rm -f *
