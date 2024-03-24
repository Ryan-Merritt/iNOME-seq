cdir -d Col0
####-####
cd Col0/data/1_raw/
prefetch -v SRR12082882
prefetch -v SRR12082883
fastq-dump --split-files --gzip --outdir . /home/u1896702/ncbi/public/sra/SRR12082882.sra
fastq-dump --split-files --gzip --outdir . /home/u1896702/ncbi/public/sra/SRR12082883.sra
mv SRR12082882_1.fastq.gz Col0_1_1.fastq.gz
mv SRR12082882_2.fastq.gz Col0_1_2.fastq.gz
mv SRR12082883_1.fastq.gz Col0_2_1.fastq.gz
mv SRR12082883_2.fastq.gz Col0_2_2.fastq.gz

conda activate fastqc
for i in *
do
  fastqc -q -o ../../qa/1_raw/ $i &
done
cd ../../qa/1_raw/
conda activate multiqc
multiqc .

cd ../../data/1_raw/
for i in *_1.*
do
  trimmomatic PE -threads 50 $i ${i/_1./_2.} ../2_trimmed/$i \
  ../2_trimmed/delete ../2_trimmed/${i/_1./_2.} ../2_trimmed/delete \
  ILLUMINACLIP:/home/u1896702/adapters/pe.fa:2:30:10 LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 MINLEN:50 2> ../../qa/2_trimmed/${i/.fastq.gz/.log} &
done

cd ../2_trimmed/
rm -f delete
conda activate fastqc
for i in *
do
  fastqc -q -o ../../qa/2_trimmed/ $i &
done

cd ../../qa/2_trimmed/
conda activate multiqc
multiqc .

cd ../../data/2_trimmed/
conda activate star
(mkdir ~/Col0/reference/Col0_STAR ;
STAR --runThreadN 10 --runMode genomeGenerate \
  --genomeDir ../../reference/Col0_STAR \
  --genomeFastaFiles ../../reference/TAIR10_chr_all.fasta \
  --sjdbGTFfile ../../reference/TAIR10_genes_transposons.gtf \
  --sjdbOverhang 75 --genomeSAindexNbases 12) &

faCount ../../reference/TAIR10_chr_all.fasta | cut -f 2,7 | awk 'END{print}' > ../../reference/Col0_effective_genome.txt
genomeSize=$(expr $(cut -f1 ../../reference/Col0_effective_genome.txt) - $(cut -f2 ../../reference/Col0_effective_genome.txt))

conda activate star
for i in *_1.*
do
  STAR --genomeDir ~/Col0/reference/Col0_STAR \
  --readFilesIn $i ${i/_1./_2.} \
  --runThreadN 10 -alignIntronMax INTRONMAX --outSAMstrandField intronMotif \
  --sjdbGTFfile ../../reference/TAIR10_genes_transposons.gtf --sjdbGTFtagExonParentGene gene_name \
  --readFilesCommand zcat --outFileNamePrefix ../3_mapped/${i%_*}_ \
  --outReadsUnmapped Fastx --quantMode GeneCounts \
  --outSAMtype BAM Unsorted &
done

conda activate samtools
cd ../3_mapped/
mv *ReadsPerGene* ../4_stats/
for i in *bam
do
  genome=${i%%_*}
  output=${i%_*}
  (samtools sort $i -o ${i/Aligned/sorted} &&
  samtools index ${i/Aligned/sorted} &&
  bamCoverage --numberOfProcessors 10 --effectiveGenomeSize ${genomeSize[$genome]} \
  --normalizeUsing CPM --binSize 10 -b ${i/Aligned/sorted} -o ../4_stats/${output}_10.bw) &
done

cd ../4_stats/
for i in *ReadsPerGene.out.tab
do
  awk 'NR>4 {print $1 "\t" $2}' $i > ${i/ReadsPerGene.out.tab/DEseq.input} &
done

ecotypes='Col0'
mkdir $ecotypes
for i in $ecotypes
do
  (cp ~/Col0/data/3_mapped/${i}*bam ./$i/
  cp ~/Col0/data/4_stats/${i}* ./$i/) &
done
