cdir -d new_rna
cd new_rna/data/1_raw/
cp ~/mnt/Shared210/Read_database/2020_10_23_Novogene_tmp/X204SC20091837-Z01-F002_1/raw_data/NH*/*fq.gz .
cp ~/mnt/Shared210/Read_database/2020_10_23_Novogene_tmp/X204SC20091837-Z01-F002_2/raw_data/NH*/*fq.gz .

conda activate fastqc
for i in *
do
  fastqc -q -o ../../qa/1_raw/ $i &
done
cd ../../qa/1_raw/
conda activate multiqc
multiqc .

cd ../../data/1_raw/
for i in *_1*
do
  trimmomatic PE -threads 50 $i ${i/_1/_2} ../2_trimmed/$i \
  ../2_trimmed/delete ../2_trimmed/${i/_1/_2} ../2_trimmed/delete \
  ILLUMINACLIP:/home/u1896702/adapters/pe.fa:2:30:10 LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 MINLEN:50 2> ../../qa/2_trimmed/${i/.fq.gz/.log} &
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
declare -A ecotypes1
ecotypes1=(["Ler"]="NH01" ["Cvi"]="NH02" ["WS2"]="NH03" ["C24"]="NH04")
conda activate star
for i in "${!ecotypes1[@]}"
do
  (mkdir ~/new_rna/reference/${i}_STAR ;
  STAR --runThreadN 10 --runMode genomeGenerate \
  --genomeDir ../../reference/${i}_STAR \
  --genomeFastaFiles ../../reference/${i}.fasta \
  --sjdbGTFfile ../../reference/${i}.gff3 \
  --sjdbOverhang 147 --genomeSAindexNbases 12) &
done

declare -A genomeSize
for i in "${!ecotypes1[@]}"
do
  faCount ../../reference/${i}.fasta | cut -f 2,7 | awk 'END{print}' > ../../reference/${i}_effective_genome.txt
  genomeSize[$i]=$(expr $(cut -f1 ../../reference/${i}_effective_genome.txt) - $(cut -f2 ../../reference/${i}_effective_genome.txt))
done

conda activate star
for i in "${!ecotypes1[@]}"
do
  STAR --genomeDir ~/new_rna/reference/${i}_STAR \
  --readFilesIn ${ecotypes1[$i]}_1.fq.gz ${ecotypes1[$i]}_2.fq.gz \
  --runThreadN 10 -alignIntronMax INTRONMAX --outSAMstrandField intronMotif \
  --sjdbGTFfile ../../reference/${i}.gff3 --sjdbGTFtagExonParentGene gene_id \
  --readFilesCommand zcat --outFileNamePrefix ../3_mapped/${i}_1_ \
  --outReadsUnmapped Fastx --quantMode GeneCounts \
  --outSAMtype BAM Unsorted &
done

declare -A ecotypes2
ecotypes2=(["Ler"]="NH05" ["Cvi"]="NH06" ["WS2"]="NH07" ["C24"]="NH08")
conda activate star
for i in "${!ecotypes2[@]}"
do
  STAR --genomeDir ~/new_rna/reference/${i}_STAR \
  --readFilesIn ${ecotypes2[$i]}_1.fq.gz ${ecotypes2[$i]}_2.fq.gz \
  --runThreadN 10 -alignIntronMax INTRONMAX --outSAMstrandField intronMotif \
  --sjdbGTFfile ../../reference/${i}.gff3 --sjdbGTFtagExonParentGene gene_id \
  --readFilesCommand zcat --outFileNamePrefix ../3_mapped/${i}_2_ \
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

ecotypes='Cvi C24 Ler WS2'
mkdir $ecotypes
for i in $ecotypes
do
  (cp ~/new_rna/data/3_mapped/${i}*bam ./$i/
  cp ~/new_rna/data/4_stats/${i}* ./$i/) &
done
