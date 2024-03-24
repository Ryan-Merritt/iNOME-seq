directory='C24'
cdir -d $directory

cd ./$directory/data/1_raw/
cp ~/mnt/Shared210/Read_database/Arabidopsis/2019_10_14_MCViP_ecotypes_BS_seq/*C24* .
cp ~/mnt/Shared210/Read_database/Arabidopsis/2019_11_11_MCViP_ecotypes_BS_seq_resequencing/*C24* .
rm -rf *_S.*

for i in *-*
do
  x=${i//-/_}
  x=${x/0159/0155}
  x=${x/L6/L8}
  cat $x >> $i
done

rm -rf *0155*

for file in *R1*
do
  (forward=$file
  reverse=${file/R1/R2}
  forwardPaired=${forward%R1*.fastq.gz}R1_paired.fastq.gz
  reversePaired=${reverse%R2*.fastq.gz}R2_paired.fastq.gz
  forwardUnpaired=${forward%R1*.fastq.gz}R1_unpaired.fastq.gz
  reverseUnpaired=${reverse%R2*.fastq.gz}R2_unpaired.fastq.gz
  log=${file%R1*}.log.txt
  trimmomatic PE -threads 10 $forward $reverse ../2_trimmed/$forwardPaired \
  ../2_trimmed/$forwardUnpaired ../2_trimmed/$reversePaired \
  ../2_trimmed/$reverseUnpaired \
  ILLUMINACLIP:/home/u1896702/adapters/pe.fa:2:30:10 LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 MINLEN:35 2> ../../qa/2_trimmed/$log) &
done

cd ../../reference/
bismark_genome_preparation --verbose --bowtie2 .

cd ../data/2_trimmed/
for file in *R1*
do
  (forward=$file
  reverse=${file/R1/R2}
  bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 \
  ../../reference/ -1 $forward -2 $reverse -o ../3_mapped/) &
done

cd ../3_mapped/
for i in *bam
do
  deduplicate_bismark --bam $i
done

for i in *deduplicated.bam
do
  bismark_methylation_extractor --bedGraph --CX_context \
  --cytosine_report --genome_folder ../../reference/ \
  -o ../4_stats/ $i
done
