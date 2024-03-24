# Building blast library
cd ~/ecotypes/RNA/reference/blast/
makeblastdb -in Arabidopsis_thaliana.TAIR10.cdna.all.fa -dbtype nucl \
-parse_seqids

# Generating genes fasta (uses R as well)
cd ~/ecotypes/RNA/reference/
samples='LER CVI C24 WS2'
for i in $samples
do
  bedtools getfasta -name -fi ./$i/${i}.fasta -bed ./$i/${i}_genes.bed \
  -fo ./$i/${i}_genes.fasta
done

# Blasting genes against Col0
cd ~/ecotypes/RNA/reference/
samples='LER CVI C24 WS2'
for i in $samples
do
  blastn -query ~/ecotypes/RNA/reference/$i/${i}_genes.fasta \
  -db ~/ecotypes/RNA/reference/blast/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  -out ./blast/${i}_blast.tab -evalue 1e-20 -outfmt 6 -max_target_seqs 1
done
