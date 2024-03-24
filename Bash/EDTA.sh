perl ~/EDTA/EDTA.pl --genome ../fasta/C24.fasta \
--cds ../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa \
--exclude ../gene_position/C24_EDTA.bed \
--overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30

perl ~/EDTA/EDTA.pl --genome ../fasta/Cvi.fasta \
--cds ../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa \
--exclude ../gene_position/Cvi_EDTA.bed \
--overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30

perl ~/EDTA/EDTA.pl --genome ../fasta/Ler.fasta \
--cds ../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa \
--exclude ../gene_position/Ler_EDTA.bed \
--overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30

perl ~/EDTA/EDTA.pl --genome ../fasta/WS2.fasta \
--cds ../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa \
--exclude ../gene_position/WS2_EDTA.bed \
--overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30
