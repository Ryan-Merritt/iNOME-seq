############################
#### Create directories ####
############################
cd ~
if [ ! -d ~/ecotypes/ ]
then
	x=~/ecotypes/
	mkdir ${x}
	mkdir ${x}data
	mkdir ${x}data/1_raw/
	mkdir ${x}data/2_trimmed/
	mkdir ${x}data/3_mapped/
	mkdir ${x}data/4_peaks/
	mkdir ${x}qa
	mkdir ${x}qa/1_raw/
	mkdir ${x}qa/2_trimmed/
	mkdir ${x}qa/3_mapped/
fi

cd ~/ecotypes/data/1_raw/
x='C24 Cvi LER WS2 met1'
for i in $x
do
	mkdir $i
	mv *_${i}_* $i
done

#######################
#### Trimming data ####
#######################
cd ~/ecotypes/data/1_raw/
for i in ls *R1*
do
  if echo $i | grep -q 0154
  then
    x=${i%R1.RUN0154.L3.fastq.gz}R2.RUN0154.L3.fastq.gz
    fp=${i%R1.RUN0154.L3.fastq.gz}R1.forward.paired.fastq.gz
    fu=${i%R1.RUN0154.L3.fastq.gz}R1.forward.unpaired.fastq.gz
    rp=${i%R1.RUN0154.L3.fastq.gz}R2.reverse.paired.fastq.gz
    ru=${i%R1.RUN0154.L3.fastq.gz}R2.reverse.unpaired.fastq.gz
    t=${i%R1.RUN0154.L3.fastq.gz}.log.txt
  elif echo $i | grep -q 0155
  then
    x=${i%R1.RUN0155.L8.fastq.gz}R2.RUN0155.L8.fastq.gz
    fp=${i%R1.RUN0155.L8.fastq.gz}R1.forward.paired.fastq.gz
    fu=${i%R1.RUN0155.L8.fastq.gz}R1.forward.unpaired.fastq.gz
    rp=${i%R1.RUN0155.L8.fastq.gz}R2.reverse.paired.fastq.gz
    ru=${i%R1.RUN0155.L8.fastq.gz}R2.reverse.unpaired.fastq.gz
    t=${i%R1.RUN0155.L8.fastq.gz}.log.txt
  elif echo $i | grep -q 0158
  then
    x=${i%R1.RUN0158.L6.fastq.gz}R2.RUN0158.L6.fastq.gz
    fp=${i%R1.RUN0158.L6.fastq.gz}R1.forward.paired.fastq.gz
    fu=${i%R1.RUN0158.L6.fastq.gz}R1.forward.unpaired.fastq.gz
    rp=${i%R1.RUN0158.L6.fastq.gz}R2.reverse.paired.fastq.gz
    ru=${i%R1.RUN0158.L6.fastq.gz}R2.reverse.unpaired.fastq.gz
    t=${i%R1.RUN0158.L6.fastq.gz}.log.txt
  fi
  trimmomatic PE -threads 10 $i $x ../2_trimmed/$fp ../2_trimmed/$fu ../2_trimmed/$rp ../2_trimmed/$ru ILLUMINACLIP:/home/u1896702/adapters/pe.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 2> ../../qa/2_trimmed/$t
done

cd ~/ecotypes/data/2_trimmed/
x='C24 Cvi LER WS2 met1'
for i in $x
do
	mkdir $i
	mv *_${i}_* $i
done
###########################
#### Preparing genomes ####
###########################
cd ~/reference/
for i in ls *
do
	bismark_genome_preparation --verbose --bowtie2 $i
done

#######################
#### Aligning data ####
#######################
cd ~/ecotypes/data/2_trimmed/paired/
for i in $(ls)
do
	for j in $(ls ./$i/*R1*)
	do
	x=${j%1.forward.paired.fastq.gz}2.reverse.paired.fastq.gz
	bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 ~/reference/$i/ \
	-1 $j -2 $x -o ../../${i}3_mapped
	done
done


for i in *output*; do awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $i > linear_${i}; done
for i in *linear*; do tr "\t" "\n" < $i | fold -w 50 > fold_${i}; done

cd ~/ecotypes/annotations/
for i in $(ls)
do
	bismark_genome_preparation --verbose --bowtie2 $i/
done

x='0h 6h 12h'
y='C24'
cd ~/ecotypes/data/2_trimmed/$y
for i in $x
do
	forward=$(ls | grep ${i} | grep R1)
	reverse=$(ls | grep ${i} | grep R2)
	bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 ~/ecotypes/annotations/$y/ -1 $forward -2 $reverse -o ../../3_mapped/$y/
done

x='0h 6h 12h'
y='WS2'
cd ~/ecotypes/data/2_trimmed/$y
for i in $x
do
	forward=$(ls | grep ${i} | grep R1)
	reverse=$(ls | grep ${i} | grep R2)
	bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 ~/ecotypes/annotations/$y/ -1 $forward -2 $reverse -o ../../3_mapped/$y/
done

x='0h 6h 12h'
y='LER'
cd ~/ecotypes/data/2_trimmed/$y
for i in $x
do
	forward=$(ls | grep ${i} | grep R1)
	reverse=$(ls | grep ${i} | grep R2)
	bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 ~/ecotypes/annotations/$y/ -1 $forward -2 $reverse -o ../../3_mapped/$y/
done

x='0h 6h 12h'
y='CVI'
cd ~/ecotypes/data/2_trimmed/$y
for i in $x
do
	forward=$(ls | grep ${i} | grep R1)
	reverse=$(ls | grep ${i} | grep R2)
	bismark --bowtie2 -N 1 -L 20 -p 4 -score_min L,0,-0.6 ~/ecotypes/annotations/$y/ -1 $forward -2 $reverse -o ../../3_mapped/$y/
done

for i in *bam
do
	x=${i%.forward.paired_bismark_bt2_pe.bam}_sorted
	samtools sort $i ./sorted/$x
done

cd ./sorted/
for i in *sorted.bam
do
	 samtools index $i
	 samtools flagstat $i > ./stats/${i%.bam}_flagstats.txt
	 samtools idxstats $i > ./stats/${i%.bam}_idxstats.txt
done
