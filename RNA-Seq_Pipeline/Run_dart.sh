#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load samtools/x86_64/1.6

#This file was made to test the RNA-Seq mapping tool Dart on a single sample.

./dart -i $HOME/Thesis/Genomes/Ptr -f SRR952886.trimmed.fastq.gz -unique -j junctions_SRR952886.tab -o SRR952886.sam.gz
