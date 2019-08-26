#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load python
module load samtools/x86_64/1.6

#This script was made to first filter the alignment file to keep only uniquely mapping reads. For STAR, this equates to a score
#of 255 (https://groups.google.com/forum/#!topic/rna-star/pQbqeQd0lNY).

for i in {886..903}; do
	samtools view -q 255 -Sb Aligned.SRR952${i}.trimmed.sam > Unique.SRR952${i}.bam
	samtools sort -o Unique.SRR952${i}.sorted Unique.SRR952${i}.bam
done


