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

for i in {292..303}; do
        samtools view -q 255 -Sb Aligned.SRR1121${i}.out.sam > Unique.SRR1121${i}.bam
        samtools sort -o Unique.SRR1121${i}.sorted Unique.SRR1121${i}.bam
done

for i in {2991..3014}; do
        samtools view -q 255 -Sb Aligned.SRR347${i}.sam > Unique.SRR347${i}.bam
        samtools sort -o Unique.SRR347${i}.sorted Unique.SRR347${i}.bam
done
