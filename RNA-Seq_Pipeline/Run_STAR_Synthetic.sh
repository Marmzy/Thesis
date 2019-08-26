#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load STAR/x86_64/2.5.2b

#Mapping the GSE54153 dataset reads (due to some samples undergone 'quality trimming' the code is messy, as the correct
#files need to be collected).

for i in {886..903}
do
	gunzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE49911/SRR952${i}.trimmed.fastq.gz
	STAR --genomeDir $HOME/Thesis/Genomes/Synthetic_Region/Index --sjdbGTFfile $HOME/Thesis/Annotations.gtf --readFilesIn /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE49911/SRR952${i}.trimmed.fastq
	
	mv Aligned.out.sam Aligned.SRR952${i}.trimmed.sam
	mv Aligned.SRR952${i}.trimmed.sam Mapped/GSE49911_Synthetic

	mv Log.final.out Log.SRR952${i}.final.out
	mv Log.SRR952${i}.final.out Mapped/GSE49911_Synthetic

	mv _STARgenome _STARgenome_SRR952${i}
	mv _STARgenome_SRR952${i} Mapped/GSE49911_Synthetic

	gzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE49911/SRR952${i}.trimmed.fastq
done

#--------------------------------------------------------------------------------------------------------------------------------#

declare quality=("292" "296" "298" "301" "302" "303")

for q in "${quality[@]}"
do
	gunzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${q}.quality.fastq.gz
	STAR --genomeDir $HOME/Thesis/Genomes/Synthetic_Region/Index --sjdbGTFfile $HOME/Thesis/Annotations.gtf --readFilesIn /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${q}.quality.fastq

	mv Aligned.out.sam Aligned.SRR1121${q}.out.sam
	mv Aligned.SRR1121${q}.out.sam Mapped/GSE54153_Synthetic

	mv Log.final.out Log.SRR1121$q.final.out
	mv Log.SRR1121${q}.final.out Mapped/GSE54153_Synthetic

	mv SJ.out.tab SJ.SRR1121${q}.out.tab
	mv SJ.SRR1121${q}.out.tab Mapped/GSE54153_Synthetic

	mv _STARgenome _STARgenome_SRR1121${q}
	mv _STARgenome_SRR1121${q} Mapped/GSE54153_Synthetic

	gzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${q}.quality.fastq
done


declare trimmed=("293" "294" "295" "297" "299" "300")

for t in "${trimmed[@]}"
do
        gunzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${t}.trimmed.fastq.gz
	STAR --genomeDir $HOME/Thesis/Genomes/Index --sjdbGTFfile $HOME/Thesis/Annotations.gtf --readFilesIn /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${t}.trimmed.fastq
	
	mv Aligned.out.sam Aligned.SRR1121${t}.out.sam
	mv Aligned.SRR1121${t}.out.sam Mapped/GSE54153_Synthetic

	mv Log.final.out Log.SRR1121${t}.final.out
	mv Log.SRR1121${t}.final.out Mapped/GSE54153_Synthetic

	mv SJ.out.tab SJ.SRR1121${t}.out.tab
	mv SJ.SRR1121${t}.out.tab Mapped/GSE54153_Synthetic

	mv _STARgenome _STARgenome_SRR1121${t}_Synthetic
	mv _STARgenome_SRR1121${t}_Synthetic Mapped/GSE54153_Synthetic

	gzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121${t}.trimmed.fastq
done

#--------------------------------------------------------------------------------------------------------------------------------#

for i in {2991..3014}
do
	gunzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE81077/SRR347${i}.trimmed.fastq.gz
        STAR --genomeDir $HOME/Thesis/Genomes/Synthetic_Region/Index --sjdbGTFfile $HOME/Thesis/Annotations.gtf --readFilesIn /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE81077/SRR347${i}.trimmed.fastq

        mv Aligned.out.sam Aligned.SRR347${i}.trimmed.sam
        mv Aligned.SRR347${i}.trimmed.sam Mapped/GSE81077_Synthetic

        mv Log.final.out Log.SRR347${i}.final.out
        mv Log.SRR347${i}.final.out Mapped/GSE49911_Synthetic

        mv _STARgenome _STARgenome_SRR347${i}
        mv _STARgenome_SRR347${i} Mapped/GSE81077_Synthetic

        gzip /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE81077/SRR347${i}.trimmed.fastq
done
