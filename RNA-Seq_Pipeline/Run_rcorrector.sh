#$ -S /bin/sh
#$ -pe serial 4
#$ -l h_vmem=8G

module load jellyfish/x86_64/2.2.6
module load perl/x86_64/5.14.1

#This script was made to correct the Illumina reads with Rcorrector. 

perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121292.quality.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153

for i in {293..295}
do
	perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121$i.trimmed.fastq.gz  -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
done

perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121296.quality.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121297.trimmed.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121298.quality.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121299.trimmed.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121300.trimmed.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153

for i in {301..303}
do
        perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR1121$i.quality.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE54153
done

for i in {886..903}; do
	perl run_rcorrector.pl -s /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE49911/SRR952${i}.trimmed.fastq.gz -od /group/biocomp/users/cadav/Thesis/rcorrector/Corrected/GSE49911
done
