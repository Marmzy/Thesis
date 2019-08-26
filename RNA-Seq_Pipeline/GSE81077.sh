#$ -S /bin/sh
#$ -pe serial 2
#$ -l h_vmem=4G

module load FastQC

#This script was made to download the GEO accession GSE81077 fastq files from EBI.
#This script was made to perform FastQC analysis a first time on all the .fastq.gz files to trim the adapter sequences, if
#this hadn't been done yet. Then, if 'Per base sequence quality' FastQC analysis reports a 'fail' for any .fastq.gz file, a
#quality trimming of the reads will be performed. SLIDINGWINDOW:5:28 will cut the read if the average quality over 5 bases
#falls below a phred score of 28 (good reads). MAXINFO:90:0.1 will trim the read as long as they have a minimum length of 90
#and it will favour longer reads, rather than correct reads.

for i in {3472991..3473014}
do
	num="${i: -1}"
        wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR347/00$num/SRR$i/SRR$i.fastq.gz
done

for num in {3472991..3473014}
do
	java -jar /software/shared/apps/x86_64/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 SRR$num.fastq.gz SRR$num.trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:35
	fastqc -o TrimmedQC --extract -f fastq SRR$num.trimmed.fastq.gz

	base=$(grep -E 'Per base sequence quality' /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE81077/TrimmedQC/SRR$num.trimmed_fastqc/fastqc_data.txt)
        result=$(echo $base | tail -c 5)

        if [ $result = 'fail' ]
        then
                java -jar /software/shared/apps/x86_64/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE81077/SRR$num.trimmed.fastq.gz /group/biocomp/users/cadav/Thesis/Transcriptome/RNA-Seq/GSE54153/SRR$num.quality.fastq.gz SLIDINGWINDOW:5:28 MAXINFO:90:0.1
		fastqc -o QualityQC --extract -f fastq SRR$num.quality.fastq.gz
        fi
done


