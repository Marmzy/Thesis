#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load EMBOSS
module load hmmer
module load perl
module load PfamScan

#This script was made to perform PfamScan on the translated pseudogene sequence.

hmmpress Pfam-A.hmm

for tool in BLAST PLAST Plant; do

	#Translating the nucleotide sequences into proteins and getting ORFs for each of the 6 frames
	for f in $HOME/Thesis/Transcriptome/Populus/Pfam/Proteins/${tool}/*.fa; do
		id=$(echo "$f" | awk -F'\/' '{print $10}' | cut -d'.' -f1)
		tail -n +2 $f > temp.txt
		sixpack -sequence temp.txt -outfile $HOME/Thesis/Transcriptome/Populus/Pfam/Sixpack/${tool}/${id}.sixpack -outseq $HOME/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}/${id}.txt
	done

	#Only keeping ORFs with a length >= 20 amino acids
	for file in $HOME/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}/*.txt; do
		id=$(echo "$file" | awk -F'\/' '{print $10}' | cut -d'.' -f1)
		tr ' ', '-' < $file > temp.txt
		
		ORFs=($(grep "^>" temp.txt))

		for orf in ${ORFs[*]}; do
			orf_id=$(echo $orf | awk -F'--' '{print $1}' | tr -d '>')
			length=$(echo $orf | awk -F'--' '{print $6}' | cut -d'a' -f1)

			if [ "$length" -ge 20 ]; then
				grep -w -A1 "$orf_id" $file >> ${id}.fa
			fi
		done

		mv ${id}.fa $HOME/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}
	done
		
	#Changing the fasta header entries, so they reflect from which sequence they originate and concatenating all the fasta files
	for fasta in $HOME/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}/*.fa; do
		entry_id=$(echo "$fasta" | awk -F'\/' '{print $10}' | cut -d'.' -f1)
		sed -i -e "s/^>/>${entry_id}/" ~/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}/${entry_id}.fa
	done

	cat $HOME/Thesis/Transcriptome/Populus/Pfam/ORF/${tool}/*.fa > ${tool}_PSgene_ORFs.fa

	#Run PfamScan on the predicted pseudogene ORFs
	pfam_scan.pl -fasta ${tool}_PSgene_ORFs.fa -dir $HOME/Thesis/Transcriptome/Populus/Pfam -outfile ${tool}_Pfam.txt -align
	
done
