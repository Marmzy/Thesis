#$ -S /bin/sh
#$ -pe serial 2
#$ -l h_vmem=4G

module load python

#This script was made to get the Reciprocal Best Hits between the predicted pseudogene regions and the proteins (blastx vs.
#tblastn). This, to make a prediction on the predicted pseudogene parent genes.

#Preparing and filtering the files
grep -v "^#" tblastnBLAST.blasted | awk '$11 ~ /e/ {print $0}' > tblastnBLAST.txt
cat tblastnPLAST.blasted | awk '$11 ~ /e/ {print $0}' > tblastnPLAST.txt
grep -v "^#" tblastnPlant.blasted | awk '$11 ~ /e/ {print $0}' > tblastnPlant.txt

grep -v "^#" BlastxBLAST.blasted | awk '$11 ~ /e/ {print $0}' > blastxBLAST.txt
cat BlastxPLAST.blasted | awk '$11 ~ /e/ {print $0}' > blastxPLAST.txt
grep -v "^#" BlastxPlant.blasted | awk '$11 ~ /e/ {print $0}' > blastxPlant.txt

#The tblastn results are not sorted by e-value correctly (by target and then e-value), so we sort by e-value in order to correctly find the RBH
sort -k1,1 -k11,11g -o tblastnBLAST.txt tblastnBLAST.txt
sort -k1,1 -k11,11g -o tblastnTLAST.txt tblastnPLAST.txt
sort -k1,1 -k11,11g -o tblastnPlant.txt tblastnPlant.txt

#Getting all the RBHs
#NOTE: Manually change the output file name in the python file!!
python3 RBH.py tblastnBLAST.txt blastxBLAST.txt
python3 RBH.py tblastnPLAST.txt blastxPLAST.txt
python3 RBH.py tblastnPlant.txt blastxPlant.txt


#Some interanchor regions have a score of 99999, which indicates they did not have a best protein. In this case, the protein with which they were discovered (using tblastn) will stand in as the parent gene
grep '99999' BLAST_RBH.txt | while read line; do
	region=$(echo "$line" | awk '{print $1}')
	prot=$(grep "${region::-2}" ~/Thesis/VanDePeerPSgenes/Shiu_Pipeline/BLAST_output/G500/PSregions_adjusted.q_G500_I500_pseudogenes.disable_count | awk '{print $1}' | cut -d'_' -f1)
	sed -i "0,/\t\t99999/{s/\t\t99999/\t${prot:1}|${prot:1:-2}\t99999/}" BLAST_RBH.txt
done

grep '99999' PLAST_RBH.txt | while read line; do
        region=$(echo "$line" | awk '{print $1}')
        prot=$(grep "${region::-2}" ~/Thesis/VanDePeerPSgenes/Shiu_Pipeline/PLAST_output/G500/query_cov_filtered.q_G500_I500_pseudogenes.disable_count | awk '{print $1}' | cut -d'_' -f1)
        sed -i "0,/\t\t99999/{s/\t\t99999/\t${prot:1}|${prot:1:-2}\t99999/}" PLAST_RBH.txt
done

grep '99999' Plant_RBH.txt | while read line; do
        region=$(echo "$line" | awk '{print $1}')
        prot=$(grep "${region::-2}" ~/Thesis/VanDePeerPSgenes/PlantPseudo/Results/final.edit.pg.xls | awk '{print $16}' | cut -d'_' -f1)
        sed -i "0,/\t\t99999/{s/\t\t99999/\t$prot|${prot::-2}\t99999/}" Plant_RBH.txt
done


#Some interanchor regions are not present (presumably filtered out, because they had an e-value that was too low). These too will also be assigned the protein with which they were discovered.
grep "^>" ~/Thesis/Transcriptome/Populus/Pseudogene_Sequences/temp.txt | while read line; do
	if ! grep -q "${line:1}" ~/Thesis/Transcriptome/Populus/RBH/BLAST_RBH.txt; then
		prot=$(grep "${line:1:-2}" ~/Thesis/VanDePeerPSgenes/Shiu_Pipeline/BLAST_output/G500/PSregions_adjusted.q_G500_I500_pseudogenes.disable_count | head -1 | awk '{print $1}' | cut -d'_' -f1)
		echo -e "${line:1}\t${prot:1}|${prot:1:-2}\t99999" >> BLAST_RBH.txt
	fi
done

grep "^>" ~/Thesis/Transcriptome/Populus/Pseudogene_Sequences/Extracted_PLAST_PSgene_Seqs.fa | while read line; do
	if ! grep -q "${line:1}" ~/Thesis/Transcriptome/Populus/RBH/PLAST_RBH.txt; then
                prot=$(grep "${line:1}" ~/Thesis/VanDePeerPSgenes/Shiu_Pipeline/PLAST_output/G500/query_cov_filtered.q_G500_I500_pseudogenes.disable_count | head -1 | awk '{print $1}' | cut -d'_' -f1)
                echo -e "${line:1}\t${prot:1}|${prot:1:-2}\t99999" >> PLAST_RBH.txt
        fi	
done	

grep "^>" ~/Thesis/Transcriptome/Populus/Pseudogene_Sequences/Extracted_Plant_PSgene_Seqs.fa | while read line; do
	if ! grep -q "${line:1}" ~/Thesis/Transcriptome/Populus/RBH/Plant_RBH.txt; then
		prot=$(grep "${line:1}" ~/Thesis/VanDePeerPSgenes/PlantPseudo/Results/final.edit.pg.xls | head -1 | awk '{print $16}' | cut -d'_' -f1)
		echo -e "${line:1}\t$prot|${prot::-2}\t99999" >> Plant_RBH.txt
	fi
done
