#$ -S /bin/sh
#$ -pe serial 2
#$ -l h_vmem=4G

module load plast
module load python

#This script was created as the PLAST output is not compatible with Shiu's pipeline's filtering scripts. This script
#therefore replaces it. The following filtering steps are used:
#This script was created as Shiu's pipeline pseudo_wrap.py script doesn't work as intended with the BLAST+ output, because it was designed for legacy BLAST.
#pseudo_wrap.py performs all the steps needed to detect pseudogenes.
#Here I made a script that basically does the same thing, but it explicetely calls the different scripts for all steps.

plast -a 2 -F T -force-query-order 1000 -p tplastn -d ~/Thesis/Genomes/Masked_seqs.fa -i ~/Thesis/Proteome/Ptrichocarpa_adjusted.fa -o Shiu_e3_outfmt1_adjusted.tab -outfmt 1 -verbose

# -E 5    = e-value < 1e-5
# -I 40   = identity > 40%
# -L 30   = match length > 30aa
# -P 0.05 = match length > 5% of the query sequence

#Adjusting the PLAST output file, as it seems to miss a few characters
awk -v PRE='=v3.0' '{if (substr($1,length($1)-4,5)!="=v3.0") {$1=$1PRE; print} else {$1; print}}' OFS="\t" Shiu_e3_outfmt1.tab > Shiu_e3_outfmt1_adjusted.tab

#Parsing the PLAST output to ensure an e-value of 1e-5 or lower
while read line; do
        array=($line)
        evalue=$(echo "${array[10]}")

        if [ "$evalue" != "0.001" ]; then
                evalue2=$(echo "${array[10]}" | cut -f2  -d "e")
                evalue1=$(echo "${array[10]}" | cut -f1  -d "e")

                if [ "$evalue2" -lt "-05" ]; then
                        echo "$line" >> evalue_filtered.txt
                elif [ "$evalue2" -eq "-05" ] && [ $evalue1 -eq "1" ]; then
                        echo "$line" >> evalue_filtered.txt
                fi
        fi

done < Shiu_e3_outfmt1_adjusted.tab

#Parsing the PLAST output to ensure an identity of 40% or greater
awk '$3>="40.00"' evalue_filtered.txt > identity_filtered.txt

#Parsing the PLAST output to ensure a length of 30 amino acids or longer of the query matching the target
awk '$4>="30"' identity_filtered.txt > match_length_filtered.txt

#Parsing the PLAST output to ensure that more than 5% of the query sequence is covered
#The output is written to query_cov_filtered.txt
python3 Query_seq_percent.py id.txt $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa match_length_filtered.txt

Gsize=500
Isize=1000

#Set max gap size for linking pseudoexons
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step2e.py query_cov_filtered.q $Gsize		

#Set intron length
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step3b.py query_cov_filtered.q_G${Gsize}.PE $Isize	

#Get a pairs file and a subj coordinate file for the phase 1 pseudogenes
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step3.5.py query_cov_filtered.q_G${Gsize}.PE_I${Isize}.PS1

#Get the putative pseudogene regions out of the genome sequence using the coordinate file just generated
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/FastaManager.py -f get_stretch4 -fasta $HOME/Thesis/Genomes/Masked_seqs.fa -coords query_cov_filtered.q_G${Gsize}.PE_I${Isize}.PS1.subj_coord

#Run tfasty
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/BlastUtility.py -f batch_sw -p tfasty36 -g query_cov_filtered.q_G${Gsize}.PE_I${Isize}.PS1.pairs -i $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -j query_cov_filtered.q_G${Gsize}.PE_I${Isize}.PS1.subj_coord.fa -bdir $HOME/Thesis/fasta-36.3.8e/bin -d 1

#Get the final output (.disable_count) file which contains the pseudogenes and info on frameshift and stop
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step6.py $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/blosum50.matrix query_cov_filtered.q_G${Gsize}.PE_I${Isize}.PS1.pairs.sw.out
