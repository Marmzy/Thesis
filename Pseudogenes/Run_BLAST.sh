#$ -S /bin/sh
#$ -pe serial 4
#$ -l h_vmem=8G

module load blast+
module load python

#This script was made to run BLAST (specifically tblastn) on the psb cluster.
#This script was made to parse the BLAST+ results file. It is the first step for finding pseudogenes following the
#Shiu's pipeline method. The other steps (files) can be found in the different "GXXX" folders, according to the parameters.
#This script was created as Shiu's pipeline pseudo_wrap.py script doesn't work as intended with the BLAST+ output, because it was designed for legacy BLAST.
#pseudo_wrap.py performs all the steps needed to detect pseudogenes.
#Here I made a script that basically does the same thing, but it explicetely calls the different scripts for all steps.

# -E 5    = e-value < 1e-5
# -I 40   = identity > 40%
# -L 30   = match length > 30aa
# -P 0.05 = match length > 5% of the query sequence

tblastn -db PSregions_DB -query Ptrichocarpa_adjusted.fa -out Ptr_PSregions.blasted -evalue 0.0001 -outfmt 7 -matrix BLOSUM50 -lcase_masking
python ../_pseudo_pkg/ParseBlast.py -f get_qualified4 -blast Ptr_PSregions_adjusted.blasted -fasta Ptrichocarpa_adjusted.fa -E 5 -I 40 -L 30 -P 0.05 -Q 1

Gsize=500
Isize=0

#Set max gap size for linking pseudoexons
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step2e.py ../PSregions_adjusted.q $Gsize		

#Set intron length
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step3b.py PSregions_adjusted.q_G${Gsize}.PE $Isize	

#Get a pairs file and a subj coordinate file for the phase 1 pseudogenes
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step3.5.py PSregions_adjusted.q_G${Gsize}.PE_I${Isize}.PS1

#Get the putative pseudogene regions out of the genome sequence using the coordinate file just generated
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/FastaManager.py -f get_stretch4 -fasta ../Masked_seqs.fa -coords PSregions_adjusted.q_G${Gsize}.PE_I${Isize}.PS1.subj_coord

#Run tfasty
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/BlastUtility.py -f batch_sw -p tfasty36 -g PSregions_adjusted.q_G${Gsize}.PE_I${Isize}.PS1.pairs -i ../Ptrichocarpa_adjusted.fa -j PSregions_adjusted.q_G${Gsize}.PE_I${Isize}.PS1.subj_coord.fa -bdir $HOME/Thesis/fasta-36.3.8e/bin -d 1

#Get the final output (.disable_count) file which contains the pseudogenes and info on frameshift and stop
python $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/script_step6.py $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/_pseudo_pkg/blosum50.matrix PSregions_adjusted.q_G${Gsize}.PE_I${Isize}.PS1.pairs.sw.out
