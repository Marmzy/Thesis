#$ -S /bin/sh
#$ -pe serial 2
#$ -l h_vmem=4G

module load python

#This script was made to filter all 'real' genes from our pseudogene results file. Actual genes are still found as the
#regions between anchorpoints still sometimes contained tandem duplicated genes. These can be recognised, because the
#genome was soft-masked, so sequences contain non-capital letters are actually genes instead of pseudogenes.

python3 PSgene_vs_gene.py PSregions_adjusted.q_G500.PE_I500.PS1.subj_coord.fa PSregions_adjusted.q_G500.PE_I500.PS1.pairs.sw.out.disable_count 500 500
