#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load gffcompare
module load python2
module load stringtie

#This module was made to run StringTie and prepare the files necessary for DESeq2 in R.

#Running StringTie on the STAR mapped reads
for i in {886..903}; do
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted -G Annotations_BLAST.gff3 -o SRR952${i}.BLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted -G Annotations_PLAST.gff3 -o SRR952${i}.PLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted -G Annotations_Plant.gff3 -o SRR952${i}.Plant.gtf -p 3
done

for i in {292..303}; do
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted -G Annotations_BLAST.gff3 -o SRR1121${i}.BLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted -G Annotations_PLAST.gff3 -o SRR1121${i}.PLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted -G Annotations_Plant.gff3 -o SRR1121${i}.Plant.gtf -p 3
done

for i in {2991..3014}; do
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted -G Annotations_BLAST.gff3 -o SRR347${i}.BLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted -G Annotations_PLAST.gff3 -o SRR347${i}.PLAST.gtf -p 3
	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted -G Annotations_Plant.gff3 -o SRR347${i}.Plant.gtf -p 3
done

#Move the gtf files to the corresponding folders (which also need to be made first)

ls -1 GSE49911/BLAST/*.gtf > GSE49911/BLAST/GSE49911_BLAST_list.txt
ls -1 GSE49911/PLAST/*.gtf > GSE49911/PLAST/GSE49911_PLAST_list.txt
ls -1 GSE49911/Plant/*.gtf > GSE49911/Plant/GSE49911_Plant_list.txt
ls -1 GSE54153/BLAST/*.gtf > GSE54153/BLAST/GSE54153_BLAST_list.txt
ls -1 GSE54153/PLAST/*.gtf > GSE54153/PLAST/GSE54153_PLAST_list.txt
ls -1 GSE54153/Plant/*.gtf > GSE54153/Plant/GSE54153_Plant_list.txt
ls -1 GSE81077/BLAST/*.gtf > GSE81077/BLAST/GSE81077_BLAST_list.txt
ls -1 GSE81077/PLAST/*.gtf > GSE81077/PLAST/GSE81077_PLAST_list.txt
ls -1 GSE81077/Plant/*.gtf > GSE81077/Plant/GSE81077_Plant_list.txt

stringtie --merge -p 3 -o GSE49911/BLAST/GSE49911_BLAST_merged.gtf -G Annotations_BLAST.gff3 GSE49911/BLAST/GSE49911_BLAST_list.txt
stringtie --merge -p 3 -o GSE49911/PLAST/GSE49911_PLAST_merged.gtf -G Annotations_PLAST.gff3 GSE49911/PLAST/GSE49911_PLAST_list.txt
stringtie --merge -p 3 -o GSE49911/Plant/GSE49911_Plant_merged.gtf -G Annotations_Plant.gff3 GSE49911/Plant/GSE49911_Plant_list.txt
stringtie --merge -p 3 -o GSE54153/BLAST/GSE54153_BLAST_merged.gtf -G Annotations_BLAST.gff3 GSE54153/BLAST/GSE54153_BLAST_list.txt
stringtie --merge -p 3 -o GSE54153/PLAST/GSE54153_PLAST_merged.gtf -G Annotations_PLAST.gff3 GSE54153/PLAST/GSE54153_PLAST_list.txt
stringtie --merge -p 3 -o GSE54153/Plant/GSE54153_Plant_merged.gtf -G Annotations_Plant.gff3 GSE54153/Plant/GSE54153_Plant_list.txt
stringtie --merge -p 3 -o GSE81077/BLAST/GSE81077_BLAST_merged.gtf -G Annotations_BLAST.gff3 GSE81077/BLAST/GSE81077_BLAST_list.txt
stringtie --merge -p 3 -o GSE81077/PLAST/GSE81077_PLAST_merged.gtf -G Annotations_PLAST.gff3 GSE81077/PLAST/GSE81077_PLAST_list.txt
stringtie --merge -p 3 -o GSE81077/Plant/GSE81077_Plant_merged.gtf -G Annotations_Plant.gff3 GSE81077/Plant/GSE81077_Plant_list.txt

gffcompare -r Annotations_BLAST.gff3 -G -o gffcompare_GSE49911_BLAST GSE49911/BLAST/GSE49911_BLAST_merged.gtf
gffcompare -r Annotations_PLAST.gff3 -G -o gffcompare_GSE49911_PLAST GSE49911/PLAST/GSE49911_PLAST_merged.gtf
gffcompare -r Annotations_Plant.gff3 -G -o gffcompare_GSE49911_Plant GSE49911/Plant/GSE49911_Plant_merged.gtf
gffcompare -r Annotations_BLAST.gff3 -G -o gffcompare_GSE54153_BLAST GSE54153/BLAST/GSE54153_BLAST_merged.gtf
gffcompare -r Annotations_PLAST.gff3 -G -o gffcompare_GSE54153_PLAST GSE54153/PLAST/GSE54153_PLAST_merged.gtf
gffcompare -r Annotations_Plant.gff3 -G -o gffcompare_GSE54153_Plant GSE54153/Plant/GSE54153_Plant_merged.gtf
gffcompare -r Annotations_BLAST.gff3 -G -o gffcompare_GSE81077_BLAST GSE81077/BLAST/GSE81077_BLAST_merged.gtf
gffcompare -r Annotations_PLAST.gff3 -G -o gffcompare_GSE81077_PLAST GSE81077/PLAST/GSE81077_PLAST_merged.gtf
gffcompare -r Annotations_Plant.gff3 -G -o gffcompare_GSE81077_Plant GSE81077/Plant/GSE81077_Plant_merged.gtf

for i in {886..903}; do
	stringtie -e -B -p 3 -G GSE49911/BLAST/GSE49911_BLAST_merged.gtf -o GSE49911/BLAST/Ballgown/SRR952${i}/SRR952${i}_BLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted
	stringtie -e -B -p 3 -G GSE49911/PLAST/GSE49911_PLAST_merged.gtf -o GSE49911/PLAST/Ballgown/SRR952${i}/SRR952${i}_PLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted
	stringtie -e -B -p 3 -G GSE49911/Plant/GSE49911_Plant_merged.gtf -o GSE49911/Plant/Ballgown/SRR952${i}/SRR952${i}_Plant.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Unique.SRR952${i}.sorted
done

for i in {292..303}; do
	stringtie -e -B -p 3 -G GSE54153/BLAST/GSE54153_BLAST_merged.gtf -o GSE54153/BLAST/Ballgown/SRR1121${i}/SRR1121${i}_BLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted
	stringtie -e -B -p 3 -G GSE54153/PLAST/GSE54153_PLAST_merged.gtf -o GSE54153/PLAST/Ballgown/SRR1121${i}/SRR1121${i}_PLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted
	stringtie -e -B -p 3 -G GSE54153/Plant/GSE54153_Plant_merged.gtf -o GSE54153/Plant/Ballgown/SRR1121${i}/SRR1121${i}_Plant.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE54153_Synthetic/Unique.SRR1121${i}.sorted
done

for i in {2991..3014}; do
	stringtie -e -B -p 3 -G GSE81077/BLAST/GSE81077_BLAST_merged.gtf -o GSE81077/BLAST/Ballgown/SRR347${i}/SRR347${i}_BLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted
	stringtie -e -B -p 3 -G GSE81077/PLAST/GSE81077_PLAST_merged.gtf -o GSE81077/PLAST/Ballgown/SRR347${i}/SRR347${i}_PLAST.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted
	stringtie -e -B -p 3 -G GSE81077/Plant/GSE81077_Plant_merged.gtf -o GSE81077/Plant/Ballgown/SRR347${i}/SRR347${i}_Plant.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE81077_Synthetic/Unique.SRR347${i}.sorted
done

#Download the prepDY.py file from: http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
#prepDE.py needs to be run from within the directory with the following commands and rename the output files to correctly reflect which dataset the reads come from and how the pseudogenes were predicted.
python2 prepDE.py -i Ballgown -p SRR952
python2 prepDE.py -i Ballgown -p SRR1121
python2 prepDE.py -i Ballgown -p SRR347


#Running StringTie on the Welch et al. filtered STAR mapped reads
for tool in BLAST PLAST Plant; do
	for i in {886..903}; do
        	stringtie /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Mappable.SRR952${i}.sorted -G Annotations_${tool}.gff3 -o SRR952${i}.${tool}.gtf -p 3
	done

	mv *.${tool}.gtf GSE49911/${tool}_Welch
	ls -1 GSE49911/${tool}_Welch/*.gtf > GSE49911/${tool}_Welch/GSE49911_${tool}_list.txt
	stringtie --merge -p 3 -o GSE49911/${tool}_Welch/GSE49911_${tool}_merged.gtf -G Annotations_${tool}.gff3 GSE49911/${tool}_Welch/GSE49911_${tool}_list.txt
	gffcompare -r Annotations_${tool}.gff3 -G -o gffcompare_GSE49911_${tool} GSE49911/${tool}_Welch/GSE49911_${tool}_merged.gtf
	mv gffcompare_GSE49911_${tool}.* GSE49911/${tool}_Welch/
	
	for i in {886..903}; do
        	stringtie -e -B -p 3 -G GSE49911/${tool}_Welch/GSE49911_${tool}_merged.gtf -o GSE49911/${tool}_Welch/Ballgown/SRR952${i}/SRR952${i}_${tool}.gtf /group/biocomp/users/cadav/Thesis/STAR/Mapped/GSE49911_Synthetic/Mappable.SRR952${i}.sorted
	done
done

#prepDE.py needs to be run from within the directory with the following commands and rename the output files to correctly reflect which dataset the reads come from and how the pseudogenes were predicted.
python2 prepDE.py -i Ballgown -p SRR952
