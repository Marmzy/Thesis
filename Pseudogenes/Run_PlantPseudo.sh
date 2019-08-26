#$ -S /bin/sh
#$ -pe serial 6
#$ -l h_vmem=12G

module load blast/x86_64/2.2.25+
module load exonerate/x86_64/2.2.0
module load fasta
module load MCScanX/x86_64/20141028
module load perl/x86_64/5.28.0
module load python2

##This script was made to run the PlantPseudo pipeline, which predicts plant pseudogenes.

##Trying to run the pipeline in 1 go produces many errors
#perl ./bin/pipeline.pl --scriptDir ./script --gff $HOME/Thesis/Annotations.gff3 --pep $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa --rawFa $HOME/Thesis/Genomes/Ptrichocarpa_210_v3.0.fa --repeatMaskedFa $HOME/Thesis/Genomes/Ptrichocarpa_210_v3.0.hardmasked.fa --eValueE -5 --idenThresh 95 --lenThresh 90 --proThresh 0.9 --qs 0 --mLenPse 50 --mLenIntron 50 --dirFile dirFile.txt --outDir ./Results

##Performing the pipeline step by step

#Step 1: Gff2Genepos.py
python2 ./script/Gff2Genepos.py $HOME/Thesis/Annotations.gff3 ./Results/gene.pos
cut -f1,2,3 ./Results/gene.pos > ./Results/masked.regions

#Step 2: fa-mask.py (not executed)
python2 ./script/fa-mask.py --region ./Results/masked.region $HOME/Thesis/Genomes/Ptrichocarpa_210_v3.0.hardmasked.fa > temp.rep.gene.masked.fa
mv temp.rep.gene.masked.fa ./Results

#Step 3: exonerate
exonerate --model protein2genome --showquerygff no --showtargetgff yes --maxintron 5000 --showvulgar yes --ryo \"%ti\\t%qi\\t%tS\\t%qS\\t%tl\\t%ql\\t%tab\\t%tae\\t%tal\\t%qab\\t%qae\\\\t%qal\\t%pi\\n\" $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa $HOME/Thesis/Genomes/Masked_seqs.fa > ./Results/exonerate2.out

#Step 4: ExtractExonerateOut.py
python2 ./script/ExtractExonerateOutv2.py ./Results/exonerate2.out ./Results/Best_Alignment2.out

#Step 5: ParseBlast.py
chmod +755 ./script/ParseBlast.py
python2 ./script/ParseBlast.py -f get_qualified4 -blast ./Results/Best_Alignment2.out -fasta $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -E 5 -I 40 -L 30 -P 0.05 -Q 1

#Step 6: Pseudo_step1.py
chmod +755 ./script/Pseudo_step1.py
python2 ./script/Pseudo_step1.py ./Results/Best_Alignment2.out?E5I40L30P5Q1.qlines 500

#Step 7: Pseudo_step2.py
chmod +755 ./script/Pseudo_step2.py
python2 ./script/Pseudo_step2.py ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE 500

#Step 8: Pseudo_step3.py
chmod +755 ./script/Pseudo_step3.py
python2 ./script/Pseudo_step3.py ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1

#Step 9: FastaManager.py	#NOTE: Instead of using the unmasked Salicoid WGD-derived interanchorpoint candidate regions fasta file, the soft-masked version is used, so "real" genes can be distinguished from pseudogenes
chmod +755 ./script/FastaManager.py
python2 ./script/FastaManager.py -f get_stretch4 -fasta ~/Thesis/Genomes/Masked_seqs.fa -coords ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord

#Step 10: BlastUtilityv2.py
python2 ./script/BlastUtilityv3.py -f batch_sw -p tfasty36 -g ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs -i $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -j ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.fa -bdir /software/shared/apps/x86_64/fasta/36.3.8d/bin -d 1

#Step 11: Pseudo_step4.py
chmod +755 ./script/Pseudo_step4.py
python2 ./script/Pseudo_step4.py ./script/blosum50.matrix ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.sw.out

#Step 12: CheckStrand.py
chmod +755 ./script/CheckStrand.py
python2 ./script/CheckStrand.py ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.sw.raw ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.Strd

#Step 13: PolyACheck.py
chmod +755 ./script/PolyACheck.py
python2 ./script/PolyACheck.py ~/Thesis/Genomes/Masked_seqs.fa ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.Strd ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.PolyA

#Step 14: CheckIntron.py
chmod +755 ./script/CheckIntron.py
python2 ./script/CheckIntron.py ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.fa ./script/env2.sh

#Step 15: SumTablev2.py
#NOTE: There were 2 irregular entries in "Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.ex.out", so I edited these manually so they resemble the other entries.
python2 ./script/SumTablev3.py ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.PolyA ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.sw.out.disable_count ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.ex.out ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.pairs.sw.out
mv pseudogene.phase2 ./Results

#Step 16: GetIntronfracv2.py
python2 ./script/GetIntronfracv3.py ./Results/gene.pos ./Results/pseudogene.phase2 $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa ./Results/Pg.xls

#Step 17: PgClassification.py
chmod +755 ./script/PgClassificationv2.py
python2 ./script/PgClassificationv2.py ./Results/Pg.xls ./Results/Pg.add.xls

#Step 18: Pggff.py,mcscanformatv2.py,Mcscan2Pglstv2.py
chmod +755 ./script/Pggff.py ./script/mcscanformatv2.py ./script/Mcscan2Pglstv2.py
./software/blast-2.2.25/bin/formatdb -i $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -p T
./software/blast-2.2.25/bin/blastall -i $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -d $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa -p blastp -e 1e-10 -b 10 -v 10 -m 8 -o blast -a 18
python2 ./script/Pggff.py ./Results/Pg.xls

python2 ./script/mcscanformatv3.py $HOME/Thesis/Proteome/Ptrichocarpa_adjusted.fa ./Results/gene.pos ./Results/spe2.gff
mv pg* ./Results
mv blast ./Results
cat ./Results/spe2.gff ./Results/pg.gff > ./Results/xyz.gff
cat ./Results/blast ./Results/pg.blast > ./Results/xyz.blast

#Step 19: MCScanX
MCScanX ./Results/xyz
python2 ./script/Mcscan2Pglstv3.py ./Results/xyz.collinearity ./Results/xyz.tandem ./Results/wgdlist ./Results/tandemlst

#Step 20:  FinalPglst.py
chmod +755 ./script/FinalPglst.py
python2 ./script/FinalPglst.py ./Results/wgdlist ./Results/tandemlst ./Results/Pg.add.xls ./Results/final.pg.xls

#Step 20.5: Remove "real" genes from final list
head -n1 ./Results/final.pg.xls > ./Results/final.edit.pg.xls

while read line; do
	pgID=$(echo "$line" | cut -f1)
	seq=$(grep -A1 "$pgID" ./Results/Best_Alignment2.out_E5I40L30P5Q1.qlines_G500.PE_I500.PS1.subj_coord.fa | tail -n1)
	
	l=$(echo "$seq" | grep -o "[a-z]" | wc -l)
	float=$(bc -l <<< $l/${#seq})

	if (( $(echo "$float < 0.3" | bc -l) )); then
		echo "$line" >> ./Results/final.edit.pg.xls
	fi
done < ./Results/final.pg.xls

#Step 21: DistanceComparev5.1.py (not executed)
