#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

#This script was made to alter the .gff3 file, so that it also contains "pseudogene" entries, so StringTie will be able work with it.
#NOTE: This script was modified to make Annotations.gff3 files for both the BLAST and PLAST output. Modify it yourself to get the respective files.

awk -F'\t' '{ if($3=="gene" || $3=="exon") print $0 }' $HOME/Thesis/Annotations.gff3 > test.gff3
mv test.gff3 Annotations.gff3

grep "^#" $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/BLAST_output/G500/PSregions_adjusted.q_G500_I500_pseudogenes.disable_count > test.txt

count=0

while read line; do
	begin=$(echo "$line" | cut -d' ' -f2 | cut -d':' -f2 | cut -d'-' -f1)
	num1=$(echo "$line" | cut -d' ' -f2 | cut -d'|' -f2 | cut -d'-' -f1)
	num2=$(echo "$line" | cut -d' ' -f2 | cut -d'|' -f2 | cut -d'-' -f2)
	start=$(($begin+$num1-1))
	stop=$(($begin+$num2-1))

	((count++))
	parent=$(echo "$line" | cut -d'_' -f1)
	parent="${parent}:exon:${count}"
	pseudo=$(echo "$line" | cut -d' ' -f2)
	chr=$(echo "$line" | cut -d' ' -f2 | cut -d':' -f1)
	
	if [[ $chr != *"scaffold"* ]]; then
		chr="Chr"$chr
	fi

	num3=$(echo "$line" | cut -d' ' -f3 | cut -d':' -f2 | cut -d'-' -f1)
        num4=$(echo "$line" | cut -d' ' -f3 | cut -d'-' -f3)
	
	if (( num3 < num4 )); then
		strand="+"
	else
		strand="-"
	fi

	echo -e "$chr\tJGI v3.1\tpseudogene\t$start\t$stop\t.\t$strand\t.\tID=${parent:1};pseudogene_id=$pseudo" >> appendix.txt

done < test.txt

rm test.txt
cp Annotations.gff3 test.txt
sort -u appendix.txt > apx.txt
cat test.txt apx.txt > merged.gff3
sort -k1,1V -k5,5n -k4,4r  merged.gff3 > Annotations_BLAST.gff3
rm apx.txt appendix.txt merged.gff3 test.txt Annotations.gff3

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

awk -F'\t' '{ if($3=="gene" || $3=="exon") print $0 }' $HOME/Thesis/Annotations.gff3 > test.gff3
mv test.gff3 Annotations.gff3

grep "^#" $HOME/Thesis/VanDePeerPSgenes/Shiu_Pipeline/PLAST_output/G500/query_cov_filtered.q_G500_I500_pseudogenes.disable_count > test.txt

count=0

while read line; do
        begin=$(echo "$line" | cut -d' ' -f2 | cut -d':' -f2 | cut -d'-' -f1)
        num1=$(echo "$line" | cut -d' ' -f2 | cut -d'|' -f2 | cut -d'-' -f1)
        num2=$(echo "$line" | cut -d' ' -f2 | cut -d'|' -f2 | cut -d'-' -f2)
        start=$(($begin+$num1-1))
        stop=$(($begin+$num2-1))

        ((count++))
        parent=$(echo "$line" | cut -d'_' -f1)
        parent="${parent}:exon:${count}"
        pseudo=$(echo "$line" | cut -d' ' -f2)
        chr=$(echo "$line" | cut -d' ' -f2 | cut -d':' -f1)

        if [[ $chr != *"scaffold"* ]]; then
                chr="Chr"$chr
        fi

        num3=$(echo "$line" | cut -d' ' -f3 | cut -d':' -f2 | cut -d'-' -f1)
        num4=$(echo "$line" | cut -d' ' -f3 | cut -d'-' -f3)

        if (( num3 < num4 )); then
                strand="+"
        else
                strand="-"
        fi

        echo -e "$chr\tJGI v3.1\tpseudogene\t$start\t$stop\t.\t$strand\t.\tID=${parent:1};pseudogene_id=$pseudo" >> appendix.txt

done < test.txt

rm test.txt
cp Annotations.gff3 test.txt
sort -u appendix.txt > apx.txt
cat test.txt apx.txt > merged.gff3
sort -k1,1V -k5,5n -k4,4r  merged.gff3 > Annotations_PLAST.gff3
rm apx.txt appendix.txt merged.gff3 test.txt Annotations.gff3

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

awk -F'\t' '{ if($3=="gene" || $3=="exon") print $0 }' $HOME/Thesis/Annotations.gff3 > test.gff3
mv test.gff3 Annotations.gff3

grep "^[0-9]" ~/Thesis/VanDePeerPSgenes/PlantPseudo/Results/final.edit.pg.xls > test.txt

count=0

while read line; do
        chr=$(echo "$line" | awk '{print $1}' | cut -d':' -f1)

        begin=$(echo "$line" | cut -d':' -f2 | cut -d'-' -f1)
        num1=$(echo "$line" | awk '{print $3}')
        num2=$(echo "$line" | awk '{print $4}')
        start=$(($begin+$num1-1))
        stop=$(($begin+$num2-1))

        ((count++))

        strand=$(echo "$line" | awk '{print $5}')
        parent=$(echo "$line" | awk '{print $16}' | cut -d'_' -f1)
        parent="${parent}:exon:${count}"
        pseudo=$(echo "$line" | awk '{print $1}')

        echo -e "Chr$chr\tJGI v3.1\tpseudogene\t$start\t$stop\t.\t$strand\t.\tID=${parent};pseudogene_id=$pseudo" >> appendix.txt

done < test.txt

rm test.txt
cp Annotations.gff3 test.txt
sort -u appendix.txt > apx.txt
cat test.txt apx.txt > merged.gff3
sort -k1,1V -k5,5n -k4,4r  merged.gff3 > Annotations_Plant.gff3
rm apx.txt appendix.txt merged.gff3 test.txt Annotations.gff3
