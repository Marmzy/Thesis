#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

#This script was made to extract all the TPM values for pseudogenes from the StringTie output.
#NOTE: It throws an error when trying to pass it through the server, so run locally!
for d in GSE49911 GSE54153 GSE81077; do
	for b in BLAST PLAST Plant; do
		if [[ "$d" == "GSE49911" ]]; then

			#Extract the TPM values and store them individual files
			for i in {886..903}; do
				paste <(cat $d/$b/Ballgown/SRR952${i}/SRR952${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f1) <(cat $d/$b/Ballgown/SRR952${i}/SRR952${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f2) > ${d}_${b}_SRR952${i}_Pseudogene_TPMs.txt
				sort ${d}_${b}_SRR952${i}_Pseudogene_TPMs.txt > ${i}_${b}.txt
			done
		
			#Combine the TPM values into a single files
			arr=(*$b.txt)
			file="${arr[0]}"
			
			for f in "${arr[@]:1}"; do
				paste "$file" <(awk '{print $2}' "$f") > _temp.txt && mv _temp.txt temp.txt
				file=temp.txt
			done

			#Add the "real" names of the pseudogenes to the end of the file
			arr=($(awk '{print $1}' temp.txt))

			for r in "${arr[@]}"; do
        			line=$(grep "${r};pseudo" Annotations_${b}.gff3 | cut -f9 | cut -d'=' -f3)
				echo "$line" >> column.txt
			done
			
			paste temp.txt column.txt > temp2.txt
			column -t temp2.txt > ${d}_${b}.table
			rm *.txt 

		
		elif [[ "$d" == "GSE54153" ]]; then

			#Extract the TPM values and store them individual files
			for i in {292..303}; do
				paste <(cat $d/$b/Ballgown/SRR1121${i}/SRR1121${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f1) <(cat $d/$b/Ballgown/SRR1121${i}/SRR1121${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f2) > ${d}_${b}_SRR1121${i}_Pseudogene_TPMs.txt
				sort ${d}_${b}_SRR1121${i}_Pseudogene_TPMs.txt > ${i}_${b}.txt
			done

			#Combine the TPM values into a single files
			arr=(*$b.txt)
                        file="${arr[1]}"

                        for f in "${arr[@]:1}"; do
                                paste "$file" <(awk '{print $2}' "$f") > _temp.txt && mv _temp.txt temp.txt
                                file=temp.txt
                        done

			#Add the "real" names of the pseudogenes to the end of the file
			arr=($(awk '{print $1}' temp.txt))

                        for r in "${arr[@]}"; do
                                line=$(grep "${r};pseudo" Annotations_${b}.gff3 | cut -f9 | cut -d'=' -f3)
				echo "$line" >> column.txt
                        done

                        paste temp.txt column.txt > temp2.txt
                        column -t temp2.txt > ${d}_${b}.table
                        rm *.txt

		else
			#Extract the TPM values and store them individual files
			for i in {2991..3014}; do
				paste <(cat $d/$b/Ballgown/SRR347${i}/SRR347${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f1) <(cat $d/$b/Ballgown/SRR347${i}/SRR347${i}_${b}.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 | cut -d'"' -f2) > ${d}_${b}_SRR347${i}_Pseudogene_TPMs.txt
				sort ${d}_${b}_SRR347${i}_Pseudogene_TPMs.txt > ${i}_${b}.txt
			done

			#Combine the TPM values into a single files
			arr=(*$b.txt)
                        file="${arr[1]}"

                        for f in "${arr[@]:1}"; do
                                paste "$file" <(awk '{print $2}' "$f") > _temp.txt && mv _temp.txt temp.txt
                                file=temp.txt
                        done

			#Add the "real" names of the pseudogenes to the end of the file
			arr=($(awk '{print $1}' temp.txt))

                        for r in "${arr[@]}"; do
                                line=$(grep "${r};pseudo" Annotations_${b}.gff3 | cut -f9 | cut -d'=' -f3)
				echo "$line" >> column.txt
                        done

                        paste temp.txt column.txt > temp2.txt
                        column -t temp2.txt > ${d}_${b}.table
                        rm *.txt
		
		fi
	done
done

#Add the sum of all TPMs to the end of the file
for f in GSE49911_BLAST GSE49911_PLAST GSE49911_Plant GSE54153_BLAST GSE54153_PLAST GSE54153_Plant GSE81077_BLAST GSE81077_PLAST GSE81077_Plant; do
        awk '{for(x=2;x<NF;x++) sum[$1]+=$x} END {for(i in sum) print i, sum[i]}' $f.table > test.txt
        sort test.txt > test2.txt
        paste $f.table <(cut -d' ' -f2 test2.txt) > test3.txt
        column -t test3.txt > ${f}_test.table
	mv ${f}_test.table ${f}.table
done

rm *.txt

#Create new files based on the total TPM count
for data in GSE49911 GSE54153 GSE81077; do
        for tool in BLAST PLAST Plant; do
                if [ "$data" == "GSE49911" ]; then
                        while read line; do
                                tpm=$(echo "$line" | awk '{print $21}')

                                if (( $(echo "$tpm >= 50" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_50.table
                                elif (( $(echo "$tpm >= 10" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_10.table
                                elif (( $(echo "$tpm >= 1" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_1.table
                                fi

                        done < ${data}_${tool}.table
                elif [ "$data" == "GSE54153" ]; then
                         while read line; do
                                tpm=$(echo "$line" | awk '{print $15}')

                                if (( $(echo "$tpm >= 50" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_50.table
                                elif (( $(echo "$tpm >= 10" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_10.table
                                elif (( $(echo "$tpm >= 1" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_1.table
                                fi

                        done < ${data}_${tool}.table
                else
                        while read line; do
                                tpm=$(echo "$line" | awk '{print $27}')

                                if (( $(echo "$tpm >= 50" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_50.table
                                elif (( $(echo "$tpm >= 10" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_10.table
                                elif (( $(echo "$tpm >= 1" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_1.table
                                fi
                        done < ${data}_${tool}.table
                fi
        done
done

########################################################################################################################################################

#This part of the script was made to extract all the TPM values for pseudogenes that were filtered usign the Welch et al. method.
for d in GSE49911; do
        for b in BLAST PLAST Plant; do
                if [[ "$d" == "GSE49911" ]]; then

                        #Extract the TPM values and store them individual files
                        for i in {886..903}; do
                                paste <(cat $d/${b}_Welch/Ballgown/SRR952${i}/SRR952${i}_BLAST.gtf | grep "exon:" | awk -F'\t' '{if($3=="transcript") print $0}' | cut -f9 | cut -d';' -f2,5 | cut -d'"' -f2,4 -d'"' -f2) > ${d}_${b}_SRR952${i}_Pseudogene_TPMs.txt
                                sort ${d}_${b}_SRR952${i}_Pseudogene_TPMs.txt > ${i}_${b}.txt
                        done

                        #Combine the TPM values into a single files
                        arr=(*$b.txt)
                        file="${arr[0]}"

                        for f in "${arr[@]:1}"; do
                                paste "$file" <(awk '{print $2}' "$f") > _temp.txt && mv _temp.txt temp.txt
                                file=temp.txt
                        done

                        #Add the "real" names of the pseudogenes to the end of the file
                        arr=($(awk '{print $1}' temp.txt))

                        for r in "${arr[@]}"; do
                                line=$(grep "${r};pseudo" Annotations_${b}.gff3 | cut -f9 | cut -d'=' -f3)
                                echo "$line" >> column.txt
                        done

                        paste temp.txt column.txt > temp2.txt
                        column -t temp2.txt > ${d}_${b}.table
                        rm *.txt
                fi
        done
done

#Add the sum of all TPMs to the end of the file
for f in GSE49911_BLAST GSE49911_PLAST GSE49911_Plant; do
        awk '{for(x=2;x<NF;x++) sum[$1]+=$x} END {for(i in sum) print i, sum[i]}' $f.table > test.txt
        sort test.txt > test2.txt
        paste $f.table <(cut -d' ' -f2 test2.txt) > test3.txt
        column -t test3.txt > ${f}_test.table
        mv ${f}_test.table ${f}.table
done

rm *.txt

#Create new files based on the total TPM count
for data in GSE49911; do
        for tool in BLAST PLAST Plant; do
                if [ "$data" == "GSE49911" ]; then
                        while read line; do
                                tpm=$(echo "$line" | awk '{print $21}')

                                if (( $(echo "$tpm >= 50" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_50.table
                                elif (( $(echo "$tpm >= 10" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_10.table
                                elif (( $(echo "$tpm >= 1" | bc -l) )); then
                                        echo "$line" >> ${data}_${tool}_1.table
                                fi

                        done < ${data}_${tool}.table
                fi
        done
done
