#$ -S /bin/sh
#$ -pe serial 3
#$ -l h_vmem=6G

module load samtools/x86_64/1.6

#This script was made to get uniquely mappable reads as determined by Welch et al., 2015, instead of simply having unique reads.
#The sam files were altered individually, so as to not create any space or memory issues and to be able to modify the script if changes are necessary.


cat Ptr_67.bedGraph | awk '$4==1 {print $0}' > MapUniq67.bedGraph
mv MapUniq67.bedGraph Mappable
cd Mappable
cat MapUniq67.bedGraph | awk '{print>$1}'
cd ..

##Freeing up some space...
gzip -9 Aligned.*


#Sort the sam file
cat Unique.SRR952897.sam | sort -k3,3 -k4,4n > temp.sam

#Initialising variables that will be used later
prev_chr=""
prev_pos=0

#Loop over all the reads (sam file)
while read line; do
        chr=$(echo "$line" | awk '{print $3}')
        pos=$(echo "$line" | awk '{print $4}')

        #If the read hails from a Synthetic region it needs it's name shortened to ensure the appropriate mappability file can be opened
        if [[ $chr == Synt* ]]; then
                chr=$(echo "$chr" | cut -d':' -f1)
        fi

        #Create a temporary file from the mappability file, containing a select amount of entries, unless the chromosome and length of the starting position of the read don't change
        if ! { [ "${#pos}" -eq "$prev_pos" ] && [ "$prev_chr" == "$chr" ]; }; then
                awk -v num="$pos" 'length($3)==length(num) {print $0}' ./Mappable/$chr > temp.txt
        fi

        #Loop over the temporary file until the end position is larger than the starting position of the read
        while read line2; do
                stop=$(echo "$line2" | cut -f3)
                stop=$((stop++))

                if [ "$stop" -ge "$pos" ]; then
                        start=$(echo "$line2" | cut -f2)
                        start=$((start++))
                        idx=$(awk -v line="$stop" '$3~line {print NR; exit}' temp.txt)

                        #If the start and end position of the mappability track are larger than the starting position of the read: break, because the read isn't uniquely mappable
                        #Else, add the read to a new file
                        #Also, update variables and remove lines from the temporary file that have already been run over
                        if [ "$start" -ge "$pos" ]; then
                                prev_chr=$chr
                                prev_pos=${#pos}
                                tail -n +"$idx" temp.txt > temp.out && mv temp.out temp.txt
                                break
                        else
                                prev_chr=$chr
                                prev_pos=${#pos}
                                tail -n +"$idx" temp.txt > temp.out && mv temp.out temp.txt
                                echo "$line" >> Mappable.SRR952897.sam
                                break
                        fi
                fi
        done < temp.txt
done < temp.sam

#Removing the original file
#rm Unique.SRR952886.sam

#Convert the SAM file to a BAM file (need to manually add the removed SAM header sequences)
for i in {886..903}; do
	len=$(zgrep "^@" Aligned.SRR952${i}.trimmed.sam.gz | wc -l)
	zcat Aligned.SRR952${i}.trimmed.sam.gz | head -n $len > Header.SRR952${i}.sam
	cat Header.SRR952${i}.sam Mappable.SRR952${i}.sam > temp.sam
	samtools view -q 255 -Sb temp.sam > temp.bam
	samtools sort -o Mappable.SRR952${i}.sorted temp.bam
done
