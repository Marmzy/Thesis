#!/usr/bin/env python
import sys,re

IN=open(sys.argv[1],'r')
OUT=open(sys.argv[2],'w')

#Query id,Subject id,% identity,alignment length,mismatches,gap openings,q. start,q. end,s. start,s. end,e-value,bit score
#%ti\\t%qi\\t%tS\\t%qS\\t%tl\\t%ql\\t%tab\\t%tae\\t%tal\\t%qab\\t%qae\\t%qal\\t%pi\\n\
#Chr01   Potri.004G209000        +       .       50495391        790     14228391        14228704        313     289     392     103     39.81

for eachline in IN:
	if re.match(r'".+', eachline):
		split=eachline.rstrip().split("\t")

		if split[2] == "+":
			start=int(split[6])
			end=int(split[7])
		else:
			start=int(split[7])
			end=int(split[6])
		length=abs(start-end)

		if length<5000:
			eval=(1/float(split[11]))*0.000001
			#query_id=split[1].split("_")[0]
			query_length=split[10].split("\\t")[1]
			query_end=split[10].split("\\t")[0]
			OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(split[1],split[0][1:],split[11],query_length,"2","2",split[9],query_end,start,end,eval,"200"))


IN.close()	
OUT.close()
