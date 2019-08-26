#!/usr/bin/env python3

from Bio import SeqIO
import sys

#Checking if a file passed through the command line
def Arguments_check():
    if len(sys.argv) <= 3:
        print("Usage:", sys.argv[0], "<id file>, <PLAST query sequence file> and <PLAST output file>")
    else:
        entry_dict = SeqLength(sys.argv[1], sys.argv[2])
        OverlapCalc(entry_dict, sys.argv[3])
        
#Calculate the length of the query sequences by id        
def SeqLength(id, query):
    
    #Create a dictionary of all query file entries (proteins)
    record_dict = SeqIO.to_dict(SeqIO.parse(query, "fasta"))
    
    #Calculating the length of the entries for each entry in the file
    with open(id) as f:
        content = f.readlines()
        entry_dict = {}
        
        #Looping over each entry in the id file and altering the entry, so it matches the record_dict key 
        for entry in content:
            entry = entry.rstrip('\n').split('annot')[0]
            entry += 'annot-version=v3.0'
            
            if entry not in entry_dict:
                entry_dict[entry] = len(record_dict[entry].seq)
    
    return entry_dict


#Calculate the query sequence coverage and filter appropiately
def OverlapCalc(dict, file):
    
    with open(file) as f:
        content = f.readlines()
        
        #For each line in the PLAST output file, calculate if the line needs to be kept (>= 5% query sequence) or not
        for line in content:
            entry = line.split('\t')[0].rstrip('\n').split('annot')[0]
            entry += 'annot-version=v3.0'
            query_length = int(dict[entry])
            M = int(line.split('\t')[3])
            
            #If the PLAST entry passes the filter, it is written to a new file
            if M/query_length >= 0.05:
                with open("query_cov_filtered.qlines", "a") as out_file:
                    out_file.write(line)
        
        
if __name__ == '__main__':
    Arguments_check()           
