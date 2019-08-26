#!/usr/bin/env python3

import sys

#Checking if a file passed through the command line
def Arguments_check():
    if len(sys.argv) <= 4:
        print("Usage:", sys.argv[0], "<subj_coord.fa file>, <.disable_count file>, <max gap size for linking pseudoexons> and <intron length>")
    else:
        GetPseudogenes(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

#Separating the actual pseudogenes from the genes
def GetPseudogenes(file, dc_file, Gsize, Isize):
    
    genes = []
    pseudogenes = []
    
    #Reading the subj_coord.fa file
    with open (file) as f:
        content = f.readlines()
        
        for line in content:
            if ">" in line:
                header = line
            else:
                header += line
                
                #If the putative pseudogene region contains only lowercase letters, it means a protein got paired against a tandem
                #duplicated gene
                if line.islower():
                    genes.append(header.split("\n")[0])
                
                #If the putative pseudogene region doesn't contain a lowercase letter, it means a protein got paired against a region
                #of DNA, presumably a pseudogene    
                elif line.isupper():
                    pseudogenes.append(header.split("\n")[0])
                
                #If the putative pseudogene region contains both uppercase and lowercase letters, it means a protein got paired against
                #a part of a tandem duplicated gene. If the sequence overlaps with no more than 30% against a tandem duplicate, it will
                #be considered a potential pseudogene
                else:
                    upper_count = sum(map(str.isupper, header))
                    
                    if (upper_count/len(header))*100 >= 70:
                        pseudogenes.append(header.split("\n")[0])
                    else:
                        genes.append(header.split("\n")[0])
    
    #Reading the .discable_count file                    
    with open (dc_file) as f3:
        dc_content = f3.readlines()
    
    entry = ""
    filename = dc_file.split(".")[0]
    firstTime = True
    
    #Looping over .discable_count file content
    for l in dc_content:
        
        if firstTime:
            firstTime = False
            h = ">" + l.split(" ")[1]
            entry += l
        
        #When a new entry is encountered, the old one is written to the appropriate file
        elif '#' in l:
            h = ">" + entry.split(" ")[1]
            
            #Write "gene" entry output to file    
            if h in genes:
                with open(filename + ".q_G" + Gsize + "_I" + Isize + "_genes.disable_count", "a") as f:
                    f.write(entry)
            
            #Write "pseudogene" entry output to file        
            elif h in pseudogenes:
                with open(filename + ".q_G" + Gsize + "_I" + Isize + "_pseudogenes.disable_count", "a") as f:
                    f.write(entry)
            
            entry = l
        
        #Add the remaining lines of the entry to the variable 'entry'    
        else:
            entry += l

    h = ">" + entry.split(" ")[1]
    
    #Write the last entry also to a file!
    if h in genes:
         with open(filename + ".q_G" + Gsize + "_I" + Isize + "_genes.disable_count", "a") as f:
             f.write(entry)
             
    elif h in pseudogenes:
         with open(filename + ".q_G" + Gsize + "_I" + Isize + "_pseudogenes.disable_count", "a") as f:
                     f.write(entry)
    

if __name__ == '__main__':
    Arguments_check()
