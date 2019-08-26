#!/usr/bin/env python3

import os
import pandas as pd
import re
import sys

#Checking if a file passed through the command line
def Arguments_check():
    if len(sys.argv) <= 4:
        print("Usage:", sys.argv[0], "<Collinear regions file>, <General feature format file>, <i-ADHoRe output file> and <Gene Family file>")
    else:
        regions = GetMultiplicon(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


def GetMultiplicon(file, annotationfile, i_adhore_file, gf_file):
    
    #Reading the file as a dataframe
    df = pd.read_csv(file, sep='\t', header=None, usecols=range(7))
    
    counter = 0
    
    filename = file.split(".")[0]
    fullpath = filename
    filename = os.path.relpath(fullpath, '/home/cadav/Thesis/MichielData/PSregions/')
    
    #Getting the multiplicon dataframes
    for multiplicon in set(df[0]):
        df_new = df.loc[df[0] == multiplicon]
        InterAnchorRegion(df_new, annotationfile, i_adhore_file, gf_file, counter, filename)


def InterAnchorRegion(df, annotationfile, i_adhore_file, gf_file, counter, filename):
    
    #Calculating the relative anchorpoint location difference between consecutive rows
    df[7] = df[3].diff() #Always negative (unless inversions)
    df[8] = df[4].diff() #Can be negative or positive

    #Checking how the genes of chromosomeB are compared to chromsomeA (increasing or decreasing)
    count = 0
    for i in df[8]:
        if i > 0:
            count += 1
        elif i < 0:
            count -= 1
    
    #Creating a list with all anchorpoints
    anchorpoints = [x for x in df[1]]
    [anchorpoints.append(y) for y in df[2]]

    #The dataframe always starts with the inversed genes, if there is an inversed sequenced
    #Here we determine how the first few relative anchorpoints are represented in the dataframe
    if df[8].iloc[1] > 0:
        isPositive = True
    else:
        isPositive = False
    
    rows_list = []

    #We loop over the dataframe and determine how the anchorpoints on chromosomeB are situated compared to chromosomeA (normal vs inversed)
    for iter, row in df[1:].iterrows():
        counter += 1
        
        if not isPositive and df[8][iter] < 0:
            isPositive = False
                
            #If there's no gene between 2 anchorpoints -> get info from chromosomeA
            #It's a potential region if there is exactly 1 gene present between the anchorpoints on chromosomeB
            if abs(row[7]) == 1 and abs(row[8]) == 2:
                ifgenes = False
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            elif abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)
                
        elif not isPositive and df[8][iter] > 0:
            isPositive = True
                
            #If there's no gene between 2 anchorpoints -> get info from chromosomeA
            #It's a potential region if there is exactly 1 gene present between the anchorpoints on chromosomeB
            if abs(row[7]) == 1 and abs(row[8]) == 2:
                ifgenes = False
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            elif abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)
                
        elif isPositive and df[8][iter] > 0:
            isPositive = True
                
            #If there's no gene between 2 anchorpoints -> get info from chromosomeA
            #It's a potential region if there is exactly 1 gene present between the anchorpoints on chromosomeB
            if abs(row[7]) == 1 and abs(row[8]) == 2:
                ifgenes = False
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            elif abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
        
                if info is not None:
                    rows_list.append(info)
                
        elif isPositive and df[8][iter] < 0:
            isPositive = False
                
            #If there's no gene between 2 anchorpoints -> get info from chromosomeA
            #It's a potential region if there is exactly 1 gene present between the anchorpoints on chromosomeB
            if abs(row[7]) == 1 and abs(row[8]) == 2:
                ifgenes = False
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            elif abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

    return rows_list
          

def InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorlist, filename):
    
    #The 2_ptr_ptr_2_data_filtered.txt file contains the relative positions (as determined by i-ADHoRe) of the genes
    #Reading the annotation file as a dataframe
    dfA = pd.read_csv(annotationfile, sep='\t', header=None, usecols=[0,3,4,8])
    dfGF = pd.read_csv(gf_file, sep='\t', header=None, usecols=[0,2])
    
    #Getting the genes from the annotation 
    if not ifgenes:
        chr = row[6]
    else:
        chr = row[5]
    
    rownum = df.loc[df[2]==row[2]].index[0]
    
    if rownum == 0:
        df = df.drop(df.loc[df[2]==row[2]].index[0])
        rownum = df.loc[df[2]==row[2]].index[0]
    
    #If the column with the relative i-ADHoRe anchorpoints positions indicates that the anchorpoint has increased to its neighbour
    if isPositive:
        anchorB1 = dfA[dfA[8].str.contains(df.loc[rownum-1][2])]
        anchorB2 = dfA[dfA[8].str.contains(row[2])]
        
    #If the column with the relative i-ADHoRe anchorpoints positions indicates that the anchorpoint has decreased to its neighbour
    else:
        anchorB1 = dfA[dfA[8].str.contains(row[2])]
        anchorB2 = dfA[dfA[8].str.contains(df.loc[rownum-1][2])]
    
    #Getting the start and end coordinates of the interanchor region
    start = int(anchorB1[4])
    end = int(anchorB2[3])
    
    isNormal = False
    isAP = False
    isHomolog = False
    isNested = False
    isOverlap = False
    isSpecial = False
    
    if ifgenes:
        anchorB1 = dfA[dfA[8].str.contains(row[1])]
        anchorB2 = dfA[dfA[8].str.contains(df.loc[rownum-1][1])]
        start = int(anchorB1[4])
        end = int(anchorB2[3])
    
    #We loop through the gff file and try to classify the gene that lies between the anchorpoints on chromosomeB
    #(for ChromosomeA, see file: "Possible_PSgenes_ChrA")
    for pos, (i,j) in enumerate(zip(dfA[3], dfA[4])):
        if chr == str(dfA[0][pos]) and start < i and j < end:

            #Making sure that the gene between the 2 anchor points is not a tandem duplicate
            tandem = Adhore(dfA[8][pos].split(";")[0][3:], i_adhore_file)
                
            #Finding out to which gene families the genes belong
            anchorB1num = dfGF.loc[dfGF[2] == str(anchorB1[8]).split(";")[0][3:].split("=")[1]].index[0]
            anchorB2num = dfGF.loc[dfGF[2] == str(anchorB2[8]).split(";")[0][3:].split("=")[1]].index[0]
                
            GFanchor1 = dfGF.loc[anchorB1num][0]
            GFanchor2 = dfGF.loc[anchorB2num][0]
            gene_family = dfGF.loc[dfGF.loc[dfGF[2] == dfA[8][pos].split(";")[0][3:]].index[0]][0]
                
            #Checking if the gene that lies between the 2 anchorpoints is another anchorpoint or not
            if tandem != -1:
                if str(dfA[8][pos].split(";")[0][3:]) in anchorlist:
                    isAP = True
                    anchorpoint_gene = dfA[8][pos].split(";")[0][3:]
                elif GFanchor1 == gene_family or GFanchor2 == gene_family:
                    isHomolog = True
                    homologous_gene = dfA[8][pos].split(";")[0][3:]
                else:
                    isNormal = True
                    normal_gene = dfA[8][pos].split(";")[0][3:]
                            
        elif chr == str(dfA[0][pos]):
                
            #Detect if there are genes lying within other genes
            if (int(anchorB1[3]) < i and j < int(anchorB1[4])) or (int(anchorB2[3]) < i and j < int(anchorB2[4])):  
                tandem = Adhore(dfA[8][pos].split(";")[0][3:], i_adhore_file)
                
                #Getting the name of the nested gene
                if tandem != -1:
                    isNested = True
                    nested_gene = dfA[8][pos].split(";")[0][3:]
                
            #Detect if there are overlapping genes, which thus only partly fulfill the requirements
            elif chr == str(dfA[0][pos]) and int(anchorB1[3]) < i and j < int(anchorB2[4]):
                tandem = Adhore(dfA[8][pos].split(";")[0][3:], i_adhore_file)
                    
                #Getting the name of the overlapping gene
                if tandem != -1:
                    isOverlap = True
                    overlapping_gene = dfA[8][pos].split(";")[0][3:]
                
            #If none of the previous categories sufficied, the gene must be special
            #Eg: An anchorpoint is nested in the gene between the 2 anchorpoints (1969 - Potri.006G043166)
            else:
                isSpecial = True
    
    #If there are no genes between the anchorpoints on chromsomeA, we need to get those coordinates
    if not ifgenes:
        chr = row[5]
        anchorB2 = dfA[dfA[8].str.contains(df.loc[rownum-1][1])]
        anchorB1 = dfA[dfA[8].str.contains(row[1])]
        start = int(anchorB1[4])
        end = int(anchorB2[3])
        
    #If there is exactly 1 gene between the anchorpoints on chromsomeA, we need to get the coordinates of the anchorpoints
    #on chromosomeB
    else:
        chr = row[6]
        anchorB1 = dfA[dfA[8].str.contains(row[2])]
        anchorB2 = dfA[dfA[8].str.contains(df.loc[rownum-1][2])]
        start = int(anchorB1[4])
        end = int(anchorB2[3])
    
    #Making sure that inversions aren't picked up
    if start < end:
        string = chr + "\t" + str(start) + "\t" + str(end)
            
        #If the gene between the 2 anchorpoints on chromosomeB is not "special"
        if isNormal:
            string = string + "\t" + normal_gene + "\n"

            with open("PSregions_" + filename + "_Normal.txt", "a") as f:
               f.write(string)
            
        #If there was only gene between the anchorpoints of chromosomeB and it was a nested gene, it is written to another file (special case)
        elif isNested:
            string = string + "\t" + nested_gene + "\n"

            with open("PSregions_" + filename + "_Nested.txt", "a") as f:
               f.write(string)
            
        #If there was some overlap with the gene between the anchorpoints of chromosomeB and the anchorpoints of chromosomeB itself, it is written to another file (special case)
        elif isOverlap:
            string = string + "\t" + overlapping_gene + "\n"

            with open("PSregions_" + filename + "_Overlap.txt", "a") as f:
                f.write(string)
            
        #If the gene between the 2 anchorpoints on chromosomeB was in fact a tandem duplicated gene, it is written to another file (special case)
        elif isHomolog:
            string = string + "\t" + homologous_gene + "\n"

            with open("PSregions_" + filename + "_Homolog.txt", "a") as f:
                f.write(string)
        
        #If the gene between the 2 anchorpoints on chromosomeB was in fact another anchorpoint, it is written to another file (special case)
        elif isAP:
            string = string + "\t" + anchorpoint_gene + "\n"

            with open("PSregions_" + filename + "_AnchorPoint.txt", "a") as f:
               f.write(string)
            
        #If the gene can't be placed in one of the previous 4 categories
        elif isSpecial:
            string = string + "\n"

            with open("PSregions_" + filename + "_Special.txt", "a") as f:
                f.write(string)


def Adhore(gene, i_adhore_file, remapped=True):
    
    #Opening the i-ADHoRe file to check if the gene is an anchorpoint
    f = open(i_adhore_file, "r")
    file_list = f.readlines()
    
    #We limit the search of the anchorpoints by the name
    r = re.compile(gene)
    newlist = filter(r.match, file_list)
    
    return int(list(newlist)[0].split("\t")[9])


if __name__ == '__main__':
    Arguments_check()
