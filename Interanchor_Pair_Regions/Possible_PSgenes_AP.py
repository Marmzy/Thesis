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
        InterAnchorRegion(df_new.sort_values(by=[4]), annotationfile, i_adhore_file, gf_file, counter, filename)


def InterAnchorRegion(df, annotationfile, i_adhore_file, gf_file, counter, filename):
    
    #Calculating the relative anchorpoint location difference between consecutive rows
    df[7] = df[3].diff() #Always negative (unless inversions)
    df[8] = df[4].diff() #Can be negative or positive

    #Checking how the genes of chromosomeB are compared to chromsomeA (increasing or decreasing)
    count = 0
    for i in df[7]:
        if i > 0:
            count += 1
        elif i < 0:
            count -= 1
    
    #Creating a list with all anchorpoints
    anchorpoints = [x for x in df[1]]
    [anchorpoints.append(y) for y in df[2]]

    #Here we determine how the first few relative anchorpoints are represented in the dataframe
    if df[7].iloc[1] > 0:
        isPositive = True
    else:
        isPositive = False
    
    rows_list = []
    
    #We loop over the dataframe and determine how the anchorpoints on chromosomeB are situated compared to chromosomeA (normal vs inversed)
    for iter, row in df[1:].iterrows():
        counter += 1
        
        #print(row, df[7][iter])
        
        if not isPositive and df[7][iter] < 0:
            isPositive = False

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            if abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)
                
        elif not isPositive and df[7][iter] > 0:
            isPositive = True

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            if abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)
                
        elif isPositive and df[7][iter] > 0:
            isPositive = True

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            if abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
        
                if info is not None:
                    rows_list.append(info)
                
        elif isPositive and df[7][iter] < 0:
            isPositive = False

            #If there's exactly 1 gene between 2 anchorpoints -> get info from chromosomeB
            #It's a potential region if there are no genes present between the anchorpoints on chromosomeB
            if abs(row[7]) == 2 and abs(row[8]) == 1:
                ifgenes = True
                info = InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorpoints, filename)
            
                if info is not None:
                    rows_list.append(info)

    return rows_list
          

def InterAnchorGenes(df, annotationfile, row, isPositive, count, ifgenes, i_adhore_file, gf_file, anchorlist, filename):
    
    #The 2_ptr_ptr_2_data_filtered.txt file contains the relative positions (as determined by i-ADHoRe) of the genes
    #Reading the annotation file as a dataframe
    dfA = pd.read_csv(annotationfile, sep='\t', header=None, usecols=[0,3,4,8])
    rownums = list(df.index.values)
    
    #Getting the genes from the annotation 
    chr = row[5]
    rownum = df.loc[df[1]==row[1]].index[0]
    previous_rownum = rownums[rownums.index(rownum)-1]
    
    #If the column with the relative i-ADHoRe anchorpoints positions indicates that the anchorpoint has increased to its neighbour
    if isPositive:
        anchorB1 = dfA[dfA[8].str.contains(df.loc[previous_rownum][1])]
        anchorB2 = dfA[dfA[8].str.contains(row[1])]
        
    #If the column with the relative i-ADHoRe anchorpoints positions indicates that the anchorpoint has decreased to its neighbour
    else:
        anchorB1 = dfA[dfA[8].str.contains(row[1])]
        anchorB2 = dfA[dfA[8].str.contains(df.loc[previous_rownum][1])]
    
    #Getting the start and end coordinates of the interanchor region
    start = int(anchorB1[4])
    end = int(anchorB2[3])
    
    isAP = False
    
    #We loop through the gff file and try to classify the gene that lies between the anchorpoints on chromosomeB
    #(for ChromosomeA, see file: "Possible_PSgenes_ChrA")
    for pos, (i,j) in enumerate(zip(dfA[3], dfA[4])):
        if chr == str(dfA[0][pos]) and start < i and j < end:

            #Making sure that the gene between the 2 anchor points is not a tandem duplicate
            tandem = Adhore(dfA[8][pos].split(";")[0][3:], i_adhore_file)
                
            #Checking if the gene that lies between the 2 anchorpoints is another anchorpoint or not
            if tandem != -1:
                if str(dfA[8][pos].split(";")[0][3:]) in anchorlist:
                    isAP = True
                    AP_gene = dfA[8][pos].split(";")[0][3:]
                    #print(start, end, AP_gene)
                    #print(row)
        
    #If there is exactly 1 gene between the anchorpoints on chromsomeA, we need to get the coordinates of the anchorpoints
    #on chromosomeB
    chr = row[6]
    anchorB1 = dfA[dfA[8].str.contains(row[2])]
    anchorB2 = dfA[dfA[8].str.contains(df.loc[previous_rownum][2])]
    end = int(anchorB1[3])
    start = int(anchorB2[4])
    
    #Making sure that inversions aren't picked up
    if start < end:
        string = chr + "\t" + str(start) + "\t" + str(end)
            
        #If the gene between the 2 anchorpoints on chromosomeB is not "special"
        if isAP:
            string = string + "\t" + AP_gene + "\n"
            
            with open("PSregions_" + filename + "_AnchorPoint.txt", "a") as f:
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