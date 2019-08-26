#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys

#Checking if a file passed through the command line
def Arguments_check():
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "<tblastn Output File> and <blastx Output File>")
    else:
        RBH(sys.argv[1], sys.argv[2])

#Get the Reciprocal Best Hits for our predicted pseudogenes and the proteins of Populus trichocarpa. This, to get a better estimate of the pseudogene parent genes.
def RBH(tblastn_file, blastx_file):
    
    #Reading the tblastn_file and blastx_file as a dataframe
    tdfn = pd.read_csv(tblastn_file, sep='\t', header=None, usecols=range(2))
    dfx = pd.read_csv(blastx_file, sep='\t', header=None, usecols=range(2))

    #Making 2 dictionaries containing all query IDs as keys and target IDs as values.
    tblastn = {}
    blastx = {}

    for q, t in zip(tdfn[0].values, tdfn[1].values):
        if q not in tblastn.keys():
            tblastn[q] = [t]
        else:
            tblastn[q].append(t)

    for q2, t2 in zip(dfx[0].values, dfx[1].values):
        if q2 not in blastx.keys():
            blastx[q2] = [t2]
        else:
            blastx[q2].append(t2)

    #Loop over all interanchor pair regions
    for key in blastx.keys():
        best_score = 99999
        best_prot = ''
        
        #Loop over all significant blastx protein hits for said interanchor pair region
        for pos, value in enumerate(blastx[key]):

            #If the blastx protein hit exists, check if it has the interanchor pair region in the tblastn result.
            #If it has the region, calculate its score
            if value in tblastn.keys():
                if key in tblastn[value]:
                    idx = tblastn[value].index(key)
                    score = pos + idx

                    #Keep the region-prot pair with the lowest RBH score
                    if score < best_score:
                        best_score = score
                        best_prot = value

        #Write the interanchor pair region and the protein with the best score to an output file
        with open("PLAST_RBH.txt", "a") as f:
            f.write(key + "\t" + best_prot + "\t" + str(best_score) + "\n")
 

if __name__ == '__main__':
    Arguments_check()
