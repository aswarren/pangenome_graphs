#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import glob
import pandas as pd

##script input:
## 1. table of filenames to organism strings
## 2. n plain text RAST/SEED tables

##1 following format e.g.
#6666666.53128.txt	Brucella sp BO2 174
#6666666.53117.txt	Brucella ovis ATTCC

#morph rast data by adding filename to columns
#adding lower base coordinate and stripping away uncessary columns
#returns pandas table
def morph_rast(input_file):
    file_path, file_extension= os.path.splitext(input_file)
    file_name=os.path.basename(input_file).replace(file_extension,'')
    initial_table=pd.read_table(input_file)
    result=pd.DataFrame({'org_id':[file_name]*len(initial_table)}) #give file name as organism ID
    result['contig_id']=initial_table['contig_id']
    result['locus_id']=initial_table['feature_id']
    result['start']=initial_table[['start','stop']].min(axis=1) #take lower base coordinate
    result['fam_id']=initial_table['figfam']#figfam ID
    result['fam_description']=initial_table['function']
    result=result[pd.notnull(result['fam_id'])] #remove things without figfams
    result=result.sort(columns=['contig_id','start'], ascending=[1,1])
    return result
    
    
    
    

def main(init_args):
    if(len(init_args)<2):
        sys.stderr.write("Usage: transform_rast.py output_file.txt [rast_table1.txt ...]\n")
        sys.exit()
    counter =0
    files=[]
    for i in init_args[1:]:
        files=files+glob.glob(i)
    if len(files)==0:
        print "no files"
        return 0
    for i in files:
        if counter == 0:
            result=morph_rast(i)
        else:
            result.append(morph_rast(i))
        counter+=1
    result.to_csv(init_args[0], sep="\t", index=False)

if __name__ == "__main__":
    main(sys.argv[1:])
