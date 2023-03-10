#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 10:12:24 2019

@author: beiki
"""
import collections
import pandas as pd
from multiprocessing import Pool
import time
import sys

""" NOTE: this script will filter out transcripts with exon length less
    tha 20nt. this script should be used for early cleaning of transcripts """
    
def nested_dict():
    return collections.defaultdict(nested_dict)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        

def get_DNAbased_bed(input,tr):
    minimum_exon_length=20 # this is to filter out transcripts with errornous exnons, to remove this filter simply set it to 0
    introns_bed={}
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    if df['ex_length'].min() >=minimum_exon_length:
        df=df.sort_values(by='ex_start', ascending=True)
        for i in range(df.shape[0]-1):
            introns_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0], df.iloc[i,2]+1, df.iloc[i+1,1]-1, tr+'intron' + str(i+1), 99, df.iloc[0,3], tr+'intron' + str(i+1))
    return(introns_bed)


def get_RNAbased_bed(input,tr):
    minimum_exon_length=20
    introns_bed={}
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    if df['ex_length'].min() >=minimum_exon_length and df.iloc[0,3]=="+":
        df=df.sort_values(by='ex_start', ascending=True)
        exonic_length=0
        for i in range(df.shape[0]-1):
            exonic_length=exonic_length + df.iloc[i,2] - df.iloc[i,1] + 1
            introns_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(tr, exonic_length-4, exonic_length+5, tr+'intron' + str(i+1), 99, df.iloc[0,3], tr+'intron' + str(i+1))# exonic_length-4, exonic_length+5 used to set a window of 10nt for intron location inorder to get better intersect with RNA-seq reads   
    elif df['ex_length'].min() >=minimum_exon_length and df.iloc[0,3]=="-":
        df=df.sort_values(by='ex_start', ascending=False)
        exonic_length=0
        for i in range(df.shape[0]-1):
            exonic_length=exonic_length + df.iloc[i,2] - df.iloc[i,1] + 1
            introns_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(tr, exonic_length-4, exonic_length+5, tr+'intron' + str(i+1), 99, df.iloc[0,3], tr+'intron' + str(i+1))# exonic_length-4, exonic_length+5 used to set a window of 10nt for intron location inorder to get better intersect with RNA-seq reads       
    return(introns_bed)

def paralle_chunk_input(chunk):
    bed=pd.DataFrame([])
    for tr in range(len(chunk)):
        info=get_DNAbased_bed(chunk[tr][1],chunk[tr][0])#RNAbased or DNAbased methods can be selected
        if len(info)>0:# for ttranscripts that has exons with lengh <20nt this item will be 0. please see get_RNAbased_bed or get_DNAbased_bed functions
            tr_bed=pd.DataFrame.from_dict(info,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'intron_start', 2: 'intron_end', 3: 'intron_id', 4: 'value', 5: 'strand', 6: 'intron_id'},axis=1)
            bed=bed.append(tr_bed)
    return(bed)


def ilen(it):
    """Yield number of items in generator."""
    return len(list(it))

input_gff=sys.argv[1]#your gff file, like "SubQ_Fat_fmlrc_proovread_cupcake_corrected_reads.collapsed.gff" or "B8198RT.collapsed.gff"
if input_gff.find('fmlrc') !=-1:
    tissue=input_gff.split("_fmlrc")[0] #tissue=input_gff.split(".collapsed")[0]
    outprefix='_fmlrc_proovread_cupcake_corrected_introns.bed'
elif input_gff.find('Trinity') !=-1:
    tissue=input_gff.split(".collapsed")[0]
    outprefix='_introns.bed'

tr_ex_info=[]
for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        info=raw[8].split("\"")
        tr_ex_info.append(info[3])

tr_ex_info_df=pd.DataFrame(tr_ex_info,columns =['tr_id'])
tr_ex_info_df=pd.DataFrame(tr_ex_info_df['tr_id'].value_counts().values, index=tr_ex_info_df['tr_id'].value_counts().index, columns=['Count']).reset_index().rename({'index':'tr_id'},axis=1)## counts the numnber of exons per transcript

tr_ex_info_dict={}
for i in range(tr_ex_info_df.shape[0]):
    tr_ex_info_dict[tr_ex_info_df.iloc[i,0]]=tr_ex_info_df.iloc[i,1]

ex_coords = nested_dict()

for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        tid=raw[-1].split('; ')[1].split()[1][1:-2]
        if tr_ex_info_dict[tid] > 1:# to get spliced transcripts
            exnumber=raw[3]
            ex_coords[tid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
     


items = list(ex_coords.items())
p = Pool() 
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, 100))
end = time.time()
print('total time (s)= ' + str(end-start)) 
introns_bed=pd.concat(processed_values)
introns_bed.to_csv(tissue + outprefix ,index = None, header=False,sep="\t")

""" single CPU check 
bed=pd.DataFrame([])
for i in range(len(items)):
    info=get_RNAbased_bed(items[i][1],items[i][0])
    tr_bed=pd.DataFrame.from_dict(info,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'intron_start', 2: 'intron_end', 3: 'intron_id', 4: 'value', 5: 'strand', 6: 'intron_id'},axis=1)
"""
    