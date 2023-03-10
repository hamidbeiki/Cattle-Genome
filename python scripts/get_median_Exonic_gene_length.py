#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:35:25 2019

@author: beiki
"""
## This program works with python3
import pandas as pd
import sys
from multiprocessing import Pool

#intervals = [[3,9],[2,6],[8,10],[15,18]]
#intervals = [[1,2],[4,5],[7,8],[3,4],[7,10],[1,2],[4,5],[6,10]]

def merge_intervals(arr):
    # Sorting based on the increasing order  
    # of the start intervals
    arr.sort(key = lambda x: x[0]) 
    # array to hold the merged intervals 
    m = []
    s = -10000
    max = -100000
    for i in range(len(arr)):
        a = arr[i]
        if a[0] > max:
            if i != 0:
                m.append([s,max])
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]
    #'max' value gives the last point of
    # that particular interval
    # 's' gives the starting point of that interval
    # 'm' array contains the list of all merged intervals
    if max != -100000 and [s, max] not in m:
        m.append([s, max])
    return(m)

def get_geneLength(gene_structure):
    df=pd.DataFrame(gene_structure)
    df['length']=df[1]-df[0]
    return(df['length'].sum())

''' quicker 
def get_geneLength(gene_structure):
    df=pd.DataFrame(gene_structure,columns=["start","end"])
    return(df.apply(lambda row: row.end - row.start, axis=1).sum())
'''

''' slow 
def get_geneLength(gene_structure):
    exon_length=pd.DataFrame([])
    for exon in range(len(gene_structure)):
        exon_length=exon_length.append(pd.DataFrame({'length':gene_structure[exon][1]-gene_structure[exon][0]},index=[0]), ignore_index=True)
    return(exon_length['length'].sum())
'''
def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def paralle_length_calculator(chunk):
    gene_length=pd.DataFrame([])
    for i in range(len(chunk)):
        gene=chunk[i][0]
        gene_structure=merge_intervals(chunk[i][1][0])
        length=get_geneLength(gene_structure)
        gene_length=gene_length.append(pd.DataFrame({'gene':gene,'exonic_length':length},index=[0]), ignore_index=True)
    return(gene_length)        

def parallle_format_df(chunk):
    """Yield successive n-sized chunks from grouped_df."""
    df=pd.DataFrame([])
    for i in range(len(chunk)):
        df=df.append(chunk[i][1])
    a=df.iloc[:,1:3].apply(list,axis=1)
    df['list']=a
    grouped=df.groupby('gene').agg({'list':lambda x: list(x)})
    return(grouped)
    
''' run parallel '''
## Step1 get formatted df
df=pd.read_csv(sys.argv[1],sep=' ',names=["gene","ex_start","ex_end"])
grouped_df=df.groupby('gene') 
p = Pool() 
dflist = []
for group in grouped_df:
    dflist.append(group)
processed_values= p.map(parallle_format_df, get_chunks(dflist, 100))
grouped=pd.concat(processed_values)
## Step2 get exonic gene length
p = Pool() 
items = list(grouped.iterrows())
processed_values= p.map(paralle_length_calculator, get_chunks(items, 100))
gene_length=pd.concat(processed_values)

print(gene_length['exonic_length'].median())

''' single cpu run 
df=pd.read_csv(sys.argv[1],sep=' ',names=["gene","ex_start","ex_end"])
a=df.iloc[:,1:3].apply(list,axis=1)
df['list']=a
grouped=df.groupby('gene').agg({'list':lambda x: list(x)})
gene_length=pd.DataFrame([])
for gene in grouped.itertuples():
    gene_structure=merge_intervals(gene[1])
    length=get_geneLength(gene_structure)
    gene_length=gene_length.append(pd.DataFrame({'gene':gene[0],'exonic_length':length},index=[0]), ignore_index=True) 

print(gene_length['exonic_length'].median())   


#for gene in grouped.itertuples():
#    if gene[0]=="PB.10002":
#        intervals=gene[1]
#        gene_structure=merge_intervals(gene[1])  
'''

''' to get median gene length included intronic region, use:
df=pd.read_csv('out2',sep=' ',names=["gene","tr_start","tr_end"])      
gene_ends=df.groupby(['gene'], sort=False)['tr_end'].max()
genes_start=df.groupby(['gene'], sort=False)['tr_start'].min()
genes_borders=pd.concat([genes_start,gene_ends],axis=1)
genes_borders['length']=genes_borders['tr_end']-genes_borders['tr_start']
genes_borders.to_csv(r'genes_borders',index = True, header=True,sep="\t")
print(genes_borders['length'].median())
'''
""" to get average exonic length of transcripts
df=pd.read_csv('delete_combined_info',sep='\t',names=["transcript","ex_start","ex_end"])
df['length']=df['ex_end']-df['ex_start']
length=df.groupby(['transcript'], sort=False)['length'].sum().reset_index()
length['length'].median()
"""