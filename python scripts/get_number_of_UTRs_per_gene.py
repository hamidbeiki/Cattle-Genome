#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 14:14:23 2021

@author: beiki
"""
import pandas as pd

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

def get_list(d):
    a=pd.DataFrame.to_numpy(d)
    l=a.tolist()
    return(l)
    
utr5=pd.read_csv('combined_5UTRs.bed',header=0,sep='\t')
pos=utr5.loc[utr5['strand']=="+"]
neg=utr5.loc[utr5['strand']=="-"]
start=pos.groupby(['tr_id'], sort=False)['exon_start'].min().reset_index()
end=pos.groupby(['tr_id'], sort=False)['exon_end'].max().reset_index()
pos=start.merge(end)
start=pos.groupby(['tr_id'], sort=False)['exon_start'].max().reset_index()
end=pos.groupby(['tr_id'], sort=False)['exon_end'].min().reset_index()
neg=start.merge(end)
df=pd.concat([pos,neg])
g=df['tr_id'].str.split('.',n = 2, expand = True)[[0,1]]
genes=(g[0] + "." + g[1]).to_frame().rename({0:'gene_id'},axis=1)
df['gene_id']=genes
genes=list(df['gene_id'].drop_duplicates())
   
utr5={}
groups = df.groupby('gene_id')
for gene, group in groups:
    d=group[['exon_start', 'exon_end']]
    l=get_list(d)
    res=merge_intervals(l)
    utr5[gene]=[len(res),res]    

utr3=pd.read_csv('combined_3UTRs.bed',header=0,sep='\t')
pos=utr3.loc[utr3['strand']=="+"]
neg=utr3.loc[utr3['strand']=="-"]
start=pos.groupby(['tr_id'], sort=False)['exon_start'].min().reset_index()
end=pos.groupby(['tr_id'], sort=False)['exon_end'].max().reset_index()
pos=start.merge(end)
start=pos.groupby(['tr_id'], sort=False)['exon_start'].max().reset_index()
end=pos.groupby(['tr_id'], sort=False)['exon_end'].min().reset_index()
neg=start.merge(end)
df=pd.concat([pos,neg])
g=df['tr_id'].str.split('.',n = 2, expand = True)[[0,1]]
genes=(g[0] + "." + g[1]).to_frame().rename({0:'gene_id'},axis=1)
df['gene_id']=genes
genes=list(df['gene_id'].drop_duplicates())

utr3={}
groups = df.groupby('gene_id')
for gene, group in groups:
    d=group[['exon_start', 'exon_end']]
    l=get_list(d)
    res=merge_intervals(l)
    utr3[gene]=[len(res),res]

utr5_df=pd.DataFrame.from_dict(utr5, orient='index').reset_index().rename({'index':'tr_id',0:'number_of_utr5',1:'UTR5'},axis=1)
utr3_df=pd.DataFrame.from_dict(utr3, orient='index').reset_index().rename({'index':'tr_id',0:'number_of_utr3',1:'UTR3'},axis=1)

utr5_df.to_csv("number_of_5UTRs_per_gene",index = False, header=True,sep="\t") 
utr3_df.to_csv("number_of_3UTRs_per_gene",index = False, header=True,sep="\t") 