#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:59:28 2020

@author: beiki
"""

import pandas as pd

df=pd.read_csv('out',sep=' ',names=["chr","tr_start","tr_end","gene_id","strand"])
gene_ends=df.groupby(['gene_id'], sort=False)['tr_end'].max()
genes_start=df.groupby(['gene_id'], sort=False)['tr_start'].min()
genes_borders=pd.concat([genes_start,gene_ends],axis=1)
genes_borders['length']=genes_borders['tr_end']-genes_borders['tr_start']
genes_borders=genes_borders.reset_index().rename({'tr_start':'gene_start', 'tr_end':'gene_end'},axis=1)
strand=df[['gene_id','chr','strand']]
strand=strand.drop_duplicates()
m=genes_borders.merge(strand)
#m.to_csv('combined_gene_borders_with_intron',index = None, header=True,sep="\t")
m=m.set_index('gene_id')
d=m.to_dict(orient="index")

bed_dict={}
for k, value in d.items():
    s=value['strand']
    if s=="+":
        ns=value['gene_start']-500
        if ns <=0:
            ns=0
        ne=value['gene_start']+100
        bed_dict[k]=[value['chr'],ns,ne,k,'99',s,k]
    else:
        ns=value['gene_end']+500
        ne=value['gene_end']-100
        if ne<=0:
            ne=0
        bed_dict[k]=[value['chr'],ne,ns,k,'99',s,k]

bed_df=pd.DataFrame.from_dict(bed_dict,orient='index')
bed_df.to_csv('combined_gene_promoters.bed',index = None, header=False,sep="\t")       