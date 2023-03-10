#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:21:01 2021

@author: beiki
"""

import pandas as pd

df=pd.read_csv('out',sep='\t',names=["chr","tr_start","tr_end","gene_id","strand"])## awk '$3=="transcript"{print $1, $4, $5,$10,$7}' combined.gff | sed 's/\"//g;s/\;//g; s/ /\t/g' >out
info=df[['gene_id','strand','chr']].drop_duplicates()
info['num']=99
pos=df.loc[df['strand']=="+"]
neg=df.loc[df['strand']=="-"]
pos_gene_ends=pos.groupby(['gene_id'], sort=False)['tr_end'].max().reset_index().rename({'tr_end':"3_end"},axis=1)
neg_gene_ends=neg.groupby(['gene_id'], sort=False)['tr_start'].min().reset_index().rename({'tr_start':"3_end"},axis=1)
pos_gene_ends['s']=pos_gene_ends['3_end']-100
pos_gene_ends['e']=pos_gene_ends['3_end']+1000
neg_gene_ends['s']=neg_gene_ends['3_end']-1000
neg_gene_ends['e']=neg_gene_ends['3_end']+100

all_3ends=pd.concat([pos_gene_ends,neg_gene_ends])
m=all_3ends.merge(info)
out=m[['chr','s','e','gene_id','num','strand','gene_id']]
out.to_csv('combined_gene_3end.bed',index = None, header=False,sep="\t")