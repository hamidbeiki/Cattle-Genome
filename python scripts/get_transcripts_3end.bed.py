#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:21:01 2021

@author: beiki
"""

import pandas as pd


df=pd.read_csv('out',sep='\t',names=["chr","ex_start","ex_end","tr_id","strand"])## awk '$3=="exon"{print $1, $4, $5,$12,$7}' combined.gff | sed 's/\"//g;s/\;//g; s/ /\t/g' >out
info=df[['tr_id','strand','chr']].drop_duplicates()
info['num']=99
pos=df.loc[df['strand']=="+"]
neg=df.loc[df['strand']=="-"]
pos_tr_ends=pos.groupby(['tr_id'], sort=False)['ex_end'].max().reset_index().rename({'ex_end':"3_end"},axis=1)
neg_tr_ends=neg.groupby(['tr_id'], sort=False)['ex_start'].min().reset_index().rename({'ex_start':"3_end"},axis=1)

pos_tr_ends['s']=pos_tr_ends['3_end']-100
pos_tr_ends['e']=pos_tr_ends['3_end']+1000
neg_tr_ends['s']=neg_tr_ends['3_end']-1000
neg_tr_ends['e']=neg_tr_ends['3_end']+100

all_3ends=pd.concat([pos_tr_ends,neg_tr_ends])
m=all_3ends.merge(info)

out=m[['chr','s','e','tr_id','num','strand','tr_id']]
out.to_csv('combined_transcript_3end.bed',index = None, header=False,sep="\t")