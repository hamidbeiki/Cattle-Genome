#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:59:28 2020

@author: beiki
"""

import pandas as pd
import sys

df=pd.read_csv(sys.argv[1],sep=' ',names=["chr","ex_start","ex_end","tr_id","strand"])
transcript_ends=df.groupby(['tr_id'], sort=False)['ex_end'].max()
transcripts_start=df.groupby(['tr_id'], sort=False)['ex_start'].min()
transcripts_borders=pd.concat([transcripts_start,transcript_ends],axis=1)
transcripts_borders['length']=transcripts_borders['ex_end']-transcripts_borders['ex_start']
transcripts_borders=transcripts_borders.reset_index().rename({'ex_start':'tr_start', 'ex_end':'tr_end'},axis=1)
strand=df[['tr_id','chr','strand']]
strand=strand.drop_duplicates()
m=transcripts_borders.merge(strand)
#m.to_csv('combined_transcript_borders_with_intron',index = None, header=True,sep="\t")
m=m.set_index('tr_id')
d=m.to_dict(orient="index")

bed_dict={}
for k, value in d.items():
    s=value['strand']
    if s=="+":
        ns=value['tr_start']-500
        if ns <=0:
            ns=0
        ne=value['tr_start']+100
        bed_dict[k]=[value['chr'],ns,ne,k,'99',s,k]
    else:
        ns=value['tr_end']+500
        ne=value['tr_end']-100
        if ne<=0:
            ne=0
        bed_dict[k]=[value['chr'],ne,ns,k,'99',s,k]

bed_df=pd.DataFrame.from_dict(bed_dict,orient='index')
bed_df.to_csv('combined_transcript_start_sites.bed',index = None, header=False,sep="\t")     



"""
window=0### selected window
bed_dict={}
for k, value in d.items():
    s=value['strand']
    if s=="+":
        ns=value['tr_start']-window
        if ns <=0:
            ns=0
        ne=value['tr_start']+window
        bed_dict[k]=[value['chr'],ns,ne,k,'99',s,k]
    else:
        ns=value['tr_end']+window
        ne=value['tr_end']-window
        if ne<=0:
            ne=0
        bed_dict[k]=[value['chr'],ne,ns,k,'99',s,k]
"""