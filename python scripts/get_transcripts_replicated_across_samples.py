#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:35:20 2020

@author: beiki
"""

import pandas as pd
import sys

input_file=sys.argv[1]#i.e. gffcmp.tracking
sample_file=sys.argv[2]#i.e tmp_files (file include list of samples gff files)
tracks=pd.read_csv(input_file,sep="\t")
samples=pd.read_csv(sample_file,sep="\t",names=["sample"])
samples['id'] = samples['sample'].str.split('\_final', 1).str[0]
info=tracks.iloc[0:,4:]
results=[]
for i in range(info.shape[0]):
    count=0
    info2=info.iloc[i,0:]
    transcripts=[]
    for j in range(len(info2)):
        if info2[j].find('-')==-1:
            tr=info2[j].split('|')[1]
            sample=samples.iloc[j,1]
            out=[tr,sample]
            transcripts.append(out)
            count=count+1
    results.append([transcripts,count])


df=pd.DataFrame(results)
replicated=df.loc[df[1]>1]
valid_trs={}
for i in range(samples.shape[0]):
    valid_trs[samples.iloc[i,1]]=[]
    

for i in range(replicated.shape[0]):
    my_list=replicated.iloc[i,0]
    for j in range(len(my_list)):
        info=my_list[j]
        valid_trs[info[1]].append(info[0])

for i in range(len(samples)):
    sample=samples.iloc[i,1]
    my_list=valid_trs[sample]
    valid_trs[sample]=pd.DataFrame(my_list).rename({0:'tr_id'},axis=1)
        
for i in range(len(samples)):
    sample=samples.iloc[i,1]
    df=valid_trs[sample]
    df.to_csv(sample + '_replicated_transcripts',index = None, header=False,sep="\t")