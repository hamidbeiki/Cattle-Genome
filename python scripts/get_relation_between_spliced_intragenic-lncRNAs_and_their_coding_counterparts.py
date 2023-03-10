#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 09:38:32 2021

@author: beiki
"""

import pandas as pd
import sys
import os
import re
import functools
import operator

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])

protein_coding=pd.read_csv('protein_coding_transcripts',names=['tr_id'])
protein_coding=set(list(protein_coding['tr_id']))
intra_lncRNAs=pd.read_csv('long_intragenic_lncRNAs',names=['tr_id'])
intra_lncRNAs=set(list(intra_lncRNAs['tr_id']))
spliced_intra_lncRNAs=pd.read_csv('long_intragenic_lncRNAs_spliced',names=['tr_id'])

files=[]
directory='/work/ABG/Hamid/Bos-taurus/FAANG_Aug/combined_tissues/redundanct_removed/noncoding-isoforms-of-coding-transcripts'
with os.scandir(directory) as entries:# or sys.argv[1]
    for entry in entries:
        if re.search('strict.ioe',entry.name):
            files.append(entry.name)
AS_event_summary={}
AS_tr_summary={}
AS_trs={}# intragenic lncRNAs resulted from AS events on refrence protein-coding transcript
for file in files: 
    AS_type=file.split('_')[5]
    out_l=[] 
    df=pd.read_csv(file,header=0,sep='\t')
    l=df[['alternative_transcripts','total_transcripts']].values.tolist()
    count=0
    for sub in l:
        info=[i.split(',', 1) for i in sub]
        ref=info[1]
        alternate=info[0]
        ref=set(list(dict.fromkeys(ref)))
        alternate=set(list(dict.fromkeys(alternate)))
        nc=list(alternate & intra_lncRNAs)
        p=list(ref & protein_coding)
        if len(nc)>0 and len(p)>0:
            count=count+1
            out_l.append(nc)
    AS_event_summary[AS_type]=count
    out_l=flatten_nested_list(out_l)
    out_l=list(dict.fromkeys(out_l))
    AS_trs[AS_type]=out_l
    AS_tr_summary[AS_type]=len(out_l)

AS_summary_df=pd.DataFrame.from_dict(AS_event_summary,orient='index').reset_index().rename({'index':'AS_event',0:'number_of_events'},axis=1)

trs=[]
for k,l in AS_trs.items():
    trs.append(l)

trs=flatten_nested_list(trs)
trs=list(dict.fromkeys(trs))
trs=pd.DataFrame(trs).rename({0:'tr_id'},axis=1)
USTs=spliced_intra_lncRNAs.loc[~ spliced_intra_lncRNAs.tr_id.isin(trs.tr_id)]#unique splice site transcripts
USTs.to_csv('intragenic_lncRNAs_resulted_from_UST_events_on_their_coding_counterparts',index = None, header=False,sep="\t")

for k,l in AS_trs.items():
    df=pd.DataFrame(l).rename({0:'tr_id'},axis=1)
    name="intragenic_lncRNAs_resulted_from_" + k + "_events_on_their_coding_counterparts"
    df.to_csv(name,index = None, header=False,sep="\t") 
    
"""        
AS_summary={}
for file in files: 
    AS_type=file.split('_')[5]    
    df=pd.read_csv(file,header=0,sep='\t')
    l=df[['alternative_transcripts','total_transcripts']].values.tolist()
    count=0
    for sub in l:
        info=[i.split(',', 1) for i in sub]
        ref=info[1]
        alternate=info[0]
        ref=set(list(dict.fromkeys(ref)))
        alternate=set(list(dict.fromkeys(alternate)))
        p=list(alternate & protein_coding)
        nc=list(ref & intra_lncRNAs)
        if len(nc)>0 and len(p)>0:
            count=count+1
    AS_summary[AS_type]=count
"""
