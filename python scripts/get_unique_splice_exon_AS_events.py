#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:28:10 2021

@author: beiki
"""
import pandas as pd

def get_list(d):
    a=d.values# instead of "pd.DataFrame.to_numpy(d)" in python 2.7
    l=a.tolist()
    return(l) 

def get_uniq_splice_exons(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]
        if a_right < b_right:
            i += 1
        else:
            j += 1
        if a_right >= b_left and b_right >= a_left and a_right != b_right and a_left != b_left:
            ranges.append([a_right,a_left]) 
    return(ranges)

EIs_df=pd.read_csv('trs_with_exons_covered_by_other_trs_introns',names=['transcript_id'])    
EIs=set(list(EIs_df['transcript_id']))
df=pd.read_csv('out',sep=' ',names=['chr','ex_start','ex_end','strand','gene_id','transcript_id']) #awk '$3=="exon"{print $1, $4, $5, $7,$10, $12}' combined_final.gff | sed 's/\"//g;s/\;//g' >out    
df2=df[['gene_id','transcript_id']].drop_duplicates()
groups=df2.groupby(['gene_id']) 
gene_info={}
for gene, group in groups:
    d=group[['transcript_id']]
    l=list(d['transcript_id'])
    gene_info[gene]=l

groups=df.groupby(['transcript_id'])
tr_info={}
for tr, group in groups:
    d=group[['ex_start','ex_end']]
    l=get_list(d)
    tr_info[tr]=l



USEs_list=[]
USTs_list=[]
for gene in gene_info:
    trs=gene_info[gene]
    s=set(trs) - EIs
    trs=list(s)
    for i in range(len(trs)-1):
        for j in range(i+1,len(trs)):
            trA=trs[i]
            trB=trs[j]
            trA_info=tr_info[trA]
            trB_info=tr_info[trB]
            if len(trA_info)>2 and len(trB_info)>2:
                A=trA_info[1:-1]
                B=trB_info[1:-1]
                USEs=get_uniq_splice_exons(A,B)
                USEs_=get_uniq_splice_exons(B,A)
                if len(USEs)>0:
                    USEs_list.append([trA,trB,len(USEs),USEs])
                if len(USEs)==len(A):
                    USTs_list.append([trA,trB])
                elif len(USEs_)==len(B):
                    USTs_list.append([trB,trA])
                    
                    
        
USEs_df=pd.DataFrame(USEs_list).rename({0:'transcript_A',1:'transcript_B',2:'number_of_event_exons',3:'event_exon'},axis=1)
USTs_df=pd.DataFrame(USTs_list).rename({0:'event_transcript',1:'refrence_transcript'},axis=1)
USEs_trs=list(USEs_df['transcript_A']) + list(USEs_df['transcript_B'])
USEs_trs = list(dict.fromkeys(USEs_trs))
USEs_trs_df=pd.DataFrame(USEs_trs).rename({0:'transcript_id'},axis=1)


USEs_df.to_csv('combined_USE_events',index = None, header=True,sep="\t")
USEs_trs_df.to_csv('USEs_transcripts',index = None, header=True,sep="\t")
USTs_df.to_csv('UST_events',index = None, header=True,sep="\t")



