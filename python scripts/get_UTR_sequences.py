#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:32:22 2020

@author: beiki
"""

'''
module load python/3.6.3-u4oaxsb
module load py-pandas
module load py-biopython/1.70-py3-wos466g
'''

from Bio import SeqIO
import pandas as pd
import sys
#import time

def DNA_based_process_groups(df):#df is dataframe for each group
    info1=[df.iloc[0,3],df.iloc[0,0],df.iloc[0,5]]#chromosome and strand
    if df.shape[0]>1:
        info2=[]#coordinates
        a=list(df['exon_start'])
        b=list(df['exon_end'])
        for i in range(len(a)):
            info2.append([a[i],b[i]])
    elif df.shape[0]==1:
        info2=[df.iloc[0,1],df.iloc[0,2]]    
    return([info1,info2])

def UTR5_RNA_based_process_groups(df):#df is dataframe for each group
    df=df.copy()
    df['ex_length']=df['exon_end'] - df['exon_start'] +1
    strand=df.iloc[0,5]
    if strand=='+':
        df=df.sort_values(by='exon_start', ascending=True)
        coords=[]
        flag=0
        for i in df['ex_length']:
            if flag==0:
                coords.append([flag,i+flag-1])#-1 is to start coords from 0
                flag=i+flag-1
            else:
                coords.append([flag+1,i+flag+1-1])#-1 is to start coords from 0
                flag=i+flag+1-1 
    elif strand=='-':
        df=df.sort_values(by='exon_start', ascending=False)
        coords=[]
        flag=0
        for i in df['ex_length']:
            if flag==0:
                coords.append([flag,i+flag-1])#-1 is to start coords from 0
                flag=i+flag-1
            else:
                coords.append([flag+1,i+flag+1-1])#-1 is to start coords from 0
                flag=i+flag+1-1 
    return(df.iloc[0,3],coords)

def UTR3_RNA_based_process_groups(df):#df is dataframe for each group
    df=df.copy()
    df['ex_length']=df['exon_end'] - df['exon_start'] +1
    strand=df.iloc[0,5]
    if strand=='+':
        df=df.sort_values(by='exon_start', ascending=False)
        coords=[]
        flag=0
        for i in df['ex_length']:
            if flag==0:
                coords.append([flag,i+flag-1])#-1 is to start coords from 0
                flag=i+flag-1
            else:
                coords.append([flag+1,i+flag+1-1])#-1 is to start coords from 0
                flag=i+flag+1-1 
    if strand=='-':
        df=df.sort_values(by='exon_start', ascending=True)
        coords=[]
        flag=0
        for i in df['ex_length']:
            if flag==0:
                coords.append([flag,i+flag-1])#-1 is to start coords from 0
                flag=i+flag-1
            else:
                coords.append([flag+1,i+flag+1-1])#-1 is to start coords from 0
                flag=i+flag+1-1 
    return(df.iloc[0,3],coords)

bed=pd.read_csv(sys.argv[1],sep='\t',header=0)#combined_5UTRs.bed
out=sys.argv[1].split('_')[1].split('.')[0]#outfile name
#a = pd.to_numeric(bed.exon_start, errors='coerce')
#bed['exon_start']=a
#a = pd.to_numeric(bed.exon_end, errors='coerce')
#bed['exon_end']=a
#start = time.time()
if out=='5UTRs':
    bed_info=bed.groupby('tr_id').apply(UTR5_RNA_based_process_groups).to_dict()# if your fasta file is genome, select DNA_based_process_groups
elif out=='3UTRs':
     bed_info=bed.groupby('tr_id').apply(UTR3_RNA_based_process_groups).to_dict()# if your fasta file is genome, select DNA_based_process_groups   
#end = time.time()
#print(end-start)



ofile = open(out+'.fasta', "w")
if out=='5UTRs':
    for rec in SeqIO.parse(sys.argv[2], "fasta"):#your transcript sequences(combined-transcripts.fa) or your genome file
        if rec.id in bed_info:
            info=bed_info[rec.id]
            seq=str(rec.seq)
            utr_seq=str()
            for c in info[1]:
                s=seq[c[0]:c[1]+1]
                utr_seq=utr_seq+s
            ofile.write(">" + rec.id + "\n" +utr_seq + "\n")
elif out=='3UTRs':
    for rec in SeqIO.parse(sys.argv[2], "fasta"):#your transcript sequences(combined-transcripts.fa) or your genome file
        if rec.id in bed_info:
            l=bed_info[rec.id][1][::-1]#reverse the order of lists in the nested list
            seq=str(rec.seq)[::-1]#reverse seq to start from 3' end of transcript
            utr_seq=str()
            for c in l:
                s=seq[c[0]:c[1]+1][::-1]#to be matched with the correct direction of transcript
                utr_seq=utr_seq+s
            ofile.write(">" + rec.id + "\n" +utr_seq + "\n")            
            

'''
coords=[]
flag=0
for i in df['ex_length']:
    if flag==0:
        coords.append([flag,i+flag-1])#-1 is to start coords from 0
        flag=i+flag-1
    else:
        coords.append([flag+1,i+flag+1-1])#-1 is to start coords from 0
        flag=i+flag+1-1   

trs=bed['tr_id'].drop_duplicates()
flag=1
for tr in trs:
    df=bed.loc[bed['tr_id']==tr]
    out=RNA_based_process_groups(df)
    print(flag)
    flag=flag+1


'''