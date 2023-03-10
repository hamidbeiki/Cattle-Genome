#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:12:46 2020

@author: beiki
"""

import pandas as pd
import sys
import collections
from multiprocessing import Pool

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
  
def get_exons_bed(input,tr):
    exons_bed={}
#    info=tr.split('.')
#    gene_id=info[0] + '.' + info[1]
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    strand=df.iloc[0,3]
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    tr_length=df.iloc[-1,2] - df.iloc[0,1]+1
    if strand=='+':
        df['%']=((df['ex_start'] - df.iloc[0,1])/float(tr_length))*100#distance to exon start to transcript start as % of tr length
    elif strand=='-':
        df['%']=((df.iloc[-1,2] - df['ex_end'])/float(tr_length))*100
    for i in range(df.shape[0]):
        exons_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0], df.iloc[i,1], df.iloc[i,2], tr, 99, strand, tr+'--' + str(df.iloc[i,4]))
    return(exons_bed)
    

def paralle_chunk_input1(chunk):
    bed=pd.DataFrame([])
    for tr in range(len(chunk)):
        info=get_exons_bed(chunk[tr][1],chunk[tr][0])#RNAbased or DNAbased methods can be selected
        tr_bed=pd.DataFrame.from_dict(info,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'exon_start', 2: 'exon_end', 3: 'tr_id', 4: 'value', 5: 'strand', 6: 'exon_id'},axis=1)
        bed=bed.append(tr_bed)
    return(bed)


def get_exon_number(ex_df,pos,strand,pos_type):
    count=0
    exon=-1
    pos=pos+1## this is because positions starts from 0 in python
    if pos_type=="end":
        pos=pos+2## because pos just show the start of stop codon relative to transcript start point (without introns). This has nothig to do with strands
    if strand=="-":
        ex_df=ex_df.sort_values(by='exon_start', ascending=False)
    for ex in range(ex_df.shape[0]):
        if count<pos and ex_df.iloc[0,7] < pos:
            count=count+ex_df.iloc[ex,7]
            dist_to_exon_end=count-pos # no need to add +1 here because this is distance not position
            exon=exon+1
    if exon==-1:
        dist_to_exon_end=ex_df.iloc[0,7] - pos# no need to add +1 here because this is distance not position
    return(exon,dist_to_exon_end)
    

def get_5UTRs(ex_df,exon,dist_to_exon_end,strand):
    if exon==-1:
        ex=1
    else:
        ex=exon+1
    if strand=="+":
        df=ex_df.iloc[0:ex,0:]
        s=df.iloc[-1,1]
        e=df.iloc[-1,2]
        info=df.iloc[-1,0:]
        e1=e-dist_to_exon_end-1
        if s<=e1:# the following conditions are amed to remove transcripts without 5' utr
            df.loc[(df['exon_start']==s)] = [[info[0],s,e1,info[3],99,info[5],info[6],info[7]]]
        elif s>e1 and df.shape[0]>1:
            df.drop(df.tail(1).index,inplace=True)
        elif s>e1 and df.shape[0]<2:
            df=pd.DataFrame([])
    elif strand=="-":
        ex_df=ex_df.sort_values(by='exon_start', ascending=False)
        df=ex_df.iloc[0:ex,0:]
        s=df.iloc[-1,1]
        e=df.iloc[-1,2]
        info=df.iloc[-1,0:]
        e1=e-dist_to_exon_end+1
        if e1<=e:# the following conditions are amed to remove transcripts without 5' utr
            df.loc[(df['exon_end']==e)] = [[info[0],e1,e,info[3],99,info[5],info[6],info[7]]]
        elif e1>e and df.shape[0]>1:
            df.drop(df.tail(1).index,inplace=True)
        elif e1>e and df.shape[0]<2:
            df=pd.DataFrame([])
    return(df)

def get_3UTRs(ex_df,exon,dist_to_exon_end,strand):
    if exon==-1:
        exon=0
    if strand=="+":
        df=ex_df.iloc[exon:,0:]
        s=df.iloc[-1,1]
        e=df.iloc[-1,2]
        info=df.iloc[-1,0:]
        s1=e-dist_to_exon_end+1
        if s1<=e:# the following conditions are amed to remove transcripts without 3' utr
            df.loc[(df['exon_start']==s)] = [[info[0],s1,e,info[3],99,info[5],info[6],info[7]]]
        elif s1>e and df.shape[0]>1:
            df.drop(df.tail(1).index,inplace=True)
        elif s1>e and df.shape[0]<2:
            df=pd.DataFrame([])    
    elif strand=="-":
        ex_df=ex_df.sort_values(by='exon_start', ascending=False)
        df=ex_df.iloc[exon:,0:]
        s=df.iloc[-1,1]
        e=df.iloc[-1,2]
        info=df.iloc[-1,0:]
        e1=e-dist_to_exon_end-1
        if s<=e1:# the following conditions are amed to remove transcripts without 3' utr
            df.loc[(df['exon_end']==e)] = [[info[0],s,e1,info[3],99,info[5],info[6],info[7]]]
        elif s>e1 and df.shape[0]>1:
            df.drop(df.tail(1).index,inplace=True)
        elif s>e1 and df.shape[0]<2:
            df=pd.DataFrame([]) 
    return(df)

def get_CDS(ex_df,exon1,dist_to_exon_end1,exon2,dist_to_exon_end2,strand):
    if exon1==-1:
        exon1=0
    if exon2==-1:
        exon2=0
    if strand=="+":
        df=ex_df.iloc[exon1:exon2+1,0:]
        info1=df.iloc[0,0:]
        info2=df.iloc[-1,0:]
        s1=info1[1]
        e1=info1[2]
        s2=info2[1]
        e2=info2[2]
        s_1=e1-dist_to_exon_end1
        e_2=e2-dist_to_exon_end2
        if s1!=s2:
            df.loc[(df['exon_start']==s1)] = [[info1[0],s_1,e1,info1[3],99,info1[5],info1[6],info1[7]]]
            df.loc[(df['exon_start']==s2)] = [[info2[0],s2,e_2,info2[3],99,info1[5],info2[6],info2[7]]]
        elif s1==s2:
            df.loc[(df['exon_start']==s1)] = [[info1[0],s_1,e_2,info1[3],99,info1[5],info1[6],info1[7]]]
    elif strand=="-":
        ex_df=ex_df.sort_values(by='exon_start', ascending=False)  
        df=ex_df.iloc[exon1:exon2+1,0:]
        info1=df.iloc[0,0:]
        info2=df.iloc[-1,0:]
        s1=info1[1]
        e1=info1[2]
        s2=info2[1]
        e2=info2[2]
        s_1=e2-dist_to_exon_end2
        e_1=e2
        e_2=e1-dist_to_exon_end1
        s_2=s1
        if s1!=s2:
            df.loc[(df['exon_end']==e1)] = [[info1[0],s_2,e_2,info1[3],99,info1[5],info1[6],info1[7]]]
            df.loc[(df['exon_end']==e2)] = [[info2[0],s_1,e_1,info2[3],99,info1[5],info2[6],info2[7]]]
        elif s1==s2:
            df.loc[(df['exon_end']==e1)] = [[info1[0],s_2,e_1,info1[3],99,info1[5],info1[6],info1[7]]]
    return(df)


def paralle_chunk_input2(chunk):
    global request_type
    res_list=[]   
    for i in range(len(chunk)):
       tr=chunk[i][0]
       strand=chunk[i][1]
       ex_df=tissue_ex_bed.loc[tissue_ex_bed['tr_id']==tr]
       cds_start=chunk[i][2]
       cds_end=chunk[i][3]
       if request_type =="5UTRs":
           pos=cds_start
           pos_type='start'
           info=get_exon_number(ex_df,pos,strand,pos_type)
           exon=info[0]
           dist=info[1]
           res=get_5UTRs(ex_df,exon,dist,strand)
           if res.shape[0]>0:# this is to remove those empty df related to transcripts with NO 5' utr
               res_list.append(res)
       elif request_type =="3UTRs":
           pos=cds_end
           pos_type='end'
           info=get_exon_number(ex_df,pos,strand,pos_type)
           exon=info[0]
           dist=info[1]
           res=get_3UTRs(ex_df,exon,dist,strand)
           if res.shape[0]>0:# this is to remove those empty df related to transcripts with NO 3' utr
               res_list.append(res)
           res_list.append(res)
       elif request_type =="CDS":
           pos=cds_start
           pos_type='start'
           info=get_exon_number(ex_df,pos,strand,pos_type)
           exon1=info[0]
           overlap1=info[1]
           pos=cds_end
           pos_type='end'
           info=get_exon_number(ex_df,pos,strand,pos_type)
           exon2=info[0]
           overlap2=info[1]
           res=get_CDS(ex_df,exon1,overlap1,exon2,overlap2,strand)
           res_list.append(res)
    out=pd.concat(res_list)        
    return(out)


def get_introns_bed(input,tr):
    introns_bed={}
#    info=tr.split('.')
#    gene_id=info[0] + '.' + info[1]
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    strand=df.iloc[0,3]
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    tr_length=df.iloc[-1,2] - df.iloc[0,1]+1
    if strand=="+":
        for i in range(df.shape[0]-1):
            percentage=((df.iloc[i,2]+1 - df.iloc[0,1])/float(tr_length))*100#distance to intron start to transcript start as % of tr length
            introns_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0], df.iloc[i,2]+1, df.iloc[i+1,1]-1, tr, 99, strand, tr+'--' + str(percentage))
    elif strand=='-':
        for i in range(df.shape[0]-1):
            percentage=((df.iloc[-1,2] - df.iloc[i+1,1]-1)/float(tr_length))*100
            introns_bed[i]= "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0], df.iloc[i,2]+1, df.iloc[i+1,1]-1, tr, 99, strand, tr+'--' + str(percentage))       
    return(introns_bed)
    
def paralle_chunk_input3(chunk):
    bed=pd.DataFrame([])
    for tr in range(len(chunk)):
        info=get_introns_bed(chunk[tr][1],chunk[tr][0])#RNAbased or DNAbased methods can be selected
        if len(info)>0:
            tr_bed=pd.DataFrame.from_dict(info,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'intron_start', 2: 'intron_end', 3: 'tr_id', 4: 'value', 5: 'strand', 6: 'intron_id'},axis=1)
            bed=bed.append(tr_bed)
    return(bed)
    
    
input_gff=sys.argv[1]# your gff file, e.g:  input_gff="Thalamus_final.collapsed.gff"
tissue=input_gff.split("_final")[0]
ex_coords = nested_dict()

for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        tid=raw[-1].split('; ')[1].split()[1][1:-2]
        exnumber=raw[3]
        ex_coords[tid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])


p = Pool()         
items = list(ex_coords.items())
processed_values= p.map(paralle_chunk_input1, get_chunks(items, 33))
tissue_ex_bed=pd.concat(processed_values)
    
p = Pool()  
processed_values= p.map(paralle_chunk_input3, get_chunks(items, 33))
tissue_int_bed=pd.concat(processed_values)


coding=pd.read_csv(tissue + '_coding_nonNMD_orfs',header=0,sep='\t')
coding=coding.rename({'ID':'tr_id'},axis=1)
idx = coding.groupby(['tr_id'])['length'].transform(max) == coding['length']# to just get two columns not entire rows: coding.groupby(['tr_id'], sort=False)['length'].max().reset_index()
coding=coding[idx]# this is to get longest orf per transcript, this way we can have a more clear deffinition ot UTRs
coding=coding.iloc[0:,0:3]
tnc=pd.read_csv(tissue + '_tnc_nonNMD_orfs',header=0,sep='\t')
tnc=tnc.rename({'ID':'tr_id'},axis=1)
idx = tnc.groupby(['tr_id'])['length'].transform(max) == tnc['length']
tnc=tnc[idx]
tnc=tnc.iloc[0:,0:3]

trs=pd.concat([coding,tnc])
info=tissue_ex_bed.merge(trs)
cds_info=info[['tr_id','strand','start','stop']].drop_duplicates()
tissue_ex_bed['exon_start']=tissue_ex_bed['exon_start'].astype(int)
tissue_ex_bed['exon_end']=tissue_ex_bed['exon_end'].astype(int)
tissue_ex_bed['length']= tissue_ex_bed.exon_end - tissue_ex_bed.exon_start + 1

items=[]
for i in range(cds_info.shape[0]):
    items.append(cds_info.iloc[i,0:])

request_type='5UTRs'
p = Pool()## global variable need to be assigned before pool()
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 33))
utr5_bed=pd.concat(processed_values)
utr5_bed=utr5_bed.iloc[0:,0:7]

del(request_type)
request_type='3UTRs'
p = Pool()
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 33))
utr3_bed=pd.concat(processed_values)
utr3_bed=utr3_bed.iloc[0:,0:7]

del(request_type)
request_type='CDS'
p = Pool()
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 33))
cds_bed=pd.concat(processed_values)
cds_bed=cds_bed.iloc[0:,0:7]

utr5_bed.to_csv(tissue + '_5UTRs.bed',index = None, header=True,sep="\t")
utr3_bed.to_csv(tissue + '_3UTRs.bed',index = None, header=True,sep="\t")
cds_bed.to_csv(tissue + '_CDS.bed',index = None, header=True,sep="\t")
tissue_ex_bed.to_csv(tissue + '_exons.bed',index = None, header=True,sep="\t")
tissue_int_bed.to_csv(tissue + '_introns.bed',index = None, header=True,sep="\t")


""" single CPU
utr5=[]
for i in range(cds_info.shape[0]):
    tr=cds_info.iloc[i,0]
    strand=cds_info.iloc[i,1]
    ex_df=tissue_ex_bed.loc[tissue_ex_bed['tr_id']==tr]
    cds_start=cds_info.iloc[i,2]
    cds_end=cds_info.iloc[i,3]
    if strand=="+":
        pos=cds_start
        info=get_exon_number(ex_df,pos)
        intr=info[0]
        overlap=info[1]
        res=get_5UTRs(ex_df,intr,overlap)
        utr5.append(res)
    elif strand=="-":
        pos=cds_end
        info=get_exon_number(ex_df,pos)
        intr=info[0]
        overlap=info[1]
        res=get_3UTRs(ex_df,intr,overlap)
        utr5.append(res)
"""    
  
