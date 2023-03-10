#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:56:24 2021

@author: beiki
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def get_biotype_change(AS,AS_type):
    global coding
    global noncoding
    coding_to_coding=0
    coding_to_noncoding=0
    noncoding_to_coding=0
    noncoding_to_noncoding=0
    events=list(AS['alternative_transcripts'])
    total=list(AS['total_transcripts'])
    for i in range(len(events)):
        ev=events[i].split(',')
        r=total[i].split(',')
        ref=list(set(r) - set(ev))
        ev_coding=list(set(ev) & set(coding))
        ev_noncoding=list(set(ev) - set(ev_coding))
        ref_coding=list(set(ref) & set(coding))
        ref_noncoding=list(set(ref) - set(ref_coding))
        if len(ref_coding)>0 and len(ev_coding)>0:
            coding_to_coding=coding_to_coding+(len(ref_coding)*len(ev_coding))
        elif len(ref_coding)>0 and len(ev_noncoding)>0:
            coding_to_noncoding=coding_to_noncoding+(len(ref_coding)*len(ev_noncoding))
        elif len(ref_noncoding)>0 and len(ev_coding)>0:
            noncoding_to_coding=noncoding_to_coding+(len(ref_noncoding)*len(ev_coding))
        elif len(ref_noncoding)>0 and len(ev_noncoding)>0:    
            noncoding_to_noncoding=noncoding_to_noncoding+(len(ref_noncoding)*len(ev_noncoding))
    return([AS_type,coding_to_coding,coding_to_noncoding,noncoding_to_coding,noncoding_to_noncoding])

def get_tissue_biotype_change(AS,AS_type,expr):
    global coding
    global noncoding
    coding_to_coding=0
    coding_to_noncoding=0
    noncoding_to_coding=0
    noncoding_to_noncoding=0
    expr=set(expr)
    coding_expr=set(expr & set(coding))
    events=list(AS['alternative_transcripts'])
    total=list(AS['total_transcripts'])
    for i in range(len(events)):
        ev=events[i].split(',')
        ev=list(set(ev) & expr)
        r=total[i].split(',')
        r=list(set(r) & expr)
        if len(r)>0 and len(ev)>0:
            ref=list(set(r) - set(ev))
            ev_coding=list(set(ev) & coding_expr)
            ev_noncoding=list(set(ev) - set(ev_coding))
            ref_coding=list(set(ref) & coding_expr)
            ref_noncoding=list(set(ref) - set(ref_coding))
            if len(ref_coding)>0 and len(ev_coding)>0:
                coding_to_coding=coding_to_coding+(len(ref_coding)*len(ev_coding))
            elif len(ref_coding)>0 and len(ev_noncoding)>0:
                coding_to_noncoding=coding_to_noncoding+(len(ref_coding)*len(ev_noncoding))
            elif len(ref_noncoding)>0 and len(ev_coding)>0:
                noncoding_to_coding=noncoding_to_coding+(len(ref_noncoding)*len(ev_coding))
            elif len(ref_noncoding)>0 and len(ev_noncoding)>0:    
                noncoding_to_noncoding=noncoding_to_noncoding+(len(ref_noncoding)*len(ev_noncoding))
    return([AS_type,coding_to_coding,coding_to_noncoding,noncoding_to_coding,noncoding_to_noncoding])

def get_event_per_tr(AS,AS_type):
    res={}
    events=list(AS['alternative_transcripts'])
    total=list(AS['total_transcripts'])
    for i in range(len(events)):
        ev=events[i].split(',')
        r=total[i].split(',')
        ref=list(set(r) - set(ev))
        num=len(ev)*len(ref)
        al=ref + ev
        for i in al:
            if i not in res:
                res[i]=num
            else:
                res[i]=res[i]+num
    df=pd.DataFrame.from_dict(res, orient='index').reset_index().rename({'index':'transcript_id',0:'event'},axis=1)
    return(df)

def get_tissue_event_per_tr(AS,AS_type,expr):
    res={}
    expr=set(expr)
    events=list(AS['alternative_transcripts'])
    total=list(AS['total_transcripts'])
    for i in range(len(events)):
        ev=events[i].split(',')
        ev=list(set(ev) & expr)
        r=total[i].split(',')
        r=list(set(r) & expr)
        if len(r)>0 and len(ev)>0:
            ref=list(set(r) - set(ev))
            num=len(ev)*len(ref)
            al=ref + ev
            for i in al:
                if i not in res:
                    res[i]=num
                else:
                    res[i]=res[i]+num
    df=pd.DataFrame.from_dict(res, orient='index').reset_index().rename({'index':'transcript_id',0:'event'},axis=1)
    df=df.loc[df['event']>0]
    return(df)
            
def get_info(ex,biotype):
    ex = ex.replace(0, np.NaN)
    ex_num=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy
    row_medians = np.nanmedian(ex_num, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
    row_counts=np.count_nonzero(~np.isnan(ex_num),axis=1)#https://stackoverflow.com/questions/21778118/counting-the-number-of-non-nan-elements-in-a-numpy-ndarray-in-python
    df1=pd.DataFrame(row_counts).rename({0:'detected_tissues'},axis=1)
    df2=pd.DataFrame(row_medians).rename({0:'median_expression'},axis=1)
    df1[biotype]=ex.index
    df2[biotype]=ex.index
    m=df1.merge(df2)
    m=m[[biotype,'detected_tissues','median_expression']]
    return(m)

def get_other_A5_A3_biotype_changes(events,event_type):
    global coding
    global noncoding
    global biotypes
    events=events.rename({'event_transcript':'transcript_id'},axis=1)
    m1=events.merge(biotypes)
    m1=m1.rename({'transcript_id':'event_transcript','biotype':'event_tr_biotype','refrence_transcript':'transcript_id'},axis=1)
    m2=m1.merge(biotypes)
    events=m2.rename({'transcript_id':'refrence_transcript','biotype':'refrence_tr_biotype'},axis=1)
    coding_coding=events.loc[(events['refrence_tr_biotype']=='protein_coding_transcripts') & (events['event_tr_biotype']=='protein_coding_transcripts')]
    coding_noncoding=events.loc[(events['refrence_tr_biotype']=='protein_coding_transcripts') & (events['event_tr_biotype']!='protein_coding_transcripts')]
    noncoding_coding=events.loc[(events['refrence_tr_biotype']!='protein_coding_transcripts') & (events['event_tr_biotype']=='protein_coding_transcripts')]
    noncoding_noncoding=events.loc[(events['refrence_tr_biotype']!='protein_coding_transcripts') & (events['event_tr_biotype']!='protein_coding_transcripts')]
    results=[event_type,coding_coding.shape[0],coding_noncoding.shape[0],noncoding_coding.shape[0],noncoding_noncoding.shape[0]]
    return(results)

def get_tissue_other_A5_A3_biotype_changes(events,event_type,expr_d):
    global coding
    global noncoding
    global biotypes
    events_expr=events.loc[(events.event_transcript.isin(expr_d.transcript_id)) & (events.refrence_transcript.isin(expr_d.transcript_id))]
    events=events_expr
    events=events.rename({'event_transcript':'transcript_id'},axis=1)
    m1=events.merge(biotypes)
    m1=m1.rename({'transcript_id':'event_transcript','biotype':'event_tr_biotype','refrence_transcript':'transcript_id'},axis=1)
    m2=m1.merge(biotypes)
    events=m2.rename({'transcript_id':'refrence_transcript','biotype':'refrence_tr_biotype'},axis=1)
    coding_coding=events.loc[(events['refrence_tr_biotype']=='protein_coding_transcripts') & (events['event_tr_biotype']=='protein_coding_transcripts')]
    coding_noncoding=events.loc[(events['refrence_tr_biotype']=='protein_coding_transcripts') & (events['event_tr_biotype']!='protein_coding_transcripts')]
    noncoding_coding=events.loc[(events['refrence_tr_biotype']!='protein_coding_transcripts') & (events['event_tr_biotype']=='protein_coding_transcripts')]
    noncoding_noncoding=events.loc[(events['refrence_tr_biotype']!='protein_coding_transcripts') & (events['event_tr_biotype']!='protein_coding_transcripts')]
    results=[event_type,coding_coding.shape[0],coding_noncoding.shape[0],noncoding_coding.shape[0],noncoding_noncoding.shape[0]]
    return(results)
        
def get_tissue_other_A5_A3_event_per_tr(events,expr_d):
    events_expr=events.loc[(events.event_transcript.isin(expr_d.transcript_id)) & (events.refrence_transcript.isin(expr_d.transcript_id))]
    other_events_per_trs=events_expr.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)
    return(other_events_per_trs)


A3=pd.read_csv('combined_SUPPA_output_A3_strict.ioe',header=0, sep='\t')
A5=pd.read_csv('combined_SUPPA_output_A5_strict.ioe',header=0, sep='\t')
AF=pd.read_csv('combined_SUPPA_output_AF_strict.ioe',header=0, sep='\t')
AL=pd.read_csv('combined_SUPPA_output_AL_strict.ioe',header=0, sep='\t')
MX=pd.read_csv('combined_SUPPA_output_MX_strict.ioe',header=0, sep='\t')
RI=pd.read_csv('combined_SUPPA_output_RI_strict.ioe',header=0, sep='\t')
SE=pd.read_csv('combined_SUPPA_output_SE_strict.ioe',header=0, sep='\t')
USE=pd.read_csv('combined_USE_events',header=0,sep='\t')
EIs_df=pd.read_csv('trs_with_exons_covered_by_other_trs_introns',names=['transcript_id'])    
intersect=pd.read_csv('intersect',sep='\t',names=['ids',1,2,3,'transcript_id',5,6,7,'p1',9,10,11,12,13,14,15,16,'p2',18])
biotypes=pd.read_csv("../final_transcript_biotypes",names=['transcript_id','biotype'],sep='\t')
gene_biotypes=pd.read_csv("../final_gene_biotypes",names=['gene_id','biotype'],sep='\t')
expressions=pd.read_csv("../quantification/RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
df=expressions.iloc[:,1:]
df=df.loc[(df!=0).any(1)]
ids=expressions['transcript_id']
expressions=df
expressions['transcript_id']=ids
gene_expressions=pd.read_csv("../quantification/RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t")
df=gene_expressions.iloc[:,1:]
df=df.loc[(df!=0).any(1)]
gene_ids=gene_expressions['gene_id']
gene_expressions=df
gene_expressions['gene_id']=gene_ids
AS_candidate_genes=pd.read_csv('AS_candidate_genes',names=['gene_id'])
ensembl=pd.read_csv('../annotation/combined_genes_ensembl_equivalent',header=0,names=['gene_id','ensembl'],sep=' ')
spliced_genes=pd.read_csv('spliced_genes',names=['gene_id'])
unspliced_trs=pd.read_csv('unspliced_transcripts_2',names=['transcript_id'])
coding=list(biotypes.loc[biotypes['biotype']=='protein_coding_transcripts']['transcript_id'])
noncoding=list(biotypes.loc[biotypes['biotype']!='protein_coding_transcripts']['transcript_id'])

m=EIs_df.merge(intersect)
groups=m.groupby(['ids'])['p1'].min().reset_index()
other_A5=groups.loc[groups['p1']==0]
other_A3=groups.loc[groups['p1']!=0]
other_A5_events = other_A5['ids'].str.split("--", n = 1, expand = True).rename({0:'event_transcript',1:'refrence_transcript'},axis=1)
other_A3_events = other_A3['ids'].str.split("--", n = 1, expand = True).rename({0:'event_transcript',1:'refrence_transcript'},axis=1)
other_A3_events_per_trs=other_A3_events.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)
other_A5_events_per_trs=other_A5_events.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)

results=[]
res=get_biotype_change(A3,'A3')
res2=get_other_A5_A3_biotype_changes(other_A3_events,'A3')
res3=[res[0],res[1]+res2[1],res[2]+res2[2],res[3]+res2[3],res[4]+res2[4]]
results.append(res3)
res=get_biotype_change(A5,'A5')
res2=get_other_A5_A3_biotype_changes(other_A5_events,'A5')
res3=[res[0],res[1]+res2[1],res[2]+res2[2],res[3]+res2[3],res[4]+res2[4]]
results.append(res3)
res=get_biotype_change(AF,'AF')
results.append(res)
res=get_biotype_change(AL,'AL')
results.append(res)
res=get_biotype_change(MX,'MX')
results.append(res)
res=get_biotype_change(RI,'RI')
results.append(res)
res=get_biotype_change(SE,'SE')
results.append(res)
results_df=pd.DataFrame(results).rename({0:'AS_event',1:"coding_to_coding",2:"coding_to_noncoding",3:"noncoding_to_coding",4:"noncoding_to_noncoding"},axis=1)   
results_df.to_csv('transcript_biotype_changes_due_to_AS_events',index = None, header=True,sep="\t")


A3_event_per_tr=get_event_per_tr(A3,'A3')
df=pd.concat([A3_event_per_tr,other_A3_events_per_trs])
A3_event_per_tr=df.groupby(['transcript_id'])['event'].sum().reset_index().rename({0:'event'},axis=1)
A5_event_per_tr=get_event_per_tr(A5,'A5')
df=pd.concat([A5_event_per_tr,other_A5_events_per_trs])
A5_event_per_tr=df.groupby(['transcript_id'])['event'].sum().reset_index().rename({0:'event'},axis=1)
AF_event_per_tr=get_event_per_tr(AF,'AF')
AL_event_per_tr=get_event_per_tr(AL,'AL')
MX_event_per_tr=get_event_per_tr(MX,'MX')
RI_event_per_tr=get_event_per_tr(RI,'RI')
SE_event_per_tr=get_event_per_tr(SE,'SE')
a=USE.groupby(['transcript_A'])['number_of_event_exons'].sum().reset_index(name='event').rename({'transcript_A':'transcript_id'},axis=1)
b=USE.groupby(['transcript_B'])['number_of_event_exons'].sum().reset_index(name='event').rename({'transcript_B':'transcript_id'},axis=1)
c=pd.concat([a,b])
USE_event_per_tr=c.groupby(['transcript_id'])['event'].sum().reset_index(name='event')

df=pd.concat([A3_event_per_tr,A5_event_per_tr,AF_event_per_tr,AL_event_per_tr,MX_event_per_tr,RI_event_per_tr,SE_event_per_tr,USE_event_per_tr])
event_per_tr=df.groupby(['transcript_id'])['event'].size().reset_index(name='total_AS_event')

l=event_per_tr['transcript_id'].str.split('.').to_list()
genes=[]
for i in l:
    g=str(i[0])+'.'+str(i[1])
    genes.append(g)

gene_id=pd.DataFrame(genes).rename({0:'gene_id'},axis=1)
df2=event_per_tr.copy()
df2['gene_id']=gene_id
event_per_gene=df2.groupby(['gene_id'])['total_AS_event'].sum().reset_index(name='total_AS_event')

A3_event_per_tr.to_csv('A3_event_per_transcript',index = None, header=True,sep="\t")
A5_event_per_tr.to_csv('A5_event_per_transcript',index = None, header=True,sep="\t")
AF_event_per_tr.to_csv('AF_event_per_transcript',index = None, header=True,sep="\t")
AL_event_per_tr.to_csv('AL_event_per_transcript',index = None, header=True,sep="\t")
MX_event_per_tr.to_csv('MX_event_per_transcript',index = None, header=True,sep="\t")
RI_event_per_tr.to_csv('RI_event_per_transcript',index = None, header=True,sep="\t")
SE_event_per_tr.to_csv('SE_event_per_transcript',index = None, header=True,sep="\t")
USE_event_per_tr.to_csv('USE_event_per_transcript',index = None, header=True,sep="\t")
event_per_tr.to_csv('total_AS_event_per_transcript',index = None, header=True,sep="\t")
event_per_gene.to_csv('total_AS_event_per_gene',index = None, header=True,sep="\t")
df=A3_event_per_tr['transcript_id'].to_frame()
df.to_csv('A3_transcripts',index = None, header=True,sep="\t")
df=A5_event_per_tr['transcript_id'].to_frame()
df.to_csv('A5_transcripts',index = None, header=True,sep="\t")

event_per_gene['total_AS_event'].quantile([0.0, .5, .90, .95])
h=event_per_gene.loc[event_per_gene['total_AS_event']>24]#genes with highest AS events
h2=h.merge(ensembl)['ensembl'].drop_duplicates().to_frame()
h2.to_csv('genes_with_highest_AS_event',index = None, header=False,sep="\t")
plt.hist(event_per_gene['total_AS_event'], bins=50)
plt.xlabel('Number AS event per gene')
plt.ylabel('Number genes')
plt.axvline(24, linestyle='--',color='red')
plt.savefig('distribution_of_genes_AS_events.png',dpi=300)
plt.close()
#plt.show()

m=event_per_gene.merge(gene_biotypes)
m.groupby(['biotype']).size().reset_index(name='count')


tissues=list(expressions.columns)[0:-1]
tissues_info={}
for tissue in tissues:
    results=[]
    df=expressions[['transcript_id',tissue]]
    expr=df.loc[df[tissue]>0]
    expr=list(expr['transcript_id'])
    expr_d=pd.DataFrame(expr).rename({0:'transcript_id'},axis=1)
    df=gene_expressions[['gene_id',tissue]]
    expr2=df.loc[df[tissue]>0]
    expr2=list(expr2['gene_id'])
    res=get_tissue_biotype_change(A3,'A3',expr)
    res2=get_tissue_other_A5_A3_biotype_changes(other_A3_events,'A3',expr_d)
    res3=[res[0],res[1]+res2[1],res[2]+res2[2],res[3]+res2[3],res[4]+res2[4]]
    results.append(res3)
    res=get_tissue_biotype_change(A5,'A5',expr)
    res2=get_tissue_other_A5_A3_biotype_changes(other_A5_events,'A5',expr_d)
    res3=[res[0],res[1]+res2[1],res[2]+res2[2],res[3]+res2[3],res[4]+res2[4]]
    results.append(res3)
    res=get_tissue_biotype_change(AF,'AF',expr)
    results.append(res)
    res=get_tissue_biotype_change(AL,'AL',expr)
    results.append(res)
    res=get_tissue_biotype_change(MX,'MX',expr)
    results.append(res)
    res=get_tissue_biotype_change(RI,'RI',expr)
    results.append(res)
    res=get_tissue_biotype_change(SE,'SE',expr)
    results.append(res)    
    biotype_change_summary=pd.DataFrame(results).rename({0:'AS_event',1:"coding_to_coding",2:"coding_to_noncoding",3:"noncoding_to_coding",4:"noncoding_to_noncoding"},axis=1)
    A3_event_per_tr=get_tissue_event_per_tr(A3,'A3',expr)
    other_A3_events_per_trs=get_tissue_other_A5_A3_event_per_tr(other_A3_events,expr_d)
    df=pd.concat([A3_event_per_tr,other_A3_events_per_trs])
    A3_event_per_tr=df.groupby(['transcript_id'])['event'].sum().reset_index().rename({0:'event'},axis=1)
    A5_event_per_tr=get_tissue_event_per_tr(A5,'A5',expr)
    other_A5_events_per_trs=get_tissue_other_A5_A3_event_per_tr(other_A5_events,expr_d)
    df=pd.concat([A5_event_per_tr,other_A5_events_per_trs])
    A5_event_per_tr=df.groupby(['transcript_id'])['event'].sum().reset_index().rename({0:'event'},axis=1)     
    AF_event_per_tr=get_tissue_event_per_tr(AF,'AF',expr)
    AL_event_per_tr=get_tissue_event_per_tr(AL,'AL',expr)
    MX_event_per_tr=get_tissue_event_per_tr(MX,'MX',expr)
    RI_event_per_tr=get_tissue_event_per_tr(RI,'RI',expr)
    SE_event_per_tr=get_tissue_event_per_tr(SE,'SE',expr)
    USE_expr=USE.loc[(USE.transcript_A.isin(expr_d.transcript_id)) & (USE.transcript_B.isin(expr_d.transcript_id))]
    a=USE_expr.groupby(['transcript_A'])['number_of_event_exons'].sum().reset_index(name='event').rename({'transcript_A':'transcript_id'},axis=1)
    b=USE_expr.groupby(['transcript_B'])['number_of_event_exons'].sum().reset_index(name='event').rename({'transcript_B':'transcript_id'},axis=1)
    c=pd.concat([a,b])
    USE_event_per_tr=c.groupby(['transcript_id'])['event'].sum().reset_index(name='event')   
    df=pd.concat([A3_event_per_tr,A5_event_per_tr,AF_event_per_tr,AL_event_per_tr,MX_event_per_tr,RI_event_per_tr,SE_event_per_tr,USE_event_per_tr])
    event_per_tr=df.groupby(['transcript_id'])['event'].size().reset_index(name='total_AS_event')
    l=event_per_tr['transcript_id'].str.split('.').to_list()
    genes=[]
    for i in l:
        g=str(i[0])+'.'+str(i[1])
        genes.append(g)
    gene_id=pd.DataFrame(genes).rename({0:'gene_id'},axis=1)
    df2=event_per_tr.copy()
    df2['gene_id']=gene_id
    event_per_gene=df2.groupby(['gene_id'])['total_AS_event'].sum().reset_index(name='total_AS_event')
    expr2_d=pd.DataFrame(expr2).rename({0:'gene_id'},axis=1)
    gene_percent=(event_per_gene.shape[0]/expr2_d.shape[0])*100
    AS_transcript_percent=(event_per_tr.shape[0]/expr_d.shape[0])*100
    unspliced_transcript_percent=(expr_d.merge(unspliced_trs).shape[0]/expr_d.shape[0])*100
    spliced_trs_from_single_spliceTrGenes=100-AS_transcript_percent-unspliced_transcript_percent
    AS_event_summary=[tissue,A3_event_per_tr['event'].sum(),A5_event_per_tr['event'].sum(),AF_event_per_tr['event'].sum(),AL_event_per_tr['event'].sum(),MX_event_per_tr['event'].sum(),RI_event_per_tr['event'].sum(),SE_event_per_tr['event'].sum(),USE_event_per_tr['event'].sum()]
    AS_event_summary2=[tissue,A3_event_per_tr['event'].median(),A5_event_per_tr['event'].median(),AF_event_per_tr['event'].median(),AL_event_per_tr['event'].median(),MX_event_per_tr['event'].median(),RI_event_per_tr['event'].median(),SE_event_per_tr['event'].median(),USE_event_per_tr['event'].median()]
    AS_tr_summary=[tissue,A3_event_per_tr.shape[0],A5_event_per_tr.shape[0],AF_event_per_tr.shape[0],AL_event_per_tr.shape[0],MX_event_per_tr.shape[0],RI_event_per_tr.shape[0],SE_event_per_tr.shape[0],USE_event_per_tr.shape[0]]
    spliced_tr_percent=(event_per_tr.shape[0]/(expr_d.shape[0]-expr_d.merge(unspliced_trs).shape[0]))*100
    unspliced_trs_expr=expr_d.merge(unspliced_trs)
    spliced_trs_expr=expr_d.loc[~ expr_d.transcript_id.isin(unspliced_trs_expr.transcript_id)]
    spliced_gene_percent=(event_per_gene.shape[0]/expr2_d.merge(spliced_genes).shape[0])*100
    m=expr2_d.merge(gene_biotypes)
    coding_gene_expr=m.loc[m['biotype']=='protein-coding'].shape[0]
    m=event_per_gene.merge(gene_biotypes)
    coding_AS_gene_expr=m.loc[m['biotype']=='protein-coding'].shape[0]
    coding_gene_percent=(coding_AS_gene_expr/coding_gene_expr)*100
    m=expr2_d.merge(gene_biotypes)
    noncoding_gene_expr=m.loc[m['biotype']!='protein-coding'].shape[0]    
    m=event_per_gene.merge(gene_biotypes)
    noncoding_AS_gene_expr=m.loc[m['biotype']!='protein-coding'].shape[0]
    noncoding_gene_percent=(noncoding_AS_gene_expr/noncoding_gene_expr)*100
    tissues_info[tissue]=[biotype_change_summary,event_per_tr,event_per_gene,AS_event_summary,AS_tr_summary,spliced_gene_percent,spliced_tr_percent,spliced_trs_from_single_spliceTrGenes,unspliced_transcript_percent,AS_event_summary2,coding_gene_percent,noncoding_gene_percent]
    print(tissue)
    
np.save('tissue_based_transcript_biotype_changes_due_to_AS_events.npy', tissues_info)
#tissues_info = np.load('tissue_based_transcript_biotype_changes_due_to_AS_events.npy',allow_pickle='TRUE').item()

for tissue in tissues:
    df=tissues_info[tissue][1]
    df.to_csv(tissue + '_AS_event_per_transcript',index = None, header=True,sep="\t")
    df=tissues_info[tissue][2]
    df.to_csv(tissue + '_AS_event_per_gene',index = None, header=True,sep="\t")

summary_df=pd.DataFrame()
for key,value in tissues_info.items(): 
    df=value[0].copy()
    df["sum"] = df.loc[:,"coding_to_coding":"noncoding_to_noncoding"].sum(axis=1)
    df_new = df.loc[:,"coding_to_coding":"noncoding_to_noncoding"].div(df["sum"], axis=0)
    df_new['tissue']=key
    df_new=df_new*100
    df_new['AS_event']=df['AS_event']
    summary_df=pd.concat([summary_df,df_new])

biotype_change_summary=summary_df
events=list(summary_df['AS_event'].drop_duplicates())
flag=1
for event in events:
    globals() ['a' + str(flag)]=summary_df.loc[summary_df['AS_event']==event]['coding_to_coding'].to_frame().rename({'coding_to_coding':"info"},axis=1)
    globals() ['a' + str(flag)]['type1']=event
    globals() ['a' + str(flag)]['type2']='coding_to_coding'
    flag=flag+1
    globals() ['a' + str(flag)]=summary_df.loc[summary_df['AS_event']==event]['coding_to_noncoding'].to_frame().rename({'coding_to_noncoding':"info"},axis=1)
    globals() ['a' + str(flag)]['type1']=event
    globals() ['a' + str(flag)]['type2']='coding_to_noncoding'
    flag=flag+1    
    globals() ['a' + str(flag)]=summary_df.loc[summary_df['AS_event']==event]['noncoding_to_coding'].to_frame().rename({'noncoding_to_coding':"info"},axis=1)
    globals() ['a' + str(flag)]['type1']=event
    globals() ['a' + str(flag)]['type2']='noncoding_to_coding'
    flag=flag+1
    globals() ['a' + str(flag)]=summary_df.loc[summary_df['AS_event']==event]['noncoding_to_noncoding'].to_frame().rename({'noncoding_to_noncoding':"info"},axis=1)
    globals() ['a' + str(flag)]['type1']=event
    globals() ['a' + str(flag)]['type2']='noncoding_to_noncoding'
    print('a' + str(flag))
    flag=flag+1
    

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_transcript_biotype_change_due_to_AS_events.png',dpi=300)
plt.close()

flag=1
for key,value in tissues_info.items(): 
    df=value[2].copy()
    globals() ['a' + str(flag)]=df.rename({'total_AS_event':'info'},axis=1)
    globals() ['a' + str(flag)]['type']=flag
    print('a' + str(flag))
    flag=flag+1


df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.rcParams["figure.figsize"] = (18,6)
plt.savefig('tissues_AS_events_per_gene.png',dpi=300)
plt.close()

AS_event_summary_l=[]
for key,value in tissues_info.items(): 
    l=value[3].copy()
    AS_event_summary_l.append(l)

AS_event_summary=pd.DataFrame(AS_event_summary_l).rename({0:'Tissue',1:'A3',2:'A5',3:'AF',4:'AL',5:'MX',6:'RI',7:'SE',8:'USE'},axis=1)
events=['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE','USE']
flag=1
for event in events:
    globals() ['a' + str(flag)]=AS_event_summary[event].to_frame().rename({event:'info'},axis=1)
    globals() ['a' + str(flag)]['type']=event
    print('a' + str(flag))
    flag=flag+1

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_AS_events_distribution.png',dpi=300)
plt.close()
    

flag=1
for key,value in tissues_info.items(): 
    df=value[2].copy()
    globals() ['a' + str(flag)]=df.rename({'total_AS_event':'info'},axis=1)
    globals() ['a' + str(flag)]['type']=flag
    print('a' + str(flag))
    flag=flag+1


df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.rcParams["figure.figsize"] = (18,6)
plt.savefig('tissues_AS_events_per_gene.png',dpi=300)
plt.close()


AS_tr_summary_l=[]
for key,value in tissues_info.items(): 
    l=value[4].copy()
    AS_tr_summary_l.append(l)

AS_tr_summary=pd.DataFrame(AS_tr_summary_l).rename({0:'Tissue',1:'A3',2:'A5',3:'AF',4:'AL',5:'MX',6:'RI',7:'SE',8:'USE'},axis=1)
events=['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE', 'USE']
flag=1
for event in events:
    globals() ['a' + str(flag)]=AS_tr_summary[event].to_frame().rename({event:'info'},axis=1)
    globals() ['a' + str(flag)]['type']=event
    print('a' + str(flag))
    flag=flag+1

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_AS_transcripts_distribution.png',dpi=300)
plt.close()

df=AS_event_summary    
df["sum"] = df.loc[:,"A3":"USE"].sum(axis=1)    
df_new = df.loc[:,"A3":"USE"].div(df["sum"], axis=0)
df_new=df_new*100
df_new['Tissue']=df['Tissue']
AS_event_summary=df_new
events=['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE','USE']
flag=1
for event in events:
    globals() ['a' + str(flag)]=AS_event_summary[event].to_frame().rename({event:'info'},axis=1)
    globals() ['a' + str(flag)]['type']=event
    print('a' + str(flag))
    flag=flag+1

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_AS_events%_distribution.png',dpi=300)
plt.close()

df=AS_tr_summary    
df["sum"] = df.loc[:,"A3":"USE"].sum(axis=1)    
df_new = df.loc[:,"A3":"USE"].div(df["sum"], axis=0)
df_new=df_new*100
df_new['Tissue']=df['Tissue']
AS_tr_summary=df_new
events=['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE','USE']
flag=1
for event in events:
    globals() ['a' + str(flag)]=AS_tr_summary[event].to_frame().rename({event:'info'},axis=1)
    globals() ['a' + str(flag)]['type']=event
    print('a' + str(flag))
    flag=flag+1

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_AS_transcripts%_distribution.png',dpi=300)
plt.close()
    
summary=[]
for key,value in tissues_info.items(): 
    summary.append([key,value[1].median()])

summary=[]
for key,value in tissues_info.items(): 
    df=value[1]
    p=(df.loc[df['total_AS_event']>1].shape[0]/df.shape[0])*100
    summary.append([key,p])
pd.DataFrame(summary)[1].median()

summary=[]
for key,value in tissues_info.items(): 
    summary.append([key,value[6],value[7],value[8]])    

summary_df=pd.DataFrame(summary).rename({0:'tissue',1:'AS_trs',2:'splice_singleSpliceGene',3:'unsplice'},axis=1)
a=summary_df['AS_trs'].to_frame().rename({'AS_trs':"info"},axis=1)
a['type']='AS_trs'
b=summary_df['splice_singleSpliceGene'].to_frame().rename({'splice_singleSpliceGene':"info"},axis=1)
b['type']='splice_singleSpliceGene'
c=summary_df['unsplice'].to_frame().rename({'unsplice':"info"},axis=1)
c['type']='unsplice'
df=pd.concat([a,b,c,])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('Further_transcript_classification.jpg',dpi=300)
plt.close()


summary=[]
for key,value in tissues_info.items(): 
    summary.append(value[9])    
summary_df=pd.DataFrame(summary).rename({0:'tissue',1:"A3",2:"A5",3:"AF",4:"AL",5:"MX",6:"RI",7:"SE",8:'USE'},axis=1)
events=['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE', 'USE']
flag=1
for event in events:
    globals() ['a' + str(flag)]=summary_df[event].to_frame().rename({event:'info'},axis=1)
    globals() ['a' + str(flag)]['type']=event
    print('a' + str(flag))
    flag=flag+1

df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissues_AS_event_per_transcripts_distribution.png',dpi=300)
plt.close()
      

summary=[]
for key,value in tissues_info.items(): 
    summary.append([key,value[5]]) 
summary_df=pd.DataFrame(summary).rename({0:'tissue',1:'info'},axis=1)
a=summary_df['info'].to_frame()
a['type']='AS_genes'
df=a
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('distribution_of_the_percentage_of_AS_genes_across_tissues.jpg',dpi=300)
plt.close()

summary=[]
for key,value in tissues_info.items(): 
    summary.append([key,value[10]]) 

pd.DataFrame(summary)
    
adult_tissues=['Brain__cortex__Frontal','Kidney','Longissimus_dorsi_ribeye_loin','Spleen']
feta_tissues=['Fetal-Brain','Fetal-Kidney','Fetal-Muscle','Fetal-Spleen']

def get_adult_fetal_shared_genes(adult_t,fetal_t):
    global tissues_info
    ad=tissues_info[adult_t][2].copy().rename({'total_AS_event':'adult_AS_event'},axis=1) 
    fet=tissues_info[fetal_t][2].copy().rename({'total_AS_event':'fetal_AS_event'},axis=1)  
    m=ad.merge(fet)
    return(m)

df=get_adult_fetal_shared_genes('Brain__cortex__Frontal','Fetal-Brain')
spearmanr(df['adult_AS_event'],df['fetal_AS_event'])
b=df['fetal_AS_event'].to_frame().rename({'fetal_AS_event':"info"},axis=1)
b['type1']='brain'
b['type2']='fetal'
a=df['adult_AS_event'].to_frame().rename({'adult_AS_event':"info"},axis=1)
a['type1']='brain'
a['type2']='adult'
df=get_adult_fetal_shared_genes('Kidney','Fetal-Kidney')
spearmanr(df['adult_AS_event'],df['fetal_AS_event'])
d=df['fetal_AS_event'].to_frame().rename({'fetal_AS_event':"info"},axis=1)
d['type1']='kidney'
d['type2']='fetal'
c=df['adult_AS_event'].to_frame().rename({'adult_AS_event':"info"},axis=1)
c['type1']='kidney'
c['type2']='adult'
df=get_adult_fetal_shared_genes('Longissimus_dorsi_ribeye_loin','Fetal-Muscle')
spearmanr(df['adult_AS_event'],df['fetal_AS_event'])
f=df['fetal_AS_event'].to_frame().rename({'fetal_AS_event':"info"},axis=1)
f['type1']='LD muscle'
f['type2']='fetal'
e=df['adult_AS_event'].to_frame().rename({'adult_AS_event':"info"},axis=1)
e['type1']='LD muscle'
e['type2']='adult'
df=get_adult_fetal_shared_genes('Spleen','Fetal-Spleen')
spearmanr(df['adult_AS_event'],df['fetal_AS_event'])
h=df['fetal_AS_event'].to_frame().rename({'fetal_AS_event':"info"},axis=1)
h['type1']='Spleen'
h['type2']='fetal'
g=df['adult_AS_event'].to_frame().rename({'adult_AS_event':"info"},axis=1)
g['type1']='Spleen'
g['type2']='adult'
df=pd.concat([b,a,d,c,f,e,h,g])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('adult_vs_fetal_tissues_AS_events_per_gene.png',dpi=300)
plt.close()


A3_trs=pd.read_csv('A3_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
A5_trs=pd.read_csv('A5_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
AF_trs=pd.read_csv('AF_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
AL_trs=pd.read_csv('AL_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
MX_trs=pd.read_csv('MX_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
RI_trs=pd.read_csv('RI_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
SE_trs=pd.read_csv('SE_transcripts',header=0,sep='\t').rename({'alternative_transcripts':'transcript_id'},axis=1)
USE_trs=pd.read_csv('USEs_transcripts',header=0,names=['transcript_id'],sep='\t')

expressions=expressions.set_index('transcript_id')
info=get_info(expressions,'transcript_id')
A3_info=A3_trs.merge(info)
A5_info=A5_trs.merge(info)
AL_info=AL_trs.merge(info)
AF_info=AF_trs.merge(info)
MX_info=MX_trs.merge(info)
RI_info=RI_trs.merge(info)
SE_info=SE_trs.merge(info)
USE_info=USE_trs.merge(info)

AS_tr_info=info.loc[info.transcript_id.isin(event_per_tr.transcript_id)]
non_AS_tr_info=info.loc[~info.transcript_id.isin(event_per_tr.transcript_id)]

AS_tr_info['detected_tissues'].median()
AS_tr_info['median_expression'].median()
non_AS_tr_info['detected_tissues'].median()
non_AS_tr_info['median_expression'].median()

A3_info.to_csv('A3_transcripts_detection_expression_info',index = None, header=True,sep="\t")
A5_info.to_csv('A5_transcripts_detection_expression_info',index = None, header=True,sep="\t")
AL_info.to_csv('AL_transcripts_detection_expression_info',index = None, header=True,sep="\t")
AF_info.to_csv('AF_transcripts_detection_expression_info',index = None, header=True,sep="\t")
MX_info.to_csv('MX_transcripts_detection_expression_info',index = None, header=True,sep="\t")
RI_info.to_csv('RI_transcripts_detection_expression_info',index = None, header=True,sep="\t")
SE_info.to_csv('SE_transcripts_detection_expression_info',index = None, header=True,sep="\t")
USE_info.to_csv('USE_transcripts_detection_expression_info',index = None, header=True,sep="\t")


a=A3_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='A3'
a['type2']='detected_tissues'
b=A3_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
b['type1']='A3'
b['type2']='median_expression'
c=A5_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='A5'
c['type2']='detected_tissues'
d=A5_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
d['type1']='A5'
d['type2']='median_expression'
e=AF_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
e['type1']='AF'
e['type2']='detected_tissues'
f=AF_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
f['type1']='AF'
f['type2']='median_expression'
g=AL_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
g['type1']='AL'
g['type2']='detected_tissues'
h=AL_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
h['type1']='AL'
h['type2']='median_expression'
j=MX_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
j['type1']='MX'
j['type2']='detected_tissues'
k=MX_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
k['type1']='MX'
k['type2']='median_expression'
l=RI_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
l['type1']='RI'
l['type2']='detected_tissues'
m=RI_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
m['type1']='RI'
m['type2']='median_expression'
n=SE_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
n['type1']='SE'
n['type2']='detected_tissues'
o=SE_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
o['type1']='SE'
o['type2']='median_expression'
p=USE_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
p['type1']='USE'
p['type2']='detected_tissues'
q=USE_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
q['type1']='USE'
q['type2']='median_expression'

df=pd.concat([a,c,e,g,j,l,n,p])
ax = sns.boxplot(x="type1", y="info", data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('AS_events_transcripts_detected_tissues.png',dpi=300)
plt.close()

df=pd.concat([b,d,f,h,k,m,o,q])
df['log']=np.log2(df['info'])
ax = sns.boxplot(x="type1", y="log", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('AS_events_transcripts_expression.png',dpi=300)
plt.close()




'''
summary=[]
for key,value in tissues_info.items(): 
    df=value[0].copy()
    df["sum"] = df.loc[:,"coding_to_coding":"noncoding_to_noncoding"].sum(axis=1)
    l=list(df["sum"])
    l.append(key)
    summary.append(l)
summary_df=pd.DataFrame(summary)
summary_df=pd.DataFrame(summary).rename({0:'A3',1:'A5',2:'AF',3:'AL',4:'MX',5:'RI',6:'SE',7:'EIs',8:'Tissue'},axis=1)

'''






