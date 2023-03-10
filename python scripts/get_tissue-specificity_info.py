#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 21:18:32 2021

@author: beiki
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from scipy.stats import gaussian_kde

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
            
def get_tissue_other_A5_A3_event_per_tr(events,expr_d):
    events_expr=events.loc[(events.event_transcript.isin(expr_d.transcript_id)) & (events.refrence_transcript.isin(expr_d.transcript_id))]
    other_events_per_trs=events_expr.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)
    return(other_events_per_trs)

    
expressions=pd.read_csv("RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
gene_expressions=pd.read_csv("RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t")
tr_biotypes=pd.read_csv("../final_transcript_biotypes",names=['transcript_id','biotype'],sep='\t')
gene_biotypes=pd.read_csv("../final_gene_biotypes",names=['gene_id','biotype'],sep='\t')
remove_genes=pd.read_csv('../remove_genes',names=['gene_id'])
remove_trs=pd.read_csv('../remove_transcripts',names=['transcript_id'])
genes_TSI=pd.read_csv("Combined_genes-tissue-specificity-scores.txt",header=0,sep="\t")
tr_TSI=pd.read_csv("Combined_transcripts-tissue-specificity-scores.txt",header=0,sep="\t")
ensembl=pd.read_csv('../annotation/combined_genes_ensembl_equivalent',names=['gene_id','ens'],sep=' ')
known_genes=pd.read_csv('../annotation/all-known-genes-main',names=['gene_id'])
known_trs=pd.read_csv('../annotation/all-known-transcripts-main',names=['transcript_id'])
AS_event_trs=pd.read_csv('../AS-events/total_AS_event_per_transcript',header=0,sep='\t')
AS_event_genes=pd.read_csv('../AS-events/total_AS_event_per_gene',header=0,sep='\t')
A3=pd.read_csv('../AS-events/combined_SUPPA_output_A3_strict.ioe',header=0, sep='\t')
A5=pd.read_csv('../AS-events/combined_SUPPA_output_A5_strict.ioe',header=0, sep='\t')
AF=pd.read_csv('../AS-events/combined_SUPPA_output_AF_strict.ioe',header=0, sep='\t')
AL=pd.read_csv('../AS-events/combined_SUPPA_output_AL_strict.ioe',header=0, sep='\t')
MX=pd.read_csv('../AS-events/combined_SUPPA_output_MX_strict.ioe',header=0, sep='\t')
RI=pd.read_csv('../AS-events/combined_SUPPA_output_RI_strict.ioe',header=0, sep='\t')
SE=pd.read_csv('../AS-events/combined_SUPPA_output_SE_strict.ioe',header=0, sep='\t')
USE=pd.read_csv('../AS-events/combined_USE_events',header=0,sep='\t')
EIs_df=pd.read_csv('../AS-events/trs_with_exons_covered_by_other_trs_introns',names=['transcript_id'])    
intersect=pd.read_csv('../AS-events/intersect',sep='\t',names=['ids',1,2,3,'transcript_id',5,6,7,'p1',9,10,11,12,13,14,15,16,'p2',18])# bedtools intersect -a combined_exons.bed -b combined_introns.bed -wo -s -f 0.9999999 | tr ' ' '\t'>intersect
utr5=pd.read_csv("../number_of_5UTRs_per_gene",sep='\t',header=0)
utr3=pd.read_csv("../number_of_3UTRs_per_gene",sep='\t',header=0)
               
gene_expressions=gene_expressions.loc[~gene_expressions.gene_id.isin(remove_genes.gene_id)]
expressions=expressions.loc[~expressions.transcript_id.isin(remove_trs.transcript_id)]
genes_TSI=genes_TSI.loc[~genes_TSI.gene_id.isin(remove_genes.gene_id)]
tr_TSI=tr_TSI[~tr_TSI.transcript_id.isin(remove_trs.transcript_id)]
utr5=utr5.loc[~utr5.gene_id.isin(remove_genes.gene_id)]
utr3=utr3.loc[~utr3.gene_id.isin(remove_genes.gene_id)]


df=expressions.iloc[:,1:]
df=df.loc[(df!=0).any(1)]
ids=expressions['transcript_id']
expressions=df
expressions['transcript_id']=ids
df=gene_expressions.iloc[:,1:]
df=df.loc[(df!=0).any(1)]
gene_ids=gene_expressions['gene_id']
gene_expressions=df
gene_expressions['gene_id']=gene_ids

genes_TSI=genes_TSI.loc[genes_TSI.gene_id.isin(gene_expressions.gene_id)]

gene_expressions=gene_expressions.set_index('gene_id')
gene_info=get_info(gene_expressions,'gene_id')
expressions=expressions.set_index('transcript_id')
tr_info=get_info(expressions,'transcript_id')

plt.hist(gene_info['detected_tissues'], bins=50,color='blue')
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of genes')
plt.savefig('hist_of_genes_detected_tissues.png')
plt.close()

plt.hist(tr_info['detected_tissues'], bins=50,color='blue')
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of transcripts')
plt.savefig('hist_of_transcripts_detected_tissues.png')
plt.close()

TS_genes=gene_info.loc[gene_info['detected_tissues']==1]
TS_genes=TS_genes.merge(gene_biotypes)
TS_transcripts=tr_info.loc[tr_info['detected_tissues']==1]
TS_transcripts=TS_transcripts.merge(tr_biotypes)

TS_genes.groupby(['biotype']).size().reset_index()
TS_transcripts.groupby(['biotype']).size().reset_index()

stats.spearmanr(tr_info['detected_tissues'],tr_info['median_expression'])
stats.spearmanr(gene_info['detected_tissues'],gene_info['median_expression'])

gene_expressions_2=gene_expressions.reset_index()
expressions_2=expressions.reset_index()
tissues=list(expressions.columns)[0:-1]
tissues_info=[]
for tissue in tissues:
    df=expressions_2[['transcript_id',tissue]]
    expr_tr=df.loc[df[tissue]>0]
    ts_tr=expr_tr.merge(TS_transcripts)
    df=gene_expressions_2[['gene_id',tissue]]
    expr_g=df.loc[df[tissue]>0]
    ts_g=expr_g.merge(TS_genes)
    tissues_info.append([tissue,ts_g.shape[0],ts_tr.shape[0]])

df=pd.DataFrame(tissues_info).rename({0:'Tissue',1:'TS_genes',2:'TS_transcripts'},axis=1)   

a=df['TS_genes'].to_frame().rename({'TS_genes':'info'},axis=1)
a['type']='tissue-specific genes'
b=df['TS_transcripts'].to_frame().rename({'TS_transcripts':'info'},axis=1)
b['type']='tissue-specific transcripts'
df2=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df2, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.legend([],[], frameon=False)
plt.savefig('tissue-specific_gene_transcript.png',dpi=300)
plt.close()

tissues_info=[]
for tissue in tissues:
    df=expressions_2[['transcript_id',tissue]]
    expr_tr=df.loc[df[tissue]>0]
    ts_tr=expr_tr.merge(TS_transcripts)
    rest_tr=expr_tr.loc[~expr_tr.transcript_id.isin(ts_tr.transcript_id)]
    df=gene_expressions_2[['gene_id',tissue]]
    expr_g=df.loc[df[tissue]>0]
    ts_g=expr_g.merge(TS_genes)
    rest_g=expr_g.loc[~expr_g.gene_id.isin(ts_g.gene_id)]
    tissues_info.append([tissue,rest_g[tissue].median(),ts_g[tissue].median(),rest_tr[tissue].median(),ts_tr[tissue].median()])

df=pd.DataFrame(tissues_info).rename({0:'Tissue',1:'non_TS_genes',2:'TS_genes',3:'non_TS_transcripts',4:'TS_transcripts'},axis=1)   
df=df.loc[df['Tissue']!='Gall_Bladder']# because it has no TS gene
df['non_TS_genes_log']=np.log2(df['non_TS_genes'])
df['TS_genes_log']=np.log2(df['TS_genes'])
df['non_TS_transcripts_log']=np.log2(df['non_TS_transcripts'])
df['TS_transcripts_log']=np.log2(df['TS_transcripts'])

a=df['non_TS_genes_log'].to_frame().rename({'non_TS_genes_log':'info'},axis=1)
a['type']='non_TS_genes'
b=df['TS_genes_log'].to_frame().rename({'TS_genes_log':'info'},axis=1)
b['type']='TS_genes'
c=df['non_TS_transcripts_log'].to_frame().rename({'non_TS_transcripts_log':'info'},axis=1)
c['type']='non_TS_transcripts'
d=df['TS_transcripts_log'].to_frame().rename({'TS_transcripts_log':'info'},axis=1)
d['type']='TS_transcripts'
df2=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type", y="info", data=df2, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue-specific_gene_transcript_expression.png',dpi=300)
plt.close()

l=TS_transcripts['transcript_id'].str.split('.').to_list()
genes=[]
for i in l:
    g=str(i[0])+'.'+str(i[1])
    genes.append(g)
gene_ids=pd.DataFrame(genes).rename({0:'gene_id'},axis=1)#genes with TS transcripts
TS_transcripts['gene_id']=gene_ids
TS_transcripts=TS_transcripts.rename({'biotype':'TS_tr_biotype'},axis=1)
TS_transcripts=TS_transcripts.merge(gene_biotypes)
TS_transcripts=TS_transcripts.rename({'biotype':'gene_biotype'},axis=1)
genes_with_TS_trs=gene_ids.drop_duplicates()
genes_with_TS_trs.merge(gene_biotypes).groupby(['biotype']).size().reset_index()
TS_transcripts.loc[(TS_transcripts['gene_biotype']=="protein-coding")]
TS_transcripts.loc[(TS_transcripts['TS_tr_biotype']!="protein_coding_transcripts") & (TS_transcripts['gene_biotype']=="protein-coding")]
d=TS_transcripts.groupby(['TS_tr_biotype','gene_biotype']).size().reset_index()
d.loc[(d['TS_tr_biotype']!='protein_coding_transcripts') & (d['gene_biotype']=='protein-coding') ][0].sum()

TS_genes.to_csv('tissue-specific_genes_info',index = None, header=True,sep="\t")
TS_transcripts.to_csv('tissue-specific_transcripts_info',index = None, header=True,sep="\t")
genes_with_TS_trs.to_csv('genes_transcribed_tissue-specific_transcripts',index = None, header=False,sep="\t")

m=EIs_df.merge(intersect)
groups=m.groupby(['ids'])['p1'].min().reset_index()
other_A5=groups.loc[groups['p1']==0]
other_A3=groups.loc[groups['p1']!=0]
other_A5_events = other_A5['ids'].str.split("--", n = 1, expand = True).rename({0:'event_transcript',1:'refrence_transcript'},axis=1)
other_A3_events = other_A3['ids'].str.split("--", n = 1, expand = True).rename({0:'event_transcript',1:'refrence_transcript'},axis=1)
other_A3_events_per_trs=other_A3_events.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)
other_A5_events_per_trs=other_A5_events.groupby(['event_transcript']).size().reset_index().rename({0:'event','event_transcript':'transcript_id'},axis=1)

tissues=list(expressions.columns)[0:-1]
tissues_info=[]
tissues_info2=[]#percentages
for tissue in tissues:
    if tissue!='Gall_Bladder':# as it has 0 TS gene and only 7 TS transcripts
        df=expressions_2[['transcript_id',tissue]]
        expr=df.loc[df[tissue]>0]
        expr=list(expr['transcript_id'])
        expr_d=pd.DataFrame(expr).rename({0:'transcript_id'},axis=1)
        df=gene_expressions_2[['gene_id',tissue]]
        expr2_d=df.loc[df[tissue]>0]
        expr2=list(expr2_d['gene_id'])
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
        TS_g_event_per_g=TS_genes.merge(event_per_gene)
        m1=TS_genes.merge(expr2_d)
        TS_t_event_per_t=TS_transcripts.merge(event_per_tr)
        m2=TS_transcripts.merge(expr_d)
        tissues_info.append([tissue,m1.shape[0],TS_g_event_per_g.shape[0],m2.shape[0],TS_t_event_per_t.shape[0]])
        p1=(TS_g_event_per_g.shape[0]/m1.shape[0])*100
        p2=(TS_t_event_per_t.shape[0]/m2.shape[0])*100
        tissues_info2.append([tissue,p1,p2])
        print(tissue)
tissue_TS_AS=pd.DataFrame(tissues_info2).rename({0:'Tissue',1:'%of_TS_genes_involve_AS',2:'%of_TS_transcripts_involve_AS'},axis=1)    
a=tissue_TS_AS['%of_TS_genes_involve_AS'].to_frame().rename({'%of_TS_genes_involve_AS':'info'},axis=1)
a['type']='TS_genes'
b=tissue_TS_AS['%of_TS_transcripts_involve_AS'].to_frame().rename({'%of_TS_transcripts_involve_AS':'info'},axis=1)
b['type']='TS_transcripts'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue-specific_genes_transcripts_AS_events.png',dpi=300)
plt.close()
'''
TSI based analysis
'''
a=genes_TSI.loc[(genes_TSI['TSI']>0.9) & (genes_TSI['TSI']<1)]
a.merge(gene_info)['detected_tissues'].median()
gene_info.loc[~gene_info.gene_id.isin(a.gene_id)]['detected_tissues'].median()
b=tr_TSI.loc[(tr_TSI['TSI']>0.9) & (tr_TSI['TSI']<1)]
b.merge(tr_info)['detected_tissues'].median()


tr_TSI_2=tr_TSI.loc[tr_TSI['TSI']<1]
tr_TSI_2['TSI'].quantile([0.0, .5, .90, .95])
genes_TSI_2=genes_TSI[genes_TSI['TSI']<1]
genes_TSI_2.quantile([0.0, .5, .90, .95])


plt.hist(genes_TSI_2['TSI'], bins=460,color='black')
plt.xlabel('Tissue Specificity Score')
plt.ylabel('Number of multi-tissue detected genes')
plt.axvline( 0.986882, linestyle='--',color='red')
plt.savefig('hist_of_multi-tissue_genes_tissue-specificity_scores.png',dpi=300)
plt.close()

plt.hist(tr_TSI_2['TSI'], bins=460,color='black')
plt.xlabel('Tissue Specificity Score')
plt.ylabel('Number of multi-tissue detected transcripts')
plt.axvline(0.982086, linestyle='--',color='red')
plt.savefig('hist_of_multi-tissue_transcripts_tissue-specificity_scores.png',dpi=300)
plt.close()


TSI_based_TS_genes=genes_TSI_2.loc[(genes_TSI_2['TSI']>0.986882) & (genes_TSI_2['TSI']<1)]
TSI_based_TS_genes=TSI_based_TS_genes.merge(gene_info)
out=TSI_based_TS_genes.merge(ensembl)['ens'].drop_duplicates().to_frame()
out.to_csv('genes_with_highest_expression_variability_across_tissues',index = None, header=False,sep="\t")

df=genes_TSI_2.merge(gene_info)
x=df["detected_tissues"]
y=df["TSI"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Number of detected tissues")
plt.ylabel("Tissue Specificity Score of multi-tissue detected genes")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_detected_tissues_in_multi-tissue_genes.png',dpi=300)
plt.close()


df=tr_TSI_2.merge(tr_info)
x=df["detected_tissues"]
y=df["TSI"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Number of detected tissues")
plt.ylabel("Tissue Specificity Score of multi-tissue detected transcripts")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_detected_tissues_in_multi-tissue_transcripts.png',dpi=300)
plt.close()





TSI_based_TS_trs=tr_TSI_2.loc[(tr_TSI_2['TSI']>0.982086) & (tr_TSI_2['TSI']<1)]
TSI_based_TS_trs=TSI_based_TS_trs.merge(tr_info)

TSI_based_TS_genes.merge(gene_biotypes).groupby(['biotype']).size().reset_index()
TSI_based_TS_trs.merge(tr_biotypes).groupby(['biotype']).size().reset_index()


TSI_based_TS_genes.merge(AS_event_genes)['total_AS_event'].median()
AS_event_genes.loc[~AS_event_genes.gene_id.isin(TSI_based_TS_genes.gene_id)]['total_AS_event'].median()

TSI_based_TS_trs.merge(AS_event_trs)['total_AS_event'].median()
AS_event_trs.loc[~AS_event_trs.transcript_id.isin(TSI_based_TS_trs.transcript_id)]['total_AS_event'].median()





#tr_TSI=pd.read_csv("Combined_transcripts-tissue-specificity-scores.txt",header=0,sep="\t")
#tr_TSI=tr_TSI[~tr_TSI.transcript_id.isin(remove_trs.transcript_id)]
df=tr_TSI.merge(AS_event_trs)
df = df[df['TSI'].notna()]
#df=df.loc[df['TSI']<1]
x=df["TSI"]
y=df["total_AS_event"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Tissue specificity score")
plt.ylabel("Number of AS events")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_AS-events_in_transcripts.png',dpi=300)
plt.close()

#tr_TSI=pd.read_csv("Combined_transcripts-tissue-specificity-scores.txt",header=0,sep="\t")
#tr_TSI=tr_TSI[~tr_TSI.transcript_id.isin(remove_trs.transcript_id)]
df=genes_TSI_2.merge(AS_event_genes)
df = df[df['TSI'].notna()]
df = df.reset_index()
del df['index']
x=df["TSI"]
y=df["total_AS_event"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Tissue specificity score")
plt.ylabel("Number of AS events")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_AS-events_in_multi-tissue_genes.png',dpi=300)
plt.close()


#tr_TSI=pd.read_csv("Combined_transcripts-tissue-specificity-scores.txt",header=0,sep="\t")
#tr_TSI=tr_TSI[~tr_TSI.transcript_id.isin(remove_trs.transcript_id)]
df=tr_TSI_2.merge(AS_event_trs)
df = df[df['TSI'].notna()]
df = df.reset_index()
del df['index']
x=df["TSI"]
y=df["total_AS_event"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Tissue specificity score")
plt.ylabel("Number of AS events")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_AS-events_in_multi-tissue_transcripts.png',dpi=300)
plt.close()
stats.spearmanr(df['TSI'],tr_info['total_AS_event'])

from scipy.stats import spearmanr

df=utr5.merge(genes_TSI)
spearmanr(df['number_of_utr5'],df['TSI'])
df.loc[df['TSI']>0.9]['number_of_utr5'].median()
df.loc[df['TSI']<=0.9]['number_of_utr5'].median()

df2=utr3.merge(genes_TSI)
"""
df=genes_TSI_2.merge(gene_info)
df['median_expression_log']=np.log2(df['median_expression'])
ll = plt.scatter(df["detected_tissues"], df["TSI"], marker='o', c=df["median_expression_log"],cmap=plt.get_cmap("jet"),s=0.5,alpha=0.6)#, s=df["-log2(p_value)"]
plt.colorbar(label="median gene expression (log2(RPKM))")
plt.xlabel("Number of detected tissues")
plt.ylabel("Tissue-specificity score")
plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_detected_tissues_in_multi-tissue_genes.png',dpi=300)
plt.close()
"""

"""
    gene_percentage=(TS_g_event_per_g.shape[0]/(TS_genes.merge(expr2_d).shape[0]))*100# % of TS genes that involved in AS
    tr_percentage=(TS_t_event_per_t.shape[0]/(TS_transcripts.merge(expr_d).shape[0]))*100# % of TS transcripts that involved in AS
    rest_g_event_per_g=event_per_gene.loc[~event_per_gene.gene_id.isin(TS_genes.gene_id)]    
    rest_t_event_per_t=event_per_tr.loc[~event_per_tr.transcript_id.isin(TS_transcripts.transcript_id)]
    TS_t_A5_event_per_t=TS_transcripts.merge(A5_event_per_tr)
    TS_t_A3_event_per_t=TS_transcripts.merge(A3_event_per_tr)
    TS_t_AF_event_per_t=TS_transcripts.merge(AF_event_per_tr)
    TS_t_AL_event_per_t=TS_transcripts.merge(AL_event_per_tr)
    TS_t_MX_event_per_t=TS_transcripts.merge(MX_event_per_tr)
    TS_t_SE_event_per_t=TS_transcripts.merge(SE_event_per_tr)
    TS_t_RI_event_per_t=TS_transcripts.merge(RI_event_per_tr)
    TS_t_USE_event_per_t=TS_transcripts.merge(USE_event_per_tr)
    rest_t_A5_event_per_t=A5_event_per_tr.loc[~A5_event_per_tr.transcript_id.isin(TS_t_A5_event_per_t.transcript_id)]
    rest_t_A3_event_per_t=A3_event_per_tr.loc[~A3_event_per_tr.transcript_id.isin(TS_t_A3_event_per_t.transcript_id)]
    rest_t_AF_event_per_t=AF_event_per_tr.loc[~AF_event_per_tr.transcript_id.isin(TS_t_AF_event_per_t.transcript_id)]
    rest_t_AL_event_per_t=AL_event_per_tr.loc[~AL_event_per_tr.transcript_id.isin(TS_t_AL_event_per_t.transcript_id)]
    rest_t_MX_event_per_t=MX_event_per_tr.loc[~MX_event_per_tr.transcript_id.isin(TS_t_MX_event_per_t.transcript_id)]
    rest_t_RI_event_per_t=RI_event_per_tr.loc[~RI_event_per_tr.transcript_id.isin(TS_t_RI_event_per_t.transcript_id)]
    rest_t_SE_event_per_t=RI_event_per_tr.loc[~RI_event_per_tr.transcript_id.isin(TS_t_SE_event_per_t.transcript_id)]
    rest_t_USE_event_per_t=RI_event_per_tr.loc[~RI_event_per_tr.transcript_id.isin(TS_t_USE_event_per_t.transcript_id)]
    tissues_info=[tissue,gene_percentage,tr_percentage,TS_g_event_per_g.shape[0],TS_t_event_per_t.shape[0]]
"""