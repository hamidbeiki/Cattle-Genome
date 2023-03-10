#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 12:05:32 2021

@author: beiki
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
#from scipy.stats import median_test

def get_info(ex,biotype):
    ex = ex.replace(0, np.NaN)
    ex_num=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy on python 3
    row_medians = np.nanmedian(ex_num, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
    row_counts=np.count_nonzero(~np.isnan(ex_num),axis=1)#https://stackoverflow.com/questions/21778118/counting-the-number-of-non-nan-elements-in-a-numpy-ndarray-in-python
    df1=pd.DataFrame(row_counts).rename({0:'detected_samples'},axis=1)
    df2=pd.DataFrame(row_medians).rename({0:'median_expression'},axis=1)
    df1[biotype]=ex.index
    df2[biotype]=ex.index
    m=df1.merge(df2)
    m=m[[biotype,'detected_samples','median_expression']]
    return(m)
    
gene_expressions=pd.read_csv("RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t")#,index_col='gene_id'
known=pd.read_csv('all-known-genes-main',names=['gene_id'],sep='\t')
gene_biotypes=pd.read_csv('final_gene_biotypes',names=['gene_id','gene_biotype'],sep='\t')
orphans=pd.read_csv('orphan-genes',names=['gene_id'],sep='\t')
tissues=list(gene_expressions.columns)[1:]
remove=pd.read_csv('remove_genes',header=0,names=['gene_id'])
gene_expressions=gene_expressions.loc[~ gene_expressions.gene_id.isin(remove.gene_id)]
gene_expressions=gene_expressions.set_index('gene_id')
gene_expressions=gene_expressions.loc[(gene_expressions!=0).any(1)]#remove genes that couldn't be quantified
known=known.loc[~ known.gene_id.isin(remove.gene_id)]
novel=gene_expressions.loc[~gene_expressions.index.isin(known.gene_id)].reset_index()['gene_id'].to_frame()

tissues_info=[]
for tissue in tissues:
    df=gene_expressions[[tissue]].reset_index()
    df=df.loc[(df[tissue]>0)]
    k=df.merge(known)
    n=df.merge(novel)
    tissues_info.append([tissue,k.shape[0],n.shape[0]])

tissues_info=pd.DataFrame(tissues_info).rename({0:'Tissue',1:'known',2:'novel'},axis=1)
a=tissues_info['known'].to_frame().rename({'known':'info'},axis=1)
a['type']='known'
b=tissues_info['novel'].to_frame().rename({'novel':'info'},axis=1)
b['type']='novel'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_known-novel_gene_number.png',dpi=300)
plt.close()    


coding=gene_biotypes.loc[gene_biotypes['gene_biotype']=='protein-coding']['gene_id'].to_frame()
noncoding=gene_biotypes.loc[gene_biotypes['gene_biotype']=='non-coding']['gene_id'].to_frame()
pseudogenes=gene_biotypes.loc[gene_biotypes['gene_biotype']=='pseudogenes']['gene_id'].to_frame()

coding_ex=gene_expressions.loc[gene_expressions.index.isin(coding.gene_id)]
noncoding_ex=gene_expressions.loc[gene_expressions.index.isin(noncoding.gene_id)]
pseudogenes_ex=gene_expressions.loc[gene_expressions.index.isin(pseudogenes.gene_id)]

coding_info=get_info(coding_ex,'genes')
noncoding_info=get_info(noncoding_ex,'genes')
pseudogenes_info=get_info(pseudogenes_ex,'genes')

coding_info['log']=np.log2(coding_info['median_expression'])
noncoding_info['log']=np.log2(noncoding_info['median_expression'])
pseudogenes_info['log']=np.log2(pseudogenes_info['median_expression'])

pseudogenes_info['detected_tissues'].median()
coding_info['detected_tissues'].median()
noncoding_info['detected_tissues'].median()

pseudogenes_info['median_expression'].median()
coding_info['median_expression'].median()
noncoding_info['median_expression'].median()

types=[coding_info,noncoding_info,pseudogenes_info]
res1=[]
res2=[]
for i in types:
    r1=[]
    r2=[]
    for j in types:
        stat,pval1=stats.mannwhitneyu(i['detected_tissues'],j['detected_tissues'])#Mann-Whitney test that is equivalent to Wilcoxon Rank Sum test with sample sizes are uequivalent: https://stats.stackexchange.com/questions/368881/wilcoxon-rank-sum-test-unequal-sample-sizes
        stat,pval2=stats.mannwhitneyu(i['median_expression'],j['median_expression'])
        if pval1<0.05 and i['detected_tissues'].median()>j['detected_tissues'].median():
            r1.append(1)
        elif pval1<0.05 and i['detected_tissues'].median()<j['detected_tissues'].median():
            r1.append(-1)
        else:
            r1.append(0)
        if pval2<0.05 and i['median_expression'].median()>j['median_expression'].median():
            r2.append(1)
        elif pval2<0.05 and i['median_expression'].median()<j['median_expression'].median():
            r2.append(-1)
        else:
            r2.append(0)
    res1.append(r1)
    res2.append(r2)

types=['coding','noncoding','pseudogenes']        
df1=pd.DataFrame(res1)
df2=pd.DataFrame(res2)
df1.columns=types
df1.index=types   
df2.columns=types
df2.index=types           
df1.to_csv("results1.txt",index = True, header=True,sep="\t") 
df2.to_csv("results2.txt",index = True, header=True,sep="\t") 


a=coding_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='coding'
a['type2']='detected_tissues'
b=coding_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='coding'
b['type2']='median_expression'
c=noncoding_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='noncoding'
c['type2']='detected_tissues'
d=noncoding_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='noncoding'
d['type2']='median_expression'
e=pseudogenes_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
e['type1']='pseudogenes'
e['type2']='detected_tissues'
f=pseudogenes_info['log'].to_frame().rename({'log':"info"},axis=1)
f['type1']='pseudogenes'
f['type2']='median_expression'
df=pd.concat([a,c,e])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
df=pd.concat([b,d,f])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('test.png',dpi=300)
plt.close() 

'''
    pseudo_lncRNAs
'''
tr_expressions=pd.read_csv("RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
remove=pd.read_csv('remove_transcripts',header=0,names=['transcript_id'])
lncRNAs=pd.read_csv('lncRNAs',names=['transcript_id'])
pseudo_lncRNAs=pd.read_csv('pseudogene_derived_lncRNAs',names=['transcript_id'])
lncRNAs=lncRNAs.loc[~lncRNAs.transcript_id.isin(pseudo_lncRNAs.transcript_id)]

tr_expressions=tr_expressions.loc[~ tr_expressions.transcript_id.isin(remove.transcript_id)]
tr_expressions=tr_expressions.set_index('transcript_id')
tr_expressions=tr_expressions.loc[(tr_expressions!=0).any(1)]#remove transcripts that couldn't be quantified

lncRNAs_ex=tr_expressions.loc[tr_expressions.index.isin(lncRNAs.transcript_id)]
pseudo_lncRNAs_ex=tr_expressions.loc[tr_expressions.index.isin(pseudo_lncRNAs.transcript_id)]

lncRNAs_info=get_info(lncRNAs_ex,'transcripts')
pseudo_lncRNAs_info=get_info(pseudo_lncRNAs_ex,'transcripts')

lncRNAs_info['median_expression'].median()
pseudo_lncRNAs_info['median_expression'].median()
lncRNAs_info['detected_tissues'].median()
pseudo_lncRNAs_info['detected_tissues'].median()

lncRNAs_info['log']=np.log2(lncRNAs_info['median_expression'])
pseudo_lncRNAs_info['log']=np.log2(pseudo_lncRNAs_info['median_expression'])


a=lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='rest_lncRNAs'
a['type2']='detected_tissues'
b=lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='rest_lncRNAs'
b['type2']='median_expression'
c=pseudo_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='pseudo_lncRNAs'
c['type2']='detected_tissues'
d=pseudo_lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='pseudo_lncRNAs'
d['type2']='median_expression'
df=pd.concat([a,c])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
df=pd.concat([b,d])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('test.png',dpi=300)
plt.close() 

expressions=pseudo_lncRNAs_ex.reset_index()
results=[]
for tissue in tissues:
    df=expressions[['transcript_id',tissue]]
    df=df.loc[(df[tissue]>0)]
    results.append([tissue,df.shape[0]])
results_df=pd.DataFrame(results).rename({0:'tissue',1:"pseudo_lncRNAs"},axis=1)     

expressions=gene_expressions.reset_index()  
results=[]
for tissue in tissues:
    df=expressions[['gene_id',tissue]]
    m=df.merge(gene_biotypes)
    df=m
    df=df.loc[(df[tissue]>0)]
    results.append([tissue,df.loc[df['gene_biotype']=="protein-coding"].shape[0],df.loc[df['gene_biotype']=="non-coding"].shape[0],df.loc[df['gene_biotype']=="pseudogenes"].shape[0]])
results_df=pd.DataFrame(results).rename({0:'tissue',1:"protein-coding",2:"non-coding",3:"pseudogene"},axis=1)     

a=results_df['protein-coding'].to_frame().rename({'protein-coding':'info'},axis=1)
a['type']='protein-coding'
b=results_df['non-coding'].to_frame().rename({'non-coding':'info'},axis=1)
b['type']='non-coding'
c=results_df['pseudogene'].to_frame().rename({'pseudogene':'info'},axis=1)
c['type']='pseudogene'
df=pd.concat([a,b,c])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_gene_number_biotype.png',dpi=300)
plt.close()
results_df['sum']=results_df.loc[:,'protein-coding':'pseudogene'].sum(axis=1)
results_df.loc[:,'protein-coding':'pseudogene'].div(results_df["sum"], axis=0)
results_df2=results_df.loc[:,'protein-coding':'pseudogene'].div(results_df["sum"], axis=0)
results_df2['tissue']=results_df['tissue']
a=results_df2['protein-coding'].to_frame().rename({'protein-coding':'info'},axis=1)
a['type']='protein-coding'
b=results_df2['non-coding'].to_frame().rename({'non-coding':'info'},axis=1)
b['type']='non-coding'
c=results_df2['pseudogene'].to_frame().rename({'pseudogene':'info'},axis=1)
c['type']='pseudogene'
df=pd.concat([a,b,c])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_gene_percentage_biotype.png',dpi=300)
plt.close()
df=results_df['sum'].to_frame().rename({'sum':'info'},axis=1)
df['type']='genes'
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_gene_number.png',dpi=300)
plt.close()

results=[]
for tissue in tissues:
    df=expressions[['gene_id',tissue]]
    df=df.loc[(df[tissue]>0)]
    known_g=known.merge(df)
    novel_g=df.loc[~ df.gene_id.isin(known_g.gene_id)]
#    stat, pval, med, tbl = median_test(known_g[tissue],novel_g[tissue]) 
    stat,pval=stats.mannwhitneyu(known_g[tissue],novel_g[tissue],alternative='greater') 
    results.append([tissue,known_g[tissue].median(),novel_g[tissue].median(),pval])
df=pd.DataFrame(results).rename({0:'tissue',1:"known",2:"novel",3:"p-value"},axis=1)     
df['-log2(p_value)']=(-np.log2(df["p-value"]+1e-306)) 
df=df.reset_index()
df.plot(kind="scatter", x="known", y="novel", alpha=0.4,
    figsize=(10,7),
    c="index", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,s=df["-log2(p_value)"])#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.show()  
   
ll = plt.scatter(df["known"], df["novel"], marker='o', c=df["-log2(p_value)"],cmap=plt.get_cmap("jet"))#, s=df["-log2(p_value)"]
plt.colorbar(label="-log2(p-value) for the comparision of expression values" "\n" "between known and novel genes in each tissue" "\n")
plt.plot(df["known"], df["known"],linestyle = 'dashed',color = 'grey')
plt.xlabel("median expression of known genes (RPKM)")
plt.ylabel("median expression of novel genes (RPKM)")
#plt.show()
plt.savefig('known_novel_gene_expression_comparision.png')
plt.close()    
    
known_ex=gene_expressions.loc[gene_expressions.index.isin(known.gene_id)]
novel_ex=gene_expressions.loc[~ gene_expressions.index.isin(known.gene_id)]    

known_info=get_info(known_ex,'genes')    
novel_info=get_info(novel_ex,'genes')     

known_info['log']=np.log2(known_info['median_expression']) 
novel_info['log']=np.log2(novel_info['median_expression'])    
   
    
a=known_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='known'
a['type2']='detected_tissues'
b=known_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='known'
b['type2']='median_expression'
c=novel_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='novel'
c['type2']='detected_tissues'
d=novel_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='novel'
d['type2']='median_expression'
df=pd.concat([a,c])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
df=pd.concat([b,d])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
#plt.show()    
plt.savefig('test.png',dpi=300)
plt.close()    
    
stat, p = stats.mannwhitneyu(known_info['median_expression'],novel_info['median_expression'],alternative='greater')
stat, p = stats.mannwhitneyu(known_info['detected_tissues'],novel_info['detected_tissues'],alternative='greater')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    