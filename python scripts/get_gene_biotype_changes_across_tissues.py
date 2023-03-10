#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 13:23:22 2021

@author: beiki
"""

#from scipy import stats
import pandas as pd
import numpy as np
from scipy.stats import median_test
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from functools import partial, reduce 
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

expressions=pd.read_csv("quantification/RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
gene_expressions=pd.read_csv("quantification/RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t",index_col='gene_id')
tissues=list(expressions.columns)[1:]
g=expressions['transcript_id'].str.split('.',n = 2, expand = True)[[0,1]]
genes=(g[0] + "." + g[1]).to_frame().rename({0:'gene_id'},axis=1)
expressions['gene_id']=genes
tr_biotypes=pd.read_csv("final_transcript_biotypes",names=['transcript_id','tr_biotype'],sep='\t')
gene_biotypes=pd.read_csv('final_gene_biotypes',names=['gene_id','gene_biotype'],sep='\t')
genes=genes.drop_duplicates()
gene_splice_info=pd.read_csv('number_of_uniq_introns_per_spliced_genes',names=['gene_id','number_of_splice'],sep='\t')
m=gene_splice_info.merge(genes,how='outer')
gene_splice_info=m.replace(np.NaN,0)
gene_splice_info.to_csv("number_of_uniq_introns_per_genes",index = False, header=True,sep="\t") 
utr3=pd.read_csv('number_of_3UTRs_per_gene',header=0,sep='\t')
utr5=pd.read_csv('number_of_5UTRs_per_gene',header=0,sep='\t')
ensembl=pd.read_csv('annotation/combined_genes_ensembl_equivalent',names=['gene_id','ENS'],sep=' ')
all_known=pd.read_csv('annotation/all-known-genes-main',names=['gene_id'])
all_known['novelty']='known'

results=[]
results2={}#switch_genes
results3={}#protein-coding genes encodes just coding transcripts
results4={}#protein-coding genes encodes mixture of coding and non-coding transcript
results5={}#non-coding genes
for tissue in tissues:
    df=expressions[['transcript_id','gene_id',tissue]]
    df=df.merge(tr_biotypes)
    df=df.merge(gene_biotypes)
    df=df.loc[(df[tissue]>0)]# & (df[tissue]<1)
    genes=set(list(df['gene_id'].drop_duplicates()))
    p=set(list(df.loc[df['gene_biotype']=="protein-coding"]['gene_id']))# protein coding genes, perhaps some of gene have protein coding transcripts detected in other tissue and just have non-coding transcripts in this tissue
    real_p=set(list(df.loc[df['tr_biotype']=="protein_coding_transcripts"]['gene_id']))# actual protein-coding genes in this tissue
    switched_genes=list(p - real_p) # genes switched from coding to non-coding in this tissue
    swithch_genes_df=pd.DataFrame(switched_genes).rename({0:'gene_id'},axis=1)
    swithch_genes_df[tissue]=1
    percentage=(len(switched_genes)/float(len(p)))*100
    real_p_df=pd.DataFrame(real_p).rename({0:"gene_id"},axis=1)
    m=real_p_df.merge(df)
    mix=set(list(m.loc[m['tr_biotype']!='protein_coding_transcripts']['gene_id'].drop_duplicates()))#genes with mixture of coding non-coding transcripts
    just_p=real_p - mix
    non_coding=genes - just_p - mix
    mix_df=pd.DataFrame(mix).rename({0:"gene_id"},axis=1)
    just_p_df=pd.DataFrame(just_p).rename({0:"gene_id"},axis=1)
    non_coding_df=pd.DataFrame(non_coding).rename({0:"gene_id"},axis=1)
    mix_df[tissue]=1
    just_p_df[tissue]=1
    non_coding_df[tissue]=1
    percentage2=(len(mix)/float(len(real_p)))*100
    percentage3=(len(just_p)/float(len(real_p)))*100
    results.append([tissue,percentage,len(switched_genes),percentage2,len(mix),percentage3,len(just_p)])
    results2[tissue]=swithch_genes_df
    results3[tissue]=just_p_df
    results4[tissue]=mix_df
    results5[tissue]=non_coding_df

results_df=pd.DataFrame(results).rename({0:'tissue',1:"%_of_coding-genes_that_were_swichched",2:"number_of_switched_genes",3:"%_of_mixed_genes",4:"#_of_mixed_genes",5:"%_of_just-p_genes",6:"#_of_just-p_genes"},axis=1)   
results_df.to_csv('tissues_coding-gene_biotype_details',index = None, header=True,sep="\t")
my_reduce = partial(pd.merge, on='gene_id', how='outer')
switched_genes=reduce(my_reduce, results2.values())
switched_genes_m=switched_genes.set_index('gene_id')
a=pd.DataFrame.to_numpy(switched_genes_m)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))#count number of switches from coding to non-coding per gene
d.index=switched_genes_m.index
switched_genes=d.reset_index().rename({'index':'gene_id',0:'number_of_tissues_gene_switched_from_coding_to_noncoding'},axis=1)


my_reduce = partial(pd.merge, on='gene_id', how='outer')
just_p=reduce(my_reduce, results3.values())
just_p_m=just_p.set_index('gene_id')
a=pd.DataFrame.to_numpy(just_p_m)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))#count number of switches from coding to non-coding per gene
d.index=just_p_m.index
just_p=d.reset_index().rename({'index':'gene_id',0:'number_of_tissues_gene_encode_just_coding_transcripts'},axis=1)
just_p['number_of_tissues_gene_encode_just_coding_transcripts'].median()
just_p['number_of_tissues_gene_encode_just_coding_transcripts'].quantile([0.0, .5, .90, .95])

my_reduce = partial(pd.merge, on='gene_id', how='outer')
mix=reduce(my_reduce, results4.values())
mix_m=mix.set_index('gene_id')
a=pd.DataFrame.to_numpy(mix_m)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))#count number of switches from coding to non-coding per gene
d.index=mix_m.index
mix=d.reset_index().rename({'index':'gene_id',0:'number_of_tissues_gene_encode_mix_transcripts'},axis=1)
mix['number_of_tissues_gene_encode_mix_transcripts'].median()
mix['number_of_tissues_gene_encode_mix_transcripts'].quantile([0.0, .5, .90, .95])

genes_info=get_info(gene_expressions,'gene_id')
mix=mix.merge(genes_info)
mix['%']=(mix['number_of_tissues_gene_encode_mix_transcripts']/mix['detected_tissues'])*100
x=mix["detected_tissues"]
y=mix["number_of_tissues_gene_encode_mix_transcripts"]
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last #https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
fig, ax = plt.subplots()
ax.scatter(x, y, marker='o', c=z, s=0.8,cmap=plt.get_cmap("jet"))
plt.xlabel("Number of detected tissues")
plt.ylabel("Number of tissues that a given gene transcribed both coding and non-coding transcripts")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('bifunctional_genes_info.png',dpi=300)
plt.close()


my_reduce = partial(pd.merge, on='gene_id', how='outer')
noncoding=reduce(my_reduce, results5.values())
noncoding_m=noncoding.set_index('gene_id')
a=pd.DataFrame.to_numpy(noncoding_m)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))#count number of switches from coding to non-coding per gene
d.index=noncoding_m.index
noncoding=d.reset_index().rename({'index':'gene_id',0:'number_of_tissues_gene_encode_just_noncoding_transcripts'},axis=1)
noncoding['number_of_tissues_gene_encode_just_noncoding_transcripts'].median()
noncoding['number_of_tissues_gene_encode_just_noncoding_transcripts'].quantile([0.0, .5, .90, .95])

data_frames = [just_p, mix, noncoding]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'], how='outer'), data_frames)
df_merged=df_merged.replace(np.NaN,0)#remaining genes coulnd'nt be quantified in any tissue
m=gene_biotypes.merge(df_merged)
m2=all_known.merge(m,on=['gene_id'],how='outer')
m2=m2.replace(np.NaN,'novel')
n=m['gene_id'].to_frame()
out=n.merge(m2,on=['gene_id'])##remaining genes coulnd'nt be quantified in any tissue

out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_mix_transcripts'].median()
out['total']=out['number_of_tissues_gene_encode_just_coding_transcripts']+out['number_of_tissues_gene_encode_mix_transcripts']+out['number_of_tissues_gene_encode_just_noncoding_transcripts']
out['%just-coding']=(out['number_of_tissues_gene_encode_just_coding_transcripts']/out['total'])*100
out['%mix']=(out['number_of_tissues_gene_encode_mix_transcripts']/out['total'])*100
out['%just-noncoding']=(out['number_of_tissues_gene_encode_just_noncoding_transcripts']/out['total'])*100
out=out.loc[out['gene_biotype']=='protein-coding']
out=out.loc[out['%just-noncoding']<100]
out.to_csv('coding-gene_biotype_changes_across_tissues',index = None, header=True,sep="\t")

switched_genes=switched_genes.loc[switched_genes.gene_id.isin(out.gene_id)]# because sume of transcripts couldn't be quantified in any tissue
switched_genes['number_of_tissues_gene_switched_from_coding_to_noncoding'].median()
switched_genes['number_of_tissues_gene_switched_from_coding_to_noncoding'].quantile([0.0, .5, .90, .95])
switched_genes.shape
###
# OR
###
switched_genes=out.loc[out['number_of_tissues_gene_encode_just_noncoding_transcripts']>0]

a=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['%just-coding'].to_frame().rename({'%just-coding':'percentage'},axis=1)
a['type1']='known'
a['type2']='just-coding'
b=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['%mix'].to_frame().rename({'%mix':'percentage'},axis=1)
b['type1']='known'
b['type2']='mix'
c=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['%just-noncoding'].to_frame().rename({'%just-noncoding':'percentage'},axis=1)
c['type1']='known'
c['type2']='just-noncoding'
d=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['%just-coding'].to_frame().rename({'%just-coding':'percentage'},axis=1)
d['type1']='novel'
d['type2']='just-coding'
e=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['%mix'].to_frame().rename({'%mix':'percentage'},axis=1)
e['type1']='novel'
e['type2']='mix'
f=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['%just-noncoding'].to_frame().rename({'%just-noncoding':'percentage'},axis=1)
f['type1']='novel'
f['type2']='just-noncoding'
df=pd.concat([a,b,c,d,e,f])
ax = sns.boxplot(x="type1", y="percentage", hue="type2", data=df, linewidth=2.5,showfliers=True)
plt.show()

a=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_just_coding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_coding_transcripts':'numbers'},axis=1)
a['type1']='known'
a['type2']='just-coding'
b=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_mix_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_mix_transcripts':'numbers'},axis=1)
b['type1']='known'
b['type2']='mix'
c=out.loc[(out['novelty']=="known") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_just_noncoding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_noncoding_transcripts':'numbers'},axis=1)
c['type1']='known'
c['type2']='just-noncoding'
d=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_just_coding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_coding_transcripts':'numbers'},axis=1)
d['type1']='novel'
d['type2']='just-coding'
e=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_mix_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_mix_transcripts':'numbers'},axis=1)
e['type1']='novel'
e['type2']='mix'
f=out.loc[(out['novelty']=="novel") & (out['gene_biotype']=="protein-coding")]['number_of_tissues_gene_encode_just_noncoding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_noncoding_transcripts':'numbers'},axis=1)
f['type1']='novel'
f['type2']='just-noncoding'
df=pd.concat([a,b,c,d,e,f])
ax = sns.boxplot(x="type1", y="numbers", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=1)
plt.show()
#plt.savefig('distribution_of_coding-gene_biotypes_A.png')

known=out.loc[out['novelty']=='known']
novel=out.loc[out['novelty']=='novel']
#lll = plt.scatter(switched_genes["%mix"], switched_genes["%just-coding"], marker='o', facecolors='none',edgecolors='r',s=14)#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
ll = plt.scatter(known["%mix"], known["%just-coding"], marker='o', c=known["%just-noncoding"],cmap=plt.get_cmap("jet"),s=2,alpha=0.4)#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
l = plt.scatter(novel["%mix"], novel["%just-coding"], marker='x', c=novel["%just-noncoding"],cmap=plt.get_cmap("jet"),s=12,alpha=0.4)#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
plt.colorbar(label="Percentage of detected tissues that a gene" "\n" "just transcribed non-coding transcripts")
plt.xlabel("Percentage of detected tissues that a gene" "\n" "transcribed mixture of coding and non-coding transcripts")
plt.ylabel("Percentage of detected tissues that a gene" "\n" "just transcribed coding transcripts")
plt.show()
##plt.savefig('distribution_of_coding-gene_biotypes_B.png')

"""
for ontology analysis
"""
groupA=out.loc[(out['%mix']<50) & (out['%just-coding']==0)]#gene that in most of their detected tissues are non-coding and never had mixture of coding and non-coding transcripts
groupB=out.loc[(out['%mix']==0) & (out['%just-coding']>50) & (out['%just-coding']<100)]#gene that in most of their detected tissues are coding and never had mixture of coding and non-coding transcripts
groupC=out.loc[(out['%mix']==100)]
m=groupB.merge(ensembl)['ENS'].to_frame()
m.to_csv('results.txt',index = None, header=False,sep="\t")

#tissues=list(switched_genes_m.columns)
#fetal_t=[x for x in tissues if x.find('Fetal')==0]
#adult_t=[x for x in tissues if x.find('Fetal')==-1]
fetal_t=['Fetal-Brain','Fetal-Kidney','Fetal-Muscle','Fetal-Spleen']
adult_t=['Brain__cortex__Frontal','Kidney','Longissimus_dorsi_ribeye_loin','Spleen']
just_p_m_fetal = just_p_m[fetal_t].rename({'Fetal-Brain':'Fetal-Brain_justP','Fetal-Kidney':'Fetal-Kidney_justP','Fetal-Muscle':'Fetal-Muscle_justP','Fetal-Spleen':'Fetal-Spleen_justP'},axis=1)
just_p_m_adult = just_p_m[adult_t].rename({'Brain__cortex__Frontal':'Adult-Brain_justP','Kidney':'Adult-Kidney_justP','Longissimus_dorsi_ribeye_loin':'Adult-Muscle_justP','Spleen':'Adult-Spleen_justP'},axis=1)
mix_m_fetal = mix_m[fetal_t].rename({'Fetal-Brain':'Fetal-Brain_M','Fetal-Kidney':'Fetal-Kidney_M','Fetal-Muscle':'Fetal-Muscle_M','Fetal-Spleen':'Fetal-Spleen_M'},axis=1)
mix_m_adult = mix_m[adult_t].rename({'Brain__cortex__Frontal':'Adult-Brain_M','Kidney':'Adult-Kidney_M','Longissimus_dorsi_ribeye_loin':'Adult-Muscle_M','Spleen':'Adult-Spleen_M'},axis=1)
noncoding_m_fetal = noncoding_m[fetal_t].rename({'Fetal-Brain':'Fetal-Brain_justN','Fetal-Kidney':'Fetal-Kidney_justN','Fetal-Muscle':'Fetal-Muscle_justN','Fetal-Spleen':'Fetal-Spleen_justN'},axis=1)
noncoding_m_adult = noncoding_m[adult_t].rename({'Brain__cortex__Frontal':'Adult-Brain_justN','Kidney':'Adult-Kidney_justN','Longissimus_dorsi_ribeye_loin':'Adult-Muscle_justN','Spleen':'Adult-Spleen_justN'},axis=1)

data_frames=[just_p_m_fetal,just_p_m_adult,mix_m_fetal,mix_m_adult,noncoding_m_fetal,noncoding_m_adult]
df_merged = reduce(lambda  left,right: pd.merge(left,right,how='outer',left_index=True, right_index=True), data_frames)
df_merged=df_merged.replace(np.NaN,0)

df_merged['fetal_sum']=df_merged.iloc[:, 0:4].sum(axis=1) + df_merged.iloc[:, 8:12].sum(axis=1) + df_merged.iloc[:, 16:20].sum(axis=1)
df_merged['adult_sum']=df_merged.iloc[:, 4:8].sum(axis=1) + df_merged.iloc[:, 12:16].sum(axis=1) + df_merged.iloc[:, 20:24].sum(axis=1)
df_merged['fetal_justN_sum']=df_merged.iloc[:, 16:20].sum(axis=1)
df_merged['adult_justP_sum']=df_merged.iloc[:, 4:8].sum(axis=1)
df_merged['adult_M_sum']=df_merged.iloc[:, 12:16].sum(axis=1)
#df_merged.loc[(df_merged['fetal_sum']>0) & (df_merged['adult_sum']==0)]#
brain=df_merged.loc[(df_merged['Fetal-Brain_justN']==1) & ((df_merged['Adult-Brain_justP']==1)  | (df_merged['Adult-Brain_M']==1))].reset_index()['gene_id'].to_frame()
kidney=df_merged.loc[(df_merged['Fetal-Kidney_justN']==1) & ((df_merged['Adult-Kidney_justP']==1)  | (df_merged['Adult-Kidney_M']==1))].reset_index()['gene_id'].to_frame()
muscle=df_merged.loc[(df_merged['Fetal-Muscle_justN']==1) & ((df_merged['Adult-Muscle_justP']==1)  | (df_merged['Adult-Muscle_M']==1))].reset_index()['gene_id'].to_frame()
spleen=df_merged.loc[(df_merged['Fetal-Spleen_justN']==1) & ((df_merged['Adult-Spleen_justP']==1)  | (df_merged['Adult-Spleen_M']==1))].reset_index()['gene_id'].to_frame()

brain=df_merged.loc[(df_merged['Fetal-Brain_justN']==1) & ((df_merged['Adult-Brain_justP']==1)  & (df_merged['Adult-Brain_M']==0))].reset_index()['gene_id'].to_frame()
kidney=df_merged.loc[(df_merged['Fetal-Kidney_justN']==1) & ((df_merged['Adult-Kidney_justP']==1)  & (df_merged['Adult-Kidney_M']==0))].reset_index()['gene_id'].to_frame()
muscle=df_merged.loc[(df_merged['Fetal-Muscle_justN']==1) & ((df_merged['Adult-Muscle_justP']==1)  & (df_merged['Adult-Muscle_M']==0))].reset_index()['gene_id'].to_frame()
spleen=df_merged.loc[(df_merged['Fetal-Spleen_justN']==1) & ((df_merged['Adult-Spleen_justP']==1)  & (df_merged['Adult-Spleen_M']==0))].reset_index()['gene_id'].to_frame()

genes=pd.concat([brain,kidney,muscle,spleen]).drop_duplicates()
m=genes.merge(ensembl)['ENS'].to_frame()
m.to_csv('results.txt',index = None, header=False,sep="\t")

df_merged.loc[(df_merged['fetal_justN_sum']==4) & ((df_merged['adult_justP_sum']==4) | (df_merged['adult_M_sum']==4))]






a=pd.DataFrame.to_numpy(just_p_m_fetal)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=just_p_m_fetal.index
just_p_fetal=d.reset_index().rename({'index':'gene_id',0:'number_of_just_just_p_fetal_tissues'},axis=1)
a=pd.DataFrame.to_numpy(just_p_m_adult)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=just_p_m_adult.index
just_p_adult=d.reset_index().rename({'index':'gene_id',0:'number_of_just_just_p_adult_tissues'},axis=1)
a=pd.DataFrame.to_numpy(noncoding_m_fetal)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=noncoding_m_fetal.index
noncoding_fetal=d.reset_index().rename({'index':'gene_id',0:'number_of_just_noncoding_fetal_tissues'},axis=1)
a=pd.DataFrame.to_numpy(noncoding_m_adult)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=noncoding_m_adult.index
noncoding_adult=d.reset_index().rename({'index':'gene_id',0:'number_of_just_noncoding_adult_tissues'},axis=1)
a=pd.DataFrame.to_numpy(mix_m_fetal)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=mix_m_fetal.index
mix_fetal=d.reset_index().rename({'index':'gene_id',0:'number_of_mix_fetal_tissues'},axis=1)
a=pd.DataFrame.to_numpy(mix_m_adult)#convert dataframe to numpy
d=pd.DataFrame(np.nansum(a,axis=1))
d.index=mix_m_adult.index
mix_adult=d.reset_index().rename({'index':'gene_id',0:'number_of_mix_adult_tissues'},axis=1)

data_frames=[just_p_fetal,noncoding_fetal,mix_fetal]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],how='outer'), data_frames)
df_fetal=df_merged.replace(np.NaN,0)
data_frames=[just_p_adult,noncoding_adult,mix_adult]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],how='outer'), data_frames)
df_adult=df_merged.replace(np.NaN,0)


data_frames=[just_p_fetal,just_p_adult,noncoding_fetal,noncoding_adult,mix_fetal,mix_adult]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],how='outer'), data_frames)
df_merged=df_merged.replace(np.NaN,0)

df_merged.loc[(df_merged['number_of_just_noncoding_fetal_tissues']>0) & (df_merged['number_of_mix_fetal_tissues']==0) & (df_merged['number_of_just_noncoding_adult_tissues']==0) & (df_merged['number_of_mix_adult_tissues']==0)]
m=df_merged.merge(ensembl)['ENS'].to_frame()
m.to_csv('results.txt',index = None, header=False,sep="\t")




'''
a=out['%just-coding'].to_frame().rename({'%just-coding':'percentage'},axis=1)
a['type']='just-coding'
b=out['%mix'].to_frame().rename({'%mix':'percentage'},axis=1)
b['type']='mix'
c=out['%just-noncoding'].to_frame().rename({'%just-noncoding':'percentage'},axis=1)
c['type']='noncoding'
df=pd.concat([a,b,c])
ax = sns.boxplot(x="type", y="percentage", data=df, linewidth=2.5,showfliers=True)
plt.show()

a=out['number_of_tissues_gene_encode_just_coding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_coding_transcripts':'numbers'},axis=1)
a['type']='just-coding'
b=out['number_of_tissues_gene_encode_mix_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_mix_transcripts':'numbers'},axis=1)
b['type']='mix'
c=out['number_of_tissues_gene_encode_just_noncoding_transcripts'].to_frame().rename({'number_of_tissues_gene_encode_just_noncoding_transcripts':'numbers'},axis=1)
c['type']='noncoding'
df=pd.concat([a,b,c])
ax = sns.boxplot(x="type", y="numbers", data=df, linewidth=2.5,showfliers=True)
plt.show()

ll = plt.scatter(out["%mix"], out["%just-coding"], marker='o', c=out["%just-noncoding"],cmap=plt.get_cmap("jet"),s=2)#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
plt.colorbar(label="Percentage of detected tissues that gene" "\n" "just transcribed non-coding transcripts")
plt.xlabel("Percentage of detected tissues that gene" "\n" "transcribed mixture of coding and non-coding transcripts")
plt.ylabel("Percentage of detected tissues that gene" "\n" "just transcribed coding transcripts")
plt.show()
'''


'''
test=out.rename({'number_of_tissues_gene_encode_just_coding_transcripts':'just_coding','number_of_tissues_gene_encode_mix_transcripts':"mix",'number_of_tissues_gene_encode_just_noncoding_transcripts':'noncoding'},axis=1)
test.loc[ (test['just_coding']>0) & (test['mix']>0) & (test['noncoding']>0) ]




'''