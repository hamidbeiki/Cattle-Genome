#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 15:57:35 2021

@author: beiki
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import collections
import csv
from scipy.stats import gaussian_kde
from scipy import stats

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

def format_mat(mat_dict):
    global tissues
    mat=pd.DataFrame.from_dict(mat_dict,orient='index')
    mat.columns=tissues
    mat=mat.reset_index()
    return(mat)

    
data=pd.read_csv('miRNA_normalised_counts-combined_samples-duplicates_removed.txt',header=0,sep='\t')
counts=pd.read_csv('miRNA_counts-combined_samples-duplicated_removed.txt',header=0,sep='\t')
novel_miRNAs=pd.read_csv('novel_miRNAs-combined_samples',names=['miRNA'],sep='\t')
tsi=pd.read_csv('miRNAs-tissue-specificity-scores.txt',names=['miRNA','TSI'],sep='\t')
data=data.set_index('miRNA')
data=data.loc[(data!=0).any(1)]
remove_miRNAs=counts.loc[~counts.miRNA.isin(data.index)]['miRNA'].to_frame()#miRNAs not expressed in any sample
counts=counts.set_index('miRNA')
counts=counts.loc[(counts!=0).any(1)]

counts_sum=counts.loc[:,'Abomasum_whole':'uterine_myometrium'].sum(axis=1).to_frame().reset_index().rename({0:'sum_counts'},axis=1)
counts_sum.loc[counts_sum['sum_counts']<=10].shape[0]

known=data.loc[~data.index.isin(novel_miRNAs.miRNA)]
novel=data.loc[data.index.isin(novel_miRNAs.miRNA)]

known_info=get_info(known,'miRNA')
novel_info=get_info(novel,'miRNA')

known_info['detected_tissues'].median()
novel_info['detected_tissues'].median()
known_info['median_expression'].median()
novel_info['median_expression'].median()

stat, p = stats.mannwhitneyu(known_info['median_expression'],novel_info['median_expression'],alternative='greater')
stat, p = stats.mannwhitneyu(known_info['detected_tissues'],novel_info['detected_tissues'],alternative='greater')

tissues=list(data.columns)
tissues_info=[]
for tissue in tissues:
    df=data[[tissue]].reset_index()
    expr=df.loc[df[tissue]>0]
    expr_known=expr.loc[~expr.miRNA.isin(novel_miRNAs.miRNA)]
    expr_novel=expr.loc[expr.miRNA.isin(novel_miRNAs.miRNA)]
    tissues_info.append([tissue,expr.shape[0],expr_known.shape[0],expr_novel.shape[0]])
tissues_info=pd.DataFrame(tissues_info).rename({0:'Tissue',1:'numper_miRNAs',2:'known',3:'novel'},axis=1)
tissues_info.sort_values('novel')

df=tissues_info['numper_miRNAs'].to_frame()
df['type']='miRNA'
ax = sns.boxplot(x="type", y="numper_miRNAs", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_miRNA_number.png',dpi=300)
plt.close()


a=tissues_info['known'].to_frame().rename({'known':'info'},axis=1)
a['type']='known'
b=tissues_info['novel'].to_frame().rename({'novel':'info'},axis=1)
b['type']='novel'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_known-novel_miRNA_number.png',dpi=300)
plt.close()

known_info['log']=np.log2(known_info['median_expression'])
novel_info['log']=np.log2(novel_info['median_expression'])
a=known_info['log'].to_frame().rename({'log':'info'},axis=1)
a['type']='known'
b=novel_info['log'].to_frame().rename({'log':'info'},axis=1)
b['type']='novel'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('miRNAs_expression.png',dpi=300)
plt.close()

a=known_info['detected_tissues'].to_frame().rename({'detected_tissues':'info'},axis=1)
a['type']='known'
b=novel_info['detected_tissues'].to_frame().rename({'detected_tissues':'info'},axis=1)
b['type']='novel'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('miRNAs_detected_tissues.png',dpi=300)
plt.close()


data_dict={}
for tissue in tissues:
    df=data[tissue].reset_index()
    expr=list(df.loc[df[tissue]>0]['miRNA'])
    data_dict[tissue]=expr
    
tissue_sim=collections.OrderedDict()
tissue_sim2=collections.OrderedDict()#this is for heatmap
flag=1
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        l1=data_dict[t1]
        l2=data_dict[t2]
        intersect=list(set(l1) & set(l2))
        num=len(intersect)
        perc=(num/float(len(l1)))*100
        v.append(num)
        v2.append(perc)
    tissue_sim[tissues[i]]=v
    tissue_sim2[tissues[i]]=v2    
    print(flag) 
    flag=flag+1

tissue_sim_mat=format_mat(tissue_sim)
tissue_sim_mat2=format_mat(tissue_sim2)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)
tissue_sim_mat.to_csv('miRNA_level_tissue_comparisions.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat2.to_csv('miRNA_level_tissue_comparisions_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal) 

''' tissue specificity '''
info=pd.concat([known_info,novel_info])
info2=get_info(counts,'miRNA')

plt.hist(info['detected_tissues'], bins=50,color='blue')
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of miRNA')
plt.savefig('hist_of_miRNAs_detected_tissues_2.png',dpi=300)
plt.close()  

plt.hist(known_info['detected_tissues'], bins=50, alpha=0.5, label='Known miRNAs')
plt.hist(novel_info['detected_tissues'], bins=50, alpha=0.5, label='Novel miRNAs')
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of miRNAs')
#plt.legend(loc='upper right')
#plt.show()
plt.savefig('hist_of_miRNAs_detected_tissues_1.png',dpi=300)
plt.close()

info['log']=np.log2(info['median_expression'])
TS_mirna=info.loc[info['detected_tissues']==1]
rest_mirna=info.loc[~info.miRNA.isin(TS_mirna.miRNA)]

a=rest_mirna['log'].to_frame().rename({'log':'info'},axis=1)
a['type']='rest'
b=TS_mirna['log'].to_frame().rename({'log':'info'},axis=1)
b['type']='TS'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue-specific_vs_rest_miRNAs_expression.png',dpi=300)
plt.close()

tissues=list(data.columns)
tissues_info2=[]
for tissue in tissues:
    df=data[[tissue]].reset_index()
    expr=df.loc[df[tissue]>0]
    ts=expr.merge(TS_mirna)
    rest=info.loc[~info.miRNA.isin(ts.miRNA)]
    tissues_info2.append([tissue,ts.shape[0]])
tissues_info2_df=pd.DataFrame(tissues_info2).rename({0:'tissue',1:'TS_mirna'},axis=1)
df=tissues_info2_df['TS_mirna'].to_frame().rename({'TS_mirna':'info'},axis=1)
df['type']='miRNAs'
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue-specific_miRNAs.png',dpi=300)
plt.close()
''' TSI based tissue specificity '''

tsi2=tsi.loc[tsi['TSI']<1]
plt.hist(tsi2['TSI'], bins=50,color='black')
plt.xlabel('Tissue Specificity Score')
plt.ylabel('Number of multi-tissue detected miRNAs')
plt.savefig('hist_of_multi-tissue_miRNAs_tissue-specificity_scores.png',dpi=300)
plt.close()

df=tsi2.merge(info)
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
plt.ylabel("Tissue Specificity Score of multi-tissue detected miRNAs")
#plt.axhline(0.986882, linestyle='--',color='red')
plt.savefig('relation_between_TSI_and_number_of_detected_tissues_in_multi-tissue_miRNAs.png',dpi=300)
plt.close()

TS_mirna_counts=TS_mirna.merge(info2,on='miRNA')
TS_mirna_counts.loc[TS_mirna_counts['median_expression_y']<10].shape[0]

multi_tissue=info.loc[info['detected_tissues']>1]
multi_tissue=multi_tissue.merge(tsi)
multi_tissue.loc[multi_tissue['TSI']>0.9]
df=multi_tissue.loc[multi_tissue['TSI']>0.9]
df['detected_tissues'].median()

 
'''
data=pd.read_csv('query.txt',header=0,sep='\t')
df=data['total_miRNA_reads'].to_frame()
df['type']='miRNA'
ax = sns.boxplot(x="type", y="total_miRNA_reads", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('miRNA_input_reads.png',dpi=300)
plt.close()
'''