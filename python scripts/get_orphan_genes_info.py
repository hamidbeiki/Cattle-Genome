#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 12:00:24 2021

@author: beiki
"""

import pandas as pd
import numpy as np
#from scipy.stats import median_test
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

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
    
gene_expressions=pd.read_csv("RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t")#,index_col='gene_id'
gene_biotypes=pd.read_csv('final_gene_biotypes',names=['gene_id','gene_biotype'],sep='\t')
orphans=pd.read_csv('orphan-genes',names=['gene_id'],sep='\t')
tissues=list(gene_expressions.columns)[1:]
remove=pd.read_csv('remove_genes',header=0,names=['gene_id'])
gene_expressions=gene_expressions.loc[~ gene_expressions.gene_id.isin(remove.gene_id)]
gene_expressions=gene_expressions.set_index('gene_id')
gene_expressions=gene_expressions.loc[(gene_expressions!=0).any(1)]#remove genes that couldn't be quantified
orphans=orphans.loc[~ orphans.gene_id.isin(remove.gene_id)]

expressions=gene_expressions.reset_index()
results=[]
for tissue in tissues:
    df=expressions[['gene_id',tissue]]
    df=df.merge(gene_biotypes)
    df=df.loc[(df[tissue]>0)]# & (df[tissue]<1)
    a=df.loc[df['gene_biotype']=='protein-coding']
    b=df.loc[df['gene_biotype']!='protein-coding']
    orph=orphans.merge(df)
    p=a.loc[a.gene_id.isin(orph.gene_id)]
    p_rest=a.loc[~ a.gene_id.isin(orph.gene_id)]
    nc=b.loc[b.gene_id.isin(orph.gene_id)]
    nc_rest=b.loc[~ b.gene_id.isin(orph.gene_id)]
#    stat, pval, med, tbl = median_test(p_rest[tissue],p[tissue])
#    stat, pval2, med, tbl = median_test(nc_rest[tissue],nc[tissue])
    stat,pval=stats.mannwhitneyu(p_rest[tissue],p[tissue])#Mann-Whitney test that is equivalent to Wilcoxon Rank Sum test with sample sizes are uequivalent: https://stats.stackexchange.com/questions/368881/wilcoxon-rank-sum-test-unequal-sample-sizes
    stat,pval2=stats.mannwhitneyu(nc_rest[tissue],nc[tissue])
    results.append([tissue,p_rest[tissue].median(),p[tissue].median(),pval,nc_rest[tissue].median(),nc[tissue].median(),pval2,nc.shape[0],p.shape[0]])
res_df=pd.DataFrame(results).rename({0:'tissue',1:'coding-homologous_median_exp',2:'coding-orphan_median_exp',3:'p-value',4:'nc-homologous_median_exp',5:'nc-orphan_median_exp',6:'p-value2',7:'number_of_nc_orphans',8:'number_of_coding_orphans'},axis=1)
res_df['p-value2']=1 # because the number of expressed nc-orphan genes in each tissue is very low and it is impossible to perform statistical test
res_df['log_p-value']=-np.log2(res_df['p-value'])   
res_df['log_p-value2']=-np.log2(res_df['p-value2']) 

gene_expressions=gene_expressions.set_index('gene_id')
ex = gene_expressions.replace(0, np.NaN)
df=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy
df_row_counts = np.count_nonzero(~np.isnan(df),axis=1)
gene_expressions=gene_expressions.reset_index()
detection_df=pd.DataFrame(df_row_counts).rename({0:'detected_tissues'},axis=1)
detection_df['gene_id']=gene_expressions['gene_id']

gene_expressions=gene_expressions.set_index('gene_id')
gene_expressions=gene_expressions.loc[(gene_expressions!=0).any(1)]
ids=pd.DataFrame(gene_expressions.index)
gene_expressions = gene_expressions.replace(0, np.NaN)
array=pd.DataFrame.to_numpy(gene_expressions)#convert dataframe to numpy
gene_expressions_row_medians = np.nanmedian(array, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
median_expression=pd.DataFrame(gene_expressions_row_medians).rename({0:'median_expression_in_detected_tissue'},axis=1)#median expression in their detected tissues
expresion_info=pd.concat([ids,median_expression],axis=1)#remaining genes could'nt be quantified in any tissue

orphans=orphans.merge(detection_df)
orphans=orphans.merge(expresion_info)
orphans=gene_biotypes.merge(orphans)

res_df.to_csv('orphan-genes-info1',index = None, header=True,sep="\t")
orphans.to_csv('orphan-genes-info2',index = None, header=True,sep="\t")# remaining genes could'nt be quantified in any tissue

orphans['log']=np.log(orphans['median_expression_in_detected_tissue'])
info=expresion_info.merge(gene_biotypes)
info['log']=np.log(info['median_expression_in_detected_tissue'])
info=info.merge(detection_df)
p=info.loc[info['gene_biotype']=='protein-coding']
nc=info.loc[info['gene_biotype']!='protein-coding']
orph_c=orphans.loc[orphans['gene_biotype']=='protein-coding']
orph_nc=orphans.loc[orphans['gene_biotype']!='protein-coding']
rest_c=p.loc[~ p.gene_id.isin(orph_c.gene_id)]
rest_nc=nc.loc[~ nc.gene_id.isin(orph_nc.gene_id)]

orph_nc['detected_tissues'].median()
orph_c['detected_tissues'].median()
rest_nc['detected_tissues'].median()
rest_c['detected_tissues'].median()

stat, p1 = stats.ttest_ind(orph_nc['detected_tissues'],rest_nc['detected_tissues'],equal_var=True,alternative='greater')
stat, p2 = stats.ttest_ind(rest_c['detected_tissues'],orph_c['detected_tissues'],equal_var=True,alternative='greater')

stat,pval=stats.mannwhitneyu(orph_c['detected_tissues'],rest_c['detected_tissues'])
stat,pval=stats.mannwhitneyu(orph_nc['detected_tissues'],rest_nc['detected_tissues'])

a=orph_nc['detected_tissues'].to_frame()
a['type1']='orphan'
a['type2']='noncoding'
b=orph_c['detected_tissues'].to_frame()
b['type1']='orphan'
b['type2']='coding'
c=rest_nc['detected_tissues'].to_frame()
c['type1']='rest'
c['type2']='noncoding'
d=rest_c['detected_tissues'].to_frame()
d['type1']='rest'
d['type2']='coding'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type1", y="detected_tissues", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()

a=orph_nc['log'].to_frame()
a['type1']='orphan'
a['type2']='noncoding'
b=orph_c['log'].to_frame()
b['type1']='orphan'
b['type2']='coding'
c=rest_nc['log'].to_frame()
c['type1']='rest'
c['type2']='noncoding'
d=rest_c['log'].to_frame()
d['type1']='rest'
d['type2']='coding'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type1", y="log", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('test.png',dpi=300)
plt.close() 

ll=plt.scatter(orph_c["detected_tissues"], orph_c["median_expression_in_detected_tissue"], marker='o',s=4,alpha=0.3,color='red')#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
l=plt.scatter(orph_nc["detected_tissues"], orph_nc["median_expression_in_detected_tissue"], marker='o',s=5,alpha=0.5,color='blue')#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
plt.xlabel("Number of detected tissues")
plt.ylabel("Median expression level (RPKM) in detected tissues")
plt.show()


res_df["log_p-value"]=np.where(res_df["log_p-value"]>=50 , 50, res_df["log_p-value"])
l=plt.scatter(res_df["nc-orphan_median_exp"], res_df["nc-homologous_median_exp"], marker='o',s=4,color='black')#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
ll=plt.scatter(res_df["coding-orphan_median_exp"], res_df["coding-homologous_median_exp"], marker='o',s=4,alpha=1,c=res_df["log_p-value"],cmap=plt.get_cmap("jet"))#
#lll=plt.plot(res_df["nc-orphan_median_exp"],res_df["nc-orphan_median_exp"],'r--',color='black')
plt.colorbar(label="-log2(p-value)")
plt.xlabel("median expression of orphan genes")
plt.ylabel("median expression of homologous genes")
plt.show()


orph_nc["median_expression_in_detected_tissue"].median()
res_df["nc-orphan_median_exp"].median()
orph_c["median_expression_in_detected_tissue"].median()
res_df["coding-orphan_median_exp"].median()


orphans_ex=gene_expressions.loc[gene_expressions.index.isin(orphans.gene_id)]
homologs_ex=gene_expressions.loc[~ gene_expressions.index.isin(orphans.index)]

orphans_info=get_info(orphans_ex,'genes')
homologs_info=get_info(homologs_ex,'genes')

orphans_info['log']=np.log2(orphans_info['median_expression'])
homologs_info['log']=np.log2(homologs_info['median_expression'])

a=orphans_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='orphans'
a['type2']='detected_tissues'
b=orphans_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='orphans'
b['type2']='median_expression'
c=homologs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='homologous'
c['type2']='detected_tissues'
d=homologs_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='homologous'
d['type2']='median_expression'
df=pd.concat([a,c])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
df=pd.concat([b,d])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
























