#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:47:34 2021

@author: beiki
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

remove=pd.read_csv('remove_genes',names=['gene_id'])
known=pd.read_csv('all-known-genes-main',names=['gene_id'],sep='\t')
biotypes=pd.read_csv('final_gene_biotypes',names=['gene_id','biotype'],sep='\t')
utr5_df=pd.read_csv('number_of_5UTRs_per_gene',header=0,sep='\t')
utr5_df=utr5_df[~ utr5_df.gene_id.isin(remove.gene_id)]
utr3_df=pd.read_csv('number_of_3UTRs_per_gene',header=0,sep='\t')
utr3_df=utr3_df[~ utr3_df.gene_id.isin(remove.gene_id)]
m=utr5_df.merge(utr3_df)
d=m.loc[(m['number_of_utr5']>1) | (m['number_of_utr3']>1)]
df=d.merge(biotypes)
df.groupby(['biotype']).size().reset_index(name='count')
df=df.loc[df['biotype']=='protein-coding']
df['number_of_utr3'].quantile([0.0, .5, .90, .95])
df['number_of_utr3'].quantile([0.0, .5, .90, .95])
df=df.loc[(df['number_of_utr5']<6) | (df['number_of_utr3']<6)]# just to remove outliers
a=df['number_of_utr5'].to_frame().rename({'number_of_utr5':'number'},axis=1)
a['type']='utr5'
b=df['number_of_utr3'].to_frame().rename({'number_of_utr3':'number'},axis=1)
b['type']='utr3'
df2=pd.concat([a,b])
ax = sns.boxplot(x="type", y="number", data=df2, linewidth=2.5,showfliers=False)
plt.savefig('ditribution-of-UTRs-per_gene.png')

df['number_of_utr5'].corr(df['number_of_utr3'])
'''
df3=df.loc[(df['number_of_utr5']<6) & (df['number_of_utr3']<6)]
plt.scatter(df3["number_of_utr5"], df3["number_of_utr3"], marker='o',alpha=0.1)
x=df3['number_of_utr5']
y=df3['number_of_utr3']
plt.hist(x, bins=50, alpha=0.5, label='number_of_5UTRs_per_gene')
plt.hist(y, bins=50, alpha=0.5, label='number_of_3UTRs_per_gene')
plt.xlabel('number_of_5UTRs_per_gene')
plt.ylabel('number_of_3UTRs_per_gene')

a=df.loc[df.gene_id.isin(known.gene_id)]
a['type']='known'
b=df.loc[~ df.gene_id.isin(known.gene_id)]
b['type']='novel'
df=pd.concat([a,b])
a=df.loc[df['type']=='known']['number_of_utr5'].to_frame().rename({'number_of_utr5':'number'},axis=1)
a['type']='utr5'
a['type2']='known'
b=df.loc[df['type']=='novel']['number_of_utr5'].to_frame().rename({'number_of_utr5':'number'},axis=1)
b['type']='utr5'
b['type2']='novel'
c=df.loc[df['type']=='known']['number_of_utr3'].to_frame().rename({'number_of_utr3':'number'},axis=1)
c['type']='utr3'
c['type2']='known'
d=df.loc[df['type']=='novel']['number_of_utr3'].to_frame().rename({'number_of_utr3':'number'},axis=1)
d['type']='utr3'
d['type2']='novel'
df2=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type", y="number", hue="type2", data=df2, linewidth=2.5,showfliers=False)
plt.show()

'''
