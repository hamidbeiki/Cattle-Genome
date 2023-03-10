#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 18:08:56 2021

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
results=[]
results2={}
for tissue in tissues:
    df=expressions[['transcript_id','gene_id',tissue]]
    df=df.merge(tr_biotypes)
    df=df.merge(gene_biotypes)
    df=df.loc[df['gene_biotype']=="protein-coding"]# limit to just protein-coding genes
    df=df.loc[(df[tissue]>0)]# & (df[tissue]<1)
    trs_per_gene=df.groupby(['gene_id']).size().reset_index(name='count').rename({'count':'trs'},axis=1)
    genes_expr=df.groupby(['gene_id'], sort=False)[tissue].sum().reset_index().rename({tissue:'expression'},axis=1) #summation of transcripts expression for each gene
    nds=df.loc[df['tr_biotype']=='non_stop_decays']['gene_id'].to_frame()
    nmds=df.loc[df['tr_biotype']=='NMDs']['gene_id'].to_frame()
    info=pd.concat([nds,nmds])    
    aberrants=info.groupby(['gene_id']).size().reset_index(name='count').rename({'count':'aberrant_trs'},axis=1)
#    info=df.loc[df['tr_biotype']!='protein_coding_transcripts']['gene_id'].to_frame()
#    aberrants=info.groupby(['gene_id']).size().reset_index(name='count').rename({'count':'aberrant_trs'},axis=1)
    res=pd.merge(pd.merge(trs_per_gene,aberrants,on='gene_id',how='outer'),genes_expr,on='gene_id')
    df= res.replace(np.NaN,0)
    df=df.merge(gene_splice_info)
    df2=df.loc[df['aberrant_trs']>0]
    rest=df.loc[~ df.gene_id.isin(df2.gene_id)]
    rest2=rest.loc[rest['number_of_splice']>0]
#    df2=res.loc[(df['number_of_splice']>0) & (df['aberrant_trs']>0)]
    percent=round(df2.shape[0]/float(df.shape[0]),2)*100# % of protein coding genes with aberant trs
    coef, p = pearsonr(df2['trs'],df2['aberrant_trs'])
    stat, p2, med, tbl = median_test(df2['number_of_splice'],rest['number_of_splice'])
    results.append([tissue,percent,coef,p,df2['number_of_splice'].median(),rest['number_of_splice'].median(),p2])
    df=df[['gene_id','aberrant_trs']].rename({'aberrant_trs':tissue},axis=1)
    results2[tissue]=df
res_df=pd.DataFrame(results).rename({0:'tissue',1:'percent_of_pGenes_with_aberants',2:'corr_number_of_trs_number_of_aberrants',3:'p-value',4:'median_number_of_introns_in_Pgenes_with_aberrants',5:'median_number_of_introns_in_rest_of_Pgenes',6:'p-value'},axis=1)
res_df['percent_of_pGenes_with_aberants'].median()
res_df['corr_number_of_trs_number_of_aberrants'].median()
res_df['median_number_of_introns_in_Pgenes_with_aberrants'].median()
res_df['median_number_of_introns_in_rest_of_Pgenes'].median()

""" merge all data frames
    in a dict #https://stackoverflow.com/questions/53935848/how-to-merge-all-data-frames-in-a-dictionary-in-python
"""
my_reduce = partial(pd.merge, on='gene_id', how='outer')
df=reduce(my_reduce, results2.values())
df=df.set_index('gene_id')
a=pd.DataFrame.to_numpy(df)#convert dataframe to numpy
d=pd.DataFrame(np.count_nonzero((a),axis=1) - np.sum(np.isnan(a),axis=1))# count the number of nonzero nonNan in per row in an array
d.index=df.index
df2=d.reset_index().rename({'gene_id':'protein_coding_gene',0:'number_of_tissues_with_aberrant'},axis=1)# remaining protein-coding genes just couldn't be quantified in any tissue.
df2.to_csv('number_of_tissues_that_each_protein-gene_has_aberrant_transcript',index = None, header=True,sep="\t")

df2=df2.rename({'protein_coding_gene':'gene_id'},axis=1)
info=pd.read_csv('number_of_aberrant_transcripts_per_gene',names=['gene_id','number_of_abberant_trs'],sep='\t')
info2=pd.read_csv('genes_detection_info',header=0,sep='\t')
info3=pd.read_csv('number_of_transcript_per_gene',names=['gene_id','number_of_transcripts'],sep='\t')
m=info.merge(df2)
m2=info2.merge(m)
m3=info3.merge(m2)

gene_expressions=gene_expressions.loc[(gene_expressions!=0).any(1)]
ids=pd.DataFrame(gene_expressions.index)
gene_expressions = gene_expressions.replace(0, np.NaN)
array=pd.DataFrame.to_numpy(gene_expressions)#convert dataframe to numpy
gene_expressions_row_medians = np.nanmedian(array, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
median_expression=pd.DataFrame(gene_expressions_row_medians).rename({0:'median_expression_in_detected_tissue'},axis=1)#median expression in their detected tissues
expresion_info=pd.concat([ids,median_expression],axis=1)#remaining genes could'nt be quantified in any tissue

m3=m3.merge(expresion_info)
rest=expresion_info[~ expresion_info.gene_id.isin(m3.gene_id)]
rest=rest.merge(gene_biotypes)
rest=rest.loc[rest['gene_biotype']=="protein-coding"]

m3['median_expression_in_detected_tissue'].median()
rest['median_expression_in_detected_tissue'].median()

m3.to_csv('protein-coding_genes_with_aberrant_transcript_info',index = None, header=True,sep="\t")
m3['diff']=m3["number_of_detected_tissue"]-m3["number_of_tissues_with_aberrant"]
m3['aberrant_%']=(m3['number_of_abberant_trs']/(m3['number_of_transcripts']))*100
m3['aberrant_tissue_%']=(m3['number_of_tissues_with_aberrant']/(m3['number_of_detected_tissue']))*100

coef, p = pearsonr(m3['number_of_abberant_trs'],m3['median_expression_in_detected_tissue'])
coef, p = spearmanr(m3['aberrant_%'],m3['median_expression_in_detected_tissue'])

plt.scatter(m3["number_of_tissues_with_aberrant"], m3["aberrant_%"], marker='o',c=m3["diff"],cmap=plt.get_cmap("jet"),alpha=0.4,s=2)
plt.xlabel("Number of tissue that a gene" "\n" "transcribed PAT transcript(s)")# PAT stands for potentially aberrant transcripts
plt.ylabel("Percentage of transcript that are PAT")
plt.colorbar(label="Number of tissue that a gene" "\n" "transcribed NO PAT transcript(s)")
plt.show()
#plt.savefig('variation_of_genes_with_PAT_transcripts.png')


plt.scatter(m3["aberrant_tissue_%"], m3["aberrant_%"], marker='o',c=m3["number_of_detected_tissue"],alpha=1,s=1,cmap=plt.get_cmap("jet"))#,c=m3["aberrant_%"],cmap=plt.get_cmap("jet"),alpha=1,s=2
plt.xlabel("Percentage of detected tissue that a gene" "\n" "transcribed PAT transcript(s)")# PAT stands for potentially aberrant transcripts
plt.ylabel("Percentage of transcript that are PAT")
plt.colorbar(label="Number of detected tissues")
#plt.show()
plt.savefig('variation_of_genes_with_PAT_transcripts.png')

s=m3.loc[(m3["aberrant_tissue_%"]==100) & (m3["aberrant_%"]>50)]
s=s.merge(ensembl)
out=s['ENS'].drop_duplicates().to_frame()
out.to_csv('protein-coding_genes_mainly_transcribe_PAT_transcripts.txt',index = None, header=False,sep="\t")

remove=pd.read_csv('remove_genes',names=['gene_id'],sep='\t')
genes_detection=info2.loc[~ info2.gene_id.isin(remove.gene_id)]
rest=genes_detection.loc[~ genes_detection.gene_id.isin(m3.gene_id)]
rest=rest.merge(gene_biotypes)
rest_coding=rest.loc[rest['gene_biotype']=="protein-coding"]
nc=rest.loc[(rest['gene_biotype']!="protein-coding") & (rest['gene_biotype']!="pseudogenes")]
ps=rest.loc[rest['gene_biotype']=="pseudogenes"]

a=m3['number_of_detected_tissue'].to_frame()
a['type']='aberrant'
b=rest_coding['number_of_detected_tissue'].to_frame()
b['type']='rest_coding'
c=nc['number_of_detected_tissue'].to_frame()
c['type']='nc_genes'
d=ps['number_of_detected_tissue'].to_frame()
d['type']='ps_genes'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type", y="number_of_detected_tissue", data=df, linewidth=2.5,showfliers=False)
plt.savefig('gene_biotype_tissue_detection.png')


plt.hist(m3['number_of_tissues_with_aberrant'], bins=50, alpha=0.5)
plt.xlabel('Number of tissues')
plt.ylabel('Number of genes with aberrant transcripts')
plt.legend(loc='upper right')
plt.savefig('test.png')
plt.close()


"""    
data=pd.read_csv('number_of_aberrant_transcripts_per_gene',names=['gene_id','number'],sep='\t')
data['number'].quantile([0.0, .5, .90, .95])
df=data.loc[data['number']<=10]
df['number'].quantile([0.0, .5, .90, .95])
plt.hist(df['number'], bins=50, alpha=0.5, label='distribution of the number of\abberant transcripts per gene')
plt.xlabel('number of abberant transcripts per gene')
plt.ylabel('number genes')

data['number'].sum()
df['number'].sum()

"""