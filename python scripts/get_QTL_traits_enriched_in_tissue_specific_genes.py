#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:36:45 2021

@author: beiki
"""

import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import scipy.stats as stats

""" requiered modules:
    module load py-pandas/0.21.1-py2-326uzkn
    module load py-scipy/1.0.0-py2-flqcuxg
    module load py-matplotlib/2.0.2-py2-m755aws
"""

#########################
#                       #
# APROACH1: Gene based  #
#                       #
#########################
TS=pd.read_csv("tissue_specific_genes",names=["gene"])
tissues=list(pd.read_csv('tissues',names=['tissue'])['tissue'])
known=pd.read_csv('annotation/combined_genes_ensembl_equivalent',sep=" ",names=["gene","ensembl"])
known2=pd.read_csv('annotation/combined_genes_ncbi_equivalent',sep=" ",names=["gene","ncbi"])
known2['ncbi']=known2['ncbi'].str.replace(r'\D', '')
known_all=pd.read_csv('annotation/all-known-genes-main',sep=" ",names=["gene"])
results_dict={}
for tissue in tissues:
    gene_QTL=pd.read_csv(tissue + '_gene_closest_QTLs',header=0,sep="\t").rename({'gene_id':'gene'},axis=1)# gene_closest_QTLs or gene_significant_QTLs
    total_qtl_genes=len(gene_QTL['gene'].drop_duplicates())
    total_genes=len(pd.read_csv('quantification/'+ tissue + '_genes',names=['gene'])['gene'])
    m=TS.merge(gene_QTL)
    if m.shape[0]>0:
        rest=gene_QTL[~ gene_QTL.gene.isin(TS.gene)]
        traits=list(m['trait'].drop_duplicates())
        results=[]
        my_dict={}#stores known ensembl TS genes for each trait
        my_dict2={}#stores FAANG TS genes each trait
        my_dict3={}#stores all know TS genes
        my_dict4={}#stores known ncbi TS genes for each trait
        for i in range(len(traits)):
            df_ref=gene_QTL.loc[gene_QTL['trait']==traits[i]]
            df=m.loc[m['trait']==traits[i]]
            genes=df['gene'].drop_duplicates().to_frame()
            k2=genes.merge(known)
            k3=genes.merge(known_all)
            k4=genes.merge(known2)
            k=k2['ensembl'].drop_duplicates().to_frame()
            k2=k2['gene'].drop_duplicates().to_frame()
            k4=k4['ncbi'].drop_duplicates().to_frame()
            my_dict[traits[i]]=k
            my_dict2[traits[i]]=genes
            my_dict3[traits[i]]=k3
            my_dict4[traits[i]]=k4
            known_genes=k3.shape[0]
            genes=genes.shape[0]
            A=len(df.merge(df_ref)['gene'].drop_duplicates())
            total_t_gene=len(df_ref['gene'].drop_duplicates())
            B=total_t_gene-A
            C=len(m['gene'].drop_duplicates())-A
            D=total_genes-A-B-C
            oddsratio,pvalue = stats.fisher_exact([[A, B], [C, D]],alternative='greater')
            results.append([traits[i],pvalue,A,B,C,D,genes,known_genes])
        results_df=pd.DataFrame(results).rename({0:'trait',1:'p-value',2:'A', 3:'B',4:'C', 5:'D',6:'genes',7:'known_genes'},axis=1) 
        pvals=np.array(results_df['p-value'])
        pval_corr = smm.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)#https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/ 
        enriched_traits=[]
        flag=0
        for l in results:
            if pval_corr[1][flag]<=0.05:
                enriched_traits.append([l[0],l[1],pval_corr[1][flag],l[2],l[3],l[4],l[5],l[6],l[7]])
            flag=flag+1
        if len(enriched_traits)>0:
            enriched_traits_df=pd.DataFrame(enriched_traits).rename({0:'trait',1:'pvalue',2:'adj-pvalue',3:'A', 4:'B',5:'C', 6:'D',7:'genes',8:'known_genes'},axis=1)             
            enriched_traits_df=enriched_traits_df.sort_values(["adj-pvalue"])
            enriched_traits_df=enriched_traits_df.loc[enriched_traits_df['genes']>=3]#minimum number of 3 genes required for enrichment
            enriched_traits_df.to_csv(tissue + 'QTL_traits_enriched_in_tissue_specific_genes_approach1',index = None, header=True,sep="\t")
            results_dict[tissue]=[my_dict,my_dict2,my_dict3,my_dict4,enriched_traits_df] 
        else:
            results_dict[tissue]="NA"
        print(tissue)
np.save('QTL_traits_enriched_in_tissue_specific_genes_approach1.npy', results_dict)
#results_dict = np.load('QTL_traits_enriched_in_tissue_specific_genes_approach1.npy',allow_pickle='TRUE').item()
for k,l in results_dict.items():
    if l !='NA' and l[4].shape[0]>0:
        print(k)
#########################
#                       #
# APROACH2: QTL based   #
#                       #
#########################
TS=pd.read_csv("tissue_specific_genes",names=["gene"])
total_qtl=len(pd.read_csv('CattleQTLdb_April2020.bed',header=0,sep="\t").iloc[0:,3].drop_duplicates())    
tissues=list(pd.read_csv('tissues',names=['tissue'])['tissue'])
known=pd.read_csv('annotation/combined_genes_ensembl_equivalent',sep=" ",names=["gene","ensembl"])
known2=pd.read_csv('annotation/combined_genes_ncbi_equivalent',sep=" ",names=["gene","ncbi"])
known2['ncbi']=known2['ncbi'].str.replace(r'\D', '')
known_all=pd.read_csv('annotation/all-known-genes-main',sep=" ",names=["gene"])
results_dict={}
for tissue in tissues:
    gene_QTL=pd.read_csv(tissue + '_gene_closest_QTLs',header=0,sep="\t").rename({'gene_id':'gene'},axis=1)# gene_closest_QTLs or gene_significant_QTLs
    m=TS.merge(gene_QTL)
    if m.shape[0]>0:
        rest=gene_QTL[~ gene_QTL.gene.isin(TS.gene)]
        traits=list(m['trait'].drop_duplicates())
        results=[]
        my_dict={}#stores known ensembl TS genes for each trait
        my_dict2={}#stores FAANG TS genes each trait
        my_dict3={}#stores all know TS genes
        my_dict4={}#stores known ncbi TS genes for each trait
        for i in range(len(traits)):
            df_ref=gene_QTL.loc[gene_QTL['trait']==traits[i]]
            df=m.loc[m['trait']==traits[i]]
            genes=df['gene'].drop_duplicates().to_frame()
            k2=genes.merge(known)
            k3=genes.merge(known_all)
            k4=genes.merge(known2)
            k=k2['ensembl'].drop_duplicates().to_frame()
            k2=k2['gene'].drop_duplicates().to_frame()
            k4=k4['ncbi'].drop_duplicates().to_frame()
            my_dict[traits[i]]=k
            my_dict2[traits[i]]=genes
            my_dict3[traits[i]]=k3
            my_dict4[traits[i]]=k4
            known_genes=k3.shape[0]
            genes=genes.shape[0]
            A=len(df.merge(df_ref)['QTL'].drop_duplicates())
            total_t_qtl=len(df_ref['QTL'].drop_duplicates())
            B=total_t_qtl-A
            C=len(m['QTL'].drop_duplicates())-A
            D=total_qtl-A-B-C
            oddsratio,pvalue = stats.fisher_exact([[A, B], [C, D]],alternative='greater')
            results.append([traits[i],pvalue,A,B,C,D,genes,known_genes])
        results_df=pd.DataFrame(results).rename({0:'trait',1:'p-value',2:'A', 3:'B',4:'C', 5:'D',6:'genes',7:'known_genes'},axis=1) 
        pvals=np.array(results_df['p-value'])
        pval_corr = smm.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)#https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/ 
        enriched_traits=[]
        flag=0
        for l in results:
            if pval_corr[1][flag]<=0.05:
                enriched_traits.append([l[0],l[1],pval_corr[1][flag],l[2],l[3],l[4],l[5],l[6],l[7]])
            flag=flag+1
        if len(enriched_traits)>0:
            enriched_traits_df=pd.DataFrame(enriched_traits).rename({0:'trait',1:'pvalue',2:'adj-pvalue',3:'A', 4:'B',5:'C', 6:'D',7:'genes',8:'known_genes'},axis=1)             
            enriched_traits_df=enriched_traits_df.sort_values(["adj-pvalue"])
            enriched_traits_df.to_csv(tissue + 'QTL_traits_enriched_in_tissue_specific_genes',index = None, header=True,sep="\t")
            results_dict[tissue]=[my_dict,my_dict2,my_dict3,my_dict4,enriched_traits_df] 
        else:
            results_dict[tissue]="NA"
        print(tissue)
np.save('QTL_traits_enriched_in_tissue_specific_genes.npy', results_dict)
#results_dict = np.load('QTL_traits_enriched_in_tissue_specific_genes.npy',allow_pickle='TRUE').item()
keys=results_dict.keys()
uniqs={}# uniq trait in each tissue that not enriched in other tissues
for k in keys:
    if results_dict[k] != "NA":
        l0=[]
        for i in keys:
            if i != k and results_dict[i] != "NA":
                l1=list(results_dict[i][4]['trait'])
                for n in l1:
                    if n not in l0:
                        l0.append(n)
        l2=list(results_dict[k][4]['trait'])
        out=[]
        for j in l2:
            if j not in l0:
                out.append(j)
        if len(out)>0:
            uniqs[k]=out
        
"""
gene_QTL=pd.read_csv('earlylactating-Mammarygland_gene_closest_QTLs',header=0,sep="\t").rename({'gene_id':'gene'},axis=1)# gene_closest_QTLs or gene_significant_QTLs
m=TS.merge(gene_QTL)
rest=gene_QTL[~ gene_QTL.gene.isin(TS.gene)]
traits=list(m['trait'].drop_duplicates())
results=[]
my_dict={}#stores known ensembl TS genes for each trait
my_dict2={}#stores FAANG TS genes each trait
my_dict3={}#stores all know TS genes
my_dict4={}#stores known ncbi TS genes for each trait
for i in range(len(traits)):
    df_ref=gene_QTL.loc[gene_QTL['trait']==traits[i]]
    df=m.loc[m['trait']==traits[i]]
    genes=df['gene'].drop_duplicates().to_frame()
    k2=genes.merge(known)
    k3=genes.merge(known_all)
    k4=genes.merge(known2)
    k=k2['ensembl'].drop_duplicates().to_frame()
    k2=k2['gene'].drop_duplicates().to_frame()
    k4=k4['ncbi'].drop_duplicates().to_frame()
    my_dict[traits[i]]=k
    my_dict2[traits[i]]=genes
    my_dict3[traits[i]]=k3
    my_dict4[traits[i]]=k4
    known_genes=k3.shape[0]
    genes=genes.shape[0]
    A=len(df.merge(df_ref)['QTL'].drop_duplicates())
    total_t_qtl=len(df_ref['QTL'].drop_duplicates())
    B=total_t_qtl-A
    C=len(m['QTL'].drop_duplicates())-A
    D=total_qtl-A-B-C
    oddsratio,pvalue = stats.fisher_exact([[A, B], [C, D]],alternative='greater')
    results.append([traits[i],pvalue,A,B,C,D,genes,known_genes])
    

results_df=pd.DataFrame(results).rename({0:'trait',1:'p-value',2:'A', 3:'B',4:'C', 5:'D',6:'genes',7:'known_genes'},axis=1) 
pvals=np.array(results_df['p-value'])
pval_corr = smm.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)#https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/ 


enriched_traits=[]
flag=0
for l in results:
    if pval_corr[1][flag]<=0.05:
        enriched_traits.append([l[0],l[1],pval_corr[1][flag],l[2],l[3],l[4],l[5],l[6],l[7]])
    flag=flag+1

enriched_traits_df=pd.DataFrame(enriched_traits).rename({0:'trait',1:'pvalue',2:'adj-pvalue',3:'A', 4:'B',5:'C', 6:'D',7:'genes',8:'known_genes'},axis=1)             
enriched_traits_df=enriched_traits_df.sort_values(["adj-pvalue"])
enriched_traits_df.to_csv('QTL_traits_enriched_in_tissue_specific_genes',index = None, header=True,sep="\t")

out=[]
t=list(enriched_traits_df['trait'])
for trait in t:
    g=list(my_dict4[trait]['ncbi'])
    if len(g)>0:
        for gene in g:
            if gene not in out:
                out.append(gene)
df=pd.DataFrame(out)
df.to_csv('resutls.txt',index = None, header=False,sep="\t")    

"""












