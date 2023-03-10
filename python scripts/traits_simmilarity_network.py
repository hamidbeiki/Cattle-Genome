#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 18:27:40 2021

@author: beiki
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from random import sample

def get_subset(trait,qtls):
    count=0
    for t in qtls:
        if t==trait:
            count=count+1
    return(count)
    
def get_qtls(rand_genes):
    global gene_qtls
    res=[]
    for gene in rand_genes:
        qtls=gene_qtls[gene]
        for qtl in qtls:
            res.append(qtl)
    return(res)
    
def random_sampling(trait,sample_length,target):
    global genes
    res=[]
    for i in range(1000):
        rand_genes=sample(genes,sample_length)
        rand_genes_qtls=get_qtls(rand_genes)
        shared_qtls=get_subset(trait,rand_genes_qtls)
        res.append(shared_qtls)
    pvalue=len([x for x in res if x > target])/1000
    return(pvalue)

QTLs_closest_gene=pd.read_csv('QTLs_closest_gene',header=0,sep='\t')
gene_enriched_with_QTLs=pd.read_csv('genes_enriched_with_QTLs_without_distance_threshold',header=0,sep='\t')
gene_biotypes=pd.read_csv("final_gene_biotypes",names=['gene_id','biotype'],sep='\t')
known_genes=pd.read_csv('annotation/all-known-genes-main',names=['gene_id'])
remove_genes=pd.read_csv('remove_genes',names=['gene_id'])
AS_event_genes=pd.read_csv('AS-events/total_AS_event_per_gene',header=0,sep='\t')
qtl_bed=pd.read_csv('CattleQTLdb_April2020.bed',names=['chr','start','end','qtl','info1','strand','info2'],sep="\t")
qtl_bed[['trait','a']]=qtl_bed.qtl.str.split('\_QTL\_',expand=True)
QTLs_closest_gene=QTLs_closest_gene.loc[~QTLs_closest_gene.gene_id.isin(remove_genes.gene_id)]

QTLs_closest_gene.merge(gene_biotypes).groupby(['biotype'])['all_dist'].median()


df=QTLs_closest_gene.loc[QTLs_closest_gene['all_dist']>0]
y=np.log2(df['all_dist'])
plt.hist(y, bins=60, alpha=1)
plt.xlabel("QTL distance to closest gene (log2 based)")
plt.ylabel("Frequency")
plt.savefig('QTL_distance_to_closest_genes.png',dpi=300)
plt.close()

2**y.median()
df=QTLs_closest_gene.loc[QTLs_closest_gene['all_dist']==0]['gene_id'].drop_duplicates().to_frame()
df.merge(known_genes)

QTL_associated_genes=QTLs_closest_gene['gene_id'].drop_duplicates().to_frame()
QTL_associated_genes.merge(known_genes)

a=AS_event_genes.merge(QTL_associated_genes)['total_AS_event'].to_frame()
b=AS_event_genes.loc[~AS_event_genes.gene_id.isin(QTL_associated_genes)]['total_AS_event'].to_frame()

stat,pval=stats.mannwhitneyu(a['total_AS_event'],b['total_AS_event'],alternative='greater')


traits=list(QTLs_closest_gene['trait'].drop_duplicates())
traits_genes={}
for trait in traits:
    df=QTLs_closest_gene.loc[QTLs_closest_gene['trait']==trait]['gene_id'].drop_duplicates().to_frame()
    traits_genes[trait]=df

trait_sim_list=[]
for i in range(len(traits)):
    res=[]
    for j in range(len(traits)):
        t1=list(traits_genes[traits[i]]['gene_id'])
        t2=list(traits_genes[traits[j]]['gene_id'])
        all_genes=list(dict.fromkeys(t1+t2))
        shared_genes=list(dict.fromkeys(list(set(t1) & set(t2))))
        p=(len(shared_genes)/len(all_genes))*100
        res.append(p)
        print([i,j])
    trait_sim_list.append(res)

trait_sim_mat=pd.DataFrame(trait_sim_list)
trait_sim_mat.columns=traits
trait_sim_mat.index=traits


traits_qtls=qtl_bed.groupby(['trait'])['qtl'].size().reset_index()
a=list(traits_qtls['trait'])
b=list(traits_qtls['qtl'])
traits_qtls={}
for i in range(len(a)):
    traits_qtls[a[i]]=b[i]

gene_qtls=QTLs_closest_gene.groupby(['gene_id'])['trait'].apply(list).to_dict()
genes=gene_qtls.keys()
trait_sim_list2=[]# traitA-->traitB
for i in range(len(traits)):
    for j in range(i+1,len(traits)):
        t1=list(traits_genes[traits[i]]['gene_id'])
        t2=list(traits_genes[traits[j]]['gene_id'])
        t1_qtls=traits_qtls[traits[i]]
        t2_qtls=traits_qtls[traits[j]]
        if t1_qtls>=10 and t2_qtls>=10:
            t2_genes_all_qtls=get_qtls(t2)
            shared_qtls=get_subset(traits[i],t2_genes_all_qtls)
            if shared_qtls>=10:
                pvalue=random_sampling(traits[i],len(t2),shared_qtls)
                trait_sim_list2.append([traits[i],traits[j],pvalue,t1_qtls,shared_qtls]) 
        print([i,j])

df1=pd.DataFrame(trait_sim_list2).rename({0:'traitA',1:'traitB',2:'pvalue',3:'traitA_qtls',4:'traitA_qtls_close_to_traitB_qtl_associated_genes'},axis=1)  

trait_sim_list3=[]#traitB-->traitA
for i in range(len(traits)):
    for j in range(i+1,len(traits)):
        t1=list(traits_genes[traits[i]]['gene_id'])
        t2=list(traits_genes[traits[j]]['gene_id'])
        t1_qtls=traits_qtls[traits[i]]
        t2_qtls=traits_qtls[traits[j]]
        if t1_qtls>=10 and t2_qtls>=10:
            t1_genes_all_qtls=get_qtls(t1)
            shared_qtls=get_subset(traits[j],t1_genes_all_qtls)
            if shared_qtls>=10:
                pvalue=random_sampling(traits[j],len(t1),shared_qtls)
                trait_sim_list3.append([traits[j],traits[i],pvalue,t1_qtls,shared_qtls]) 
        print([i,j])

df2=pd.DataFrame(trait_sim_list3).rename({0:'traitA',1:'traitB',2:'pvalue',3:'traitA_qtls',4:'traitA_qtls_close_to_traitB_qtl_associated_genes'},axis=1)  

cytoscape_input=pd.concat([df1,df2])
cytoscape_input.to_csv('QTL-based_trait_similarities_cytoscape_input.txt',index = None, header=True,sep="\t")#http://manual.cytoscape.org/en/stable/Styles.html        

cytoscape_input2=cytoscape_input.loc[cytoscape_input['pvalue']<0.001]
cytoscape_input2.to_csv('QTL-based_SIGNIFICANT_trait_similarities_cytoscape_input.txt',index = None, header=True,sep="\t")#http://manual.cytoscape.org/en/stable/Styles.html        



"""
traits_qtls=qtl_bed.groupby(['trait'])['qtl'].size().reset_index()
a=list(traits_qtls['trait'])
b=list(traits_qtls['qtl'])
traits_qtls={}
for i in range(len(a)):
    traits_qtls[a[i]]=b[i]

trait_sim_list2=[]
for i in range(len(traits)):
    for j in range(i+1,len(traits)):
        t1=list(traits_genes[traits[i]]['gene_id'])
        t2=list(traits_genes[traits[j]]['gene_id'])
        t1_qtls=traits_qtls[traits[i]]
        t2_qtls=traits_qtls[traits[j]]
        if t1_qtls>=10 and t2_qtls>=10:
            all_genes=list(dict.fromkeys(t1+t2))
            shared_genes=list(dict.fromkeys(list(set(t1) & set(t2))))
            if len(shared_genes)>=10:
                df=pd.DataFrame(shared_genes).rename({0:'gene_id'},axis=1).merge(QTLs_closest_gene)
                shared_t1_qtls=df.loc[df['trait']==traits[i]].shape[0]
                shared_t2_qtls=df.loc[df['trait']==traits[j]].shape[0]
                if shared_t1_qtls>0 and shared_t2_qtls>0:
                    p=(len(shared_genes)/len(all_genes))*100
                    p1=(shared_t1_qtls/t1_qtls)*100
                    p2=(shared_t2_qtls/t2_qtls)*100
                    trait_sim_list2.append([traits[i] + '--' + traits[j],len(shared_genes),p,t1_qtls,p1,t2_qtls,p2])
        print([i,j])

trait_sim_mat2=pd.DataFrame(trait_sim_list2).rename({0:'traits',1:'#of_shared_genes',2:'%of_both_tissue_genes_shared',3:'traitA_qtls',4:'%traitA_qtls_associated_with_shared_genes',5:'traitB_qtls',6:'%traitB_qtls_associated_with_shared_genes'},axis=1)
df=trait_sim_mat2.loc[(trait_sim_mat2['%traitA_qtls_associated_with_shared_genes']>50) & (trait_sim_mat2['%traitB_qtls_associated_with_shared_genes']>50)] 
df[['traitA','traitB']]=df.traits.str.split('\-\-',expand=True)
del df['traits']
cytoscape_input=df[['traitA','traitB','#of_shared_genes', '%of_both_tissue_genes_shared', 'traitA_qtls','%traitA_qtls_associated_with_shared_genes', 'traitB_qtls','%traitB_qtls_associated_with_shared_genes']]
cytoscape_input.to_csv('QTL-based_trait_similarities_cytoscape_input.txt',index = None, header=True,sep="\t")#http://manual.cytoscape.org/en/stable/Styles.html




"""
'''
traits_qtls=qtl_bed.groupby(['trait'])['qtl'].size().reset_index()
trait_sim_list2=[]
for i in range(len(traits)):
    for j in range(i+1,len(traits)):
        t1=list(traits_genes[traits[i]]['gene_id'])
        t2=list(traits_genes[traits[j]]['gene_id'])
        all_genes=list(dict.fromkeys(t1+t2))
        shared_genes=list(dict.fromkeys(list(set(t1) & set(t2))))
        p=(len(shared_genes)/len(all_genes))*100
        trait_sim_list2.append([traits[i] + '--' + traits[j],p])
        print([i,j])
        
trait_sim_mat2=pd.DataFrame(trait_sim_list2).rename({0:'traits',1:'%_of_shared_genes'},axis=1)
summary=trait_sim_mat2.groupby(['traits']).max()
df=summary.loc[summary['%_of_shared_genes']>0].reset_index()
df[['traitA','traitB']]=df.traits.str.split('\-\-',expand=True)
del df['traits']
cytoscape_input=df[['traitA','traitB','%_of_shared_genes']]
t=list(cytoscape_input['traitA'])+list(cytoscape_input['traitB'])
t=list(dict.fromkeys(t))
res=[]
for trait in t:
    all_q=qtl_bed.loc[qtl_bed['trait']==trait].shape[0]
#    ass_q=gene_closest_QTLs[ gene_closest_QTLs['trait']==trait].shape[0]
    dist=QTLs_closest_gene[ QTLs_closest_gene['trait']==trait]['all_dist'].median()
    res.append([trait,all_q,dist])
info=pd.DataFrame(res).rename({0:'trait',1:'number_of_associated_qtls',2:'median_qtl_dist'},axis=1)
info=info.rename({'trait':'traitA'},axis=1)
m=cytoscape_input.merge(info)    
    



cytoscape_input.to_csv('QTL-based_trait_similarities_cytoscape_input.txt',index = None, header=True,sep="\t")#http://manual.cytoscape.org/en/stable/Styles.html

df=summary.loc[summary['%_of_shared_genes']>50].reset_index()
df[['traitA','traitB']]=df.traits.str.split('\-\-',expand=True)
del df['traits']
cytoscape_input=df[['traitA','traitB','%_of_shared_genes']]
cytoscape_input.to_csv('QTL-based_highly_simmilar_traits_cytoscape_input.txt',index = None, header=True,sep="\t")
'''

