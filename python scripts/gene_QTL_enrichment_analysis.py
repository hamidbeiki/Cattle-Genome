#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:20:37 2021

@author: beiki
"""

import pandas as pd
import collections
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
import math

""" requiered modules:
    module load py-pandas/0.21.1-py2-326uzkn
    module load py-scipy/1.0.0-py2-flqcuxg
    module load py-matplotlib/2.0.2-py2-m755aws
"""

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_gene_startEnd(input):
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    chromosome=df.iloc[0,0]
    gene_end=df.iloc[0,2]
    gene_start=int(gene_start)
    gene_end=int(gene_end)
    strand=df.iloc[0,3]
    return(chromosome,gene_start,gene_end,strand)
    
def get_genes_window(gene_start,gene_end,window,strand):
#if you are intrested in to put window jus on genes 5prime end, replace get_genes_window 
#function with this function 
    if strand=="+":
        if gene_start>1 and gene_start<=window:
            window_start=1
            window_end=gene_end
        elif gene_start>window:
            window_start=gene_start-window
            window_end=gene_end
        elif gene_start==1:
            window_start=1
            window_end=gene_end
    elif strand=="-":
        s=gene_end
        e=gene_start
        window_start=e
        window_end=s+window
    return(window_start,window_end)

def get_geneWindow_dict(input,gene,window_size):
    '''no need to filter for minimum_exon_length because transcripts has been already firtered for this criteria in upstream steps
        this will atually consider the last exons on genes in "-" strand as first exon that is fine because I'm lookig for divergent transcription'''
    geneBorders_dict={}
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    strand=df.iloc[0,3]
    gene_end=df.iloc[0,2]
    gene_start=int(gene_start)
    gene_end=int(gene_end)
    info=get_genes_window(gene_start,gene_end,window_size,strand)
    window_start=info[0]
    window_end=info[1]
    geneBorders_dict[gene] = "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0],window_start,window_end,gene,99,df.iloc[0,3],gene)      
    return(geneBorders_dict)

def get_geneWindow_bed(candidates_dict,window_size):
    out_dict={}
    for gene in candidates_dict:
        gene_dict=get_geneWindow_dict(candidates_dict[gene],gene,window_size)
        out_dict.update(gene_dict)
    out_bed=pd.DataFrame.from_dict(out_dict,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'promoter_start', 2: 'gene_end', 3: 'gene_id', 4: 'value', 5: 'strand', 6: 'gene_id'},axis=1)
    return(out_bed) 

def get_genes_closest_qtl(candidates_bed):
    global qtl_dict
    genes_closest_qtl={} 
    for i in range(candidates_bed.shape[0]):
        info=candidates_bed.iloc[i,0:]
        gene=info[3]
        chrom=info[0]
        start=int(info[1])
        end=int(info[2])
        targets=qtl_dict[chrom]
        flag='A'
        for key, value in targets.iteritems():#key is qtl name and value is a list containing it's coordinates
            if value[0]>=start and value[0]<=end:#all qtls in CattleQTLdb has length of 4nt, so I just checked for their start (value[0])
                distance1=start-value[0]
                distance2=value[0]-end
            else:
                distance1=abs(start-value[0])
                distance2=abs(end-value[0])
            dist=[distance1,distance2]
            dist.sort()
            dist=dist[0]
            if dist<flag:
                flag=dist
                qtl=key
        if dist<0:# if a gene locate inside a gene, dist will be negative. To get a correct estimate of dist median we conver this to 0
            dist=0
        genes_closest_qtl[gene]=[qtl,dist]#this gives the distance of the closest qtls to gene
    return(genes_closest_qtl)

def get_gene_qtl_distances(candidates_bed):
    global qtl_dict
    genes_qtl_dist={} 
    for i in range(candidates_bed.shape[0]):
        info=candidates_bed.iloc[i,0:]
        gene=info[3]
        chrom=info[0]
        start=int(info[1])
        end=int(info[2])
        targets=qtl_dict[chrom]
        for key, value in targets.iteritems():#key is qtl name and value is a list containing it's coordinates
            if value[0]>=start and value[0]<=end:#all qtls in CattleQTLdb has length of 4nt, so I just checked for their start (value[0])
                distance1=start-value[0]
                distance2=value[0]-end
            else:
                distance1=abs(start-value[0])
                distance2=abs(end-value[0])
            dist=[distance1,distance2]
            dist.sort()
            dist=dist[0]
            if dist<0:# if a gene locate inside a gene, dist will be negative. To get a correct estimate of dist median we conver this to 0
                dist=0
            ID=gene + '_' + key
            genes_qtl_dist[ID]=[gene,key,dist]#this gives the distance of the closest qtls to gene
    genes_qtl_dist_df=pd.DataFrame.from_dict(genes_qtl_dist,orient='index').rename({0:'gene_id',1: 'QTL', 2: 'all_dist'},axis=1)
    return(genes_qtl_dist_df)

""" RUN """
qtl_bed=pd.read_csv('CattleQTLdb_April2020.bed',header=0,sep="\t")
tissues=list(pd.read_csv('tissues',names=['tissue'])['tissue'])
input_gff=('combined.gff')
ex_coords = nested_dict()
for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        gid=raw[-1].split('; ')[0].split()[1]
        gid=re.sub('["]', '', gid)
        exnumber=raw[3]
        ex_coords[gid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
        

candidates_dict={}
for gene in ex_coords:
    candidates_dict[gene]=ex_coords[gene]
        
window_size=0
candidates_bed=get_geneWindow_bed(candidates_dict,window_size)
all_genes_df=candidates_bed.iloc[0:,3].drop_duplicates().to_frame()

qtl_dict=nested_dict()
for i in range(qtl_bed.shape[0]):
    info=qtl_bed.iloc[i,0:]
    qtl_dict[info[0]][info[3]]=[info[1],info[2]]
    
main_all_df=get_gene_qtl_distances(candidates_bed)

all_df=main_all_df.loc[main_all_df.groupby('QTL')['all_dist'].idxmin(), :]
subset=all_df.loc[all_df['all_dist']>0]
y=np.log(subset['all_dist'])
#threshold=y.mean()-2*y.std()#6.758886487783972 rhat is equal to threshold**2.71828=862 bases
#threshold=np.log(1)## in bellow codes, all 0 distances converted to 1, so this should be OK
plt.hist(y, bins=60, alpha=1)
plt.savefig('QTL_enrichment_threshold.png')
plt.close()

main_all_df.to_csv('gene_QTL_distances',index = None, header=True,sep="\t")
all_df['trait'] = all_df['QTL'].str.split('\_QTL\_\(', 1).str[0]
all_df.to_csv('QTLs_closest_gene',index = None, header=True,sep="\t")

all_df['log_dist']=np.log(all_df['all_dist']+1)
sign=all_df.loc[all_df['log_dist']<=threshold].copy()
sign['trait'] = sign['QTL'].str.split('\_QTL\_\(', 1).str[0]
sign.to_csv('gene_significant_QTLs',index = None, header=True,sep="\t")

for tissue in tissues:
    genes=pd.read_csv('quantification/'+ tissue + '_genes',names=['gene_id'])
    m=main_all_df.merge(genes)
    all_df=m.loc[m.groupby('QTL')['all_dist'].idxmin(), :]
    all_df['trait'] = all_df['QTL'].str.split('\_QTL\_\(', 1).str[0]
    all_df.to_csv(tissue + '_gene_closest_QTLs',index = None, header=True,sep="\t")
    sign=all_df# because Jim asked me to remove threshold
    traits=list(sign['trait'].drop_duplicates())
    trait_dict={}
    for i in traits:
        df=sign.loc[sign['trait']==i]
        trait_dict[i]=df        
    genes=list(sign['gene_id'].drop_duplicates())        
    gene_dict={}
    for i in genes:
        df=sign.loc[sign['gene_id']==i]
        l=list(df['trait'])
        l=list(dict.fromkeys(l))
        gene_dict[i]=[l,df.shape[0]]    
    total_qtils=qtl_bed.shape[0]
    results=[]
    for gene in genes:
        info=gene_dict[gene]
        traits=info[0]
        qtls=info[1]
        for trait in traits:
            df1=trait_dict[trait]
            df2=df1.loc[df1['gene_id']==gene]
            if df2.shape[0]>=3:# minimum of three genes required
                a=df2.shape[0]
                b=df1.shape[0]-a
                c=qtls-a
                d=total_qtils-a-b-c
                obs = np.array([[a, b], [c, d]])
                oddsratio, p=fisher_exact(obs, alternative='greater')
                results.append([gene,trait,p,a,b,c,d])    
    results_df=pd.DataFrame(results).rename({0:'gene',1:'trait',2:'p-value',3:'A', 4:'B',5:'C', 6:'D'},axis=1)            
    pvals=np.array(results_df['p-value'])
    pval_corr = smm.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)#https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/ 
    enriched_genes=[]
    flag=0
    for l in results:
        if pval_corr[1][flag]<=0.05:
            enriched_genes.append([l[0],l[1],l[2],pval_corr[1][flag],l[3],l[4],l[5],l[6]])
        flag=flag+1
    enriched_genes_df=pd.DataFrame(enriched_genes).rename({0:'gene',1:'trait',2:'pvalue',3:'adj-pvalue',4:'A', 5:'B',6:'C', 7:'D'},axis=1)             
    enriched_genes_df=enriched_genes_df.sort_values(["adj-pvalue"])
    enriched_genes_df.to_csv(tissue + '_genes_enriched_with_QTLs_without_distance_threshold',index = None, header=True,sep="\t")
    print(tissue)
    
"""" run this if you want to have the analysis genomewide
del(sign)
sign=all_df# because Jim asked me to remove threshold

traits=list(sign['trait'].drop_duplicates())
trait_dict={}
for i in traits:
    df=sign.loc[sign['trait']==i]
    trait_dict[i]=df
    
genes=list(sign['gene_id'].drop_duplicates())        
gene_dict={}
for i in genes:
    df=sign.loc[sign['gene_id']==i]
    l=list(df['trait'])
    l=list(dict.fromkeys(l))
    gene_dict[i]=[l,df.shape[0]]

total_qtils=qtl_bed.shape[0]
results=[]
for gene in genes:
    info=gene_dict[gene]
    traits=info[0]
    qtls=info[1]
    for trait in traits:
        df1=trait_dict[trait]
        df2=df1.loc[df1['gene_id']==gene]
        if df2.shape[0]>=3:# minimum of three genes required
            a=df2.shape[0]
            b=df1.shape[0]-a
            c=qtls-a
            d=total_qtils-a-b-c
            obs = np.array([[a, b], [c, d]])
#            g, p, dof, expctd = chi2_contingency(obs)#for description of enrichment analysis please see #https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/
            oddsratio, p=fisher_exact(obs, alternative='greater')
            results.append([gene,trait,p,a,b,c,d])
#            results=results.append(pd.DataFrame({'gene': gene,'trait':trait, 'pvalue':p, 'A':a,'B':b,'C':c,'D':d}, index=[0]), ignore_index=True)
results_df=pd.DataFrame(results).rename({0:'gene',1:'trait',2:'p-value',3:'A', 4:'B',5:'C', 6:'D'},axis=1)            
pvals=np.array(results_df['p-value'])
pval_corr = smm.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)#https://bbglab.irbbarcelona.org/courses/bco/documents/enrichment.pdf OR https://online.stat.psu.edu/stat555/node/14/ 
enriched_genes=[]
flag=0
for l in results:
    if pval_corr[1][flag]<=0.05:
        enriched_genes.append([l[0],l[1],l[2],pval_corr[1][flag],l[3],l[4],l[5],l[6]])
    flag=flag+1
#in case if you used chi-square test
#enriched_genes=[]
#flag=0
#for l in results:
#    if l[3]>0 and l[4]>0 and l[3]/float(l[4]) >l[5]/float(l[6]):
#        enriched_genes.append([l[0],l[1],l[2],pval_corr[1][flag],l[3],l[4],l[5],l[6]])
#    elif l[3]>0 and l[4]==0:
#        enriched_genes.append([l[0],l[1],l[2],pval_corr[1][flag],l[3],l[4],l[5],l[6]])
#    flag=flag+1
   
enriched_genes_df=pd.DataFrame(enriched_genes).rename({0:'gene',1:'trait',2:'pvalue',3:'adj-pvalue',4:'A', 5:'B',6:'C', 7:'D'},axis=1)             
enriched_genes_df.to_csv('genes_enriched_with_QTLs_without_distance_threshold',index = None, header=True,sep="\t")
""""            
































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


#https://online.stat.psu.edu/stat555/node/14/
#https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.chi2_contingency.html