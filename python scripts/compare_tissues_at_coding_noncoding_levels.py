#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:26:57 2020

@author: beiki
"""

# This has been written in python 3
import pandas as pd
import re
import numpy as np
import sys
import collections  
import csv

def format_mat(mat_dict):
    global tissues
    mat=pd.DataFrame.from_dict(mat_dict,orient='index')
    mat.columns=tissues
    mat=mat.reset_index()
    return(mat)
    
#### RUN ####
input_files='input_files'
files=pd.read_csv(input_files,sep="\t",names=['file'])

tissues=[]
for i in range(files.shape[0]):
    file=files.iloc[i,0]
    tissue=re.split('mRNA-|_final|.collapsed', file)[1]
    tissues.append(tissue)

file1=sys.argv[1]#combined_transcript_biotypes
file2=sys.argv[2]#combined_genes_biotypes
tissue_tr = np.load('tissue_transcript.npy',allow_pickle='TRUE').item()   
tr_tissue = np.load('transcript_tissue.npy',allow_pickle='TRUE').item()  
tissue_gene = np.load('tissue_gene.npy',allow_pickle='TRUE').item()
gene_tissue = np.load('gene_tissue.npy',allow_pickle='TRUE').item()
t_biotypes=pd.read_csv(file1,header=0,sep="\t")
g_biotypes=pd.read_csv(file2,header=0,sep="\t")


non_coding_genes=[]
protein_coding_genes=[]
for i in range(g_biotypes.shape[0]):
    g=g_biotypes.iloc[i,0]
    b=g_biotypes.iloc[i,1]
    if b=="genes_with_coding_noncoding_tr" or b=="genes_with_just_coding_tr":
        protein_coding_genes.append(g)
    else:
        non_coding_genes.append(g)
non_coding_trs=[]
protein_coding_trs=[]
for i in range(t_biotypes.shape[0]):
    t=t_biotypes.iloc[i,0]
    b=t_biotypes.iloc[i,1]
    if b=="protein_coding":
        protein_coding_trs.append(t)
    elif b!="protein_coding":
        non_coding_trs.append(t)

""" TRANSCRIPT level """       
tissue_sim_tr_coding=collections.OrderedDict()
tissue_sim_tr2_coding=collections.OrderedDict()#this is for heatmap
tissue_sim_tr_noncoding=collections.OrderedDict()
tissue_sim_tr2_noncoding=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    v3=[]
    v4=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_tr[t1]
        b=tissue_tr[t2]
        value=list(set(a) & set(b))
        value_coding=list(set(value) & set(protein_coding_trs))
        r1=list(set(a) & set(protein_coding_trs))
        v.append(len(value_coding))
        v2.append(round((len(value_coding)/float(len(r1)))*100,0))
        value_noncoding=list(set(value) & set(non_coding_trs))
        r2=list(set(a) & set(non_coding_trs))
        v3.append(len(value_noncoding))
        v4.append(round((len(value_noncoding)/float(len(r2)))*100,0))
    tissue_sim_tr_coding[tissues[i]]=v
    tissue_sim_tr2_coding[tissues[i]]=v2
    tissue_sim_tr_noncoding[tissues[i]]=v3
    tissue_sim_tr2_noncoding[tissues[i]]=v4
    
tissue_sim_mat_tr_coding=format_mat(tissue_sim_tr_coding)
tissue_sim_mat_tr2_coding=format_mat(tissue_sim_tr2_coding)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)
tissue_sim_mat_tr_noncoding=format_mat(tissue_sim_tr_noncoding)
tissue_sim_mat_tr2_noncoding=format_mat(tissue_sim_tr2_noncoding)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)
          
""" GENE level """
tissue_sim_gene_coding=collections.OrderedDict()
tissue_sim_gene2_coding=collections.OrderedDict()#this is for heatmap
tissue_sim_gene_noncoding=collections.OrderedDict()
tissue_sim_gene2_noncoding=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    v3=[]
    v4=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_gene[t1]
        b=tissue_gene[t2]
        value=list(set(a) & set(b))
        value_coding=list(set(value) & set(protein_coding_genes))
        r1=list(set(a) & set(protein_coding_genes))
        v.append(len(value_coding))
        v2.append(round((len(value_coding)/float(len(r1)))*100,0))
        value_noncoding=list(set(value) & set(non_coding_genes))
        r2=list(set(a) & set(non_coding_genes))
        v3.append(len(value_noncoding))
        v4.append(round((len(value_noncoding)/float(len(r2)))*100,0))
    tissue_sim_gene_coding[tissues[i]]=v
    tissue_sim_gene2_coding[tissues[i]]=v2
    tissue_sim_gene_noncoding[tissues[i]]=v3
    tissue_sim_gene2_noncoding[tissues[i]]=v4
    
tissue_sim_mat_gene_coding=format_mat(tissue_sim_gene_coding)
tissue_sim_mat_gene2_coding=format_mat(tissue_sim_gene2_coding)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)
tissue_sim_mat_gene_noncoding=format_mat(tissue_sim_gene_noncoding)
tissue_sim_mat_gene2_noncoding=format_mat(tissue_sim_gene2_noncoding)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)


tissue_sim_mat_tr_coding.to_csv('transcript_level_tissue_comparisions_coding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_tr2_coding.to_csv('transcript_level_tissue_comparisions_coding_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal) 
tissue_sim_mat_tr_noncoding.to_csv('transcript_level_tissue_comparisions_noncoding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_tr2_noncoding.to_csv('transcript_level_tissue_comparisions_noncoding_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal) 


tissue_sim_mat_gene_coding.to_csv('gene_level_tissue_comparisions_coding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_gene2_coding.to_csv('gene_level_tissue_comparisions_HEATMAP_coding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)  
tissue_sim_mat_gene_noncoding.to_csv('gene_level_tissue_comparisions_noncoding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_gene2_noncoding.to_csv('gene_level_tissue_comparisions_HEATMAP_noncoding.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)  
