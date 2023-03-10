#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Aug 25 10:17:57 2020

@author: beiki
"""
import re
import pandas as pd
import sys
import collections
import csv
from collections import OrderedDict

def format_mat(mat_dict):
    global tissues
    mat=pd.DataFrame.from_dict(mat_dict,orient='index')
    mat.columns=tissues
    mat=mat.reset_index()
    return(mat)
    
input_files=sys.argv[1]#a file containing list of files that contained list of expressed genes (known genes, e.g Ensembl or NCBI, in each tissue
files=pd.read_csv(input_files,sep="\t",names=['file'])


tissue_gene={}# detected genes per tissue
tissues=[]
for i in range(files.shape[0]):
    file=files.iloc[i,0]
    data=pd.read_csv(file,sep="\t",names=['genes'])
    tissue=re.split('mRNA-|_expressed|_genes.txt', file)[1]
    tissues.append(tissue)
    tissue_gene[tissue]=data

""" To re-order tissues so that tissue be in the same order as transcriptome assembly results
tissues_trinity=pd.read_csv('tissues_trinity.txt',names=['tissues'],sep="\t")
t=[]
for i in range(tissues_trinity.shape[0]):
    a=tissues_trinity.iloc[i,0]
    t.append(a)
t2=[]
for i in t:
    if i in tissues:
        t2.append(i)
tissues=t2
"""    
#get tissue pairs
count=0
tissue_pairs={}
for i in range(len(tissues)-1):
    for j in range(i+1,len(tissues)):
        tissue_pairs[count]=[tissues[i],tissues[j]]
        count=count+1

tissue_sim_gene=collections.OrderedDict()
tissue_sim_gene2=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_gene[t1]
        b=tissue_gene[t2]
        value=a.merge(b)
        v.append(value.shape[0])
        v2.append(round((value.shape[0]/float(a.shape[0]))*100,0))
    tissue_sim_gene[tissues[i]]=v
    tissue_sim_gene2[tissues[i]]=v2
    
tissue_sim_mat_gene=format_mat(tissue_sim_gene)
tissue_sim_mat_gene2=format_mat(tissue_sim_gene2)
tissue_sim_mat_gene.to_csv('gene_level_tissue_comparisions.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_gene2.to_csv('gene_level_tissue_comparisions_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)  


""" Rplot
library(gplots)
jpeg(file="Rplot-transcript_level_tissue_comparisions_HEATMAP.jpeg",width=15, height=15, units="in", res=300)
results=as.matrix(read.table(file="transcript_level_tissue_comparisions_HEATMAP.txt",header=TRUE,row.names=1))
my_palette <- colorRampPalette(c("white", "orange", "red"))(n=nrow(results)^2)## for less color details: n=100
heatmap.2(results,margins=c(15,15),trace="none",col=my_palette,srtCol=45,dendrogram='none',cexCol=0.9,cexRow=0.9,keysize=0.75,Rowv=FALSE,Colv=FALSE)
dev.off()

jpeg(file="Rplot-gene_level_tissue_comparisions_HEATMAP.jpeg",width=15, height=15, units="in", res=300)
results=as.matrix(read.table(file="gene_level_tissue_comparisions_HEATMAP.txt",header=TRUE,row.names=1))
my_palette <- colorRampPalette(c("white", "orange", "red"))(n=nrow(results)^2)## for less color details: n=100
heatmap.2(results,margins=c(15,15),trace="none",col=my_palette,srtCol=45,dendrogram='none',cexCol=0.9,cexRow=0.9,keysize=0.75,Rowv=FALSE,Colv=FALSE)
dev.off()
"""