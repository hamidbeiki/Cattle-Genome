#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:27:47 2020

@author: beiki
"""
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp


data=pd.read_csv('Cerebral_Cortex_gene_recombinationHotspots-relation',header=0,sep="\t")
windows=data['window'].to_frame()
data=pd.read_csv('Cerebral_Cortex_gene_CNVs-relation',header=0,sep="\t")#or Cerebral_Cortex_gene_recombinationHotspots-relation or Cerebral_Cortex_gene_CNVs-relation
merge1=data.merge(windows)#because the number of windows in hotspots in 10 times smaller than CNV results.
x=np.array(merge1['background%'])
y=np.array(merge1['bidirection%'])
pvalue1=ks_2samp(x, y)[1]
print(pvalue1)
y=np.array(merge1['convergent%'])
pvalue2=ks_2samp(x, y)[1]
print(pvalue2)  
x=np.array(merge1['bidirection%'])
y=np.array(merge1['convergent%'])
pvalue3=ks_2samp(x, y)[1]      
print(pvalue3) 

data1=pd.read_csv('Cerebral_Cortex_gene_recombinationHotspots-relation',header=0,sep="\t")
windows=data['window'].to_frame()
data=pd.read_csv('Cerebral_Cortex_gene_CNVs-relation',header=0,sep="\t")#or Cerebral_Cortex_gene_recombinationHotspots-relation or Cerebral_Cortex_gene_CNVs-relation
data2=data.merge(windows)#because the number of windows in hotspots in 10 times smaller than CNV results.
x=np.array(data1['bidirection%'])
y=np.array(data2['bidirection%'])
pvalue=ks_2samp(x, y)[1]
print(pvalue)
x=np.array(data1['convergent%'])
y=np.array(data2['convergent%'])
pvalue=ks_2samp(x, y)[1]
print(pvalue)