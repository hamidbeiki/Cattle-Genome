#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 12:06:06 2021

@author: beiki
"""

""" 
    NOTE: on nova use
    module load miniconda3/4.3.30-qdauveb
    source activate /home/beiki/.conda/envs/py-libs
"""
import pandas as pd
import numpy as np
import networkx
import functools
import operator


def isNaN(num):
    return num != num

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])

def remove_space(l):
    out=[]
    for i in l:
        j = i.replace(' ','')
        out.append(j)
    return(out)
    
data=pd.read_csv('Orthogroups.csv',header=0,sep="\t")
orthologous_dict={}
orthologous_list=[]
for i in range(data.shape[0]):
    group=data.iloc[i,0]
    info=data.iloc[i,1:][85]
    if not isNaN(info):
        l=remove_space(list(info.split(',')))
        orthologous_dict[group]=l
        orthologous_list.append(l)

l=flatten_nested_list(orthologous_list)
df=pd.DataFrame(l)
np.save('CattleFAANG_orthologous.npy', orthologous_dict)
df.to_csv('cattleFAANG-ortholog-genes',index = None, header=False,sep="\t")
       
        