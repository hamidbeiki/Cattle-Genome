#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 15:14:54 2020

@author: beiki
"""
"""
    This progam is based on Python 3
    """

""" 
    NOTE: on nova use
    module load miniconda3/4.3.30-qdauveb
    source activate /home/beiki/.conda/envs/py-libs
"""
import pandas as pd
import networkx
import functools
import operator

def pairs(lst):
    i = iter(lst)
    first = prev = item = i.__next__()# use i.next() in python 2.7
    for item in i:
        yield prev, item
        prev = item
    yield item, first

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])
    
data=pd.read_csv('Duplications.csv',header=0,sep="\t")# Duplications.csv is output fom orthofinder program

within_species_paralog=[]#https://useast.ensembl.org/info/genome/compara/homology_types.html
between_species_paralog=[]#the results of real duplications followed by gene losses
for i in range(data.shape[0]):
    l1=data.iloc[i,5].split(',')
    l2=data.iloc[i,6].split(',')
    res1=[x.split('_')[1] for x in l1 if 'FAANG' in x]
    res2=[x.split('_')[1] for x in l2 if 'FAANG' in x]
    if len(res1)>0 and len(res2)>0:
        within_species_paralog.append(res1+res2)
    elif len(res1)>0 and len(res2)==0:
        for j in res1:
            between_species_paralog.append(j)
    elif len(res1)==0 and len(res2)>0:
        for j in res2:
            between_species_paralog.append(j)            
        
g = networkx.Graph()#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
for sub_list in within_species_paralog:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

within_species_paralog_groups=list(networkx.connected_components(g))#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
within_species_paralog_genes=flatten_nested_list(within_species_paralog_groups)
df=pd.DataFrame(within_species_paralog_genes).rename({0:'gene_id'},axis=1)
df.to_csv('within_species_paralogus_genes',index = None, header=False,sep="\t") 

l = list(dict.fromkeys(between_species_paralog))
df2=pd.DataFrame(l).rename({0:'gene_id'},axis=1)
df2.to_csv('between_species_paralogus_genes',index = None, header=False,sep="\t") 


"""
data=pd.read_csv('cattleFAANG.pep.all__v__cattle.pep.all.csv',header=0,sep="\t")
result=[]
for i in range(data.shape[0]):
    l=data.iloc[i,1].split(',')
    for j in l:
        result.append(j)
"""        