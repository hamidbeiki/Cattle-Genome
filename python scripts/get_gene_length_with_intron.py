# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd

df=pd.read_csv('out',sep=' ',names=["chr","tr_start","tr_end","gene_id"])
gene_ends=df.groupby(['gene_id'], sort=False)['tr_end'].max()
genes_start=df.groupby(['gene_id'], sort=False)['tr_start'].min()
genes_borders=pd.concat([genes_start,gene_ends],axis=1)
genes_borders['length']=genes_borders['tr_end']-genes_borders['tr_start']
genes_borders=genes_borders.reset_index().rename({'tr_start':'gene_start', 'tr_end':'gene_end'},axis=1)
print(genes_borders['length'].mean())

'''to sum all columns'''
#genes_borders.sum(axis = 0, skipna = True)
