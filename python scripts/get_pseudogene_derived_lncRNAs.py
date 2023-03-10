#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:39:06 2021

@author: beiki
"""

import pandas as pd

lncRNAs=pd.read_csv('b',names=['gene','lncRNA'],sep='\t')#file b is a file after adding gene_id to lncRNAs file
pseudogenes=pd.read_csv('pseudogenes',names=['gene'],sep='\t')
m=lncRNAs.merge(pseudogenes)['lncRNA'].to_frame()
m.to_csv('pseudogene_derived_lncRNAs',index = None, header=True,sep="\t")