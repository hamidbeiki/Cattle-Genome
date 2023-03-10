#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 09:32:29 2021

@author: beiki
"""

import pandas as pd
import sys
import re

input_file=sys.argv[1]#Testis_bidirection_promoter_events or Abomasum_whole_convergent_gene_events
diverngent_type=sys.argv[2]
tissue=re.split('_bidirection|_convergent',input_file)[0]

data=pd.read_csv(input_file,header=0,sep='\t')
if diverngent_type=='bidirection':
    a=data.loc[(data['trA_biotype']=='protein_coding') & (data['trB_biotype']!='protein_coding')]
    a=a[['trB','trB_seq','trA','trA_seq']].rename({'trB':'A','trB_seq':'B','trA':'C','trA_seq':'D'},axis=1)
    b=data.loc[(data['trB_biotype']=='protein_coding') & (data['trA_biotype']!='protein_coding')]
    b=b[['trA','trA_seq','trB','trB_seq']].rename({'trA':'A','trA_seq':'B','trB':'C','trB_seq':'D'},axis=1)
    df=pd.concat([a,b])
elif diverngent_type=='convergent':
    a=data.loc[(data['trA_biotype']=='protein_coding_transcripts') & (data['trB_biotype']!='protein_coding_transcripts')]
    a=a[['trB','trB-seq','trA','trA-seq']].rename({'trB':'A','trB-seq':'B','trA':'C','trA-seq':'D'},axis=1)
    b=data.loc[(data['trB_biotype']=='protein_coding_transcripts') & (data['trA_biotype']!='protein_coding_transcripts')]
    b=b[['trA','trA-seq','trB','trB-seq']].rename({'trA':'A','trA-seq':'B','trB':'C','trB-seq':'D'},axis=1)
    df=pd.concat([a,b])

df.to_csv(tissue + '_' + diverngent_type + '_LncTar_input',index = None, header=False,sep="\t")