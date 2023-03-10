#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 14:43:08 2021

@author: beiki
"""
import pandas as pd
#import dask.dataframe as dd

miranda=pd.read_csv('miranda/miranda-transcripts3UTR-target-prediction-output',sep='\t',names=['mirna','transcript_id'])
pita=pd.read_csv('PITA/PITA-transcripts3UTR-target-prediction-output',sep='\t',names=['mirna','transcript_id'])

a=(miranda['mirna'] + '---' + miranda['transcript_id']).to_frame()
b=(pita['mirna'] + '---' + pita['transcript_id']).to_frame()

m=a.merge(b)
#m=dd.merge(a, b)
m[['miRNA', 'transcript_id']] = m[0].str.split('---', 1, expand=True)
out=m[['miRNA', 'transcript_id']]
out.to_csv('../miranda-PITA_mirna-transcripts3UTR-predicted_tragets',index = None, header=True,sep="\t")
