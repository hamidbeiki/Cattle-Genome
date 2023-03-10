#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:05:09 2021

@author: beiki
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

df=pd.read_csv('out',names=['read','dist_percentage'],sep=' ')#bedtools intersect -a Adipose_sub-cutaneous_WTTS_reads_unmapped_to_TTS.bed -b combined_exons.bed -wb  | awk '{print $4, $13}' | sed 's/\-\-/\t/g' | awk '{print $1, $3}' | sort | uniq >out
tissue=sys.argv[1]
df=df.groupby(['read'], sort=False)['dist_percentage'].max().reset_index()


df['groups']=pd.cut(df['dist_percentage'],
                      bins=[-1,40.,60.,80.,np.inf],
                      labels=[1,2,2.1,3])

df["groups"].replace({1: 3}, inplace=True)#group 1 is actually single exon transcripts, so it shoul be changed to last exon
df["groups"].replace({2:1,2.1:1}, inplace=True)# group 1 is middle exons, infact this group is last exons too but as I selected "max" at line "16" I filtered out "0.000%" which belong to unslpiced transcripts
df["groups"].replace({3: 2}, inplace=True)# group2 is last exons

plt.hist(df['groups'])
plt.savefig(tissue + '_WTTS_reads_unmapped_to_TTS_mapped_to_exons_hist.png')
plt.close()

df2=df['groups'].to_frame()
info=df2.groupby(df2.columns.tolist(),as_index=False).size()
info["groups"].replace({1.0: "first_exon",2.0:"middle_exon",3.0:"last_exon"}, inplace=True)
info.to_csv(tissue + '_WTTS_reads_unmapped_to_TTS_mapped_to_exons_hist.txt',index = None, header=True,sep="\t")
