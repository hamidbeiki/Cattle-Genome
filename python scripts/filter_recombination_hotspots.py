#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 13:41:19 2020

@author: beiki
"""

import pandas as pd

data=pd.read_csv('data2',header=0,sep='\t')# datais a file containing the all columns of 'cattle_recombinationHotspots_PMC5923192_Additional_file_1.xlsx' file
a=pd.concat([data["HO_f"].to_frame(),data["HO_m"].to_frame().rename({'HO_m':'HO_f'},axis=1),
             data["JE_f"].to_frame().rename({'JE_f':'HO_f'},axis=1),
             data["JE_m"].to_frame().rename({'JE_m':'HO_f'},axis=1),
             data["BS_f"].to_frame().rename({'BS_f':'HO_f'},axis=1),
             data["BS_m"].to_frame().rename({'BS_m':'HO_f'},axis=1),
             data["AY_f"].to_frame().rename({'AY_f':'HO_f'},axis=1),
             data["AY_m"].to_frame().rename({'AY_m':'HO_f'},axis=1)])

thr=a['HO_f'].mean()+2*a['HO_f'].std()

res=[]
for i in range(data.shape[0]):
    if data.iloc[i,3]>= thr or data.iloc[i,4]>= thr or data.iloc[i,5]>= thr or data.iloc[i,6]>= thr or data.iloc[i,7]>= thr or data.iloc[i,8]>= thr or data.iloc[i,9]>= thr or data.iloc[i,10]>= thr:
        res.append(data.iloc[i,0])
        
df = pd.DataFrame(res).rename({0:'spot'},axis=1)
hotspots_bed=pd.read_csv('recombination_hotspots_ARS-UCD1.2.bed',sep="\t",names=(['chr','start','end','spot','digit','strand','info']))
m=hotspots_bed.merge(df)

m.to_csv('filtered_recombination_hotspots_ARS-UCD1.2.bed',index = None, header=True,sep="\t") 
