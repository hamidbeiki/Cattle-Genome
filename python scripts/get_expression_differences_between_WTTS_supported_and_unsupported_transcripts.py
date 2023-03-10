#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:40:42 2021

@author: beiki
"""

import numpy as np
import pandas as pd
from scipy.stats import median_test
import functools
import operator
from scipy import stats

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])

transcript_tissue=np.load('transcript_tissue.npy',allow_pickle='TRUE').item()
tissue_transcript=np.load('tissue_transcript.npy',allow_pickle='TRUE').item()
wtts_tissues=pd.read_csv('WTTS/tissues',names=['tissue'])
wtts_tissues=list(wtts_tissues['tissue'])
trs=[]
for k,l in transcript_tissue.items():
    if len(set(l) & set(wtts_tissues))>1:
        trs.append(k)

valid_trs_df=pd.read_csv('WTTS/valid_transcripts_unstranded',names=['transcript_id'])
valid_trs=list(valid_trs_df['transcript_id'])

t=[]
for tr in valid_trs:
    l=transcript_tissue[tr]
    l2=list(set(l) & set(wtts_tissues))
    t.append(len(l2))

df=pd.DataFrame(t)
len(df.loc[df[0]==1])

rpkm_=pd.read_csv('quantification/RSEM_FAANG_transcripts_modified_fpkm',header=0,sep="\t",index_col=0)
s=rpkm_.loc[:, rpkm_.columns.isin(wtts_tissues)].reset_index()
rpkm=valid_trs_df.merge(s)
maxValues = rpkm.max(axis = 1).to_frame()
rpkm_=rpkm_.reset_index()

""" export expressed transcripts per tissue
cols=list(s.columns)
for i in range(1,s.shape[1]):
    t=cols[i]
    df=pd.DataFrame(s.iloc[:,i])
    df['transcript_id']=s['transcript_id']
    trs=pd.DataFrame(df.loc[df[cols[i]]>0]['transcript_id'])
    trs.to_csv('WTTS/' + t + '_transcripts',index = None, header=False,sep="\t")
"""
    
"""
rest=s[~ s.transcript_id.isin(rpkm.transcript_id)]
trs=[]
for i in wtts_tissues:
    if i in tissue_transcript.keys():
        tr=tissue_transcript[i]
        trs.append(tr)
trs_f=flatten_nested_list(trs)
trs = pd.DataFrame(list(dict.fromkeys(trs_f))).rename({0:'transcript_id'},axis=1)
m=trs.merge(rest)
rest=m
rest_maxValues = rest.max(axis = 1).to_frame()
"""
cols=list(rpkm.columns)
res={}
for i in range(1,rpkm.shape[1]):
    df=pd.read_csv("WTTS/" + cols[i] + "_valid_transcripts_unstranded", names=["transcript_id"])
    m=df.merge(s)
    df=m
    rest=pd.DataFrame(s[~ s.transcript_id.isin(df.transcript_id)][cols[i]])
    r=rest.loc[rest[cols[i]]>0]
    num=df.loc[df[cols[i]]>0].shape[0]
    a=df.loc[df[cols[i]]>0]
    perc=df.loc[df[cols[i]]>1].shape[0]/float(df.loc[df[cols[i]]>0].shape[0])
    y=r[cols[i]]
    x=a[cols[i]]
#    stat, p, med, tbl = median_test(x,y)
    stat, p=stats.ttest_ind(x, y, equal_var=True,alternative='greater')#https://stackoverflow.com/questions/25064506/what-scipy-statistical-test-do-i-use-to-compare-sample-means
    res[cols[i]]=[num,perc,x.mean(),y.mean(),p]

df=pd.DataFrame.from_dict(res, orient='index').reset_index().rename({"index":"tissue",0:'expr_trs',1:'%of >1rpkm',2:'supported_expr',3:'un-supported_expr',4:'p-value'},axis=1)
df['diff']=df['supported_expr']-df['un-supported_expr']
df['diff'].mean()


"""


#expression based
cols=list(rpkm.columns)
res={}
for i in range(1,rpkm.shape[1]):
    df=pd.DataFrame(rpkm.iloc[:,i])
    df['transcript_id']=rpkm['transcript_id']
    rest=pd.DataFrame(s[~ s.transcript_id.isin(df.transcript_id)][cols[i]])
    r=rest.loc[rest[cols[i]]>0]
    num=df.loc[df[cols[i]]>0].shape[0]
    a=df.loc[df[cols[i]]>0]
    perc=df.loc[df[cols[i]]>1].shape[0]/float(df.loc[df[cols[i]]>0].shape[0])
    y=r[cols[i]]
    x=a[cols[i]]
#    stat, p, med, tbl = median_test(x,y)
#    res[cols[i]]=[num,perc,x.median(),y.median(),p]
    stat, p=stats.ttest_ind(x, y, equal_var=True,alternative='greater')#https://stackoverflow.com/questions/25064506/what-scipy-statistical-test-do-i-use-to-compare-sample-means
    res[cols[i]]=[num,perc,x.mean(),y.mean(),p]    

#detection based
cols=list(rpkm.columns)
res={}
for i in range(1,rpkm.shape[1]):
    df=pd.DataFrame(rpkm.iloc[:,i])
    df['transcript_id']=rpkm['transcript_id']
    detected=pd.DataFrame(tissue_transcript[cols[i]]).rename({0:'transcript_id'},axis=1)
    df=detected.merge(df)
    rest=pd.DataFrame(detected[~ detected.transcript_id.isin(df.transcript_id)]['transcript_id'])
    rest2=rest.loc[:, rest.columns.isin(wtts_tissues)].reset_index()
    m=pd.DataFrame(rest.merge(rpkm_)[cols[i]])
    m2=m.loc[m[cols[i]]>0]
    num=df.loc[df[cols[i]]>0].shape[0]
    a=df.loc[df[cols[i]]>0]
    perc=df.loc[df[cols[i]]>1].shape[0]/float(df.loc[df[cols[i]]>0].shape[0])
    y=m2[cols[i]]
    x=a[cols[i]]
#    stat, p, med, tbl = median_test(x,y)
#   res[cols[i]]=[num,perc,x.median(),y.median(),p]
    stat, p=stats.ttest_ind(x, y, equal_var=True,alternative='greater')#https://stackoverflow.com/questions/25064506/what-scipy-statistical-test-do-i-use-to-compare-sample-means
    res[cols[i]]=[num,perc,x.mean(),y.mean(),p]  



    
cols=list(rpkm.columns)
res={}
for i in range(1,rpkm.shape[1]):
    df=pd.DataFrame(rpkm.iloc[:,i])
    df['transcript_id']=rpkm['transcript_id']
#    rest=pd.DataFrame(s[~ s.transcript_id.isin(df.transcript_id)][cols[i]])
#    r=rest.loc[rest[cols[i]]>0]
    detected=pd.DataFrame(tissue_transcript[cols[i]]).rename({0:'transcript_id'},axis=1)
    rest=pd.DataFrame(detected[~ detected.transcript_id.isin(df.transcript_id)]['transcript_id'])
    rest2=rest.loc[:, rest.columns.isin(wtts_tissues)].reset_index()
    m=pd.DataFrame(rest.merge(rpkm_)[cols[i]])
    m2=m.loc[m[cols[i]]>0]
    num=df.loc[df[cols[i]]>0].shape[0]
    a=df.loc[df[cols[i]]>0]
    perc=df.loc[df[cols[i]]>1].shape[0]/float(df.loc[df[cols[i]]>0].shape[0])
    y=m2[cols[i]]
    x=a[cols[i]]
    stat, p, med, tbl = median_test(x,y)
    res[cols[i]]=[num,perc]

pd.DataFrame.from_dict(res, orient='index').reset_index().rename({"index":"tissue",0:'expr_trs',1:'%of >1rpkm'},axis=1)

import matplotlib.pyplot as plt
plt.hist(maxValues[0], 25, histtype='step', color='black')
plt.savefig('hist.png')
plt.close()

detected=pd.DataFrame(tissue_transcript[cols[i]]).rename({0:'transcript_id'},axis=1)
rest=pd.DataFrame(detected[~ detected.transcript_id.isin(df.transcript_id)]['transcript_id'])
rest2=rest.loc[:, rest.columns.isin(wtts_tissues)].reset_index()
m=pd.DataFrame(rest.merge(rpkm_)[cols[i]])
m2=m.loc[m[cols[i]]>0]
y=m2[cols[i]]
x=a[cols[i]]
stat, p, med, tbl = median_test(x,y)

"""








