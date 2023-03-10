#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:19:10 2019

@author: beiki
"""

import collections
from collections import OrderedDict
import pandas as pd
import re
from multiprocessing import Pool
import time
import numpy as np
import pybedtools
import sys
#import itertools
## it is requred to rum 'module load bedtools2 before starting python ##
# on python 3, like on Ceres: change "iteritems" to "items"

input=sys.argv[1]# your gff file, like input="Cerebral_Cortex_fmlrc_proovread_cupcake_corrected_reads.collapsed.gff"
tissue=input.split("_fmlrc")[0]

def nested_dict():
    return collections.defaultdict(nested_dict)


def print_dict(dictionary, ident = '', braces=1):
    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print('%s%s%s%s' %(ident,braces*'[',key,braces*']'))
            print_dict(value, ident+'  ', braces+1)
        else:
            print(ident+'%s = %s' %(key, value))
        
""" slow
def get_df(input):
    output=pd.DataFrame([])
    for key,value in input.items(): # on python 2.7 use iteritems
        output=output.append(pd.DataFrame({'feature': key, 'coordinate': value}, index=[0]), ignore_index=True)
    return(output)

def get_tr_candidates3(df):
    candidates={}
    transcripts=df['feature'].tolist()
    comparisions=list(itertools.combinations(transcripts, 2))
    strand=re.split('\(|\)',df.iloc[0,:][0])[1]
    for pair in range(len(comparisions)):
        if compare_tr_begining_end(df.loc[df['feature']==comparisions[pair][0]].iloc[0,:][0],df.loc[df['feature']==comparisions[pair][1]].iloc[0,:][0],strand)=="TRUE" and get_RI_transcripts(comparisions[pair][0],comparisions[pair][1],strand).shape[0]>0:
            candidates[comparisions[pair][0]]=comparisions[pair][1]
    if len(candidates)>0:
        candidates=pd.DataFrame.from_dict(candidates,orient='index')[0].to_frame().reset_index().rename({'index':'transcriptA',0:'transcriptB'},axis=1)
    else:
        candidates=pd.DataFrame([])
    return(candidates)
"""    
    
def get_df(input):
    output=pd.DataFrame([])
    output=pd.DataFrame.from_dict(input,orient='index')[0].to_frame().reset_index().rename({'index':'feature',0:'coordinate'},axis=1)
    output=output[['coordinate','feature']]
    return(output)
   

def compare_tr_begining_end(trA_info,trB_info,strand):
    fuzzyTR_3=100
    fuzzyTR_5=1000
    a1=re.split('\-|\:|\(',trA_info)[1]
    a2=re.split('\-|\:|\(',trA_info)[2]
    b1=re.split('\-|\:|\(',trB_info)[1]
    b2=re.split('\-|\:|\(',trB_info)[2]
    if strand=="+" and abs(int(a1)-int(b1))<fuzzyTR_5 and abs(int(a2)-int(b2))<fuzzyTR_3:
        return("TRUE")
    elif strand=="-" and abs(int(a1)-int(b1))<fuzzyTR_3 and abs(int(a2)-int(b2))<fuzzyTR_5:
        return("TRUE")
    else:
        return("FALSE")
    

def get_tr_candidates(df):
    tr_numbers=len(df.index)
    candidates={}
    strand=re.split('\(|\)',df.iloc[0,:][0])[1]
    for i in range(tr_numbers-1):
        for j in range(i+1,tr_numbers):
            info=len(ex_coords[df.iloc[i,1]])+len(ex_coords[df.iloc[j,1]])## to prevent uspliced-unspliced pairs
            if info>2 and compare_tr_begining_end(df.iloc[i,:][0],df.iloc[j,:][0],strand)=="TRUE" and compare_tr(df.iloc[i,:],df.iloc[j,:],strand)=='TRUE':
                candidates[df.iloc[i,1]]=df.iloc[j,1]
    if len(candidates)>0:
        candidates=pd.DataFrame.from_dict(candidates,orient='index')[0].to_frame().reset_index().rename({'index':'transcriptA',0:'transcriptB'},axis=1)
    else:
        candidates=pd.DataFrame([])
    return(candidates)


def get_first_second_order(a,b,end,strand):
    result=[]
    if strand=="+":
        if end=='5' and (int(a[0])<int(b[0])):
            result.append(a)
            result.append(b)
            return(result)
        elif (end=='5') and (int(b[0])<int(a[0])):
            result.append(b)
            result.append(a)
            return(result)
        elif end=='5' and (int(a[0])==int(b[0])):
            result.append(a)
            result.append(b)
            return(result)
        elif end=='3' and (int(a[1])<int(b[1])):
            result.append(a)
            result.append(b)
            return(result)
        elif end=='3' and (int(b[1])<int(a[1])):
            result.append(b)
            result.append(a)
            return(result)
        else:
            result.append(a)
            result.append(b)
            return(result)
    if strand=="-":
        if end=='5' and (int(a[1])>int(b[1])):
            result.append(a)
            result.append(b)
            return(result)
        elif (end=='5') and (int(b[1])>int(a[1])):
            result.append(b)
            result.append(a)
            return(result)
        elif end=='5' and (int(a[1])==int(b[1])):
            result.append(a)
            result.append(b)
            return(result)
        elif end=='3' and (int(a[0])>int(b[0])):
            result.append(a)
            result.append(b)
            return(result)
        elif end=='3' and (int(b[0])>int(a[0])):
            result.append(b)
            result.append(a)
            return(result)
        else:
            result.append(a)
            result.append(b)
            return(result)
        return(result)
        
        
def get_first_last_exon(tr,strand):
    if strand=="+":
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=True)
        last_exon=tr_ex_sort.tail(1).iloc[0,0]
        last_exon_coords=re.split('\:|\-|\(|\)',last_exon)[1:3]
        first_exon=tr_ex_sort.head(1).iloc[0,0]
        first_exon_coords=re.split('\:|\-|\(|\)',first_exon)[1:3]
        return(first_exon_coords,last_exon_coords) 
    elif strand=="-":
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=False)
        last_exon=tr_ex_sort.tail(1).iloc[0,0]
        last_exon_coords=re.split('\:|\-|\(|\)',last_exon)[1:3]
        last_exon_coords=re.split('\:|\-|\(|\)',last_exon)[1:3]
        first_exon=tr_ex_sort.head(1).iloc[0,0]
        first_exon_coords=re.split('\:|\-|\(|\)',first_exon)[1:3]
        return(first_exon_coords,last_exon_coords)


def get_first2_last2_exons(tr,strand):
    if strand=="+":
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=True)
        last2_exon=tr_ex_sort.tail(2)
        first2_exon=tr_ex_sort.head(2)
        return(first2_exon,last2_exon)
    if strand=="-":
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=False)
        last2_exon=tr_ex_sort.tail(2)
        first2_exon=tr_ex_sort.head(2)
        return(first2_exon,last2_exon)
    return(first2_exon,last2_exon)


def format_df(df):
    feature=df['feature']
    output=pd.DataFrame([])
    a=df['coordinate'].str.split('\(|\)')
    strand=a.iloc[0][1]
    for key,value in a.iteritems():
        info=pd.Series(a[key][0]).str.split('\:|\-')[0]
        output=output.append(pd.DataFrame({'chr': info[0], 'start': info[1], 'end':info[2], 'name':feature[key], 'info':'90', 'strand':strand}, index=[0]), ignore_index=True)
        output = output[['chr', 'start', 'end','name', 'info', 'strand']]
    return(output)
    
        
def get_Ex_intersect_df(setA,setB):
    output=dict()
    for i in range(0,2):
        A=format_df(setA[i]).sort_values(by='start')
        B=format_df(setB[i]).sort_values(by='start')
        bedA=pybedtools.BedTool.from_dataframe(A)
        bedB=pybedtools.BedTool.from_dataframe(B)
        A_and_B = bedA.intersect(bedB,wao=True,s=True)
        output[i]=A_and_B
    return(output)
        

def parse_EX_intersect_df(inputs):
    result=dict()
    for i in range(0,2):
        intersect=inputs[i]
        df = pd.read_table(intersect.fn,names=['chromA', 'startA', 'stopA', 'exA', 'infoA', 'strandA','chromB', 'startB', 'stopB', 'exB', 'infoB', 'strandB','overlap'])
        a=df.groupby(['exA']).size().reset_index(name='count')
        b=df.groupby(['exB']).size().reset_index(name='count')
        if (len(a.loc[a['count']>1].index)>0) or (len(b.loc[b['count']>1].index)>0):
            result[i]='multicover'
        else:
            result[i]='none'
    return(result)
    
    

def compare_tr(a,b,strand):
    result=''
    fuzzyTR_3=100
    fuzzyTR_5=1000
    fuzzySplice=5
    if strand=='+':
        trA_last_exon=get_first_last_exon(a[1],strand)[1]
        trB_last_exon=get_first_last_exon(b[1],strand)[1]
        trA_first_exon=get_first_last_exon(a[1],strand)[0]
        trB_first_exon=get_first_last_exon(b[1],strand)[0]
        trA_last_exon_range=range(int(trA_last_exon[0]),int(trA_last_exon[1]))
        trB_last_exon_range=range(int(trB_last_exon[0]),int(trB_last_exon[1]))
        trA_first_exon_range=range(int(trA_first_exon[0]),int(trA_first_exon[1]))
        trB_first_exon_range=range(int(trB_first_exon[0]),int(trB_first_exon[1]))
        last_exon_overlap=list(set(trA_last_exon_range) & set(trB_last_exon_range))
        first_exon_overlap=list(set(trA_first_exon_range) & set(trB_first_exon_range))
        if (
                len(first_exon_overlap)>0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_first_exon[0])-int(trB_first_exon[0]))<=fuzzyTR_5
                           ) and (abs(int(trA_last_exon[1])-int(trB_last_exon[1]))<=fuzzyTR_3
                                 ) and (
                                         abs(int(trA_first_exon[1])-int(trB_first_exon[1]))<=fuzzySplice
                                         ) and (
                                                 abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzySplice):
            result='TRUE'
            return(result)
        elif (
                len(first_exon_overlap)>0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_first_exon[0])-int(trB_first_exon[0]))<=fuzzyTR_5
                           ) and (
                                   abs(int(trA_last_exon[1])-int(trB_last_exon[1]))<=fuzzyTR_3
                                   ) and (
                                           abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzySplice):
            trA_first2_last2_exons=get_first2_last2_exons(a[1],strand)
            trB_first2_last2_exons=get_first2_last2_exons(b[1],strand)
            Ex_intersects=get_Ex_intersect_df(trA_first2_last2_exons,trB_first2_last2_exons)
            decisions=parse_EX_intersect_df(Ex_intersects)
            if (decisions[0]=='multicover') or (decisions[1]=='multicover'):
                result='TRUE'
                return(result)
            else:
               result='FALSE'
               return(result)
            return(result)
        elif (
                len(first_exon_overlap)==0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_last_exon[1])-int(trB_last_exon[1]))<=fuzzyTR_3
                           ):
            order=get_first_second_order(trA_first_exon,trB_first_exon,'5',strand)
            a_fuzzy_end=int(order[0][1])+fuzzySplice
            b_fuzzy_end=int(order[1][1])-fuzzySplice
            if a_fuzzy_end == b_fuzzy_end:
                result='TRUE'
                return(result)
            else:
                result='FALSE'
                return(result)
            return(result)
        elif (
                len(first_exon_overlap)>0) and (len(last_exon_overlap)==0
                   ) and (
                           abs(int(trA_first_exon[0])-int(trB_first_exon[0]))<=fuzzyTR_5
                           ):
            order=get_first_second_order(trA_last_exon,trB_last_exon,'3',strand)
            b_fuzzy_end=int(order[0][1])+fuzzySplice
            a_fuzzy_end=int(order[1][1])-fuzzySplice
            if a_fuzzy_end == b_fuzzy_end:
                result='TRUE'
                return(result)
            else:
                result='FALSE'
                return(result)
            return(result)
        return(result)
    elif strand=='-':
        trA_last_exon=get_first_last_exon(a[1],strand)[1]
        trB_last_exon=get_first_last_exon(b[1],strand)[1]
        trA_first_exon=get_first_last_exon(a[1],strand)[0]
        trB_first_exon=get_first_last_exon(b[1],strand)[0]
        trA_last_exon_range=range(int(trA_last_exon[0]),int(trA_last_exon[1]))
        trB_last_exon_range=range(int(trB_last_exon[0]),int(trB_last_exon[1]))
        trA_first_exon_range=range(int(trA_first_exon[0]),int(trA_first_exon[1]))
        trB_first_exon_range=range(int(trB_first_exon[0]),int(trB_first_exon[1]))
        last_exon_overlap=list(set(trA_last_exon_range) & set(trB_last_exon_range))
        first_exon_overlap=list(set(trA_first_exon_range) & set(trB_first_exon_range))
        if (
                len(first_exon_overlap)>0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_first_exon[1])-int(trB_first_exon[1]))<=fuzzyTR_5
                           ) and (abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzyTR_3
                                 ) and (
                                         abs(int(trA_first_exon[0])-int(trB_first_exon[0]))<=fuzzySplice
                                         ) and (
                                                 abs(int(trA_last_exon[1])-int(trB_last_exon[1]))<=fuzzySplice):
            result='TRUE'
            return(result)
        elif (
                len(first_exon_overlap)==0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_last_exon[1])-int(trB_last_exon[1]))<=fuzzyTR_3
                           ) and (
                                   abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzySplice
                                   ) and (
                                           abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzyTR_3):
            order=get_first_second_order(trA_first_exon,trB_first_exon,'5',strand)
            a_fuzzy_end=int(order[0][1])+fuzzySplice
            b_fuzzy_end=int(order[1][1])-fuzzySplice
            if a_fuzzy_end == b_fuzzy_end:
                result='TRUE'
                return(result)
            else:
                result='FALSE'
                return(result)
            return(result)
        elif (
                len(first_exon_overlap)==0) and (len(last_exon_overlap)>0
                   ) and (
                           abs(int(trA_last_exon[0])-int(trB_last_exon[0]))<=fuzzyTR_3
                           ):
            order=get_first_second_order(trA_first_exon,trB_first_exon,'5',strand)
            a_fuzzy_end=int(order[1][1])+fuzzySplice
            b_fuzzy_end=int(order[0][1])-fuzzySplice
            if a_fuzzy_end == b_fuzzy_end:
                result='TRUE'
                return(result)
            else:
                result='FALSE'
                return(result)
            return(result)
        elif (
                len(first_exon_overlap)>0) and (len(last_exon_overlap)==0
                   ) and (
                           abs(int(trA_first_exon[1])-int(trB_first_exon[1]))<=fuzzyTR_5
                           ):
            order=get_first_second_order(trA_last_exon,trB_last_exon,'3',strand)
            b_fuzzy_end=int(order[1][1])+fuzzySplice
            a_fuzzy_end=int(order[0][1])-fuzzySplice
            if a_fuzzy_end == b_fuzzy_end:
                result='TRUE'
                return(result)
            else:
                result='FALSE'
                return(result)
            return(result)
        else:
         result='FALSE'
         return(result)
    return(result)

        
def compare_tr_ex_numbers(trA,trB):
    trA_ex=get_df(ex_coords.get(trA))
    trB_ex=get_df(ex_coords.get(trB))
    result=pd.DataFrame([])
    if len(trA_ex.index)<len(trB_ex.index):
        result=result.append(pd.DataFrame({'min': trA, 'max': trB}, index=[0]), ignore_index=True)
    else:
        result=result.append(pd.DataFrame({'min': trB, 'max': trA}, index=[0]), ignore_index=True)
    return(result)


def ilen(it):
    """Yield number of items in generator."""
    return len(list(it))


def get_TrEx_number(tr,strand):
    if strand=='+':
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=True)
        tr_ex_sort['exon_number']=np.arange(len(tr_ex_sort))
        tr_ex_sort.iloc[0,2]='first'
        tr_ex_sort.iloc[-1,2]='last'
        return(tr_ex_sort)
    elif strand=='-':
        tr_ex=get_df(ex_coords.get(tr))
        tr_ex_sort=tr_ex.sort_values(by='coordinate',ascending=False)
        tr_ex_sort['exon_number']=np.arange(len(tr_ex_sort))
        tr_ex_sort.iloc[0,2]='first'
        tr_ex_sort.iloc[-1,2]='last'
        return(tr_ex_sort)
    return(tr_ex_sort)
    
def get_AllEx_intersects_df(trA_ex,trB_ex):
    A=format_df(trA_ex).sort_values(by='start')
    B=format_df(trB_ex).sort_values(by='start')
    bedA=pybedtools.BedTool.from_dataframe(A)
    bedB=pybedtools.BedTool.from_dataframe(B)
    A_and_B = bedA.intersect(bedB,wo=True,s=True)
    output=A_and_B
    return(output)

 
def parse_AllEx_intersects_df(Ex_intersects,type_):
    if type_=='multi':
        multicover=nested_dict()
        df = pd.read_table(Ex_intersects.fn,names=['chromA', 'startA', 'stopA', 'exA', 'infoA', 'strandA','chromB', 'startB', 'stopB', 'exB', 'infoB', 'strandB','overlap'])
        a=df.groupby(['exA']).size().reset_index(name='count')
        b=df.groupby(['exB']).size().reset_index(name='count')
        if a.loc[a['count']>1].shape>0:
            for i in range(a.loc[a['count']>1].shape[0]):
                multicover['trA'][i]=a.loc[a['count']>1].reset_index().iloc[i,1]
        if b.loc[b['count']>1].shape>0:
            for j in range(b.loc[b['count']>1].shape[0]):
                multicover['trB'][j]=b.loc[b['count']>1].reset_index().iloc[j,1]
        return(multicover)  
    elif type_=='single':     
        singlecover=nested_dict()
        df = pd.read_table(Ex_intersects.fn,names=['chromA', 'startA', 'stopA', 'exA', 'infoA', 'strandA','chromB', 'startB', 'stopB', 'exB', 'infoB', 'strandB','overlap'])
        a=df.groupby(['exA']).size().reset_index(name='count')
        b=df.groupby(['exB']).size().reset_index(name='count')
        if a.loc[a['count']==1].shape>0:
            for i in range(a.loc[a['count']==1].shape[0]):
                singlecover['trA'][i]=a.loc[a['count']==1].reset_index().iloc[i,1]
        if b.loc[b['count']==1].shape>0:
            for j in range(b.loc[b['count']==1].shape[0]):
                singlecover['trB'][j]=b.loc[b['count']==1].reset_index().iloc[j,1]
        return(singlecover)

        
def get_multicoverEX(trA,trB,strand,df,info):
    fuzzySplice=5
    output=pd.DataFrame([])
    if strand=="+" and len(info.items())>0:
        for i in range(0,len(info.keys())):
            multicoverEX=0
            coveredEx=0
            tr=info.keys()[i]
            tr_ex_number=get_TrEx_number(eval(tr),strand)
            counts=df.groupby(['ex'+tr[2]]).size().reset_index(name='count') 
            if tr=='trA':
                tr2='trB' 
                tr2_ex_number=get_TrEx_number(eval('trB'),strand)
            else:
                tr2='trA' 
                tr2_ex_number=get_TrEx_number(eval('trA'),strand)
            for j in range(0,len(info[tr].values())): 
                ex=info[tr].values()[j] 
                if (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['stop'+tr2[2]].values[0])))<=fuzzySplice): 
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last'): 
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last'): ######## CAUTION: for both last!
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last' 
                            ):
                    if (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                    )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice 
                        ) and (
                            abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['stop'+tr2[2]].values[0])))<=fuzzySplice): 
                        multicoverEX=multicoverEX+1
                        coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
            if tr=='trA' and 'trB' not in info: 
                output=output.append(pd.DataFrame({'transcript': eval('trB'), 'multicoverEXs': 0, 'coveredEXs': 0}, index=[0]), ignore_index=True) 
            elif tr=='trB' and 'trA' not in info: 
                output=output.append(pd.DataFrame({'transcript': eval('trA'), 'multicoverEXs': 0, 'coveredEXs': 0}, index=[0]), ignore_index=True) 
            output=output.append(pd.DataFrame({'transcript': eval(tr), 'multicoverEXs': multicoverEX,'coveredEXs': coveredEx}, index=[0]), ignore_index=True) 
            output = output[['transcript', 'multicoverEXs', 'coveredEXs']]
        return(output)                      
    elif strand=="-" and len(info.items())>0: 
        for i in range(0,len(info.keys())):
            multicoverEX=0
            coveredEx=0
            tr=info.keys()[i]
            tr_ex_number=get_TrEx_number(eval(tr),strand)
            counts=df.groupby(['ex'+tr[2]]).size().reset_index(name='count') 
            if tr=='trA':
                tr2='trB'
                tr2_ex_number=get_TrEx_number(eval('trB'),strand)
            else:
                tr2='trA'
                tr2_ex_number=get_TrEx_number(eval('trA'),strand)
            for j in range(0,len(info[tr].values())): 
                ex=info[tr].values()[j] 
                if (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice): 
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last'): 
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last'): 
                    multicoverEX=multicoverEX+1
                    coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last' 
                            ):
                    if (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                    )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice 
                        ) and (
                            abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr2[2]].values[0])))<=fuzzySplice): 
                        multicoverEX=multicoverEX+1
                        coveredEx=coveredEx+counts.loc[counts['ex'+tr[2]]==ex].iloc[0,1] 
            if tr=='trA' and 'trB' not in info: 
                output=output.append(pd.DataFrame({'transcript': eval('trB'), 'multicoverEXs': 0, 'coveredEXs': 0}, index=[0]), ignore_index=True) 
            elif tr=='trB' and 'trA' not in info: 
                output=output.append(pd.DataFrame({'transcript': eval('trA'), 'multicoverEXs': 0, 'coveredEXs': 0}, index=[0]), ignore_index=True) 
            output=output.append(pd.DataFrame({'transcript': eval(tr), 'multicoverEXs': multicoverEX,'coveredEXs': coveredEx}, index=[0]), ignore_index=True)  
            output = output[['transcript', 'multicoverEXs', 'coveredEXs']]
        return(output)               
    return(output)  


def get_sharedEX(trA,trB,strand,df,info):
    fuzzySplice=5
    singlecoverEX=0
    if strand=="+" and len(info.items())>0:
        for i in range(0,len(info.keys())):
            singlecoverEX=0
            tr=info.keys()[i]
            tr_ex_number=get_TrEx_number(eval(tr),strand)
            if tr=='trA':
                tr2='trB' 
                tr2_ex_number=get_TrEx_number(eval('trB'),strand)
            else:
                tr2='trA' 
                tr2_ex_number=get_TrEx_number(eval('trA'),strand)
            for j in range(0,len(info[tr].values())): 
                ex=info[tr].values()[j] 
                if (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]=='first' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr2[2]].values[0])))<=fuzzySplice): 
                    singlecoverEX=singlecoverEX+1
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='last'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice): 
                    singlecoverEX=singlecoverEX+1                
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first' and tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='last'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]!='first' and tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[-1:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last'
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0]) - int(
                                    df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr2[2]].values[0])))<=fuzzySplice):
                                singlecoverEX=singlecoverEX+1
        return(singlecoverEX)
    elif strand=="-" and len(info.items())>0:      
        for i in range(0,len(info.keys())):
            singlecoverEX=0
            tr=info.keys()[i]
            tr_ex_number=get_TrEx_number(eval(tr),strand)
            if tr=='trA':
                tr2='trB'
                tr2_ex_number=get_TrEx_number(eval('trB'),strand)
            else:
                tr2='trA'
                tr2_ex_number=get_TrEx_number(eval('trA'),strand)
            for j in range(0,len(info[tr].values())): 
                ex=info[tr].values()[j] 
                if (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='first'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]=='first' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice): 
                    singlecoverEX=singlecoverEX+1
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]=='last'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]=='last' 
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr2[2]].values[0])))<=fuzzySplice): 
                    singlecoverEX=singlecoverEX+1
                elif (tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='first' and tr_ex_number.loc[tr_ex_number['feature']==ex].iloc[0,2]!='last'
                    ) and (
                            tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]!='first' and tr2_ex_number.loc[tr2_ex_number['feature']==df.loc[df['ex'+tr[2]]==ex].iloc[0:]['ex'+tr2[2]].values[0]].iloc[0,2]!='last'
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0]
                            )-int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['start'+tr2[2]].values[0])))<=fuzzySplice
                            ) and (abs((int(df.loc[df['ex'+tr[2]]==ex].iloc[0:]['stop'+tr[2]].values[0])))<=fuzzySplice):
                                singlecoverEX=singlecoverEX+1                                
        return(singlecoverEX)                          
    return(singlecoverEX)         
    
    
def get_RI_transcripts(trA,trB,strand):
    retained_intronsTR=pd.DataFrame([])
    trA_ex=get_df(ex_coords.get(trA))
    trB_ex=get_df(ex_coords.get(trB))
    Ex_intersects=get_AllEx_intersects_df(trA_ex,trB_ex)
    df = pd.read_table(Ex_intersects.fn,names=['chromA', 'startA', 'stopA', 'exA', 'infoA', 'strandA','chromB', 'startB', 'stopB', 'exB', 'infoB', 'strandB','overlap'])
    info=parse_AllEx_intersects_df(Ex_intersects,'multi')
    multicover_ex=get_multicoverEX(trA,trB,strand,df,info)
    info2=parse_AllEx_intersects_df(Ex_intersects,'single')
    shared_ex=get_sharedEX(trA,trB,strand,df,info2)
    if multicover_ex.shape[0]>0:
        if get_TrEx_number(trA,strand).shape[0]==multicover_ex.loc[multicover_ex['transcript']==trB].iloc[0,2]+multicover_ex.loc[multicover_ex['transcript']==trA].iloc[0,1]+shared_ex:
            retained_intronsTR=retained_intronsTR.append(pd.DataFrame({'RItranscript': trB, 'reference': trA}, index=[0]), ignore_index=True)
        elif get_TrEx_number(trB,strand).shape[0]==multicover_ex.loc[multicover_ex['transcript']==trA].iloc[0,2]+multicover_ex.loc[multicover_ex['transcript']==trB].iloc[0,1]+shared_ex:
            retained_intronsTR=retained_intronsTR.append(pd.DataFrame({'transcript': trA,'reference': trB}, index=[0]), ignore_index=True)
    return(retained_intronsTR)


def get_tr_strand(tr):
    tr_ex=get_df(ex_coords.get(tr))
    info=tr_ex.iloc[0,0]
    strand=re.split('\(|\)',info)[1]
    return(strand)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
def paralle_chunk_input(chunk):
    output=pd.DataFrame([])
    for geneset in range(len(chunk)):
        df=get_df(chunk[geneset][1])
        candidates=get_tr_candidates(df)
        if len(candidates.index)>0:
            output=output.append(candidates)
    return(output)

    
def paralle_chunk_input2(chunk):
    retained_intronsTR=pd.DataFrame([])
    for i in range(len(chunk)):
        trA=chunk[i][1][0]
        trB=chunk[i][1][1]
        strand=get_tr_strand(trA)
        RI_tr=get_RI_transcripts(trA,trB,strand)
        if RI_tr.shape[0]>0:
            retained_intronsTR=retained_intronsTR.append(pd.DataFrame({'RItranscript': RI_tr.iloc[0,0],'reference': RI_tr.iloc[0,1]}, index=[0]), ignore_index=True)
    return(retained_intronsTR)


ex_coords = collections.defaultdict(dict)
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2]=='exon':
        tid = raw[-1].split('; ')[1].split()[1][1:-2]
        exid=raw[2]+raw[3]+tid
        ex_coords[tid][exid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])       

gid="PB.1"
count=0
gene_tr={}
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2] == 'transcript':
        tid = raw[-1].split('; ')[1].split()[1][1:-2]
        gid_new = raw[-1].split('; ')[0].split()[1][1:-1]
        if gid==gid_new:
            count=count+1
            gene_tr[gid]=count
        if gid!=gid_new:
            count=0
            count=count+1
            gid=gid_new
            gene_tr[gid]=count
            
            
tr_coords = nested_dict()
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2] == 'transcript':
        tid = raw[-1].split('; ')[1].split()[1][1:-2]
        gid = raw[-1].split('; ')[0].split()[1][1:-1]
        if len(ex_coords[tid])>=1 and gene_tr[gid]>1: ## this export spliced/uspliced transcripts related to multit-ranscript genes  
            tr_coords[gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

len(tr_coords)

def get_optimum_number_of_genes_per_chunk(number_of_genes,chunk_size):
    """ this is to prevent getting chunks with only one gene
        (paralle_chunk_input function can't differntiate them
         with multi-gene sets"""
    chunk_size=chunk_size+1
    for i in reversed(range(1,chunk_size)):
        if number_of_genes%i==0:
            answer=i
            break
    return(answer)



    
""" to get gene with highest number ot ranscripts:
import operator
max(gene_tr.iteritems(), key=operator.itemgetter(1))[0]
gene_tr['PB.2146']   
"""
     
    
"""step1 Super efficient multi-threading.copleted in 10 minutes on 32 cpus!"""

p = Pool() 
optimum_number_of_genes=get_optimum_number_of_genes_per_chunk(len(tr_coords),3)   
items = list(tr_coords.items())
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets (i.e. uptimum_number_of_genes)in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, optimum_number_of_genes))
end = time.time()
print('total time (s)= ' + str(end-start)) 
candidates=pd.concat(processed_values)
candidates.to_csv(tissue + '_fmlrc_proovread_cupcake_corrected_reads_RI.candidates',index = None, header=True,sep="\t")


""" step2 multi-threading, completed in 11 minutes just on 7 cpus!"""
      
p = Pool() 
items = list(candidates.iterrows())
start = time.time()
# chunks=get_chunks(items, 30), chunk=next(chunks), ilen(chunks), len(chunk)
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 30))
end = time.time()
print('total time (s)= ' + str(end-start)) 
retained_intronsTR=pd.concat(processed_values)
RI_counts=retained_intronsTR.groupby(['RItranscript']).size().reset_index(name='count')                    
retained_intronsTR.to_csv(tissue + '_fmlrc_proovread_cupcake_corrected_reads.retained_intronsTR',index = None, header=True,sep="\t")

