#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:26:57 2020

@author: beiki
"""

""" IMPORTANT NOTE:
    This program assume that there is no redundant
    transcript in the unput gff files or all transcript
    in each tissue are unique
    """

"""
    This progam is based on Python 3
    """

""" 
    NOTE: on nova use
    module load miniconda3/4.3.30-qdauveb
    source activate /home/beiki/.conda/envs/py-libs
"""

import re
import pandas as pd
import time
import sys
import pybedtools
from multiprocessing import Pool
import multiprocessing
from collections import OrderedDict
import csv
from collections import Iterable
import collections   
import networkx#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
import dask.dataframe as dd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import copy
import swifter
import functools
import operator

def get_ex_coords(input):
    ex_coords = []
    flag=0
    tid_old='PB.1.1'
    if input.find('gtf') !=-1:
        input_type='gtf'
    elif input.find('gff') !=-1:
        input_type='gff'
    if input_type=="gff":
        for line in open(input):
            raw=line.strip().split("\t")
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            if raw[2]=='exon':
                ex_coords.append([raw[0], raw[3], raw[4], tid,99,raw[6],flag+1])
            if tid==tid_old and tid!='PB.1.1':
                flag=flag+1
            else:
                flag=0
                tid_old=tid
    elif input_type=="gtf":
        for line in open(input):
            raw=line.strip().split("\t")
            if raw[2]=='exon':
                tid = raw[-1].split('; ')[0].split()[1][1:-1]
                ex_coords.append([raw[0], raw[3], raw[4], tid,99,raw[6],tid])    
    bed=pd.DataFrame(ex_coords)
    df=pd.DataFrame({'tid':['PB.1.1']})
    a=bed[~bed[3].isin(df['tid'])]
    b=bed[bed[3].isin(df['tid'])]
    if b.shape[0]>1:
        l=[]
        flag=0
        for i in range(b.shape[0]):
            info=list(b.iloc[i,0:])
            info[6]=flag+1
            l.append(info)
            flag=flag+1
        b=pd.DataFrame(l)
        bed=pd.concat([a,b]).reset_index().iloc[0:,1:]
    return(bed)

def get_tr_num_ex(input):
    tr_ex_info=[]
    if input.find('gtf') !=-1:
        input_type='gtf'
    elif input.find('gff') !=-1:
        input_type='gff'
    if input_type=="gff":
        for line in open(input):
            raw=line.strip().split("\t")
            if raw[2] == 'exon':
                info=raw[8].split("\"")
                tr_ex_info.append(info[3])        
    if input_type=="gtf":
        for line in open(input):
            raw=line.strip().split("\t")
            if raw[2] == 'exon':
                info=raw[8].split("\"")
                tr_ex_info.append(info[1])
    tr_ex_info_df=pd.DataFrame(tr_ex_info,columns =['tr_id'])
    tr_ex_info_df=pd.DataFrame(tr_ex_info_df['tr_id'].value_counts().values, index=tr_ex_info_df['tr_id'].value_counts().index, columns=['Count']).reset_index().rename({'index':'tr_id'},axis=1)## counts the numnber of exons per transcript
    return(tr_ex_info_df)
        
def get_ex_coords_dict(input,tissue):
    ex_coords = {}
    if input.find('gtf') !=-1:
        input_type='gtf'
    elif input.find('gff') !=-1:
        input_type='gff'
    if input_type=="gff":
        for line in open(input):
            raw=line.strip().split("\t")
            if raw[2]=='exon':
                tid = raw[-1].split('; ')[1].split()[1][1:-2]
                tid=tid.replace('PB',tissue)
                if ex_coords.get(tid):
                    ex_coords[tid].append([int(raw[3]), int(raw[4])])
                else:
                    ex_coords[tid]=[[int(raw[3]), int(raw[4])]]
    elif input_type=="gtf":
        for line in open(input):
            raw=line.strip().split("\t")
            if raw[2]=='exon':
                tid = raw[-1].split('; ')[0].split()[1][1:-1]
                tid=tid.replace('PB',tissue)
                if ex_coords.get(tid):
                    ex_coords[tid].append([int(raw[3]), int(raw[4])])
                else:
                    ex_coords[tid]=[[int(raw[3]), int(raw[4])]]
    return(ex_coords)

def check_exon_similarity(exonA,exonB):
    global fuzzySplice
    intervals=[]
    intervals.append(exonA)
    intervals.append(exonB)
    if abs((int(intervals[0][0]))-(int(intervals[1][0])))<=fuzzySplice and abs((int(intervals[0][1]))-(int(intervals[1][1])))<=fuzzySplice:
        return("same")
    else:
        return("different")

def check_tr_begining_end2(A,B,strand,tr_type):
    global fuzzySplice
    fuzzyTR_3=1000
    if tr_type=='spliced':
        if strand=="+" and abs(int(A[0][1])-int(B[0][1]))<=fuzzySplice and abs(int(A[-1][0])-int(B[-1][0]))<=fuzzySplice:
            return("TRUE")
        elif strand=="-" and abs(int(A[0][0])-int(B[0][0]))<=fuzzySplice and abs(int(A[-1][1])-int(B[-1][1]))<=fuzzySplice:
            return("TRUE")
        else:
            return("FALSE") 
    elif tr_type=='unspliced':
        if strand=="+" and abs(int(A[0][1])-int(B[0][1]))<=fuzzyTR_3:
            return("TRUE")
        elif strand=="-" and abs(int(A[0][0])-int(B[0][0]))<=fuzzyTR_3:
            return("TRUE")
        else:
            return("FALSE")           

def check_tr_begining_end(A,B,strand,tr_type):
    global fuzzySplice
#    fuzzyTR_3=1000
    if tr_type=='spliced':
        if abs(int(A[0][1])-int(B[0][1]))<=fuzzySplice and abs(int(A[-1][0])-int(B[-1][0]))<=fuzzySplice:
            return("TRUE")
        else:
            return("FALSE") 
    elif tr_type=='unspliced':
#        if abs(int(A[0][1])-int(B[0][1]))<=fuzzyTR_3:
        return("TRUE")
#        else:
#            return("FALSE") 

def compare_trs(list1,list2):
    flag=0
    for i in range(len(list1)):
        if check_exon_similarity(list1[i],list2[i])=="same":
            flag=flag+1
    return(flag)
    

def get_candidates(m,required_skipp_exons):
    out=[]
    for i in range(m.shape[0]):
        info=m.iloc[i,0:]
        if info[4]==info[1] and info[4]-info[5]<=required_skipp_exons: # this can be either skipped exn or less exon on query
            out.append(info)
    df=pd.DataFrame(out)
    return(df)
    
def get_masked_exons(name):
    global intersections_dict
    out1=[]# ref masked exons
    a=pd.DataFrame(intersections_dict[name])
    b=a.pivot_table(index=[13], aggfunc='size').reset_index().rename({0:'counts'},axis=1)
    my_list=list(b.loc[b['counts']>1][13])# query retained intron exons
    for i in range(len(my_list)):
        c=a.loc[a[13]==my_list[i]][6].to_frame()
        for j in range(c.shape[0]):
            out1.append(c.iloc[j,0]-1)# -1 becaue exons starts from 1
    out2=[x - 1 for x in my_list]# -1 becaue exons starts from 1
    out3=[]# query masked exons
    b=a.pivot_table(index=[6], aggfunc='size').reset_index().rename({0:'counts'},axis=1)
    my_list=list(b.loc[b['counts']>1][6])# ref retained intron exons
    for i in range(len(my_list)):
        c=a.loc[a[6]==my_list[i]][13].to_frame()
        for j in range(c.shape[0]):
            out3.append(c.iloc[j,0]-1)# -1 becaue exons starts from 1
    out4=[x - 1 for x in my_list]# -1 becaue exons starts from 1
    return(out1,out2,out4,out3)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
    

def paralle_chunk_input_1(chunk):
    bed_dicts={}
    for item in range(len(chunk)):
        file=chunk[item]
        tissue=re.split('mRNA-|_final|.collapsed', file)[1]
        bed=get_ex_coords(file)
        bed=bed.replace(to_replace='PB',value=tissue,regex=True)
        bed_dicts[tissue]=pybedtools.BedTool.from_dataframe(bed)
    return(bed_dicts)
    
def paralle_chunk_input_2(chunk):
    intersections=pd.DataFrame([])
    for item in chunk:
        i=item[0]
        j=item[1]
        tissueA=tissues[i]
        tissueB=tissues[j]
        intersection=bed_dicts[tissueA].intersect(bed_dicts[tissueB],wa=True,wb=True,s=True).to_dataframe()
        #on my linux machine, also add:
        intersection.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]
        intersections=pd.concat([intersections,intersection],sort=True) 
        #intersections=pd.concat([intersections,intersection])
    return(intersections)

def get_dask_parallel(chunk):
    global intersections
    for item in range(len(chunk)):
        df=intersections.iloc[chunk[item][0]:chunk[item][1],0:]
        d=dd.from_pandas(df,npartitions=-1)
        return(d)

def get_items(candidates):
    items=[]
    global isize
    for i in range(0,candidates.shape[0],isize):
        if i<candidates.shape[0]:
            items.append([i,i+isize])
        else:
            items.append([i,i+(isize+1)])
    return(items)
        
def aggregate(row):
    global intersections_dict
    if row.info in intersections_dict:
        intersections_dict[row.info].append(list(row))
    else:
        intersections_dict[row.info]=[list(row)]
    return(True)

def check_transcripts(l1,l2):                
    for i in range(len(l2)-1):
        if i>0 and i <len(l2):
            return(True)

def get_lists_intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]
        if a_right < b_right:
            i += 1
        else:
            j += 1
        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)
    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges
            
def paralle_chunk_input1(chunk):
    """ exact same structures """
    global group1
    global ex_coords
    global tr_strand_chr
    same_trs_spliced={}
    same_trs_unspliced={}
    for item in range(len(chunk)):
        candidates=group1.iloc[chunk[item][0]:chunk[item][1],0:]
        for i in range(candidates.shape[0]):
            trs=candidates.iloc[i,0].split('---')
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            A=[l1[0],l1[-1]]# ref tr beginig and end exons
            B=[l2[0],l2[-1]]# query tr beginig and end exons
            if len(l1)>1 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
                if len(l1)>2:
                    list1=l1[1:-1]
                    list2=l2[1:-1]
                    same_ex=compare_trs(list1,list2)# beginig and end exons aere already simmilar
                    if (len(l1)-2)==same_ex:
                        if trs[0] not in same_trs_spliced:
                            same_trs_spliced[trs[0]]=[trs[1]]
                        else:
                            same_trs_spliced[trs[0]].append(trs[1])
                elif len(l1)==2:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
            elif len(l1)==1 and len(l2)==1 and check_tr_begining_end(l1,l2,strand,'unspliced')=='TRUE':
                if trs[0] not in same_trs_unspliced:
                    same_trs_unspliced[trs[0]]=[trs[1]]
                else:
                    same_trs_unspliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_unspliced)

def paralle_chunk_input2(chunk):
    """ same 3' end exon but extra exon(s) on the 5' end of RNA-seq transcript
        same 3' end but with reained intron exon at RNA-seq or Iso-seq transcript 
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group2
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}# same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group2.iloc[chunk[item][0]:chunk[item][1],0:]    
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            strand=tr_strand_chr[trs[0]][0][0]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            A=[l1[first_r],l1[last_r]]# ref tr beginig and end exons
            B=[l2[first_q],l2[last_q]]# query tr beginig and end exons
            if len(l1)>2 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
                if len(l1)>2:
                    list1=l1[first_r+1:last_r]
                    list2=l2[first_q+1:last_q]
                    if len(list1)==len(list2):
                        same_ex=compare_trs(list1,list2)
                        if (len(list1))==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])
                    elif len(list2)<len(list1):
                        masked_info=get_masked_exons(name)
                        if len(masked_info[0])>0:
                            flag=[x - 1 for x in masked_info[0]]# adjust exon number
                            list1_2=[ list1[j] for j in range(len(list1)) if j not in flag ]
                            flag=[x - 1 for x in masked_info[1]]# adjust exon number
                            list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2) == len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained:
                                    same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                                else:
                                    same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})
                    elif len(list2)>len(list1):#!!!!!!!!!!!!!!!!!!!!!
                        masked_info=get_masked_exons(name)
                        if len(masked_info[2])>0:
                            flag=[x - 1 for x in masked_info[2]]# adjust exon number
                            list1_2=[ list1[j] for j in range(len(list1)) if j not in flag ]
                            flag=[x - 1 for x in masked_info[3]]# adjust exon number
                            list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2) == len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained2:
                                    same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}#masked_info[2] show the retained intron exon in refrence trancript
                                else:
                                    same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]}) 
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)
    
def paralle_chunk_input3(chunk):
    """ same 5' end exon but extra exon(s) on the 3' end of RNA-seq transcript. on the reverse strand this would be same 3' but extra exon on 5' end of the rna-seq transcript
        same  as above but with retained exon on RNA-seq or Iso-seq transcript 
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group3
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}# same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group3.iloc[chunk[item][0]:chunk[item][1],0:]            
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            A=[l1[first_r],[5,10]]# ref tr beginig and end exons 
            B=[l2[first_q],[5,10]]# query tr beginig and end exons 
            if len(l1)>2 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
                if len(l1)>2:
                    list1=l1[first_r+1:(last_r)]
                    list2=l2[first_q+1:(last_q)]
                    if len(list1)==len(list2):
                        same_ex=compare_trs(list1,list2)# beginig and end exons aere already simmilar
                        if (len(list1))==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])
                    elif len(list1)<len(list2):
                        masked_info=get_masked_exons(name)
                        if len(masked_info[0])>0:
                            masked_info[0].sort()
                            list1_2=[ l1[j] for j in range(len(l1)) if j>0 and j<len(l1)-1 and j not in masked_info[0] ]
                            list2_2=[ l2[k] for k in range(len(l2)) if k>0 and k<len(l2)-1 and k not in masked_info[1] ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2)==len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained:
                                    same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                                else:
                                    same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})
                    elif len(list2)<len(list1):#!!!!!!!!!!!!!
                        masked_info=get_masked_exons(name)
                        if len(masked_info[2])>0:
                            masked_info[3].sort()
                            list1_2=[ l1[j] for j in range(len(l1)) if j>0 and j<len(l1)-1 and j not in masked_info[2] ]
                            list2_2=[ l2[k] for k in range(len(l2)) if k>0 and k<len(l2)-1 and k not in masked_info[3] ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2)==len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained2:
                                    same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}#masked_info[2] show the retained intron exon in refrence trancript
                                else:
                                    same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                                
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)

def paralle_chunk_input4(chunk):
    """same splice junctions at all ref exons but additional exons at both 5' end and 3' end of RNA-seq transcript
       same as above but with retained-intron exon at rna-seq or Iso-seq transcript
       transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group4
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group4.iloc[chunk[item][0]:chunk[item][1],0:]  
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            if len(l1)>2 and len(l2)>1:
                list1=l1[first_r:last_r+1]
                list2=l2[first_q:last_q+1]
                if len(list1)==len(list2) and len(list2)>1:
                    A=[[0,list1[0][1]],[list1[-1][0],0]]
                    B=[[0,list2[0][1]],[list2[-1][0],0]]
                    if check_tr_begining_end(A,B,strand,'spliced')=='TRUE':
                        same_ex=compare_trs(list1[1:-1],list2[1:-1])# beginig and end exons aere already simmilar
                        if (len(list1)-2)==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])
                elif len(list2)<len(list1) and len(list2)>1:
                    masked_info=get_masked_exons(name)
                    if len(masked_info[0])>0:
                        list1_2=[ list1[j] for j in range(len(list1)) if j not in masked_info[0] ]
                        flag=[x - first_q for x in masked_info[1]]# adjust exon number
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                        if len(list1_2)>0 and len(list2_2)>0:
                            A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                            B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                            if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!
                                if trs[0] not in same_trs_spliced_with_retained:
                                    same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                                else:
                                    same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})
                        else:
                            if trs[0] not in same_trs_spliced_with_retained:
                                same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                            else:
                                same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})                                
                elif len(list1)<len(list2) and len(list1)>1 and first_q==0:#!!!!!!
                    masked_info=get_masked_exons(name)
                    if len(masked_info[2])>0:
                        list1_2=[ list1[j] for j in range(len(list1)) if j not in masked_info[2] ]
                        flag=[x - first_q for x in masked_info[3]]# adjust exon number
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                        if len(list1_2)>0 and len(list2_2)>0:
                            A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                            B=[[0,list2_2[0][1]],[list2_2[-1][0],5]] 
                            if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!!
                                if trs[0] not in same_trs_spliced_with_retained2:
                                    same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}#masked_info[2] show the retained intron exon in refrence trancript
                                else:
                                    same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]}) 
                        else:
                            if trs[0] not in same_trs_spliced_with_retained2:
                                same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}#masked_info[2] show the retained intron exon in refrence trancript
                            else:
                                same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                             
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)

def paralle_chunk_input5(chunk):
    """ same 3' end exon but extra exon(s) on the 5' end of Iso-seq transcript or reverse in + strand
        same as above but with retained-intron exon at rna-seq or Iso-seq transcript
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group5
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group5.iloc[chunk[item][0]:chunk[item][1],0:] 
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            if len(l1)>2 and len(l2)>1:
                list1=l1[first_r:last_r+1]
                list2=l2[first_q:last_q+1]
                if len(list1)==len(list2) and len(list2)>1:
                    A=[[0,list1[0][1]],[list1[-1][0],0]]
                    B=[[0,list2[0][1]],[list2[-1][0],0]]
                    if check_tr_begining_end(A,B,strand,'spliced')=='TRUE':
                        same_ex=compare_trs(list1[1:-1],list2[1:-1])# beginig and end exons aere already simmilar
                        if (len(list1)-2)==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])                            
                elif len(list2)<len(list1) and len(list2)>1: #!!!!!!!!
                    masked_info=get_masked_exons(name)
                    if len(masked_info[0])>0 and (first_r==0 or first_q==0):#!!!!!!!
                        list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[0] ]
                        flag=[x - first_q for x in masked_info[1]]# adjust exon number
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                        if len(list1_2)>0 and len(list2_2)>0:
                            A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                            B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                            if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!!!!!!!!
                                if trs[0] not in same_trs_spliced_with_retained:
                                    same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                                else:
                                    same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})  
                        else:
                            if trs[0] not in same_trs_spliced_with_retained:
                                same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                            else:
                                same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})                             
                elif len(list1)<len(list2) and len(list1)>1:#!!!!!!
                    masked_info=get_masked_exons(name)
                    if len(masked_info[2])>0:#!!!!!
                        list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[2] ]#!!!!
                        flag=[x - first_q for x in masked_info[3]]# adjust exon number !!!!
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                        if len(list1_2)>0 and len(list2_2)>0:
                            A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                            B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                            if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):
                                if trs[0] not in same_trs_spliced_with_retained2:
                                    same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}
                                else:
                                    same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})  
                        else:
                            if trs[0] not in same_trs_spliced_with_retained2:
                                same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}
                            else:
                                same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                            
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)
        
def paralle_chunk_input6(chunk):
    """same 5' end exon but extra exon(s) on the 3' end of RNA-seq transcript. on the reverse strand this would be same 3' but extra exon on 5' end of the rna-seq transcript
        same  as above but with retained exon on RNA-seq or Iso-seq transcript 
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group6
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group6.iloc[chunk[item][0]:chunk[item][1],0:] 
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            A=[l1[first_r],l1[last_r]]# ref tr beginig and end exons 
            B=[l2[first_q],l2[last_q]]# query tr beginig and end exons 
            if len(l1)>1 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
                if len(l1)>2:
                    list1=l1[(first_r+1):(last_r)]
                    list2=l2[(first_q+1):(last_q)]
                    if len(list1)==len(list2):
                        same_ex=compare_trs(list1,list2)# beginig and end exons aere already simmilar
                        if (len(list1))==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])
                    elif len(list1)<len(list2):
                        masked_info=get_masked_exons(name)
                        if len(masked_info[0])>0 and (first_r==0 or first_q==0):#!!!!!!
                            masked_info[0].sort()
                            list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[0] ]
                            list2_2=[ list2[k] for k in range(len(list2)) if (k+first_q) not in masked_info[1] ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2)==len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained:
                                    same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                                else:
                                    same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]})  
                    elif len(list2)<len(list1):
                        masked_info=get_masked_exons(name)
                        if len(masked_info[2])>0:
                            masked_info[3].sort()
                            list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[2] ]
                            list2_2=[ list2[k] for k in range(len(list2)) if (k+first_q) not in masked_info[3] ]
                            if len(list1_2)==len(list2_2) and compare_trs(list1_2,list2_2)==len(list1_2):
                                if trs[0] not in same_trs_spliced_with_retained2:
                                    same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}
                                else:
                                    same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                                                              
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)
    
def paralle_chunk_input7(chunk):
    """same splice junctions at all ref exons but additional exons at both 5' end and 3' end of RNA-seq transcript
        same as above but with retained-intron exon at rna-seq or Iso-seq transcript
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group7
    global coords
    global ex_coords
    global tr_strand_chr
    global fuzzySplice
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group7.iloc[chunk[item][0]:chunk[item][1],0:]
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=candidates.iloc[i,0].split('---')
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=tr_strand_chr[trs[0]][0][0]
            l1=ex_coords[trs[0]]
            l2=ex_coords[trs[1]]
            if len(l1)>2 and len(l2)>1:
                list1=l1[first_r:last_r+1]
                list2=l2[first_q:last_q+1]
                if len(list1)==len(list2) and len(list2)>1:
                    A=[[0,list1[0][1]],[list1[-1][0],0]]
                    B=[[0,list2[0][1]],[list2[-1][0],0]]
                    if check_tr_begining_end(A,B,strand,'spliced')=='TRUE':
                        same_ex=compare_trs(list1[1:-1],list2[1:-1])# beginig and end exons aere already simmilar
                        if (len(list1)-2)==same_ex:
                            if trs[0] not in same_trs_spliced:
                                same_trs_spliced[trs[0]]=[trs[1]]
                            else:
                                same_trs_spliced[trs[0]].append(trs[1])                            
                elif len(list2)<len(list1) and len(list2)>1 and (first_r==0 or first_q==0):#!!!!
                    masked_info=get_masked_exons(name)
                    if len(masked_info[0])>0:
                        list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[0] ]
                        flag=[x - first_q for x in masked_info[1]]# adjust exon number
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ]
                        if len(list1_2)>0 and len(list2_2)==list1_2 and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):
                            if trs[0] not in same_trs_spliced_with_retained:
                                same_trs_spliced_with_retained[trs[0]]={trs[1]:masked_info[1]}
                            else:
                                same_trs_spliced_with_retained[trs[0]].update({trs[1]:masked_info[1]}) 
                elif len(list1)<len(list2) and len(list1)>1 and (first_r==0 or first_q==0): #!!!!!
                    masked_info=get_masked_exons(name)
                    if len(masked_info[2])>0:
                        list1_2=[ list1[j] for j in range(len(list1)) if (j+first_r) not in masked_info[2] ]
                        flag=[x - first_q for x in masked_info[3]]# adjust exon number
                        list2_2=[ list2[k] for k in range(len(list2)) if k not in flag ] 
                        if len(list1_2)>0 and len(list2_2)==list1_2 and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):
                            if trs[0] not in same_trs_spliced_with_retained2:
                                same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}
                            else:
                                same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                        
            if len(l1)==2 and len(l2)>1:
                info=get_lists_intersections(l1,l2)
                if len(info)==2 and abs(l1[0][1]-info[0][1])<=fuzzySplice and abs(l1[1][0]-info[1][0])<=fuzzySplice:
                    if trs[0] not in same_trs_spliced:
                        same_trs_spliced[trs[0]]=[trs[1]]
                    else:
                        same_trs_spliced[trs[0]].append(trs[1])
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)
    

    
#intervals = [[3,9],[2,6],[8,10],[15,18]]
#intervals = [[1,2],[4,5],[7,8],[3,4],[7,10],[1,2],[4,5],[6,10]]

def merge_intervals(arr):
    # Sorting based on the increasing order  
    # of the start intervals
    arr.sort(key = lambda x: x[0]) 
    # array to hold the merged intervals 
    m = []
    s = -10000
    max = -100000
    for i in range(len(arr)):
        a = arr[i]
        if a[0] > max:
            if i != 0:
                m.append([s,max])
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]
    #'max' value gives the last point of
    # that particular interval
    # 's' gives the starting point of that interval
    # 'm' array contains the list of all merged intervals
    if max != -100000 and [s, max] not in m:
        m.append([s, max])
    return(m)

def check_exon_position(exonA,exonB,ex_type):
    if ex_type=='first':
        if exonB[0]>exonA[1]:
            return("included")
        else:
            return("notincluded")
    elif ex_type=='last':
        if exonB[1]<=exonA[1]:
            return("included")
        else:
            return("notincluded")   

def pairs(lst):
    i = iter(lst)
    first = prev = item = i.__next__()# use i.next() in python 2.7
    for item in i:
        yield prev, item
        prev = item
    yield item, first

#lists = [[1,2,3],[3,5,6],[8,9,10],[11,12,13]]
#lists=[['A','B'],['B','C','D'],['D','E','F'],['O','G']]
    
def get_key(tr,my_dict):
    for key, value in my_dict.items():
        for key2 in value.keys():
            if key2==tr:
                return(key)
                
def adjust_ex(exon_list,num_list):
    out=[]
    for i in range(len(exon_list)):
        if i not in num_list:
            out.append(exon_list[i])
    return(out)
    
 
def get_tr_splicing_info(gff_input,tissue):
    tr_ex_info=[]
    chr_strand_info={}
    for line in open(gff_input):
        raw=line.strip().split("\t")
        if raw[2] == 'exon':
            info=raw[8].split("\"")
            tr=info[3].replace('PB',tissue)
            tr_ex_info.append(tr)
            chr_strand_info[tr]=[raw[0],raw[6]]
    tr_ex_info_df=pd.DataFrame(tr_ex_info,columns =['tr_id'])
    tr_ex_info_df=pd.DataFrame(tr_ex_info_df['tr_id'].value_counts().values, index=tr_ex_info_df['tr_id'].value_counts().index, columns=['Count']).reset_index().rename({'index':'tr_id'},axis=1)## counts the numnber of exons per transcript
    spliced=tr_ex_info_df.loc[tr_ex_info_df['Count']>1]
    unspliced=tr_ex_info_df.loc[tr_ex_info_df['Count']==1]
    return(spliced,unspliced,chr_strand_info)

   
def filter_intersection(chunk):
    global intersection
    df=pd.DataFrame([])
    for item in range(len(chunk)):
        sections=intersection.iloc[chunk[item][0]:chunk[item][1],0:]
        sections=sections[sections[[6,10]].nunique(axis=1) != 1]
        df=pd.concat([df,sections])
    return(df)

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])

"""def flatten_nested_list(lis):
    return functools.reduce(operator.concat, lis)
    """
""" in python 2.7 
def flatten_nested_list(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, basestring):
             for x in flatten_nested_list(item):
                 yield x
         else:        
             yield item
"""
             
def get_key2(tr,my_dict):
    for key, value in my_dict.items():
        for key2 in value:
            if key2==tr:
                return(key)

def format_mat(mat_dict):
    global tissues
    mat=pd.DataFrame.from_dict(mat_dict,orient='index')
    mat.columns=tissues
    mat=mat.reset_index()
    return(mat)

def get_comon(l1,l2):
    return(list(set(l1) & set(l2)))    

def recursive_check(mylist):
    mergable=[]
    separate_merge=[]
    for i in mylist:
        if i in same_trs_spliced:
            info=get_comon(mylist,same_trs_spliced[i])
            if len(info)==1:
                mergable.append(i)
                mergable.append(info[0])
            elif len(info)>1:
                info2=recursive_check(info)[0]
                if len(info2)>0:
                    mergable.append(i)
                    for j in info2:
                        mergable.append(j)
                else:
                    for j in info:
                        separate_merge.append([i,j])                        
    return(mergable,separate_merge) 

def get_items_fromLists(lst):
    items=[]
    global isize
    for i in range(0,len(lst),isize):
        if i<len(lst):
            items.append([i,i+isize])
        else:
            items.append([i,i+(isize+1)])
    return(items)

def get_value(i,j,d1):
    n1= i + '--' + j 
    n2= j + '--' + i 
    if n1 in d1:
        return(d1[n1])
    elif n2 in d1:
        return(d1[n2])
    else:
        return('none')
        
def separate_lists(l):
    global same_trs_spliced
    l1=[]
    l2=[]
    for i in l:
        if i in same_trs_spliced:
            l1.append(i)
        else:
            l2.append(i)
    return(l1,l2)
    
def get_adjust_dict(l1):
    global same_trs_spliced
    subdict=dict((k, same_trs_spliced[k]) for k in l1)
    adj={}
    for tr in l1:
        v1=subdict[tr]
        nv=[]
        for k,v in subdict.items():#python3
            if tr in v and k not in nv:
                nv.append(k)
            for i in v1:
                if i not in nv:
                    nv.append(i)
        adj[tr]=nv
    return(adj)

def get_comonality(l1,adjd):
    d1={}# a dictionary containing comonality between keys in l1
    for i in range(len(l1)):
        for j in range(i+1,len(l1)):
            a=l1[i]
            b=l1[j]
            va=adjd[a]
            vb=adjd[b]
            c=get_comon(va,vb)
            n1=a + '--' + b
            n2=b + '--' + a
            d1[n1]=c
            d1[n2]=c
    return(d1)

def get_comonality2(l1,l2,adjd):
    d2={}#a dictionary containing dicometized similarity between transcripts in l 
    l12=l1[:]
    universe=list(adjd.values())
    patches=[]
    for i in l2:
        l12.append(i)
    flag=len(l1)
    o2=[]
    count=0
    for i in l1:
        o=[]
        v=adjd[i]
        for j in l12:
            if j in v:
                o.append(1)
            else:
                o.append(0)
        o[count]=1
        count=count+1
        ol=o[flag:]
        o2.append(ol)
        d2[i]=o
    df=pd.DataFrame(o2)
    counter=0
    for i in range(len(l2)):
        s=[]
        t1=l2[i]
        a=list(df.iloc[0:,i])
        for j in range(len(l2)):
            t2=l2[j]
            if check_nested([t1,t2],universe):
                s.append(1)
                patches.append([t1,t2])
            else:
                s.append(0)
        s[counter]=1
        counter=counter+1
        for k in s:
            a.append(k)
        n=l2[i]
        d2[n]=a
    patches=list(set(tuple(sorted(sub)) for sub in patches))# remove exact duplicates. takes [['hamid', 'mina'], ['mina', 'hamid']] and return [('hamid', 'mina')]
    return(d2,l12,patches)
        
def get_indexes(v2,l1):
    idx=[]
    o=[ idx.append(x) for x in range(len(l1)) if l1[x] in v2 ]
    return(idx)


def get_elements(ls,idx,l12):
    out=[]
    for i in range(len(ls)):
        if i in idx and ls[i]>0:
            out.append(l12[i])
    return(out)
    
    
def check_nested(el,m2):
    if len(m2)>0:
        for i in m2:
            if len(get_comon(el,i))==len(el):
                return(True)
    else:
        return(False)

def nested_sum(L):
    total = 0  # don't use `sum` as a variable name
    for i in L:
        if isinstance(i, list):  # checks if `i` is a list
            total += nested_sum(i)
        else:
            total += i
    return total

def get_mylist(v2,l1,l12,d2):
    m2=[]
    info=[]
    idx=get_indexes(v2,l12)
    for s in v2:
        ls=d2[s]
        el=get_elements(ls,idx,l12)
        if not check_nested(el,m2):
            mylist=[]
            m2.append(el)
            for element in el:
                ls2=d2[element]
                idx2=get_indexes(el,l12)
                l=[ls2[i] for i in idx2]
                mylist.append(l)
            check=nested_sum(mylist)
            info.append([check])
    return(info,m2)
    
def get_redundants(m):
    redundants=[]
    for i in m:
        if check_nested2(i,m):
            redundants.append(i)
    return(redundants)

def remove_redundants(m,r):
    res=[]
    o=[res.append(x) for x in m if x not in r]
    return(res)

def get_ordered_trs(nestedlist):
    f=list(flatten_nested_list(nestedlist))
    f=list(dict.fromkeys(f))
    d={}
    for i in f:
        ex= ex_coords[i]
        s=ex[0][0]
        e=ex[-1][1]
        l=e-s
        d[i]=l
    d=OrderedDict(sorted(d.items(), key=lambda x: x[1],reverse=True))
    k=d.keys()
    return(k)
                    
def get_final_list(nestedlist,ordered_trs):
    final_list=[]
    p=[]
    for i in ordered_trs:
        out=[]
        for j in nestedlist:
            if i in j:
                if j in p:
                    pass
                else:
                    out.append(j)
                    p.append(j)
        if len(out)>0:
            final_list.append(out)
    return(final_list)

def get_m_set(final_list):
    m=[]
    for i in final_list:
        o=list(flatten_nested_list(i))
        o=list(dict.fromkeys(o))
        m.append(o)
    return(m)

def clean_res(nestedlist):
    ordered_trs=get_ordered_trs(nestedlist)
    final_list=get_final_list(nestedlist,ordered_trs)
    res=get_m_set(final_list)
    return(res)

def get_patch(trs,patches):
    out=[]
    res=[]
    for i in patches:
        for tr in trs:
            if tr in i:
                l=list(i)
                l.remove(tr)
                out.append(l)
    out=list(flatten_nested_list(out))
    out=list(dict.fromkeys(out))
    for i in out:
        if i not in trs:
            res.append(i)
    return(res)           

def get_mergable_same_trs_spliced_parallel(chunk):
    global same_trs_spliced
    global merges
    m=[]#list mergable transcripts
    for item in range(len(chunk)):
        chunk_k=merges[chunk[item][0]:chunk[item][1]]
        for lis in chunk_k:
            l=list(lis)
            info=separate_lists(l)
            l1=info[0]# list of transcripts that existed in same_trs_spliced as key
            l2=info[1]# list of transcripts that NOT existed in same_trs_spliced as key
            adjd=get_adjust_dict(l1)
            d1=get_comonality(l1,adjd)# a dictionary containing comonality between keys in l1
            info=get_comonality2(l1,l2,adjd)
            d2=info[0]
            l12=info[1]
            patches=info[2]
            for i,v in adjd.items():#python3!
                for j in v:
                    v2=get_value(i,j,d1)
                    if v2!='none' and len(v2)>0:
                        if len(v2)==1:
                            if i not in v2:
                                v2.append(i)
                            if j not in v2:
                                v2.append(j)
                            if not check_nested(v2,m):
                                m.append(v2)
                        elif len(v2)>1:
                            info=get_mylist(v2,l1,l12,d2)
                            checks=info[0]
                            els=info[1]
                            c=0
                            for el in els:
                                check=checks[c][0]
                                c=c+1
                                if check == (len(el)**2):
                                    if i not in el:
                                        el.append(i)
                                    if j not in el:
                                        el.append(j)
                                    if not check_nested(el,m):
                                        m.append(el)
                    else:
                        if not check_nested([i,j],m):
                            m.append([i,j])                                   
    m=list(set(tuple(sorted(sub)) for sub in m))# remove exact duplicates. takes [['hamid', 'mina'], ['mina', 'hamid']] and return [('hamid', 'mina')]
    if len(m)>1:
        r=get_redundants(m)#redundants
        nestedlist=remove_redundants(m,r)
        res1=clean_res(nestedlist)
        r=get_redundants(res1)
        res=remove_redundants(res1,r)
    elif len(m)==1:
        res=m
    return(res)
"""
def check_nested2(i,mergable_same_trs_spliced):
    for j in mergable_same_trs_spliced:
        if j!=i and (set(i).issubset(set(j))):
            return(True)
"""
def check_nested2(i,mergable_same_trs_spliced):
    for j in mergable_same_trs_spliced:
        if j!=i and (set(i) <=set(j)):
            return(True)
   
def get_redundancy_parallel(chunk):
    global mergable_same_trs_spliced
    redundants=[]
    for item in range(len(chunk)):
        chunk_k=mergable_same_trs_spliced[chunk[item][0]:chunk[item][1]]  
        for i in chunk_k:
            if check_nested2(i,mergable_same_trs_spliced):
                redundants.append(i)
    return(redundants)


def get_tissue_parallel(chunk):
    global remain_spliced
    d={}
    for item in range(len(chunk)):
        df=remain_spliced.iloc[chunk[item][0]:chunk[item][1],0:]
        for i in range(df.shape[0]):
            tr=df.iloc[i,0]
            t=tr.split('.')[0]
            d[tr]=[t]
    return(d)    
#################
#       RUN     #
#################    

#####################################
#  phase1: get simmilar transcripts #
##################################### 
    
input_files=sys.argv[1]#a file containing list of gff files to be merged. e.g,input_files
files=pd.read_csv(input_files,sep="\t",names=['file'])

items=[]
tissues=[]
for i in range(files.shape[0]):
    file=files.iloc[i,0]
    items.append(file)
    tissue=re.split('mRNA-|_final|.collapsed', file)[1]
    tissues.append(tissue)
    
cpus=multiprocessing.cpu_count()
p = Pool()      
start = time.time()
processed_values= p.map(paralle_chunk_input_1, get_chunks(items, 1))
end = time.time()
print('total time (s)= ' + str(end-start)) 

bed_dicts={}
for i in range(len(processed_values)):
    bed_dicts.update(processed_values[i])

items=[]
for i in range(len(tissues)):
    for j in range(i+1,len(tissues)):
        items.append([i,j])

p = Pool()      
start = time.time()
processed_values= p.map(paralle_chunk_input_2, get_chunks(items, 1))
end = time.time()
print('total time (s)= ' + str(end-start)) 
intersections=pd.concat(processed_values)
######################

intersections['info']=intersections[3] + '---' + intersections[10]

ex_coords={}
for i in range(files.shape[0]):
    file=files.iloc[i,0]
    tissue=re.split('mRNA-|_final|.collapsed', file)[1]
    res_dict=get_ex_coords_dict(file,tissue)
    ex_coords.update(res_dict)

tr_ex_info=[]
for i in range(files.shape[0]):
    file=files.iloc[i,0]
    tissue=re.split('mRNA-|_final|.collapsed', file)[1]
    for line in open(file):
        raw=line.strip().split("\t")
        if raw[2] == 'exon':
            info=raw[8].split("\"")
            tr_ex_info.append([info[3].replace('PB',tissue),raw[6],raw[0]])
    
tr_ex_info_df=pd.DataFrame(tr_ex_info,columns =['tr_id','strand','chr'])
info=tr_ex_info_df.drop_duplicates()
#tr_strand=info[['tr_id','strand']]
#tr_strand=tr_strand.set_index('tr_id').to_dict()
info=dd.from_pandas(info,npartitions=cpus*2)
tr_strand_chr=info.groupby('tr_id')[['strand','chr']].apply(lambda g: g.values.tolist(),meta=('int')).compute().to_dict()
del(info)

tr_ex_info_df=pd.DataFrame(tr_ex_info_df['tr_id'].value_counts().values, index=tr_ex_info_df['tr_id'].value_counts().index, columns=['Count']).reset_index().rename({'index':'tr_id'},axis=1)## counts the numnber of exons per transcript
del(tr_ex_info)
spliced=tr_ex_info_df.loc[tr_ex_info_df['Count']>1]
unspliced=tr_ex_info_df.loc[tr_ex_info_df['Count']==1]
intersections.columns=[0, 1,2,'tr_id',4,5,6,7,8,9,10,11,12,13,'info']
m=intersections.merge(tr_ex_info_df)
m.columns=[0, 1,2,3,4,5,6,7,8,9,'tr_id',11,12,13,'info','ref_count']
m2=m.merge(tr_ex_info_df)
del(m)
m2.columns=[0, 1,2,3,4,5,6,7,8,9,10,11,12,13,'info','ref_count','query_count']
intersections=m2
del(m2)
info=intersections.pivot_table(index=['info'], aggfunc='size').reset_index().rename({0:'counts'},axis=1)
info['ref_tr'],info['query_tr']=info['info'].str.split('---',1).str
info.columns = ['info','counts','tr_id','query_tr']
m1=info.merge(tr_ex_info_df)
del(info)
m1.columns = ['info','counts','ref_tr','tr_id','ref_count']
m=m1.merge(tr_ex_info_df)
del(m1)
m.columns = ['info','counts','ref_tr','query_tr','ref_count','query_count']

intersections=dd.from_pandas(intersections,npartitions=cpus*2)# use npartitions=-1 to use all cpus According to StackOverflow, it is advised to partition the Dataframe in about as many partitions as cores your computer has, or a couple times that number, as each partition will run on a different thread and communication between them will become too costly if there are too many. https://towardsdatascience.com/trying-out-dask-dataframes-in-python-for-fast-data-analysis-in-parallel-aa960c18a915
"""
#isize=1000000
#items=get_items(intersections)
#p = Pool()      
#start = time.time()
#processed_values= p.map(get_dask_parallel, get_chunks(items, 1))
#end = time.time()
#print('total time (s)= ' + str(end-start)) 
#intersections=dd.concat(processed_values)
##intersections.to_parquet('intersections.parquet', engine='pyarrow')
"""
info1=intersections.groupby('info')[13].min().compute().reset_index()
info1.columns = ['info','first_intersected_query_exon']
info2=intersections.groupby('info')[13].max().compute().reset_index() 
info2.columns = ['info','last_intersected_query_exon']
info3=intersections.groupby('info')[6].min().compute().reset_index()
info3.columns = ['info','first_intersected_ref_exon']       
info4=intersections.groupby('info')[6].max().compute().reset_index() 
info4.columns = ['info','last_intersected_ref_exon']       


m2=m.merge(info1)
del(m)
m3=m2.merge(info2)
del(m2)
m4=m3.merge(info3)
del(m3)
m5=m4.merge(info4)
del(m4)
m5=dd.from_pandas(m5,npartitions=cpus*2)
c=m5[(m5.counts==m5.ref_count)]
candidates=c[(c.query_count==c.last_intersected_query_exon)]## all candidates; last exon of query and ref intersected. this is because our Iso-seq data are 3'end selected
group1=candidates[(candidates.counts==candidates.ref_count) & (candidates.ref_count==candidates.query_count)].compute()#candidates for exact same structure in all re and query exons
candidates=candidates.compute()# compute() convert dask to pandas dataframe
group2=candidates[~candidates['info'].isin(group1['info'])]# candidates for same atructure at all Iso-seq transcripts in RNA-seq transcripts plus extra exon(s) at RNA-seq transcript 5' end or vise versa (same 5' extra exon on 3' in reverse strand). This is because Iso-seq data are not 5' end selected
candidates=c[(c.query_count > c.last_intersected_query_exon)]#all ref exons intersected with query exon but query transcript has additional exons on their 3' end
group3=candidates[(candidates.first_intersected_query_exon == 1)].compute()#same 5' end
group4=candidates[(candidates.first_intersected_query_exon != 1)].compute()#different 5' end
c=m5[(m5.counts == m5.query_count)].compute()
c2=c[~c['info'].isin(group1['info'])]
c2=dd.from_pandas(c2,npartitions=cpus*2)
a=c2[(c2.ref_count > c2.query_count)]
a['flag']=a.last_intersected_query_exon+(a.ref_count-a.query_count)## flag show the query last exon number compatible with ref exons
del(c)
del(c2)
b=a[(a.query_count > 1)]
del(a)
group5=b[(b.ref_count == b.flag)].compute()##same 3' end but different 5' end
group6=b[(b.first_intersected_ref_exon == 1)].compute()##same 5' end but different 3' end
o=pd.concat([group5['info'].to_frame(),group6['info'].to_frame()])
c=b[~b['info'].isin(o['info'])]
group7=c[(c.ref_count > c.last_intersected_ref_exon)].compute()
del(m5)
del(o)
del(c)
del(candidates)
del(info1)
del(info2)
del(info3)
del(info4)

p.close()
p=pd.concat([group1['info'].to_frame(),group2['info'].to_frame(),group3['info'].to_frame(),group4['info'].to_frame(),group5['info'].to_frame(),group6['info'].to_frame(),group7['info'].to_frame()]).drop_duplicates()
intersections=intersections.compute()
intersections=intersections[intersections['info'].isin(p['info'])]
del(p)
intersections_dict={}
a=intersections.swifter.apply(lambda row:aggregate(row),axis=1)#swifter tutorial https://gdcoder.com/speed-up-pandas-apply-function-using-dask-or-swifter-tutorial/  
del(intersections)
#intersections_dict = np.load('intersections_dict.npy',allow_pickle='TRUE').item()
#group1.to_csv('group1',index = None, header=True,sep="\t")
#group2.to_csv('group2',index = None, header=True,sep="\t")
#group3.to_csv('group3',index = None, header=True,sep="\t")
#group4.to_csv('group4',index = None, header=True,sep="\t")
#group5.to_csv('group5',index = None, header=True,sep="\t")
#group6.to_csv('group6',index = None, header=True,sep="\t")
#group7.to_csv('group7',index = None, header=True,sep="\t")
#spliced.to_csv('spliced',index = None, header=True,sep="\t")
#unspliced.to_csv('unspliced',index = None, header=True,sep="\t")
#np.save('intersections_dict.npy', intersections_dict)
#np.save('tr_strand_chr.npy', tr_strand_chr)
#np.save('ex_coords.npy', ex_coords)

fuzzySplice=20

isize=50#size of each item (number of itersections rows passed to each cpu) in items
p = Pool()
items=get_items(group1) 
start = time.time()
processed_values= p.map(paralle_chunk_input1, get_chunks(items, 30))
end = time.time()
print('total time (s)= ' + str(end-start)) 

same_trs_spliced={}
same_trs_unspliced={}
for i in range(0,2):
    if i==0:
        outputs = [result[i] for result in processed_values]
        if len(outputs)>0:
            for j in range(len(outputs)):
                d=outputs[j]
                k=list(d.keys())
                for key in k:
                    if key not in same_trs_spliced:
                        same_trs_spliced[key]=d[key]
                    else:
                        v=d[key]
                        for value in v:
                            same_trs_spliced[key].append(value)
    elif i==1:
        outputs = [result[i] for result in processed_values]
        if len(outputs)>0:
            for k in range(len(outputs)):
                d=outputs[k]
                for key,value in d.items():
                    if key not in same_trs_unspliced:
                        same_trs_unspliced[key]=value
                    else:
                        for v in value:
                           same_trs_unspliced[key].append(v) 

"""
#start=time.time()
#intersections_dict=intersections.groupby('info').apply(lambda g: g.values.tolist(),meta=('int')).compute().to_dict()#dictionary which 'info' are the keys and lest of columns fro each info are values (fprmatted as list)
#end = time.time()
#print('total time (s)= ' + str(end-start))
#del(intersections)
"""

del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool() 
items=get_items(group2)
start = time.time()
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])
                    

del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool()
items=get_items(group3) 
start = time.time()
processed_values= p.map(paralle_chunk_input3, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])

del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool() 
items=get_items(group4)
start = time.time()
processed_values= p.map(paralle_chunk_input4, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])


del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool() 
items=get_items(group5)
start = time.time()
processed_values= p.map(paralle_chunk_input5, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])

del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool() 
items=get_items(group6)
start = time.time()
processed_values= p.map(paralle_chunk_input6, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])

del(processed_values)
#p.close()
isize=50#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool() 
items=get_items(group7)
start = time.time()
processed_values= p.map(paralle_chunk_input7, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(0,3):
    if i==0:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced:
                    same_trs_spliced[key]=dict_t[key]
                else:
                    elements=dict_t[key]
                    ref_elements=same_trs_spliced[key]
                    for x in range(len(elements)):
                        if elements[x] not in ref_elements:
                            same_trs_spliced[key].append(elements[x]) 
    if i==1:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained:
                    same_trs_spliced_with_retained[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained[key].update(dict_t[key]) 
    if i==2:
        for result in processed_values:
            dict_t=result[i]
            for key in dict_t:
                if key not in same_trs_spliced_with_retained2:
                    same_trs_spliced_with_retained2[key]=dict_t[key]
                else:
                    same_trs_spliced_with_retained2[key].update(dict_t[key])

#save
#np.save('same_trs_spliced.npy', same_trs_spliced)
#np.save('same_trs_spliced_with_retained.npy', same_trs_spliced_with_retained)

#load
#same_trs_spliced = np.load('same_trs_spliced.npy',allow_pickle='TRUE').item()


##############################
#  phase2: merge transcripts #
#         (spliced)          #
############################## 
same_trs_spliced2=copy.deepcopy(same_trs_spliced)

for key, value in same_trs_spliced2.items():
    value.append(key)

lists=list(same_trs_spliced2.values())

g = networkx.Graph()#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

merges=list(networkx.connected_components(g))#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
#del(same_trs_spliced2)

isize=1#size of each item (number of itersections rows passed to each cpu) in items
p = Pool()
items=get_items_fromLists(merges) 
start = time.time()
processed_values= p.map(get_mergable_same_trs_spliced_parallel, get_chunks(items, 1))
#p.close()
end = time.time()
print('total time (s)= ' + str(end-start))             
mergable_same_trs_spliced=[]
for i in range(len(processed_values)): 
    for j in range(len(processed_values[i])):
        mergable_same_trs_spliced.append(processed_values[i][j])

#with open("mergable_same_trs_spliced.txt", "w") as file:
#    file.write(str(mergable_same_trs_spliced))
#with open("mergable_same_trs_spliced.txt", "r") as file:
#    mergable_same_trs_spliced = eval(file.readline())
        
        
#test_list = [[1, 0, -1], [-1, 0, 1], [-1, 0, 1], [1, 2, 3], [3, 4, 1],[2,1,3]]
#test_list=[['hamid','mina','mona'],['hamid','mona','mina']]
#res = list(set(map(lambda i: tuple(sorted(i)), test_list))) 
"""
isize=5#size of each item (number of itersections rows passed to each cpu) in items
#p = Pool()
items=get_items_fromLists(mergable_same_trs_spliced) 
start = time.time()
processed_values= p.map(get_redundancy_parallel, get_chunks(items, 1))
#p.close()
end = time.time()
print('total time (s)= ' + str(end-start)) 
redundants=[]
for i in range(len(processed_values)): 
    for j in range(len(processed_values[i])):
        redundants.append(processed_values[i][j])
out=[]
for i in mergable_same_trs_spliced:
    if i not in redundants:
        out.append(i)
        
mergable_same_trs_spliced=out
del(out)
"""
tr_tissues={}
transcripts_dict={}
spliced_merged_id_info={}
counter=1
for l in mergable_same_trs_spliced:
    l_tissues=[]
    intervals=[]
    ll=[]
    for i in l:
        tissue=i.split('.')[0]
        if tissue not in l_tissues:
            l_tissues.append(tissue)
    key='merge' + '.' + str(counter) + '.' + str(counter)
    counter=counter+1
    tr_tissues[key]=l_tissues
    transcripts_dict[key]=l
    for tr in l:
        ll.append(ex_coords[tr])
    s=sorted(ll, key=len,reverse=True)
    ref_ex=s[0]
    first=ref_ex[0]
    last=ref_ex[-1]
    rest=s[1:]
    for exon in ref_ex:
        intervals.append(exon)
    for ex_set in rest:
        info=ex_set[:]
        tr_first=ex_set[0]
        tr_last=ex_set[-1]
        if check_exon_position(first,tr_first,'first')=="included":
            info=info[1:]
        if check_exon_position(last,tr_last,'last')=="included":
            info=info[0:-1]
        for exon in info:
            intervals.append(exon)
    new_ref_ex=merge_intervals(intervals)
    transcripts_dict[key]=new_ref_ex
    spliced_merged_id_info[key]=l[0]

#processed=list(flatten_nested_list(mergable_same_trs_spliced))          
v=list(same_trs_spliced2.values())
processed=list(flatten_nested_list(v))     
processed=list(dict.fromkeys(processed))
processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)

remain_spliced=spliced[~spliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

#####################################
#  phase3: get remained transcripts #
#               (spliced)           #
#####################################

# spliced transcripts in both ref and query that have unique structure compared to other transcriptome.
for i in range(remain_spliced.shape[0]):
    key=remain_spliced.iloc[i,0]
    exons=ex_coords[key]
    transcripts_dict[key]=exons

#global remain_spliced
isize=50#size of each item (number of itersections rows passed to each cpu) in items
p = Pool()
items=get_items(remain_spliced) 
start = time.time()
processed_values= p.map(get_tissue_parallel, get_chunks(items, 30))
end = time.time()
print('total time (s)= ' + str(end-start)) 

for i in range(len(processed_values)):
    tr_tissues.update(processed_values[i])
    
############################################
#  phase4: get transcripts representatives #
#           (unspliced)                    #
############################################  

for key, value in same_trs_unspliced.items():
    value.append(key)

lists=same_trs_unspliced.values()

g = networkx.Graph()
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

m=list(networkx.connected_components(g))

updated_same_trs_unspliced={}
for l in m:
    l=list(l)
    key=l[0]
    val=l[1:]
    updated_same_trs_unspliced[key]=val
    

##############################
#  phase5: merge transcripts #
#         (unspliced)        #
############################## 
processed=[] 
for key, value in updated_same_trs_unspliced.items():
    t=[]
    t.append(key.split('.')[0])
    intervals=[]
    ref_ex=ex_coords[key]
    for exon in ref_ex:
        intervals.append(exon)
    for tr in value:
        processed.append(tr)
        info=ex_coords[tr]
        for exon in info:
            intervals.append(exon)
        tis=tr.split('.')[0]
        if tis not in t:
            t.append(tis)
    new_ref_ex=merge_intervals(intervals)
    transcripts_dict[key]=new_ref_ex
    tr_tissues[key]=t
    
v=list(same_trs_unspliced.values())
processed=list(flatten_nested_list(v))     
processed=list(dict.fromkeys(processed))
processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

""" I YOU WANT TO INCLUDE TISSUE SPECIFIC UNSPLICED TRANSCRIPTS 
#processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
#remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

#unspliced_umerged_final={}# unspliced transcripts in both ref and query that have unique structure compared to other transcriptome.
#for i in range(remain_unspliced.shape[0]):
#    key=remain_unspliced.iloc[i,0]
#    exons=get_exons(key)
#    unspliced_umerged_final[key]=exons
"""
""""""""""""  

##############################
#  phase6: get gene ids      # 
#          and format trids  #
############################## 

results={}
for key, value in transcripts_dict.items():
    if key.find('merge')==-1:
        info=tr_strand_chr[key][0]
        results[key]=[info[1],value[0][0],value[-1][1],key,'99',info[0],key]
    else:
        tr=spliced_merged_id_info[key]
        info=tr_strand_chr[tr][0]
        results[key]=[info[1],value[0][0],value[-1][1],key,'99',info[0],key]

bed=pd.DataFrame.from_dict(results,orient='index')##if genes defined as group of transcripts with exoni c overlap, this shoud include exons. otherwise transcripts
bed=pybedtools.BedTool.from_dataframe(bed)
intersection=bed.intersect(bed,wa=True,wb=True,s=True).to_dataframe()
del(bed)
intersection.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]
p = Pool() 
items=[]
for i in range(0,intersection.shape[0],1000):
    if i<intersection.shape[0]:
        items.append([i,i+1000])
    else:
        items.append([i,i+11000])

start = time.time()
processed_values= p.map(filter_intersection, get_chunks(items, 1))
end = time.time()
print('total time (s)= ' + str(end-start)) 
intersection=pd.concat(processed_values)

df1=intersection[[3,13]].rename({3:'A',13:'B'},axis=1) 
df2=intersection[[13,3]].rename({13:'A',3:'B'},axis=1)
df=pd.concat([df1,df2]).drop_duplicates()
df=dd.from_pandas(df,npartitions=cpus*2)
my_dict=df.groupby('A')['B'].apply(lambda g: g.values.tolist(),meta=('int')).compute().to_dict()
#my_dict=df.groupby('A')['B'].apply(lambda g: g.values.tolist()).to_dict()
del(df1)
del(df2)
del(df)

for key, value in my_dict.items():
    value.append(key)

lists=my_dict.values()
g = networkx.Graph()
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

m=list(networkx.connected_components(g))

result_dict={}
for l in m:
    l=list(l)
    k=l[0]
    result_dict[k]=l

val=list(result_dict.values())
processed=list(flatten_nested_list(val))

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
all_trs=pd.DataFrame(transcripts_dict.keys()).rename({0:'tr_id'},axis=1)
r=all_trs[~all_trs['tr_id'].isin(processed['tr_id'])]['tr_id'].unique()#

for i in r:
    result_dict[i]=[i]

""" this section is aimed to sort final transcripts
    based on their location to get sorted gene ids 
    THIS CAN BE SAFELY IGNORED """ 
""" THIS JUST WORKS IN PYTHON2.7
d={}
for key in result_dict.keys():
    if key.find('merge')==-1:
        info=tr_strand_chr[key][0]
        value=transcripts_dict[key]
        chrom=info[1].split('chr')[1]
        if chrom.find('X')==-1 and chrom.find('ref')==-1 and chrom.find('MT')==-1:
            chrom=int(chrom)
        d[key]=[chrom,value[0][0]]
    else:
        tr=spliced_merged_id_info[key]
        info=tr_strand_chr[tr][0]
        chrom=info[1].split('chr')[1]
        if chrom.find('X')==-1 and chrom.find('ref')==-1 and chrom.find('MT')==-1:
            chrom=int(chrom)
        d[key]=[chrom,value[0][0]]

d=OrderedDict(sorted(d.items(), key=lambda x: x[1],reverse=False))    

out=collections.OrderedDict()
for key in list(d.keys()):
    val=result_dict[key]
    out[key]=val 
    
result_dict=out
"""
""""""
tr_ids={}
g_counter=0
for key, value in result_dict.items():
    g_counter=g_counter+1
    gid='PB.' + str(g_counter)
    tr_counter=0
    for i in value:
        tr_counter=tr_counter+1
        tid=gid + '.' + str(tr_counter)
        tr_ids[i]=tid
        
##############################
#  phase7: create            #
#       .gff file            #
############################## 

gff_list=[]
for key, value in transcripts_dict.items():
    if key.find('merge')==-1:
        info=tr_strand_chr[key][0]
    else:
        tr=spliced_merged_id_info[key]
        info=tr_strand_chr[tr][0]  
    tid=tr_ids[key]
    a=tid.split('.')[0:2]
    gid=a[0]+'.'+a[1]
    transcript=[info[1],'CattleFAANG','transcript', value[0][0],value[-1][1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
    gff_list.append(transcript)
    for i in value:
        exon=[info[1],'CattleFAANG','exon',i[0],i[1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
        gff_list.append(exon)
    
    
gff_df=pd.DataFrame(gff_list) 



remain_unspliced_list=list(remain_unspliced['tr_id'])
gff_list2=[]#some of transcripts in this might be intersected with transcripts in 'gff_list'. As I'm going to include TSS unspliced transcripts related to known genes that were not captured by 'gff_list', I don't care about these transcripts
remain_unspliced_tissues={}
for i in remain_unspliced_list:
    t=i.split('.')[0]
    g_counter=g_counter+1
    gid='PB.' + str(g_counter)
    tid=gid + '.' + str(1)
    info=tr_strand_chr[i][0]
    value=ex_coords[i]
    transcript=[info[1],'CattleFAANG','transcript', value[0][0],value[-1][1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
    exon=[info[1],'CattleFAANG','exon',value[0][0],value[-1][1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
    gff_list2.append(transcript)
    gff_list2.append(exon)
    remain_unspliced_tissues[tid]=[t]

gff_df2=pd.DataFrame(gff_list2)
####################################
#       phase8: get                #
# transcripts/genes replication    #   
#    across tissues                #
#################################### 
#change keys to formatted keys
tr_tissue={}#detected tissue per transcript
for k,v in tr_tissues.items():
    new_key=tr_ids[k]
    tr_tissue[new_key]=v

del(tr_tissues)

gene_tissue={}#detected tissue per gene
for k,v in tr_tissue.items():
    g=k.split('.')[0] + '.' + k.split('.')[1]
    if g not in gene_tissue:
        gene_tissue[g]=v
    else:
        v2=gene_tissue[g][:]
        for i in v:
            if i not in v2:
                v2.append(i)
        gene_tissue[g]=v2
        
df=pd.DataFrame.from_dict(tr_tissue,orient='index').reset_index()
df_tr=df.sort_values(by=['index'])

df=pd.DataFrame.from_dict(gene_tissue,orient='index').reset_index()
df_g=df.sort_values(by=['index'])

df=pd.DataFrame.from_dict(remain_unspliced_tissues,orient='index').reset_index()
df_tr2=df.sort_values(by=['index'])

#########################
#   phase9: get         #
# pairwise tissue       #   
#   comparisions        #
######################### 

#get tissue pairs
count=0
tissue_pairs={}
for i in range(len(tissues)-1):
    for j in range(i+1,len(tissues)):
        tissue_pairs[count]=[tissues[i],tissues[j]]
        count=count+1

""" TRANSCRIPT level """       
tissue_tr={}# detected transcripts per tissue
for i in tissues:
    val=[]
    for k,v in tr_tissue.items():
        if list(set([i]) & set(v)):
            val.append(k)
    tissue_tr[i]=val


tissue_sim_tr=collections.OrderedDict()
tissue_sim_tr2=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_tr[t1]
        b=tissue_tr[t2]
        value=list(set(a) & set(b))
        v.append(len(value))
        v2.append(round((len(value)/float(len(a)))*100,0))
    tissue_sim_tr[tissues[i]]=v
    tissue_sim_tr2[tissues[i]]=v2

    
tissue_sim_mat_tr=format_mat(tissue_sim_tr)
tissue_sim_mat_tr2=format_mat(tissue_sim_tr2)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)


""" GENE level """

tissue_gene={}# detected genes per tissue
for i in tissues:
    val=[]
    for k,v in gene_tissue.items():
        if list(set([i]) & set(v)):
            val.append(k)
    tissue_gene[i]=val

tissue_sim_gene=collections.OrderedDict()
tissue_sim_gene2=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_gene[t1]
        b=tissue_gene[t2]
        value=list(set(a) & set(b))
        v.append(len(value))
        v2.append(round((len(value)/float(len(a)))*100,0))
    tissue_sim_gene[tissues[i]]=v
    tissue_sim_gene2[tissues[i]]=v2

tissue_sim_mat_gene=format_mat(tissue_sim_gene)
tissue_sim_mat_gene2=format_mat(tissue_sim_gene2)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)

#########################
#       pahes10:        #
#    export outputs     #
######################### 

gff_df.to_csv('combined.gff',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)
df_tr.to_csv('transcripts_tracking.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
c=df_tr.count(axis=1)-1
c=list(c) 
plt.hist(c, 25, histtype='step', color='black')
plt.savefig('transcripts_replication_hist.png') 
plt.close()  
print('number of multi-tissue transcripts:')
print('\n')
print(len([i for i in c if i > 1]))
print('number of single-tissue transcripts:')
print('\n')
print(len([i for i in c if i == 1]))
df_g.to_csv('genes_tracking.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
c=df_g.count(axis=1)-1
c=list(c) 
plt.hist(c, 25, histtype='step', color='black')
plt.savefig('genes_replication_hist.png')   
print('number of multi-tissue genes:')
print('\n')
print(len([i for i in c if i > 1]))
print('number of single-tissue genes:')
print('\n')
print(len([i for i in c if i == 1]))
tissue_sim_mat_tr.to_csv('transcript_level_tissue_comparisions.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_tr2.to_csv('transcript_level_tissue_comparisions_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal) 
tissue_sim_mat_gene.to_csv('gene_level_tissue_comparisions.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_gene2.to_csv('gene_level_tissue_comparisions_HEATMAP.txt',index = None, header=True,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)  
np.save('tissue_gene.npy', tissue_gene)
np.save('gene_tissue.npy', gene_tissue)
np.save('tissue_transcript.npy', tissue_tr)
np.save('transcript_tissue.npy', tr_tissue)
tss_unspliced_tr=remain_unspliced['tr_id'].to_frame()
tss_unspliced_tr.to_csv('tissue_specific_unsplised_transcripts.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)
gff_df2.to_csv('tissue_specific_unsplised_transcripts.gff',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)    
df_tr2.to_csv('tissue_specific_unsplised_transcripts_tracking.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)
#load
#import numpy as np
#tissue_gene = np.load('tissue_gene.npy',allow_pickle='TRUE').item()


#intersections_dict = np.load('intersections_dict.npy',allow_pickle='TRUE').item()
#tr_strand_chr = np.load('tr_strand_chr.npy',allow_pickle='TRUE').item()
#ex_coords = np.load('ex_coords.npy',allow_pickle='TRUE').item()
#group1=pd.read_csv('group1',header=0,sep="\t")
#group2=pd.read_csv('group2',header=0,sep="\t")
#group3=pd.read_csv('group3',header=0,sep="\t")
#group4=pd.read_csv('group4',header=0,sep="\t")
#group5=pd.read_csv('group5',header=0,sep="\t")
#group6=pd.read_csv('group6',header=0,sep="\t")
#group7=pd.read_csv('group7',header=0,sep="\t")
#spliced=pd.read_csv('spliced',header=0,sep="\t")
#unspliced=pd.read_csv('unspliced',header=0,sep="\t")

#with open("mergable_same_trs_spliced.txt", "w") as file:
#    file.write(str(mergable_same_trs_spliced))
#with open("mergable_same_trs_spliced.txt", "r") as file:
#    mergable_same_trs_spliced = eval(file.readline())



""" based on mering all intersected transcripts(NOT CORRECT, just for my own recoreds)   
############################################
#  phase2: get transcripts representatives #
#             (spliced)                    #
############################################  

""merge intersected lists within neseted list""
for key, value in same_trs_spliced.items():
    value.append(key)

lists=same_trs_spliced.values()

g = networkx.Graph()#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

m=list(networkx.connected_components(g))#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections    
    
updated_same_trs_spliced={}
for l in m:
    info={}
    for i in l:
        info[i]=ex_coords[i]
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))#sort dictionari by length of values (list) https://stackoverflow.com/questions/50863093/how-to-sort-a-python-dictionary-based-on-the-length-of-the-list-of-the-values
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    updated_same_trs_spliced[new_key]=new_value

updated_same_trs_spliced_with_retained={}
for key, value in same_trs_spliced_with_retained.items():
    info={}
    info[key]=ex_coords[key]
    for key2 in value:
        info[key2]=ex_coords[key2]
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    if new_key not in updated_same_trs_spliced_with_retained:
        updated_same_trs_spliced_with_retained[new_key]=new_value
    elif new_key in updated_same_trs_spliced_with_retained:
        v=updated_same_trs_spliced_with_retained[new_key][:]
        for i in new_value:
            if i not in v:
                v.append(i)
        updated_same_trs_spliced_with_retained[new_key]=v


updated_same_trs_spliced_with_retained2={}
for key, value in same_trs_spliced_with_retained2.items():
    info={}
    info[key]=ex_coords[key]
    for key2 in value:
        info[key2]=ex_coords[key2]
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    if new_key not in updated_same_trs_spliced_with_retained2:
        updated_same_trs_spliced_with_retained2[new_key]=new_value
    elif new_key in updated_same_trs_spliced_with_retained2:
        v=updated_same_trs_spliced_with_retained2[new_key][:]
        for i in new_value:
            if i not in v:
                v.append(i)
        updated_same_trs_spliced_with_retained2[new_key]=v
                   

##############################
#  phase3: merge transcripts #
#         (spliced)          #
##############################  

spliced_merged_structures={}
for key, value in updated_same_trs_spliced.items():
    intervals=[]
    ref_ex=ex_coords[key]
    first=ref_ex[0]
    last=ref_ex[-1]
    for exon in ref_ex:
        intervals.append(exon)
    for tr in value:
        info=ex_coords[tr]
        tr_first=info[0]
        tr_last=info[-1]
        if check_exon_position(first,tr_first,'first')=="included":
            info=info[1:]
        if check_exon_position(last,tr_last,'last')=="included":
            info=info[0:-1]
        for exon in info:
            intervals.append(exon)
    new_ref_ex=merge_intervals(intervals)
    spliced_merged_structures[key]=new_ref_ex


spliced_with_retained_merged_structures={}
for key, value in updated_same_trs_spliced_with_retained.items():
    intervals=[]
    if key in same_trs_spliced_with_retained:
        info_dict=same_trs_spliced_with_retained[key]
        ref_ex=ex_coords[key]
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        for key2,value2 in info_dict.items():
            query_ex=ex_coords[key2]
            l=len(query_ex)
            query_ex=adjust_ex(query_ex,[value2[0]])
#            query_ex.pop(value)
            if 0 < value2[0] < l:
                if check_exon_position(first,query_ex[0],'first')=="included":
                    query_ex=query_ex[1:]
                if check_exon_position(last,query_ex[-1],'last')=="included":
                    query_ex=query_ex[0:-1]
            elif value2[0]==0 and check_exon_position(last,query_ex[-1],'last')=="included":
                query_ex=query_ex[0:-1]
            elif value2[0]==l and check_exon_position(first,query_ex[0],'first')=="included":
                query_ex=query_ex[1:]
            for exon in query_ex:
                intervals.append(exon)
        new_ref_ex=merge_intervals(intervals)
        spliced_with_retained_merged_structures[key]=new_ref_ex
    else:
        main_key=get_key(key,same_trs_spliced_with_retained)
        info=updated_same_trs_spliced_with_retained[key]
        info_main=same_trs_spliced_with_retained[main_key]
        ref_ex=ex_coords[key]
        ref_ex=adjust_ex(ref_ex,[info_main[key][0]])
#        ref_ex.pop(info_main[key][0])
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        main_ex=ex_coords[main_key]
        for exon in main_ex:
            intervals.append(exon)
        for key2,value2 in info_main.items():
            if key2 !=key:
                query_ex=ex_coords[key]
                l=len(query_ex)
                query_ex=adjust_ex(query_ex,[value2[0]])
                if 0 < value2[0] < l:
                    if check_exon_position(first,query_ex[0],'first')=="included":
                        query_ex=query_ex[1:]
                    if check_exon_position(last,query_ex[-1],'last')=="included":
                        query_ex=query_ex[0:-1]
                elif value2[0]==0 and check_exon_position(last,query_ex[-1],'last')=="included":
                    query_ex=query_ex[0:-1]
                elif value2[0]==l and check_exon_position(first,query_ex[0],'first')=="included":
                    query_ex=query_ex[1:]
                for exon in query_ex:
                    intervals.append(exon)
            new_ref_ex=merge_intervals(intervals)
            spliced_with_retained_merged_structures[main_key]=new_ref_ex                
        
spliced_with_retained_merged_structures2={}
for key, value in updated_same_trs_spliced_with_retained2.items():
    intervals=[]
    if key in same_trs_spliced_with_retained2:
        info_dict=same_trs_spliced_with_retained2[key]
        ref_ex=ex_coords[key]
        out=[]
        for value2 in info_dict.values():
            out.append(value2[0])
        ref_ex=adjust_ex(ref_ex,out)
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        for tr in info_dict.keys():
            info=ex_coords[tr]
            tr_first=info[0]
            tr_last=info[-1]
            if check_exon_position(first,tr_first,'first')=="included":
                info=info[1:]
            if check_exon_position(last,tr_last,'last')=="included":
                info=info[0:-1]
            for exon in info:
                intervals.append(exon)
        new_ref_ex=merge_intervals(intervals)
        spliced_with_retained_merged_structures2[key]=new_ref_ex
    else:
        main_key=get_key(key,same_trs_spliced_with_retained2)
        info=updated_same_trs_spliced_with_retained2[key]
        info_main=same_trs_spliced_with_retained2[main_key]
        ref_ex=ex_coords[key]
        main_ex=ex_coords[main_key]
        for exon in ref_ex:
            intervals.append(exon)
        out=[]
        for value2 in info_main.values():
            out.append(value2[0])
        main_ex=adjust_ex(main_ex,out)
        for exon in main_ex:
            intervals.append(exon)
        for tr in info_main.keys():
            if tr !=key:
                q_ex=ex_coords[tr]
                for exon in q_ex:
                    intervals.append(exon)
        new_ref_ex=merge_intervals(intervals)
        spliced_with_retained_merged_structures2[main_key]=new_ref_ex


" merge merged trascripts "
spliced_merged_final={}
for key,value in spliced_merged_structures.items():
    intervals=[]
    info=[]
    info.append(value)
    if key in spliced_with_retained_merged_structures.keys():
        info.append(spliced_with_retained_merged_structures[key])
    if key in spliced_with_retained_merged_structures2.keys():
        info.append(spliced_with_retained_merged_structures2[key])   
    info.sort(key=len,reverse=True)#sort list of lists according to length of sublists https://stackoverflow.com/questions/30346356/how-to-sort-list-of-lists-according-to-length-of-sublists
    ref_ex=info[0]
    first=ref_ex[0]
    last=ref_ex[-1]
    for ex in ref_ex:
        intervals.append(ex)
    if len(info)>=2:
        info=info[1:]
        for exons in info:
            tr_first=exons[0]
            tr_last=exons[-1]
            if check_exon_position(first,tr_first,'first')=="included":
                exons=exons[1:]
            if check_exon_position(last,tr_last,'last')=="included":
                exons=exons[0:-1]
            for exon in exons:
                intervals.append(exon)
    new_ref_ex=merge_intervals(intervals)
    spliced_merged_final[key]=new_ref_ex
    
for key,value in spliced_with_retained_merged_structures.items():
    if key not in spliced_merged_final.keys():
        intervals=[]
        info=[]
        info.append(value)
        if key in spliced_with_retained_merged_structures2.keys():
            info.append(spliced_with_retained_merged_structures2[key])
        info.sort(key=len,reverse=True)
        ref_ex=info[0]
        first=ref_ex[0]
        last=ref_ex[-1]
        for ex in ref_ex:
            intervals.append(ex)
        if len(info)>=2:
            info=info[1:]
            for exons in info:
                tr_first=exons[0]
                tr_last=exons[-1]
                if check_exon_position(first,tr_first,'first')=="included":
                    exons=exons[1:]
                if check_exon_position(last,tr_last,'last')=="included":
                    exons=exons[0:-1]
                for exon in exons:
                    intervals.append(exon) 
        new_ref_ex=merge_intervals(intervals)
        spliced_merged_final[key]=new_ref_ex


for key,value in spliced_with_retained_merged_structures2.items():
    if key not in spliced_merged_final.keys():    
        spliced_merged_final[key]=value

############################################
#  phase4: get transcripts representatives #
#           (unspliced)                    #
############################################  

for key, value in same_trs_unspliced.items():
    value.append(key)

lists=same_trs_unspliced.values()

g = networkx.Graph()
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

m=list(networkx.connected_components(g))

updated_same_trs_unspliced={}
for l in m:
    l=list(l)
    key=l[0]
    val=l[1:]
    updated_same_trs_unspliced[key]=val
    
#####################################
#  phase5: get remained transcripts #
#               (spliced)           #
#####################################  

processed=[]
for key,value in updated_same_trs_spliced.items():
    processed.append(key)
    for i in value:
        processed.append(i)
        
for key,value in updated_same_trs_spliced_with_retained.items():
    processed.append(key)
    for i in value:
        processed.append(i)   

for key,value in updated_same_trs_spliced_with_retained2.items():
    processed.append(key)
    for i in value:
        processed.append(i) 

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_spliced=spliced[~spliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

remained_spliced_dict={}# spliced transcripts in both ref and query that have unique structure compared to other transcriptome.
for i in range(remain_spliced.shape[0]):
    key=remain_spliced.iloc[i,0]
    exons=ex_coords[key]
    remained_spliced_dict[key]=exons


transcripts_dict=spliced_merged_final.copy()

transcripts_dict.update(remained_spliced_dict)## mixure of both merged and unmerged spliced transcripts

##############################
#  phase6: merge transcripts #
#         (unspliced)        #
############################## 
 
merged_structures2={}#unspliced
processed=[]
for key, value in updated_same_trs_unspliced.items():
    processed.append(key)
    intervals=[]
    ref_ex=ex_coords[key]
    for exon in ref_ex:
        intervals.append(exon)
    for tr in value:
        processed.append(tr)
        info=ex_coords[tr]
        for exon in info:
            intervals.append(exon)
    new_ref_ex=merge_intervals(intervals)
    merged_structures2[key]=new_ref_ex

transcripts_dict.update(merged_structures2)## mixure of merged, unmerged spliced transcripts and merge unspliced transcripts

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

"" I YOU WANT TO INCLUDE TISSUE SPECIFIC UNSPLICED TRANSCRIPTS 
processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

unspliced_umerged_final={}# unspliced transcripts in both ref and query that have unique structure compared to other transcriptome.
for i in range(remain_unspliced.shape[0]):
    key=remain_unspliced.iloc[i,0]
    exons=get_exons(key)
    unspliced_umerged_final[key]=exons
""
""""""""""""     

##############################
#  phase7: get gene ids      # 
#          and format trids  #
############################## 

results={}
for key, value in transcripts_dict.items():
    info=tr_strand_chr[key][0]
    results[key]=[info[1],value[0][0],value[-1][1],key,'99',info[0],key]

bed=pd.DataFrame.from_dict(results,orient='index')##if genes defined as group of transcripts with exoni c overlap, this shoud include exons. otherwise transcripts
bed=pybedtools.BedTool.from_dataframe(bed)
intersection=bed.intersect(bed,wa=True,wb=True,s=True).to_dataframe()
del(bed)
intersection.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]
p = Pool() 
items=[]
for i in range(0,intersection.shape[0],1000):
    if i<intersection.shape[0]:
        items.append([i,i+1000])
    else:
        items.append([i,i+11000])

start = time.time()
processed_values= p.map(filter_intersection, get_chunks(items, 1))
end = time.time()
print('total time (s)= ' + str(end-start)) 
intersection=pd.concat(processed_values)

df1=intersection[[3,13]].rename({3:'A',13:'B'},axis=1) 
df2=intersection[[13,3]].rename({13:'A',3:'B'},axis=1)
df=pd.concat([df1,df2]).drop_duplicates()
df=dd.from_pandas(df,npartitions=-1)
my_dict=df.groupby('A')['B'].apply(lambda g: g.values.tolist(),meta=('int')).compute().to_dict()
#my_dict=df.groupby('A')['B'].apply(lambda g: g.values.tolist()).to_dict()
del(df1)
del(df2)
del(df)

for key, value in my_dict.items():
    value.append(key)

lists=my_dict.values()
g = networkx.Graph()
for sub_list in lists:
    for edge in pairs(sub_list):
            g.add_edge(*edge)

m=list(networkx.connected_components(g))

result_dict={}
for l in m:
    l=list(l)
    k=l[0]
    result_dict[k]=l

val=result_dict.values()[:]
processed=list(flatten_nested_list(val))

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
all_trs=pd.DataFrame(transcripts_dict.keys()).rename({0:'tr_id'},axis=1)
r=all_trs[~all_trs['tr_id'].isin(processed['tr_id'])]['tr_id'].unique()#

for i in r:
    result_dict[i]=[i]

"" this section is aimed to sort final transcripts
    based on their location to get sorted gene ids 
    THIS CAN BE SAFELY IGNORED "" 

d={}
for key in result_dict.keys():
    info=tr_strand_chr[key][0]
    value=transcripts_dict[key]
    chrom=info[1].split('chr')[1]
    if chrom.find('X')==-1 and chrom.find('ref')==-1 and chrom.find('MT')==-1:
        chrom=int(chrom)
    d[key]=[chrom,value[0][0]]
    
d=OrderedDict(sorted(d.items(), key=lambda x: x[1],reverse=False))

out=collections.OrderedDict()
for key in d.keys():
    val=result_dict[key]
    out[key]=val 
    
result_dict=out

""""""
tr_ids={}
g_counter=0
for key, value in result_dict.items():
    g_counter=g_counter+1
    gid='PB.' + str(g_counter)
    tr_counter=0
    for i in value:
        tr_counter=tr_counter+1
        tid=gid + '.' + str(tr_counter)
        tr_ids[i]=tid
        
##############################
#  phase8: create            #
#       .gff file            #
############################## 

gff_list=[]
for key, value in transcripts_dict.items():
    info=tr_strand_chr[key][0]
    tid=tr_ids[key]
    a=tid.split('.')[0:2]
    gid=a[0]+'.'+a[1]
    transcript=[info[1],'CattleFAANG','transcript', value[0][0],value[-1][1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
    gff_list.append(transcript)
    for i in value:
        exon=[info[1],'CattleFAANG','exon',i[0],i[1],'.',info[0],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
        gff_list.append(exon)
    
    
gff_df=pd.DataFrame(gff_list)    


####################################
#       phase9: get                #
# transcripts/genes replication    #   
#    across tissues                #
#################################### 

tr_tissue_rep={}
for key, value in updated_same_trs_spliced.items():
    val=value[:]
    if key in same_trs_spliced_with_retained:
        d=same_trs_spliced_with_retained[key]
        for i in d.keys():
            val.append(i)
    if key in same_trs_spliced_with_retained2:
        d=same_trs_spliced_with_retained2[key]
        for i in d.keys():
            val.append(i)    
    v=[]
    v.append(key.split('.')[0])
    for i in val:
        t=i.split('.')[0]
        if t not in v:
            v.append(t)
    tr_tissue_rep[key]=v

a=spliced_merged_final.keys()
b=tr_tissue_rep.keys()
rem=list(set(a) ^ set(b))

""
count=0
for i in rem:
    if i in updated_same_trs_spliced_with_retained.keys():
       count=count+1 
    else:
        main_key=get_key2(i,updated_same_trs_spliced_with_retained)
        count=count+1
""

for i in rem:
    v=[]
    if i in updated_same_trs_spliced_with_retained:
        val=updated_same_trs_spliced_with_retained[i]
        if i in updated_same_trs_spliced_with_retained2:
            for j in updated_same_trs_spliced_with_retained2[i]:
                val.append(j)
        elif get_key2(i,updated_same_trs_spliced_with_retained2):
            m=get_key2(i,updated_same_trs_spliced_with_retained2)
            for j in updated_same_trs_spliced_with_retained2[m]:
                val.append(j)
        for t in val:
            if t not in v:
                v.append(t.split('.')[0])
    elif get_key2(i,updated_same_trs_spliced_with_retained):
        main_key=get_key2(i,updated_same_trs_spliced_with_retained)
        val=updated_same_trs_spliced_with_retained[main_key]
        if main_key in updated_same_trs_spliced_with_retained2:
            for j in updated_same_trs_spliced_with_retained2[main_key]:
                val.append(j)
        elif get_key2(main_key,updated_same_trs_spliced_with_retained2):
            m=get_key2(main_key,updated_same_trs_spliced_with_retained2)
            for j in updated_same_trs_spliced_with_retained2[m]:
                val.append(j)
        for t in val:
            if t not in v:
                v.append(t.split('.')[0])
    tr_tissue_rep[i]=v


for key in remained_spliced_dict:
        tr_tissue_rep[key]=[key.split('.')[0]] 
        

for key in merged_structures2:
        tr_tissue_rep[key]=[key.split('.')[0]] 

"" change keys to formatted keys "" 
tr_tissue={}# detected tissue per transcripts
for key,value in tr_tissue_rep.items():
    k=tr_ids[key]
    v=value[:]
    tr_tissue[k]=list(dict.fromkeys(v))

df=pd.DataFrame.from_dict(tr_tissue,orient='index').reset_index()
df_tr=df.sort_values(by=['index'])

gene_tissue={}#detected tissue per gene
for k,v in tr_tissue.items():
    g=k.split('.')[0] + '.' + k.split('.')[1]
    if g not in gene_tissue:
        gene_tissue[g]=v
    else:
        v2=gene_tissue[g][:]
        for i in v:
            if i not in v2:
                v2.append(i)
        gene_tissue[g]=v2

df=pd.DataFrame.from_dict(gene_tissue,orient='index').reset_index()
df_g=df.sort_values(by=['index'])
        
#########################
#   phase11: get        #
# pairwise tissue       #   
#   comparisions        #
######################### 

#get tissue pairs
count=0
tissue_pairs={}
for i in range(len(tissues)-1):
    for j in range(i+1,len(tissues)):
        tissue_pairs[count]=[tissues[i],tissues[j]]
        count=count+1

"" TRANSCRIPT level ""       
tissue_tr={}# detected transcripts per tissue
for i in tissues:
    val=[]
    for k,v in tr_tissue.items():
        if list(set([i]) & set(v)):
            val.append(k)
    tissue_tr[i]=val


tissue_sim_tr=collections.OrderedDict()
tissue_sim_tr2=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_tr[t1]
        b=tissue_tr[t2]
        value=list(set(a) & set(b))
        v.append(len(value))
        v2.append(round((len(value)/float(len(a)))*100,0))
    tissue_sim_tr[tissues[i]]=v
    tissue_sim_tr2[tissues[i]]=v2

    
tissue_sim_mat_tr=format_mat(tissue_sim_tr)
tissue_sim_mat_tr2=format_mat(tissue_sim_tr2)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal)

" GENE level "

tissue_gene={}# detected genes per tissue
for i in tissues:
    val=[]
    for k,v in gene_tissue.items():
        if list(set([i]) & set(v)):
            val.append(k)
    tissue_gene[i]=val

tissue_sim_gene=collections.OrderedDict()
tissue_sim_gene2=collections.OrderedDict()#this is for heatmap
for i in range(len(tissues)):
    v=[]
    v2=[]
    for j in range(len(tissues)):
        t1=tissues[i]
        t2=tissues[j]
        a=tissue_gene[t1]
        b=tissue_gene[t2]
        value=list(set(a) & set(b))
        v.append(len(value))
        v2.append(round((len(value)/float(len(a)))*100,0))
    tissue_sim_gene[tissues[i]]=v
    tissue_sim_gene2[tissues[i]]=v2

tissue_sim_mat_gene=format_mat(tissue_sim_gene)
tissue_sim_mat_gene2=format_mat(tissue_sim_gene2)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)

#########################
#       pahes12:        #
#    export outputs     #
######################### 

gff_df.to_csv('combined.gff',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
df_tr.to_csv('transcripts_tracking.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
c=df_tr.count(axis=1)-1
c=list(c) 
plt.hist(c, 25, histtype='step', color='black')
plt.savefig('transcripts_replication_hist.png')   
print('number of multi-tissue transcripts:')
print('\n')
print(len([i for i in c if i > 1]))
print('number of single-tissue transcripts:')
print('\n')
print(len([i for i in c if i == 1]))
df_g.to_csv('genes_tracking.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
c=df_g.count(axis=1)-1
c=list(c) 
plt.hist(c, 25, histtype='step', color='black')
plt.savefig('genes_replication_hist.png')   
print('number of multi-tissue genes:')
print('\n')
print(len([i for i in c if i > 1]))
print('number of single-tissue genes:')
print('\n')
print(len([i for i in c if i == 1]))
tissue_sim_mat_tr.to_csv('transcript_level_tissue_comparisions.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_tr2.to_csv('transcript_level_tissue_comparisions_HEATMAP.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total transcript of row name (lowee diagonal) and transcript of col name (upper diagonal) 
tissue_sim_mat_gene.to_csv('gene_level_tissue_comparisions.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
tissue_sim_mat_gene2.to_csv('gene_level_tissue_comparisions_HEATMAP.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)#simmilarity % as total genes of row name (lowee diagonal) and genes of col name (upper diagonal)  
np.save('tissue_gene.npy', tissue_gene)
np.save('gene_tissue.npy', gene_tissue)
np.save('tissue_transcript.npy', tissue_tr)
np.save('transcript_tissue.npy', tr_tissue)
tss_unspliced_tr=remain_unspliced['tr_id'].to_frame()
tss_unspliced_tr.to_csv('tissue_specific_unsplised_transcripts.txt',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)

"""

"""
tissue_pair_tr={}#transcripts commom in each pair of tissues
for pair in tissue_pairs.values():
    t1=pair[0]
    t2=pair[1]
    a=tissue_tr[t1]
    b=tissue_tr[t2]
    value=list(set(a) & set(b))
    key=t1 + '---' + t2
    tissue_pair_tr[key]=value


#longest = max(len(item) for item in tissue_pair_tr.values())
tissue_pair_tr=OrderedDict(sorted(tissue_pair_tr.items(), key=lambda x: len(x[1]),reverse=True))# sort tissue_pair_tr so dat most simillar tissue pairs become firrst in dict
"""
        
"""
# program to Find maximum  
# length list in a nested list 
  
def FindMaxLength(lst): 
    maxList = max((x) for x in lst) 
    maxLength = max(len(x) for x in lst ) 
  
    return maxList, maxLength 

lst=updated_same_trs_spliced.values()
FindMaxLength(lst)
"""
                 
""""""
"""
out=[]
out2=[]
for key, value in same_trs_spliced.items():
    out.append(key)
    for i in value:
        out2.append(i)

for key, value in same_trs_spliced_with_retained.items():
    out.append(key)
    for i in value:
        out2.append(i)

for key, value in same_trs_spliced_with_retained2.items():
    out.append(key)
    for i in value:
        out2.append(i)
        
ref_spliced_replicated=list(dict.fromkeys(out))
ref_spliced_replicated=pd.DataFrame(ref_spliced_replicated).rename({0:'tr_id'},axis=1)
query_spliced_replicated=list(dict.fromkeys(out2))
query_spliced_replicated=pd.DataFrame(query_spliced_replicated).rename({0:'tr_id'},axis=1)

out=[]
out2=[]
for key, value in same_trs_unspliced.items():
    out.append(key)
    for i in value:
        out2.append(i)

ref_unspliced_replicated=list(dict.fromkeys(out))
ref_unspliced_replicated=pd.DataFrame(ref_unspliced_replicated).rename({0:'tr_id'},axis=1)
query_unspliced_replicated=list(dict.fromkeys(out2))
query_unspliced_replicated=pd.DataFrame(query_unspliced_replicated).rename({0:'tr_id'},axis=1)
"""

""" Slow method
def merge_lists(lsts):
#     this function merge intersected lists within neseted list 
    sets = [set(lst) for lst in lsts if lst]
    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return sets
#lsts=[['A','B'],['B','C','D'],['D','E','F'],['O','G']]
                   
for key, value in same_trs_spliced.items():
    value.append(key)

val=same_trs_spliced.values()
m=merge_lists(val)
"""


"""""
tr_ids={}
g_counter=0
processed=[]
for key, value in result_dict.items():
    g_counter=g_counter+1
    gid='PB.' + str(g_counter)
    tr_counter=0
    for i in value:
        processed.append(i)
        tr_counter=tr_counter+1
        tid=gid + '.' + str(tr_counter)
        tr_ids[i]=tid
        
processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
all_trs=pd.DataFrame(transcripts_dict.keys()).rename({0:'tr_id'},axis=1)
r=all_trs[~all_trs['tr_id'].isin(processed['tr_id'])]['tr_id'].unique()#

for i in r:
    g_counter=g_counter+1
    gid='PB.' + str(g_counter)
    tr_counter=1
    tid=gid + '.' + str(tr_counter)
    tr_ids[i]=tid
"""