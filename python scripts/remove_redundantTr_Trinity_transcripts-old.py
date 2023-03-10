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

import re
import pandas as pd
import time
import sys
import pybedtools
from multiprocessing import Pool
from collections import OrderedDict
import csv
from collections import Iterable
import collections   
import networkx#https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
import dask.dataframe as dd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import swifter

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
        bed_dicts[tissue]=pybedtools.BedTool.from_dataframe(bed).sort()
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
            elif len(l1)==1 and check_tr_begining_end(l1,l2,strand,'unspliced')=='TRUE':
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
            if len(l1)>1 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
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
                if len(l1)==2 and not check_transcripts(l1,l2):
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
            if len(l1)>1 and check_tr_begining_end(A,B,strand,'spliced')=='TRUE':# len(A)>0: meaning both query and ref have multiple exons
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
            elif len(l1)==2 and not check_transcripts(l1,l2):
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
            elif len(l1)==2 and not check_transcripts(l1,l2):
                if trs[0] not in same_trs_spliced:
                    same_trs_spliced[trs[0]]=[trs[1]] 
                else:
                    same_trs_spliced[trs[0]]=[trs[1]]
    return(same_trs_spliced,same_trs_spliced_with_retained,same_trs_spliced_with_retained2)

def paralle_chunk_input5(chunk):
    """ same 3' end exon but extra exon(s) on the 5' end of Iso-seq transcript or reverse in + strand
        same as above but with retained-intron exon at rna-seq or Iso-seq transcript
        transcripts with retained intron exons at both Iso-seq and RNA-seq transcripts were considered as different transcripts"""
    global group5
    global ex_coords
    global tr_strand_chr
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
            elif len(l1)==2 and not check_transcripts(l1,l2):
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
            elif len(l1)==2 and not check_transcripts(l1,l2):
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
            elif len(l1)==2 and not check_transcripts(l1,l2):
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
    first = prev = item = i.next()
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
    
def flatten_nested_list(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, basestring):
             for x in flatten_nested_list(item):
                 yield x
         else:        
             yield item
             
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
    
def get_merged_structures_parallel(chunk):  
    global mergable_same_trs_spliced
    global ex_coords
    spliced_merged_final={}
    ids_info={}
    ids_tissue={}
    for item in range(len(chunk)):
        lsts=mergable_same_trs_spliced[chunk[item][0]:chunk[item][1]]
        ids=chunk[item][0]
        count=1
        for l in lsts:
            l_tissues=[]
            intervals=[]
            ll=[]
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
            key='merge' + '.' + str(ids) + '.' + str(count)
            spliced_merged_final[key]=new_ref_ex
            ids_info[key]=l[0]
            for i in l:
                tissue=i.split('.')[0]
                if tissue not in l_tissues:
                    l_tissues.append(tissue)
            ids_tissue[key]=l_tissues
            count=count+1
    return(spliced_merged_final,ids_info,ids_tissue)  

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
    
input_file=sys.argv[1]# gff file 'mRNA-Testis_final.collapsed.gff'


tissues=[]
tissue=re.split('mRNA-|_final|.collapsed', input_file)[1]
tissues.append(tissue)
items=[input_file,input_file]
    
p = Pool()      
start = time.time()
processed_values= p.map(paralle_chunk_input_1, get_chunks(items, 1))
end = time.time()
print('total time (s)= ' + str(end-start)) 

bed_dicts={}
for i in range(len(processed_values)):
    bed_dicts.update(processed_values[i])



intersections=bed_dicts[tissue].intersect(bed_dicts[tissue],wa=True,wb=True,s=True).to_dataframe()
intersections.columns = [0, 1,2,'a',4,5,6,7,8,9,'b',11,12,13]
intersections=intersections[(intersections.a != intersections.b)]
intersections.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]

######################

intersections['info']=intersections[3] + '---' + intersections[10]

ex_coords=get_ex_coords_dict(input_file,tissue)

tr_ex_info=[]
for line in open(input_file):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        info=raw[8].split("\"")
        tr_ex_info.append([info[3].replace('PB',tissue),raw[6],raw[0]])
    
tr_ex_info_df=pd.DataFrame(tr_ex_info,columns =['tr_id','strand','chr'])
info=tr_ex_info_df.drop_duplicates()
#tr_strand=info[['tr_id','strand']]
#tr_strand=tr_strand.set_index('tr_id').to_dict()
info=dd.from_pandas(info,npartitions=-1)
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

intersections=dd.from_pandas(intersections,npartitions=36)# use npartitions=-1 to use all cpus According to StackOverflow, it is advised to partition the Dataframe in about as many partitions as cores your computer has, or a couple times that number, as each partition will run on a different thread and communication between them will become too costly if there are too many. https://towardsdatascience.com/trying-out-dask-dataframes-in-python-for-fast-data-analysis-in-parallel-aa960c18a915

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
m5=dd.from_pandas(m5,npartitions=-1)
c=m5[(m5.counts==m5.ref_count)]
candidates=c[(m5.query_count==m5.last_intersected_query_exon)]## all candidates; last exon of query and ref intersected. this is because our Iso-seq data are 3'end selected
group1=candidates[(candidates.counts==candidates.ref_count) & (candidates.ref_count==candidates.query_count)].compute()#candidates for exact same structure in all re and query exons
candidates=candidates.compute()# compute() convert dask to pandas dataframe
group2=candidates[~candidates['info'].isin(group1['info'])]# candidates for same atructure at all Iso-seq transcripts in RNA-seq transcripts plus extra exon(s) at RNA-seq transcript 5' end or vise versa (same 5' extra exon on 3' in reverse strand). This is because Iso-seq data are not 5' end selected
candidates=c[(c.query_count > c.last_intersected_query_exon)]#all ref exons intersected with query exon but query transcript has additional exons on their 3' end
group3=candidates[(candidates.first_intersected_query_exon == 1)].compute()#same 5' end
group4=candidates[(candidates.first_intersected_query_exon != 1)].compute()#different 5' end
c=m5[(m5.counts == m5.query_count)].compute()
c2=c[~c['info'].isin(group1['info'])]
c2=dd.from_pandas(c2,npartitions=-1)
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

p=pd.concat([group1['info'].to_frame(),group2['info'].to_frame(),group3['info'].to_frame(),group4['info'].to_frame(),group5['info'].to_frame(),group6['info'].to_frame(),group7['info'].to_frame()]).drop_duplicates()
intersections=intersections.compute()
intersections=intersections[intersections['info'].isin(p['info'])]
intersections_dict={}
a=intersections.swifter.apply(lambda row:aggregate(row),axis=1)#swifter tutorial https://gdcoder.com/speed-up-pandas-apply-function-using-dask-or-swifter-tutorial/  
del(intersections)
#np.save('intersections_dict.npy', intersections_dict)
#intersections_dict = np.load('intersections_dict.npy',allow_pickle='TRUE').item()


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
                same_trs_spliced.update(outputs[j])
    elif i==1:
        outputs = [result[i] for result in processed_values]
        if len(outputs)>0:
            for k in range(len(outputs)):
                same_trs_unspliced.update(outputs[k])


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
same_trs_spliced=OrderedDict(sorted(same_trs_spliced.items(), key=lambda x: len(x[1]),reverse=True))
mergable_same_trs_spliced=[]
count=0
for key, value in same_trs_spliced.items():
    count=count+1
    vals=[]
    info=recursive_check(value)
    if len(info[0])==0 and len(info[1])>0:
        for i in info[1]:
            i.append(key)
            mergable_same_trs_spliced.append(i)      
    elif len(info[0])==0 and len(info[1])==0:
        v=value[:]
        v.append(key)
        if len(v)==2:
            mergable_same_trs_spliced.append(v)
        else:
            r=len(ex_coords[key])
            v2=[]
            for j in value:
                if  len(ex_coords[j])<=r:
                    v2.append(j)
                else:
                    mergable_same_trs_spliced.append([j,key])
            v2.append(key)
            if len(v2)>1:
                mergable_same_trs_spliced.append(v2)
    elif len(info[0])>0 and len(info[1])==0:
        info[0].append(key)
        mergable_same_trs_spliced.append(info[0])

""" The following section doesn't need to be included in combine gff files program
    and just required when the gola is to remove within tissue redundancy"""
mergable_same_trs_spliced=list(map(sorted, mergable_same_trs_spliced))
fset = set(frozenset(x) for x in mergable_same_trs_spliced)
mergable_same_trs_spliced = [list(x) for x in fset]
del(fset)
""""""
#global mergable_same_trs_spliced
#global ex_coords
isize=100#size of each item (number of itersections rows passed to each cpu) in items
p = Pool()
items=get_items_fromLists(mergable_same_trs_spliced) 
start = time.time()
processed_values= p.map(get_merged_structures_parallel, get_chunks(items, 10))
end = time.time()
print('total time (s)= ' + str(end-start)) 
transcripts_dict={}
spliced_merged_id_info={}
tr_tissues={}#detected tissue per transcript
for i in range(len(processed_values)):
    transcripts_dict.update(processed_values[i][0])
    spliced_merged_id_info.update(processed_values[i][1])
    tr_tissues.update(processed_values[i][2])

processed=list(flatten_nested_list(mergable_same_trs_spliced))          
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
        tis=tr.split('.')[0]
        if tis not in t:
            t.append(tis)
    new_ref_ex=merge_intervals(intervals)
    transcripts_dict[key]=new_ref_ex
    tr_tissues[key]=t
    

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

for i in range(remain_unspliced.shape[0]):
    tr=remain_unspliced.iloc[i,0]
    t=tr.split('.')[0]
    v=ex_coords[tr]
    transcripts_dict[tr]=v
    tr_tissues[tr]=t

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

""" this section is aimed to sort final transcripts
    based on their location to get sorted gene ids 
    THIS CAN BE SAFELY IGNORED """ 

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


#########################
#       pahes8:        #
#    export outputs     #
######################### 

gff_df.to_csv('mRNA-' + tissue + '_final.collapsed' + '_redundant_removed'+ '.gff',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  
