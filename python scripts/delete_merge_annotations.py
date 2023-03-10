#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 23:47:48 2020

@author: beiki
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
    global intersection
    out1=[]# ref masked exons
    a=intersection.loc[intersection['info']==name]
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
        
def paralle_chunk_input1(chunk):
    """ exact same structures """
    global group1
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_unspliced={}
    for item in range(len(chunk)):
        candidates=group1.iloc[chunk[item][0]:chunk[item][1],0:]
        for i in range(candidates.shape[0]):
            trs=info[candidates.iloc[i,0]]
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}# same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group2.iloc[chunk[item][0]:chunk[item][1],0:]    
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
                if len(l1)==2:
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}# same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group3.iloc[chunk[item][0]:chunk[item][1],0:]            
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
            elif len(l1)==2:
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}# same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group4.iloc[chunk[item][0]:chunk[item][1],0:]  
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
                        A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                        B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                        if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!
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
                        A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                        B=[[0,list2_2[0][1]],[list2_2[-1][0],5]] 
                        if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!!
                            if trs[0] not in same_trs_spliced_with_retained2:
                                same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}#masked_info[2] show the retained intron exon in refrence trancript
                            else:
                                same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})               
                elif len(l1)==2:
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group5.iloc[chunk[item][0]:chunk[item][1],0:] 
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
                        A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                        B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                        if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):#!!!!!!!!!!!
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
                        A=[[0,list1_2[0][1]],[list1_2[-1][0],5]]
                        B=[[0,list2_2[0][1]],[list2_2[-1][0],5]]
                        if len(list1_2)==len(list2_2) and check_tr_begining_end(A,B,strand,'spliced')=='TRUE' and compare_trs(list1_2[1:-1],list2_2[1:-1]) == (len(list1_2)-2):
                            if trs[0] not in same_trs_spliced_with_retained2:
                                same_trs_spliced_with_retained2[trs[0]]={trs[1]:masked_info[2]}
                            else:
                                same_trs_spliced_with_retained2[trs[0]].update({trs[1]:masked_info[2]})                         
            elif len(l1)==2:
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group6.iloc[chunk[item][0]:chunk[item][1],0:] 
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            trs=info[candidates.iloc[i,0]]
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
            elif len(l1)==2:
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
    global ref_ex_coords
    global query_ex_coords
    global info
    same_trs_spliced={}
    same_trs_spliced_with_retained={}#same transcripts but with retained intron exon(s) at RNA-seq transcript
    same_trs_spliced_with_retained2={}#same transcripts but with retained intron exon(s) at Iso-seq transcript
    for item in range(len(chunk)):
        candidates=group7.iloc[chunk[item][0]:chunk[item][1],0:]
        for i in range(candidates.shape[0]):
            name=candidates.iloc[i,0]
            trs=info[name]
            first_q=candidates.iloc[i,6]-1#firs exon in query that intersected with ref; -1 as it exon numbers atarted from 1
            last_q=candidates.iloc[i,7]-1
            first_r=candidates.iloc[i,8]-1
            last_r=candidates.iloc[i,9]-1
            strand=ref_bed.loc[ref_bed[3]==trs[0]].iloc[0,5]
            l1=ref_ex_coords[trs[0]]
            l2=query_ex_coords[trs[1]]
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
            elif len(l1)==2:
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

def get_exons(key):
    global ref_ex_coords
    global query_ex_coords
    global tissue1
    global tissue2
    if key.find(tissue1)==0:
        ex=ref_ex_coords[key]
    else:
        ex=query_ex_coords[key]
    return(ex)

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
#################
#       RUN     #
#################    

#####################################
#  phase1: get simmilar transcripts #
#####################################  

ref=sys.argv[1]#"mRNA-Testis_final.collapsed.gff" 
query=sys.argv[2]#"mRNA-Thymus_final.collapsed.gff"
tissue1=re.split('mRNA-|_final|.collapsed', ref)[1]
tissue2=re.split('mRNA-|_final|.collapsed', query)[1]

fuzzySplice=20

ref_tr_info= get_tr_num_ex(ref).rename({'tr_id':'ref_tr','Count':'ref_count'},axis=1)          
ref_tr_info=ref_tr_info.replace(to_replace='PB',value=tissue1,regex=True)
query_tr_info= get_tr_num_ex(query).rename({'tr_id':'query_tr','Count':'query_count'},axis=1)  
query_tr_info=query_tr_info.replace(to_replace='PB',value=tissue2,regex=True)


ref_bed=get_ex_coords(ref)
ref_bed=ref_bed.replace(to_replace='PB',value=tissue1,regex=True)
query_bed=get_ex_coords(query)
query_bed=query_bed.replace(to_replace='PB',value=tissue2,regex=True)
query_bed2=pybedtools.BedTool.from_dataframe(query_bed)
ref_bed2=pybedtools.BedTool.from_dataframe(ref_bed)
intersection=ref_bed2.intersect(query_bed2,wa=True,wb=True,s=True).to_dataframe() 
# on my linux machine:
#intersection.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]
intersection['info']=intersection[3] + '---' + intersection[10]
intersection=intersection[intersection[[3,10]].nunique(axis=1) != 1]
info=intersection.pivot_table(index=['info'], aggfunc='size').reset_index().rename({0:'counts'},axis=1)
info['ref_tr'],info['query_tr']=info['info'].str.split('---',1).str
m1=info.merge(query_tr_info)
m=m1.merge(ref_tr_info)
info1=intersection.loc[intersection.groupby('info')[13].idxmin(), :][[13,'info']].rename({13:'first_intersected_query_exon'},axis=1)
info2=intersection.loc[intersection.groupby('info')[13].idxmax(), :][[13,'info']].rename({13:'last_intersected_query_exon'},axis=1)
info3=intersection.loc[intersection.groupby('info')[6].idxmin(), :][[6,'info']].rename({6:'first_intersected_ref_exon'},axis=1)
info4=intersection.loc[intersection.groupby('info')[6].idxmax(), :][[6,'info']].rename({6:'last_intersected_ref_exon'},axis=1)
m2=m.merge(info1)
m3=m2.merge(info2)
m4=m3.merge(info3)
m5=m4.merge(info4)
c=m5[m5[['counts','ref_count']].nunique(axis=1) == 1]
candidates=c[c[['query_count','last_intersected_query_exon']].nunique(axis=1) == 1]## all candidates; last exon of query and ref intersected. this is because our Iso-seq data are 3'end selected
group1=candidates[candidates[['counts','ref_count','query_count']].nunique(axis=1) == 1]#candidates for exact same structure in all re and query exons
group2=candidates[~candidates['info'].isin(group1['info'])]# candidates for same atructure at all Iso-seq transcripts in RNA-seq transcripts plus extra exon(s) at RNA-seq transcript 5' end or vise versa (same 5' extra exon on 3' in reverse strand). This is because Iso-seq data are not 5' end selected
candidates=c.query('query_count > last_intersected_query_exon')#all ref exons intersected with query exon but query transcript has additional exons on their 3' end
group3=candidates.loc[candidates['first_intersected_query_exon']==1]  #same 5' end
group4=candidates.loc[candidates['first_intersected_query_exon']!=1]  #different 5' end
c=m5[m5[['counts','query_count']].nunique(axis=1) == 1]
c2=c[~c['info'].isin(group1['info'])]
a=c2.query('ref_count > query_count').copy()
a['flag']=a['last_intersected_query_exon'].copy()+(a['ref_count'].copy()-a['query_count'].copy())## flag show the query last exon number compatible with ref exons
b=a.loc[a['query_count']>1]
group5=b[b[['ref_count','flag']].nunique(axis=1) == 1]##same 3' end but different 5' end
group6=b[b['first_intersected_ref_exon']==1]##same 5' end but different 3' end
o=pd.concat([group5['info'].to_frame(),group6['info'].to_frame()])
c=b[~b['info'].isin(o['info'])]
group7=c.query('ref_count > last_intersected_ref_exon')

candidates=group1
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]

ref_ex_coords=get_ex_coords_dict(ref,tissue1)
query_ex_coords=get_ex_coords_dict(query,tissue2)

            
p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
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

candidates=group2
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input2, get_chunks(items, 30))
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

candidates=group3
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input3, get_chunks(items, 30))
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

candidates=group4
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input4, get_chunks(items, 30))
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

candidates=group5
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input5, get_chunks(items, 30))
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

candidates=group6
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input6, get_chunks(items, 30))
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


candidates=group7
info={}
for i in range((candidates).shape[0]):
    a=candidates.iloc[i,0:]
    info[a[0]]=[a[2],a[3]]


p = Pool() 
items=[]
for i in range(0,candidates.shape[0],10):
    if i<candidates.shape[0]:
        items.append([i,i+10])
    else:
        items.append([i,i+11])
        
start = time.time()
processed_values= p.map(paralle_chunk_input7, get_chunks(items, 30))
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
#same_trs_spliced['PB.12366.11']

tr_replication_dict={}
for key, value in same_trs_spliced.items():
    v=value[:]
    if key in same_trs_spliced_with_retained.keys():
        d=same_trs_spliced_with_retained[key]
        for i in d.keys():
            v.append(i)
    if key in same_trs_spliced_with_retained2.keys():
        d=same_trs_spliced_with_retained2[key]
        for i in d.keys():
            v.append(i)
    tr_replication_dict[key]=v

for key, value in same_trs_spliced_with_retained.items():
    if key not in tr_replication_dict.keys():
        v=value.keys()
        if key in same_trs_spliced_with_retained2.keys():
            d=same_trs_spliced_with_retained2[key]
            for i in d.keys():
                v.append(i)
        tr_replication_dict[key]=v

for key, value in same_trs_spliced_with_retained2.items():
    if key not in tr_replication_dict.keys():
        tr_replication_dict[key]=value.keys()

for key, value in same_trs_unspliced.items():
    tr_replication_dict[key]=value[:]

tr_replication_info=pd.DataFrame.from_dict(tr_replication_dict,orient='index').reset_index().rename({'index':'tr_id'},axis=1)
tr_replication_info.to_csv(r'transcript_replication_info.txt',index = None, header=False,sep="\t")

############################################
#  phase2: get transcripts representatives #
#             (spliced)                    #
############################################  

"""merge intersected lists within neseted list"""
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
        info[i]=get_exons(i)
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))#sort dictionari by length of values (list) https://stackoverflow.com/questions/50863093/how-to-sort-a-python-dictionary-based-on-the-length-of-the-list-of-the-values
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    updated_same_trs_spliced[new_key]=new_value

                        
updated_same_trs_spliced_with_retained={}
for key, value in same_trs_spliced_with_retained.items():
    info={}
    info[key]=get_exons(key)
    for key2 in value.keys():
        info[key2]=get_exons(key2)
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    if new_key not in updated_same_trs_spliced_with_retained.keys():
        updated_same_trs_spliced_with_retained[new_key]=new_value
    elif new_key in updated_same_trs_spliced_with_retained.keys():
        v=updated_same_trs_spliced_with_retained[new_key][:]
        for i in new_value:
            if i not in v:
                v.append(i)
        updated_same_trs_spliced_with_retained[new_key]=v

updated_same_trs_spliced_with_retained2={}
for key, value in same_trs_spliced_with_retained2.items():
    info={}
    info[key]=get_exons(key)
    for key2 in value.keys():
        info[key2]=get_exons(key2)
    info=OrderedDict(sorted(info.items(), key=lambda x: len(x[1]),reverse=True))
    k=info.keys()
    new_key=k[0]
    new_value=k[1:]
    if new_key not in updated_same_trs_spliced_with_retained2.keys():
        updated_same_trs_spliced_with_retained2[new_key]=new_value
    elif new_key in updated_same_trs_spliced_with_retained2.keys():
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
    ref_ex=get_exons(key)
    first=ref_ex[0]
    last=ref_ex[-1]
    for exon in ref_ex:
        intervals.append(exon)
    for tr in value:
        info=get_exons(tr)
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
    if key.find(tissue1)==0:
        info_dict=same_trs_spliced_with_retained[key]
        ref_ex=get_exons(key)
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        for key2,value2 in info_dict.items():
            query_ex=get_exons(key2)
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
        ref_ex=get_exons(key)
        ref_ex=adjust_ex(ref_ex,[info_main[key][0]])
#        ref_ex.pop(info_main[key][0])
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        main_ex=get_exons(main_key)
        for exon in main_ex:
            intervals.append(exon)
        for key2,value2 in info_main.items():
            if key2 !=key:
                query_ex=get_exons(key2)
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
    if key.find(tissue1)==0:
        info_dict=same_trs_spliced_with_retained2[key]
        ref_ex=get_exons(key)
        out=[]
        for value2 in info_dict.values():
            out.append(value2[0])
        ref_ex=adjust_ex(ref_ex,out)
        first=ref_ex[0]
        last=ref_ex[-1]
        for exon in ref_ex:
            intervals.append(exon)
        for tr in info_dict.keys():
            info=get_exons(tr)
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
        ref_ex=get_exons(key)
        main_ex=get_exons(main_key)
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
                q_ex=get_exons(tr)
                for exon in q_ex:
                    intervals.append(exon)
        new_ref_ex=merge_intervals(intervals)
        spliced_with_retained_merged_structures2[main_key]=new_ref_ex
            
""" merge merged trascripts """
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

info=get_tr_splicing_info(ref,tissue1)
ref_spliced_df=info[0]
ref_unspliced_df=info[1]
chr_strand_info=info[2]
info=get_tr_splicing_info(query,tissue2)
query_spliced_df=info[0]
query_unspliced_df=info[1]
chr_strand_info.update(info[2])

spliced=pd.concat([ref_spliced_df,query_spliced_df])
unspliced=pd.concat([ref_unspliced_df,query_unspliced_df])

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
    exons=get_exons(key)
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
    ref_ex=get_exons(key)
    for exon in ref_ex:
        intervals.append(exon)
    for tr in value:
        processed.append(tr)
        info=get_exons(tr)
        for exon in info:
            intervals.append(exon)
    new_ref_ex=merge_intervals(intervals)
    merged_structures2[key]=new_ref_ex

transcripts_dict.update(merged_structures2)## mixure of merged, unmerged spliced transcripts and merge unspliced transcripts

processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

""" I YOU WANT TO INCLUDE TISSUE SPECIFIC UNSPLICED TRANSCRIPTS 
processed=pd.DataFrame(processed).rename({0:'tr_id'},axis=1)
remain_unspliced=unspliced[~unspliced['tr_id'].isin(processed['tr_id'])]#in both ref and query

unspliced_umerged_final={}# unspliced transcripts in both ref and query that have unique structure compared to other transcriptome.
for i in range(remain_unspliced.shape[0]):
    key=remain_unspliced.iloc[i,0]
    exons=get_exons(key)
    unspliced_umerged_final[key]=exons
"""
""""""""""""    

##############################
#  phase7: get gene ids      # 
#          and format trids  #
############################## 

results={}
for key, value in transcripts_dict.items():
    info=chr_strand_info[key]
    results[key]=[info[0],value[0][0],value[-1][1],key,'99',info[1],key]

bed=pd.DataFrame.from_dict(results,orient='index')##if genes defined as group of transcripts with exoni c overlap, this shoud include exons. otherwise transcripts
bed=pybedtools.BedTool.from_dataframe(bed)
intersection=bed.intersect(bed,wa=True,wb=True,s=True).to_dataframe()
#intersection.columns = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13]
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

my_dict=df.groupby('A')['B'].apply(lambda g: g.values.tolist()).to_dict()


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
    based on their location to get sorted gene ids """ 

d={}
for key in result_dict.keys():
    info=chr_strand_info[key]
    value=transcripts_dict[key]
    chrom=info[0].split('chr')[1]
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
#  phase8: create & export   #
#       .gff file            #
############################## 

gff_list=[]
for key, value in transcripts_dict.items():
    info=chr_strand_info[key]
    tid=tr_ids[key]
    a=tid.split('.')[0:2]
    gid=a[0]+'.'+a[1]
    transcript=[info[0],'CattleFAANG','transcript', value[0][0],value[-1][1],'.',info[1],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
    gff_list.append(transcript)
    for i in value:
        exon=[info[0],'CattleFAANG','exon',i[0],i[1],'.',info[1],'.','gene_id "'+gid+'";' ' transcript_id "'+tid+'";']
        gff_list.append(exon)
    
    
gff_df=pd.DataFrame(gff_list)    
gff_df.to_csv('merged.gff',index = None, header=False,sep="\t",quoting=csv.QUOTE_NONE)  

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