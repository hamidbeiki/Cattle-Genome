#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:46:30 2020

@author: beiki
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:02:10 2020

@author: beiki
"""
import collections
from collections import OrderedDict
import pandas as pd
import re
from multiprocessing import Pool
import time
import sys


input=sys.argv[1]# your gff file, like input="Cerebral_Cortex_fmlrc_proovread_cupcake_corrected_reads.collapsed.gff" or "gffcmp.annotated.gtf"
if input.find('fmlrc_proovread_cupcake') !=-1:
    tissue=input.split("_fmlrc")[0]
    out_prefix="fmlrc_proovread_cupcake_corrected"
elif input.find('final.collapsed') !=-1:
    tissue=input.split("_final")[0]
    out_prefix="final"
elif input.find('_tmp') !=-1:
    tissue=input.split("_tmp")[0]
    out_prefix="valid_splice"
elif input.find('gtf') !=-1:
    tissue='irrelevant'
    
    

gff=open(input, 'r').read()

if input.find('gtf') !=-1:
    input_type="gtf"
else:
    input_type="gff"
############## functions ##############
def nested_dict():
    return collections.defaultdict(nested_dict)

def get_premature_mRNA_events(df):
    """ premature mRNA means retained intron transcripts which except their RI exons, the rest of exons are exactly same as their refrence transcripts"""
    tr_numbers=len(df.index)
    premature_events=pd.DataFrame([])
    strand=re.split('\(|\)',df.iloc[0,:][0])[1]
    for i in range(tr_numbers-1):
        for j in range(i+1,tr_numbers):
            info=len(ex_coords[df.iloc[i,1]])+len(ex_coords[df.iloc[j,1]])## to prevent uspliced-unspliced pairs
            if info>2:
                res=check_for_premature_mRNAs(df.iloc[i,1],df.iloc[j,1],strand)
                if res!="clean":
                    premature_transcript=res[0]
                    ref_transcript=res[1]
                    premature_events=premature_events.append(pd.DataFrame({'premature_transcript': premature_transcript,'ref_transcript': ref_transcript}, index=[0]), ignore_index=True)
    return(premature_events)

def check_tr_begining_end(A,B,strand):
    fuzzyTR_3=100
    fuzzyTR_5=1000
    if strand=="+" and abs(int(A[0][0])-int(B[0][0]))<=fuzzyTR_5 and abs(int(A[-1][1])-int(B[-1][1]))<=fuzzyTR_3:
        return("TRUE")
    elif strand=="-" and abs(int(A[0][1])-int(B[0][1]))<=fuzzyTR_5 and abs(int(A[-1][0])-int(B[-1][0]))<=fuzzyTR_3:
        return("TRUE")
    else:
        return("FALSE")   

def get_tr_exon_list(tr):
    global ex_coords
    my_dict=ex_coords.get(tr)
    sorted_dict=OrderedDict(sorted(my_dict.items(), key=lambda t: t[0]))
    my_list=[]
    for key, value in sorted_dict.iteritems():
        formatted=value.replace("-",",")
        a=re.split('\(|\:',formatted)[1]
        b=list(a.split(","))
        my_list.append(b)
    return(my_list)

def check_for_premature_mRNAs(trA,trB,strand):
    global A
    global B
    A=get_tr_exon_list(trA)
    B=get_tr_exon_list(trB)
    if check_tr_begining_end(A,B,strand)=="TRUE":
        res=premature_mRNA_test(A,B,strand)
        if res!="clean" and res[0]==0:
            return(trA,trB)
        elif res!="clean" and res[0]==1:
            return(trB,trA)
        else:
            return("clean")
    else:
        return("clean")
        
def check_intersection(exonA,exonB):
    N=2
    intervals=[]
    intervals.append(exonA)
    intervals.append(exonB)
  
    # First interval 
    l = int(intervals[0][0]) 
    r = int(intervals[0][1]) 
  
    # Check rest of the intervals  
    # and find the intersection 
    for i in range(1,N): 
  
        # If no intersection exists 
        if (int(intervals[i][0]) > r or int(intervals[i][1]) < l): 
            return("notintersected") 
  
        # Else update the intersection 
        else: 
            l = max(l, int(intervals[i][0])) 
            r = min(r, int(intervals[i][1]))
            return("intersected") 
#def check_intersection(query,ref):
#    """tjis function won always work Converting a list to a set sometimes changes element order (https://stackoverflow.com/questions/9792664/converting-a-list-to-a-set-changes-element-order)
#    """
#    xs=set(query)
#    xy=set(ref)
#    if len(xs.intersection(xy))>0:
#        return("intersected")
#    else:
#        return("non_intersected")

def check_exons(strand,request_type,query_tr,ex_num,ref_tr,ref_num,start,end):
    fuzzyTR_3=100
    fuzzyTR_5=1000
    fuzzySplice=5
    ref_tr_exons=len(ref_tr)
    query_tr_exons=len(query_tr)
    if strand=="+":
        if (ref_tr_exons-1)==ref_num and (query_tr_exons-1)==ex_num:
            fuzzyTR_3=100
            fuzzyTR_5=fuzzySplice
        elif ex_num==0 and ref_num==0:
            fuzzyTR_5=1000
            fuzzyTR_3=fuzzySplice
        else:
            fuzzyTR_5=fuzzySplice
            fuzzyTR_3=fuzzySplice
        if request_type=="start" and abs(int(query_tr[ex_num][0])-int(start))<=fuzzyTR_5:
            return("True")
        elif request_type=="end" and int(query_tr[ex_num][1])<=(int(end)+fuzzyTR_3): ## because I'm intrested in exons with equal or smaller 3' end
            return("True")
        elif request_type=="equality" and abs(int(query_tr[ex_num][0])-int(start))<=fuzzyTR_5 and abs(int(query_tr[ex_num][1])-int(end))<=fuzzyTR_3:
            return("True")
        else:
            return("False")
    elif strand=="-":
        if ex_num==0 and ref_num==0:
            fuzzyTR_3=100
            fuzzyTR_5=fuzzySplice
        elif (ref_tr_exons-1)==ref_num and (query_tr_exons-1)==ex_num:
            fuzzyTR_5=1000
            fuzzyTR_3=fuzzySplice
        else:
            fuzzyTR_5=fuzzySplice
            fuzzyTR_3=fuzzySplice             
        if strand=="-" and request_type=="start" and abs(int(query_tr[ex_num][1])-int(start))<=fuzzyTR_5:
            return("True")
        elif strand=="-" and request_type=="end" and int(query_tr[ex_num][0])>=(int(end)-fuzzyTR_3):## because I'm intrested in exons with equal or smaller 3' end
            return("True")
        elif strand=="-" and request_type=="equality" and abs(int(query_tr[ex_num][1])-int(start))<=fuzzyTR_5 and abs(int(query_tr[ex_num][0])-int(end))<=fuzzyTR_3:
            return("True")
        else:
            return("False")


def judge(lenA,lenB,masked_exons,shared_exons,multi_covered_exons):
    if multi_covered_exons["A"]>0 and (shared_exons["A"]+masked_exons["A"]+multi_covered_exons["A"])==lenA and (shared_exons["A"]+masked_exons["B"]+multi_covered_exons["B"])==lenB:
        return(0,1)#code "0" means traA is premature and it's reference is trB (code "1")
    elif multi_covered_exons["B"]>0 and (shared_exons["B"]+masked_exons["B"]+multi_covered_exons["B"])==lenB and (shared_exons["B"]+masked_exons["A"]+multi_covered_exons["A"])==lenA:
        return(1,0)
    else:
        return("clean")
       
def premature_mRNA_test(A,B,strand):
    s="AB"
    listsA=list(s)
    s="BA"
    listsB=list(s)
    masked_exons={}
    shared_exons={}
    multi_covered_exons={}
    if strand=="+" and check_intersection(A[0],B[0])=="intersected" and check_intersection(A[-1],B[-1])=="intersected":
        for l in range(len(listsA)):
            query_tr=globals()[listsA[l]]
            ref_tr=globals()[listsB[l]]
            shared_exons[listsA[l]]=0
            masked_exons[listsA[l]]=0
            multi_covered_exons[listsB[l]]=0
            for exon in range(len(ref_tr)):
                start_exon=ref_tr[exon]
                start=start_exon[0]
                end=start_exon[1]
                same_starts=[]
                covered_ends=[]
                for j in range(len(query_tr)):
                    if check_intersection(query_tr[j],start_exon)=="intersected":
                        if check_exons(strand,'start',query_tr,j,ref_tr,exon,start,end)=="True":
                            same_starts.append(j)
                        if check_exons(strand,'equality',query_tr,j,ref_tr,exon,start,end)=="True":
                            shared_exons[listsA[l]]=shared_exons[listsA[l]]+1
                [ covered_ends.append(j) for j in range(len(query_tr)) if check_intersection(query_tr[j],start_exon)=="intersected" and check_exons(strand,'end',query_tr,j,ref_tr,exon,start,end)=="True" ]
                same_starts_2=[]
                covered_ends_2=[]
                [ same_starts_2.append(same_starts[x])  for x in range(len(same_starts)) if query_tr[same_starts[x]]!=ref_tr[exon] ]
                [ covered_ends_2.append(covered_ends[x])  for x in range(len(covered_ends)) if query_tr[covered_ends[x]]!=ref_tr[exon] ]
                if  len(covered_ends_2)>1 and len(same_starts_2)>0:
                    x=set(covered_ends_2)
                    if len(x.intersection(same_starts))>0:
                        masked_exons[listsA[l]]=masked_exons[listsA[l]]+len(covered_ends_2)#this shows the number of exons in this transcript that masked by multicover exons in other transcript
                        multi_covered_exons[listsB[l]]=multi_covered_exons[listsB[l]]+1
        decision=judge(len(A),len(B),masked_exons,shared_exons,multi_covered_exons)
        return(decision)
    elif strand=="-" and check_intersection(A[0],B[0])=="intersected" and check_intersection(A[-1],B[-1])=="intersected":
        for l in range(len(listsA)):
            query_tr=globals()[listsA[l]]
            ref_tr=globals()[listsB[l]]
            shared_exons[listsA[l]]=0
            masked_exons[listsA[l]]=0
            multi_covered_exons[listsB[l]]=0
            for exon in reversed(range(len(ref_tr))):
                start_exon=ref_tr[exon]
                start=start_exon[1]
                end=start_exon[0]
                same_starts=[]
                covered_ends=[]
                for j in range(len(query_tr)):
                    if check_intersection(query_tr[j],start_exon)=="intersected":
                        if check_exons(strand,'start',query_tr,j,ref_tr,exon,start,end)=="True":
                            same_starts.append(j)
                        if check_exons(strand,'equality',query_tr,j,ref_tr,exon,start,end)=="True":
                            shared_exons[listsA[l]]=shared_exons[listsA[l]]+1               
                [ covered_ends.append(j) for j in reversed(range(len(query_tr))) if check_intersection(query_tr[j],start_exon)=="intersected" and check_exons(strand,'end',query_tr,j,ref_tr,exon,start,end)=="True"  ]
                same_starts_2=[]
                covered_ends_2=[]
                [ same_starts_2.append(same_starts[x])  for x in range(len(same_starts)) if query_tr[same_starts[x]]!=ref_tr[exon] ]
                [ covered_ends_2.append(covered_ends[x])  for x in range(len(covered_ends)) if query_tr[covered_ends[x]]!=ref_tr[exon] ]
                if  len(covered_ends_2)>1 and len(same_starts_2)>0:
                    x=set(covered_ends_2)
                    if x.intersection(same_starts):
                        masked_exons[listsA[l]]=masked_exons[listsA[l]]+len(covered_ends_2)#this shows the number of exons in this transcript that masked by multicover exons in other transcript
                        multi_covered_exons[listsB[l]]=multi_covered_exons[listsB[l]]+1
        decision=judge(len(A),len(B),masked_exons,shared_exons,multi_covered_exons)    
        return(decision)
    else:
        return("clean")

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_df(input):
    output=pd.DataFrame([])
    output=pd.DataFrame.from_dict(input,orient='index')[0].to_frame().reset_index().rename({'index':'feature',0:'coordinate'},axis=1)
    output=output[['coordinate','feature']]
    return(output)
        
def paralle_chunk_input(chunk):
    output=pd.DataFrame([])
    for geneset in range(len(chunk)):
        df=get_df(chunk[geneset][1])
        premature_events=get_premature_mRNA_events(df)
        if len(premature_events.index)>0:
            output=output.append(premature_events)
    return(output)

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
    
############## run ##############

ex_coords = collections.defaultdict(dict)
if input_type=="gff":
    for line in open(input):
        raw=line.strip().split("\t")
        if raw[2]=='exon':
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            exid=raw[2]+raw[3]+tid
            ex_coords[tid][exid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6]) 
elif input_type=="gtf":
    for line in open(input):
        raw=line.strip().split("\t")
        if raw[2]=='exon':
            tid = raw[-1].split('; ')[0].split()[1][1:-1]
            exid=raw[2]+raw[3]+tid
            ex_coords[tid][exid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

gid="PB.1"
count=0
gene_tr={}
if input_type=="gff":
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
elif input_type=="gtf":
    for line in open(input):
        raw=line.strip().split("\t")
        if raw[2] == 'transcript':
            tid = raw[-1].split('; ')[0].split()[1][1:-1]
            gid_new = raw[-1].split('; ')[1].split()[1][1:-1]
            if gid==gid_new:
                count=count+1
                gene_tr[gid]=count
            if gid!=gid_new:
                count=0
                count=count+1
                gid=gid_new
                gene_tr[gid]=count    
                            
            
tr_coords = nested_dict()
if input_type=="gff":
    for line in open(input):
        raw=line.strip().split("\t")
        if raw[2] == 'transcript':
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            gid = raw[-1].split('; ')[0].split()[1][1:-1]
            if len(ex_coords[tid])>=1 and gene_tr[gid]>1: ## this export spliced/uspliced transcripts related to multit-ranscript genes  
                tr_coords[gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])
elif input_type=="gtf":
    for line in open(input):
        raw=line.strip().split("\t")
        if raw[2] == 'transcript':
            tid = raw[-1].split('; ')[0].split()[1][1:-1]
            gid = raw[-1].split('; ')[1].split()[1][1:-1]
            if len(ex_coords[tid])>=1 and gene_tr[gid]>1: ## this export spliced/uspliced transcripts related to multit-ranscript genes  
                tr_coords[gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6]) 

""" to get gene with highest number ot ranscripts:
import operator
max(gene_tr.iteritems(), key=operator.itemgetter(1))[0]
gene_tr['PB.2146']   
"""
     
    
"""efficient multi-threading.copleted in 4 minutes on 8 cpus!"""

p = Pool() 
optimum_number_of_genes=get_optimum_number_of_genes_per_chunk(len(tr_coords),3)   
items = list(tr_coords.items())
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets (i.e. uptimum_number_of_genes)in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, optimum_number_of_genes))
end = time.time()
print('total time (s)= ' + str(end-start)) 
premature_mRNA_events=pd.concat(processed_values)
a=premature_mRNA_events["premature_transcript"].to_frame()
b=premature_mRNA_events["ref_transcript"].to_frame()
premature_mRNAs=a[~a.premature_transcript.isin(b.ref_transcript)].drop_duplicates()

if tissue.find('irrelevant')==-1:
    premature_mRNA_events.to_csv(tissue + '_' + out_prefix + '_' + 'premature_mRNA_events',index = None, header=True,sep="\t")
    premature_mRNAs.to_csv(tissue + '_' + out_prefix + '_' + 'premature_mRNAs',index = None, header=True,sep="\t")
elif tissue.find('irrelevant')!=-1:
    premature_mRNA_events.to_csv('premature_mRNA_events',index = None, header=True,sep="\t")
    premature_mRNAs.to_csv('premature_mRNAs',index = None, header=True,sep="\t")   


""" to test speed of program on a gene with highest number of transcripts Cerebral_Cortex (885 transcripts!)
df=get_df(tr_coords.get("PB.2146"))
start=time.time()
premature_mRNAs=get_premature_mRNA_events(df)
end=time.time()
end-start
""" 

""" to test the accuracy of program(mainly premature_mRNA_test function)
A=[['5', '10'], ['15', '40']]
B=[['5', '10'], ['15', '20'], ['25', '30'], ['35', '40']]

A=[['5', '10'], ['15', '20'], ['25', '45'], ['50', '55'],['60', '65']]
B=[['5', '10'], ['15', '20'], ['25', '30'], ['35', '45'],['50', '65']]

## for reverse strand check the following pairs of transcripts: 
PB.16938.42	PB.16938.41## correct
PB.16938.45	PB.16938.42## incorrect
"""
  