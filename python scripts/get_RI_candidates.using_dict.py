# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import collections
import pandas as pd
import re
from multiprocessing import Pool
import time
import pybedtools

input="Sscrofa.brain.error-corrected.untrimmed.uniq.gDNAcleaned.sorted.collapsed.gff"

def nested_dict():
    return collections.defaultdict(nested_dict)


def print_dict(dictionary, ident = '', braces=1):
    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print '%s%s%s%s' %(ident,braces*'[',key,braces*']')
            print_dict(value, ident+'  ', braces+1)
        else:
            print ident+'%s = %s' %(key, value)



tr_coords = nested_dict()
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2] == 'transcript':
        tid = raw[-1].split('; ')[1].split()[1][1:-2]
        gid = raw[-1].split('; ')[0].split()[1][1:-1]
        tr_coords[gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

print_dict(tr_coords)

ex_coords = collections.defaultdict(dict)
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2]=='exon':
        tid = raw[-1].split('; ')[1].split()[1][1:-2]
        exid=raw[2]+raw[3]+tid
        ex_coords[tid][exid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])            


def get_df(input):
    output=pd.DataFrame([])
    for key,value in input.iteritems():
        output=output.append(pd.DataFrame({'feature': key, 'coordinate': value}, index=[0]), ignore_index=True)
    return(output)

def get_tr_candidates(df):
    tr_numbers=len(df.index)
    candidates=pd.DataFrame([])
    strand=re.split('\:|\-|\(|\)',df.iloc[0,:][0])[3]
    for i in range(tr_numbers-1):
        for j in range(i+1,tr_numbers):
            if compare_tr(df.iloc[i,:],df.iloc[j,:],strand)=='TRUE':
                candidates=candidates.append(pd.DataFrame({'transcriptA': df.iloc[i,1], 'transcriptB': df.iloc[j,1]}, index=[0]), ignore_index=True)
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
    a=df['coordinate'].str.split('\:|\-|\(|\)')
    for key,value in a.iteritems():
        output=output.append(pd.DataFrame({'chr': value[0], 'start': value[1], 'end':value[2], 'name':feature[key], 'info':'90', 'strand':value[3]}, index=[0]), ignore_index=True)
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
    fuzzySplice=10
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
    if strand=='-':
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


def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def ilen(it):
    """Yield number of items in generator."""
    return len(list(it))

def paralle_chunk_input(chunk):
    output=pd.DataFrame([])
    for geneset in range(len(chunk)):
        df=get_df(chunk[geneset][1])
        candidates=get_tr_candidates(df)
        if len(candidates.index)>0:
            output=output.append(candidates)
    return(output)


"""Super efficient multi-threading.copleted in 10 minutes!"""

p = Pool() 
items = list(tr_coords.items())
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, 70))
end = time.time()
print('total time (s)= ' + str(end-start)) 
candidates=pd.concat(processed_values)


"""unefficient multiprocessing considering each chromosome as a group, take an hour to finish"""
#tr_coords = nested_dict()
#for line in open(input):
#    raw=line.strip().split("\t")
#    if raw[2] == 'transcript':
#        tid = raw[-1].split('; ')[1].split()[1][1:-2]
#        gid = raw[-1].split('; ')[0].split()[1][1:-1]
#        tr_coords[raw[0]][gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

#def paralle_chunk_input(chunk):
#    output=pd.DataFrame([])
#    for i in range(len(chunk)):
#        test=chunk[i][1]
#        for gene,transcript in test.items():
#            df=get_df(transcript)
#            candidates=get_tr_candidates(df)
#            if len(candidates.index)>0:
#                output=output.append(candidates)
#    return(output)    

#p = Pool() 
#items = list(tr_coords.items())
#chunksize = 3
#chunks = [items[i:i + chunksize ] for i in range(0, len(items), chunksize)]
#start = time.time()
#processed_values= p.map(paralle_chunk_input, chunks)
#candidates=pd.concat(processed_values)
#end = time.time()
#print('total time (s)= ' + str(end-start))    

""" unefficient multiprocessing option2 """
#from multiprocessing import Pool
#tr_coords = nested_dict()
#for line in open(input):
#    raw=line.strip().split("\t")
#    if raw[2] == 'transcript':
#        tid = raw[-1].split('; ')[1].split()[1][1:-2]
#        gid = raw[-1].split('; ')[0].split()[1][1:-1]
#        tr_coords[raw[0]][gid][tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

#def paralle_input(values):
#    output=pd.DataFrame([])
#    for gene, transcripts in values.items():
#        df=get_df(transcripts)
#        candidates=get_tr_candidates(df)
#        if len(candidates.index)>0:
#            output=output.append(candidates)
#    return(output)       
#p = Pool()  
#keys, values= zip(*tr_coords.iteritems())
#processed_values= p.map(paralle_input,[values[20],values[28]])
#candidates=pd.concat(processed_values)






""" Others """
#for i in range(0,len(values)):
#    print("%d\t%d" % (i,len(values[i].keys())))


for chr,info in tr_coords.iteritems():
    for gene, transcripts in info.items():
        df=get_df(transcripts)
        candidates=get_tr_candidates(df)
        if len(candidates.index)>0:
            print(candidates)
        


    


#for gene, transcripts in tr_coords.iteritems():
#    df=get_df(transcripts)
#    candidates=get_tr_candidates(df)
#    if len(candidates.index)>0:
#        print(candidates)
#        for line in range(len(candidates.index)):
#            info=compare_tr_ex_numbers(candidates.iloc[line,0],candidates.iloc[line,1])
#            compare_tr_ex(info.iloc[0,1],info.iloc[0,0])





def compare_tr_ex(trA,trB):
    fuzzy_3=10
    fuzzy_5=10
    trA_ex=get_df(ex_coords.get(trA))
    trB_ex=get_df(ex_coords.get(trB))
    trA_ex=trA_ex.drop(columns="feature")
    trA_ex_sort=trA_ex.sort_values(by='coordinate')
    trA_ex_sort.insert(0, 'exon_number', range(len(trA_ex)))
    strand=re.split('\:|\-|\(|\)',trA_ex_sort.iloc[0,1])[3]
    trB_ex=trB_ex.drop(columns="feature")
    trB_ex_sort=trB_ex.sort_values(by='coordinate')
    trB_ex_sort.insert(0, 'exon_number', range(len(trB_ex)))
    similar_ex=pd.DataFrame([])
    different_ex=pd.DataFrame([])
    if strand=='+':
        for i in range(len(trA_ex)):
            for j in range(len(trB_ex)):
                if i==0 and j==0 and (int(re.split('\:|\-|\(|\)',trA_ex_sort.iloc[i,1])[2])-int(re.split('\:|\-|\(|\)',trA_ex_sort.iloc[j,1])[2]))<=fuzzy_3:
                    similar_ex=similar_ex.append(pd.DataFrame({'transcriptA': trA, 'exonA': trA_ex_sort.iloc[i,1], 'transcriptB': trB,'exonB': trB_ex_sort.iloc[i,1]}, index=[0]), ignore_index=True)
                    

                
 








    
        
            

    
        
                

  


                


            
            
 ############ others ########
 
for index,i in enumerate(df['transcript']):
    print(df[index])
    


int(df.loc[df['transcript']=='PB.2886.1']['coordinate'][0].split()[1])-int(df.loc[df['transcript']=='PB.2886.1']['coordinate'][0].split()[0])