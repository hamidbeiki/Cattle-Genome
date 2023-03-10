#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 09:55:47 2020

@author: beiki
"""

import collections
import pandas as pd
from multiprocessing import Pool
import time
import sys
import re

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_firstExon_bed(input,gene):
    '''no need to filter for minimum_exon_length because transcripts has been already firtered for this criteria in upstream steps
        this will atually consider the last exons on genes in "-" strand as first exon that is fine because I'm lookig for divergent transcription'''
    firstExon_bed=pd.DataFrame([])
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    df=df.sort_values(by='ex_start', ascending=True)
    firstExon_bed=firstExon_bed.append(pd.DataFrame({'chr': df.iloc[0,0], 'exon_start': df.iloc[0,1],'exon_end':df.iloc[0,2],
                                                             'gene':gene,'value':99,'strand':df.iloc[0,3],
                                                             'name':gene}, index=[0]), ignore_index=True)
    firstExon_bed = firstExon_bed[['chr','exon_start','exon_end','gene','value','strand','name']]
    return(firstExon_bed)

def get_geneBorders_bed(input,gene):
    '''no need to filter for minimum_exon_length because transcripts has been already firtered for this criteria in upstream steps
        this will atually consider the last exons on genes in "-" strand as first exon that is fine because I'm lookig for divergent transcription'''
    geneBorders_bed=pd.DataFrame([])
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    gene_end=df.iloc[0,2]
    geneBorders_bed=geneBorders_bed.append(pd.DataFrame({'chr': df.iloc[0,0], 'gene_start': gene_start,'gene_end':gene_end,
                                                             'gene':gene,'value':99,'strand':df.iloc[0,3],
                                                             'name':gene}, index=[0]), ignore_index=True)
    geneBorders_bed = geneBorders_bed[['chr','gene_start','gene_end','gene','value','strand','name']]
    return(geneBorders_bed)

def get_geneBordersAndPromoter_bed(input,gene):
    '''no need to filter for minimum_exon_length because transcripts has been already firtered for this criteria in upstream steps
        this will atually consider the last exons on genes in "-" strand as first exon that is fine because I'm lookig for divergent transcription'''
    geneBorders_bed=pd.DataFrame([])
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    gene_end=df.iloc[0,2]
    promoter_lengh=150
    if df.iloc[0,3]=="+":
        geneBorders_bed=geneBorders_bed.append(pd.DataFrame({'chr': df.iloc[0,0], 'gene_start': int(gene_start) - int(promoter_lengh),'gene_end':gene_end,
                                                                 'gene':gene,'value':99,'strand':df.iloc[0,3],
                                                                 'name':gene}, index=[0]), ignore_index=True)
    elif df.iloc[0,3]=="-":
        geneBorders_bed=geneBorders_bed.append(pd.DataFrame({'chr': df.iloc[0,0], 'gene_start': gene_start,'gene_end': int(gene_end) + int(promoter_lengh),
                                                                 'gene':gene,'value':99,'strand':df.iloc[0,3],
                                                                 'name':gene}, index=[0]), ignore_index=True)        
    geneBorders_bed = geneBorders_bed[['chr','gene_start','gene_end','gene','value','strand','name']]
    return(geneBorders_bed)
    

def paralle_chunk_input(chunk):
    bed=pd.DataFrame([])
    for gene in range(len(chunk)):
        if len(chunk[gene][1])>=1:
            gene_bed=get_geneBordersAndPromoter_bed(chunk[gene][1],chunk[gene][0])#RNAbased or DNAbased methods can be selected
            bed=bed.append(gene_bed)
    return(bed)

def ilen(it):
    """Yield number of items in generator."""
    return len(list(it))

input_file=sys.argv[1]#your gff file       
ex_coords = nested_dict()
for line in open(input_file):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        gid=raw[-1].split('; ')[0].split()[1]
        gid=re.sub('["]', '', gid)
        exnumber=raw[3]
        ex_coords[gid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
 
items = list(ex_coords.items())  
p = Pool() 
start = time.time()
processed_values= p.map(paralle_chunk_input, get_chunks(items, 400))
end = time.time()
print('total time (s)= ' + str(end-start)) 
geneBorders_bed=pd.concat(processed_values)
geneBorders_bed.to_csv(r'geneBorders_bed',index = None, header=False,sep="\t")

'''
for line in open(input_file):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        gid=raw[-1].split('; ')[0].split()[1]
        gid=re.sub('["]', '', gid)
        if gid=="PB.608":
            print(gid)'''
            
 