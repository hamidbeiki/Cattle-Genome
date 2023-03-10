#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:15:44 2020

@author: beiki
"""
import pandas as pd
import pybedtools
import collections
import re
import time
import sys
from multiprocessing import Pool
from scipy.stats import ks_2samp
import numpy as np

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_genes_window(gene_start,gene_end,window,strand):
#if you are intrested in to put window jus on genes 5prime end, replace get_genes_window 
#function with this function 
    if strand=="+":
        if gene_start>1 and gene_start<=window:
            window_start=1
            window_end=gene_end
        elif gene_start>window:
            window_start=gene_start-window
            window_end=gene_end
        elif gene_start==1:
            window_start=1
            window_end=gene_end
    elif strand=="-":
        s=gene_end
        e=gene_start
        window_start=e
        window_end=s+window
    return(window_start,window_end)

def get_geneWindow_dict(input,gene,window_size):
    '''no need to filter for minimum_exon_length because transcripts has been already firtered for this criteria in upstream steps
        this will atually consider the last exons on genes in "-" strand as first exon that is fine because I'm lookig for divergent transcription'''
    geneBorders_dict={}
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    strand=df.iloc[0,3]
    gene_end=df.iloc[0,2]
    gene_start=int(gene_start)
    gene_end=int(gene_end)
    info=get_genes_window(gene_start,gene_end,window_size,strand)
    window_start=info[0]
    window_end=info[1]
    geneBorders_dict[gene] = "{0},{1},{2},{3},{4},{5},{6}".format(df.iloc[0,0],window_start,window_end,gene,99,df.iloc[0,3],gene)      
    return(geneBorders_dict)

def get_geneWindow_bed(candidates_dict,window_size):
    out_dict={}
    for gene in candidates_dict:
        gene_dict=get_geneWindow_dict(candidates_dict[gene],gene,window_size)
        out_dict.update(gene_dict)
    out_bed=pd.DataFrame.from_dict(out_dict,orient='index')[0].str.split(',',expand=True).rename({0:'chr',1: 'promoter_start', 2: 'gene_end', 3: 'gene_id', 4: 'value', 5: 'strand', 6: 'gene_id'},axis=1)
    return(out_bed)

def get_gene_startEnd(input):
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df=df.sort_values(by='ex_start', ascending=True)
    gene_start=df.iloc[0,1]
    df=df.sort_values(by='ex_end', ascending=False)
    chromosome=df.iloc[0,0]
    gene_end=df.iloc[0,2]
    gene_start=int(gene_start)
    gene_end=int(gene_end)
    strand=df.iloc[0,3]
    return(chromosome,gene_start,gene_end,strand)
    
def get_genes_related_to_CNVs(intersections):
    global ex_coords
    df=pd.DataFrame([])
    for i in range(intersections.shape[0]):
        gene=intersections.iloc[i,3]
        info=get_gene_startEnd(ex_coords[gene])
        gene_start=info[1]
        gene_end=info[2]
        cnv=intersections.iloc[i,10]
        cnv_start=intersections.iloc[i,8]
        strand=intersections.iloc[i,5]
        if strand=="+":
            distance=abs(gene_start-cnv_start)
        elif strand=="-":
            distance=abs(gene_end-cnv_start)
        df=df.append(pd.DataFrame({'gene': gene, 'cnv': cnv, 'distance': distance}, index=[0]), ignore_index=True)## capture genes with minimum distance to CNVs (from their 5prime end) as related gene
    info=df.loc[df.groupby('cnv')['distance'].idxmin(), :]
    genes=info["gene"].to_frame().drop_duplicates()## list of genes related (closest) to CNVs
#    cnvs=info["cnv"].to_frame().drop_duplicates()## list of cnvs related to genes
    return(genes)
    
def paralle_chunk_input2(chunk):
    global bidirection_promoter_events
    global ex_coords
    global hotspots_bed
    global cnv_bed
    global all_candidates_dict
    global all_candidates
    results=pd.DataFrame([])
    for i in range(len(chunk)):
        window_size=chunk[i]#We found a ~1-Mb window size to be optimal for capturing cis effects of CNVs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4538210/
        candidates=pd.concat([bidirection_promoter_events.iloc[0:,0].to_frame().rename({'geneA':'gene_id'},axis=1),bidirection_promoter_events.iloc[0:,1].to_frame().rename({'geneB':'gene_id'},axis=1)]).drop_duplicates().values.tolist()
        candidates_dict={}
        for gene in ex_coords:
            if candidates.count([gene]):
                candidates_dict[gene]=ex_coords[gene]
        candidates_bed=get_geneWindow_bed(candidates_dict,window_size)
        candidates_bed=pybedtools.BedTool.from_dataframe(candidates_bed)
        candidates_bed=candidates_bed.sort()
        intersections=candidates_bed.intersect(hotspots_bed,wa=True,wb=True)
        if intersections.total_coverage() >0:
            intersections=intersections.to_dataframe()
            genes1=get_genes_related_to_CNVs(intersections)
        else:
            genes1=pd.DataFrame([])    
        intersections=candidates_bed.intersect(cnv_bed,wa=True,wb=True)
        if intersections.total_coverage() >0:
            intersections=intersections.to_dataframe()
            genes2=get_genes_related_to_CNVs(intersections)
        else:
            genes2=pd.DataFrame([])
        if genes1.shape[0]>0 and genes2.shape[0]>0:
            m=genes1.merge(genes2)
            hotspot_percentage=(m.shape[0]/float(genes1.shape[0]))*100
            cnv_percentage=(m.shape[0]/float(genes2.shape[0]))*100
        else:
            hotspot_percentage=0
            cnv_percentage=0  
        results=results.append(pd.DataFrame({'window': chunk[i], 'Bidirection_hotspot%_shared_with_cnv': hotspot_percentage,'Bidirection_cnv%_shared_with_hotspot':cnv_percentage}, index=[0]), ignore_index=True)
    return(results)
        
def paralle_chunk_input3(chunk):
    global ex_coords
    global hotspots_bed
    global cnv_bed
    global convergent_gene_events
    results=pd.DataFrame([])
    for i in range(len(chunk)):
        window_size=chunk[i]#We found a ~1-Mb window size to be optimal for capturing cis effects of CNVs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4538210/
        candidates=pd.concat([convergent_gene_events.iloc[0:,0].to_frame().rename({'geneA':'gene_id'},axis=1),convergent_gene_events.iloc[0:,2].to_frame().rename({'geneB':'gene_id'},axis=1)]).drop_duplicates().values.tolist()
        candidates_dict={}
        for gene in ex_coords:
            if candidates.count([gene]):
                candidates_dict[gene]=ex_coords[gene]
        candidates_bed=get_geneWindow_bed(candidates_dict,window_size)
        candidates_bed=pybedtools.BedTool.from_dataframe(candidates_bed)
        candidates_bed=candidates_bed.sort()
        intersections=candidates_bed.intersect(hotspots_bed,wa=True,wb=True)
        if intersections.total_coverage() >0:
            intersections=intersections.to_dataframe()
            genes3=get_genes_related_to_CNVs(intersections)
        else:
            genes3=pd.DataFrame([])
        intersections=candidates_bed.intersect(cnv_bed,wa=True,wb=True)
        if intersections.total_coverage() >0:
            intersections=intersections.to_dataframe()
            genes4=get_genes_related_to_CNVs(intersections)
        else:
            genes4=pd.DataFrame([])
        if genes3.shape[0]>0 and genes4.shape[0]>0:
            m=genes3.merge(genes4)
            hotspot_percentage=(m.shape[0]/float(genes3.shape[0]))*100
            cnv_percentage=(m.shape[0]/float(genes4.shape[0]))*100
        else:
            hotspot_percentage=0
            cnv_percentage=0              
        results=results.append(pd.DataFrame({'window': chunk[i], 'Convergent_hotspot%_shared_with_cnv': hotspot_percentage,'Convergent_cnv%_shared_with_hotspot':cnv_percentage}, index=[0]), ignore_index=True)
    return(results)
        
input_gff=sys.argv[1] # input_gff="Cerebral_Cortex_final.collapsed.gff"
tissue=input_gff.split("_final")[0]
hotspots_bed=pd.read_csv('filtered_recombination_hotspots_ARS-UCD1.2.bed',sep="\t",header=0)## all data are in 'recombination_hotspots_ARS-UCD1.2.bed' file
hotspots_bed=pybedtools.BedTool.from_dataframe(hotspots_bed)
hotspots_bed=hotspots_bed.sort()
cnv_bed=pd.read_csv('CNVs_ARS-UCD1.2.bed',header=0,sep="\t")
cnv_bed=pybedtools.BedTool.from_dataframe(cnv_bed)
cnv_bed=cnv_bed.sort()
bidirection_promoter_events=pd.read_csv(tissue + '_bidirection_promoter_events',header=0,sep="\t")
convergent_gene_events=pd.read_csv(tissue + '_convergent_gene_events',header=0,sep="\t")

ex_coords = nested_dict()
for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        gid=raw[-1].split('; ')[0].split()[1]
        gid=re.sub('["]', '', gid)
        exnumber=raw[3]
        ex_coords[gid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
        
all_candidates=[]
for gene in ex_coords:
    all_candidates.append(gene)

all_candidates_dict={}

for gene in ex_coords:
    if all_candidates.count(gene):
        all_candidates_dict[gene]=ex_coords[gene]
        
windows=[]
for window in range(0,1200100,10000):
      windows.append(window) 

p = Pool() 
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets (i.e. uptimum_number_of_genes)in each chunk ##
processed_values= p.map(paralle_chunk_input2, get_chunks(windows, 3))
end = time.time()
#print('total time (s)= ' + str(end-start)) 
bidirection_res=pd.concat(processed_values)
bidirection_res=bidirection_res[['window', 'Bidirection_hotspot%_shared_with_cnv', 'Bidirection_cnv%_shared_with_hotspot']]

p = Pool() 
start = time.time()
processed_values= p.map(paralle_chunk_input3, get_chunks(windows, 3))
end = time.time()
convergent_res=pd.concat(processed_values)
convergent_res=convergent_res[['window', 'Convergent_hotspot%_shared_with_cnv', 'Convergent_cnv%_shared_with_hotspot']]

merge1=bidirection_res.merge(convergent_res)
merge1.to_csv( tissue + '_gene_cnv_hotspot_comparision',index = None, header=True,sep="\t") 

x=np.array(merge1['Bidirection_hotspot%_shared_with_cnv'])
y=np.array(merge1['Bidirection_cnv%_shared_with_hotspot'])
pvalue1=ks_2samp(x, y)[1]
print("p-value for comparision of shared genes between hotspots & CNVs in Bidirection genes using Two-sample KS test")
print(pvalue1)
x=np.array(merge1['Convergent_hotspot%_shared_with_cnv'])
y=np.array(merge1['Convergent_cnv%_shared_with_hotspot'])
pvalue2=ks_2samp(x, y)[1]
print("p-value for comparision of shared genes between hotspots & CNVs in Convergent genes using Two-sample KS test")
print(pvalue2)

"""to find the point that cnv and hot spot diverged significantly
results=pd.DataFrame([])
for i in range(10,merge1.shape[0],10):# this is to find the point that cnv and hot spot diverged significantly
    data=merge1.iloc[0:i+1,0:]
    x=np.array(data['Bidirection_hotspot%_shared_with_cnv'])
    y=np.array(data['Bidirection_cnv%_shared_with_hotspot'])
    pvalue1=ks_2samp(x, y)[1]
    x=np.array(data['Convergent_hotspot%_shared_with_cnv'])
    y=np.array(data['Convergent_cnv%_shared_with_hotspot'])
    pvalue2=ks_2samp(x, y)[1]
    results=results.append(pd.DataFrame({'distance': data.iloc[-1,0], 'bidirection_logPvalue': -np.log2(pvalue1), 'convergent_logPvalue': -np.log2(pvalue2)}, index=[0]), ignore_index=True)

results=results[['distance','bidirection_logPvalue', 'convergent_logPvalue']]
"""      