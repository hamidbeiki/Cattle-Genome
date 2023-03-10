#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 22:56:11 2019

@author: beiki
"""

import os
import re
import pysam
import time
import pandas as pd
from multiprocessing import Pool
import sys

''' to get error rate in a given region, use samfile = pysam.AlignmentFile("ex1.bam", "rb" ) function, more info:
    https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
    '''

''' parallel processing '''


def get_chunks(df,n):
    """Yield a list of n-sized chunks ."""
    for i in range(0, df.shape[0], n):
        chunk=df.iloc[i:i+n,0].to_frame()
        yield chunk

def get_mismatchesANDdeletions(md):
    mismatchesANDdeletions=0
    for s in md:
        if s in "ACGT":
            mismatchesANDdeletions += 1
            
    return(mismatchesANDdeletions)
    
def get_reads_qc(chunk):
    reads_qc=pd.DataFrame([])
    for file in range(chunk.shape[0]):
        bamfile=pysam.AlignmentFile(chunk.iloc[file,0],mode='rb')
        for read in bamfile:
            if not read.is_unmapped:
                stats=read.get_cigar_stats()
                md=read.get_tag('MD')
                mismatch=get_mismatchesANDdeletions(md)-stats[0][2]
                match=stats[0][0]-mismatch
                read_length=read.infer_query_length()
                softClip=stats[0][4]
                reads_qc=reads_qc.append(pd.DataFrame({'read': read.query_name,'read_length': read_length,'MATCH': match
                , 'INSERTION': stats[0][1],'DELETION': stats[0][2],'MISMATCH': mismatch, 'softClip': softClip}, index=[0]) # stats[0][4] reflect softClip more info: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
                , ignore_index=True)
    return(reads_qc)

files=pd.DataFrame([])
with os.scandir(sys.argv[1]) as entries:
    for entry in entries:
        if re.search('ptmp',entry.name):
            files=files.append(pd.DataFrame({'file': entry.name}, index=[0]), ignore_index=True)

p = Pool()
start = time.time() 
processed_values= p.map(get_reads_qc, get_chunks(files, 5))
end = time.time()
print('total time (s)= ' + str(end-start)) 
reads_qc=pd.concat(processed_values)
reads_qc['ERROR_RATE']=(reads_qc['INSERTION']+reads_qc['DELETION']+reads_qc['MISMATCH']+reads_qc['softClip'])/(reads_qc['read_length']+reads_qc['DELETION'])

print(reads_qc['ERROR_RATE'].mean())
print(reads_qc['ERROR_RATE'].median())
reads_qc.to_csv(r'reads_qc',index = None, header=True,sep="\t")
''' single cpu

reads_qc=pd.DataFrame([])

start = time.time()
with os.scandir('/project/cattle_genome_assemblies/IsoSeq_analysis/Jreecy/Iso-seq/correction_evaluation/test2') as entries:
    for entry in entries:
        if re.search('output',entry.name):
            bamfile=pysam.AlignmentFile(entry.name,mode='rb')
            for read in bamfile:
                stats=read.get_cigar_stats()
                #    mismatch=stats[0][10]-stats[0][1]-stats[0][2]
                reads_qc=reads_qc.append(pd.DataFrame({'read': read.query_name, 'MATCH': stats[0][0]
                , 'INSERTION': stats[0][1],'DELETION': stats[0][2],'DIFFERENCES': stats[0][10]}, index=[0])
                , ignore_index=True)
        
end = time.time()
print('total time (s)= ' + str(end-start))

'''
