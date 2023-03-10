#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 22:23:40 2019

@author: beiki
"""
import os
import re
import pysam
import time
import pandas as pd
from multiprocessing import Pool
import sys

''' this program has been written in Python 3.7.5 '''
''' positional arguments:
    (A) path to splitted bam (sorted by read name) file
    (B) tissue name
    '''

def get_chunks(df,n):
    """Yield a list of n-sized chunks ."""
    for i in range(0, df.shape[0], n):
        chunk=df.iloc[i:i+n,0].to_frame()
        yield chunk
        
        
def get_reads_longest_alingment_info(chunk):
    output=pd.DataFrame([])
    for file in range(chunk.shape[0]):
        bamfile=pysam.AlignmentFile(chunk.iloc[file,0],mode='rb')
        old_seq=None
        old_seq_align_len=None
        for read in bamfile: # to remove unmapped reads, insert 'if not read.is_unmapped' as next line
            if read.query_name !=old_seq:
                output=output.append(pd.DataFrame({'read': read.query_name, 
                                                      'query_length': read.query_length, 'alignment_length':read.query_alignment_length,
                                                      'seq':read.query_sequence}, index=[0]) , ignore_index=True)
                old_seq=read.query_name
                old_seq_align_len=read.query_alignment_length
            elif read.query_name == old_seq and read.query_alignment_length > old_seq_align_len:
                row=int(output.loc[output['read']==read.query_name].index[0])
                output.iloc[row,2]=read.query_alignment_length
                output.iloc[row,3]=read.query_sequence
                old_seq_align_len=read.query_alignment_length
    return(output)
                    
                    
def get_reads_longest_alingment_fasta(chunk):
    output=[]
    for file in range(chunk.shape[0]):
        bamfile=pysam.AlignmentFile(chunk.iloc[file,0],mode='rb')
        old_seq=None
        old_seq_align_len=None
        for read in bamfile: # to remove unmapped reads, insert 'if not read.is_unmapped' as next line
            if read.query_name !=old_seq:
                output.append('>' + read.query_name)
                output.append(read.query_sequence)
                old_seq=read.query_name
                old_seq_align_len=read.query_alignment_length
            elif read.query_name == old_seq and read.query_alignment_length > old_seq_align_len:
                for i, e in enumerate(output):
                    if e == ('>' + read.query_name):
                        output[i+1]=read.query_sequence
                old_seq_align_len=read.query_alignment_length
    return(pd.DataFrame(output))

files=pd.DataFrame([])
with os.scandir(sys.argv[1]) as entries:
    for entry in entries:
        if re.search('ptmp',entry.name):
            files=files.append(pd.DataFrame({'file': entry.name}, index=[0]), ignore_index=True)

p = Pool()
start = time.time() 
processed_values= p.map(get_reads_longest_alingment_fasta, get_chunks(files, 5))
end = time.time()
print('total time (s)= ' + str(end-start)) 
df=pd.concat(processed_values)
df.to_csv((sys.argv[2]+'_fmlrc_proovread_cupcake_corrected_reads.fasta'),index = None, header=False,sep="\t")