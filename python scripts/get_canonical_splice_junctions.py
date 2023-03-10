#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 13:04:38 2020

@author: beiki
"""

from Bio import SeqIO
import collections
import pandas as pd
from multiprocessing import Pool
import time




fasta_file='Cerebral_Cortex_fmlrc_proovread_cupcake_transcripts_with_introns.fa'
input='Cerebral_Cortex_fmlrc_proovread_cupcake_corrected_reads.collapsed.gff'

def nested_dict():
    return collections.defaultdict(nested_dict)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_intron_validation(input,tr):
    global seq_dict
    minimum_exon_length=10 # this is to filter out transcripts with errornous exnons, to remove this filter simply set it to 0
    introns_info=pd.DataFrame([])
    df=pd.DataFrame([v.split(';') for k, v in input.iteritems()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 3.7, like on Ceres: change "iteritems" to "items"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    if df['ex_length'].min() >=minimum_exon_length and df.iloc[0,3]=="+":
        df=df.sort_values(by='ex_start', ascending=True)
        start=df.iloc[0,1]
        valid_introns=0
        introns=0
        tr_seq=str(seq_dict[tr])
        for i in range(df.shape[0]-1):
            introns=introns+1
            intron_start=df.iloc[i,2]+1-start-1# '-1' used because in "bedtools getfasta" locations start from '0'
            intron_end=df.iloc[i+1,1]-1-start-1# '-1' used because in "bedtools getfasta" locations start from '0'
            intron_acceptor=tr_seq[intron_start:intron_start+2]
            intron_donor=tr_seq[intron_end-1:intron_end+1]
            if (intron_acceptor=='GT' or intron_acceptor=='gt') and (intron_donor=='AG' or intron_donor=='ag'):
                valid_introns=valid_introns+1
        introns_info=introns_info.append(pd.DataFrame({'tr': tr, 'number_of_introns': introns,'number_of_valid_introns':valid_introns,'strand':df.iloc[0,3]}, index=[0]), ignore_index=True)
        introns_info = introns_info[['tr','number_of_introns','number_of_valid_introns','strand']]
    elif df['ex_length'].min() >=minimum_exon_length and df.iloc[0,3]=="-":
        df=df.sort_values(by='ex_start', ascending=False)
        start=df.iloc[0,2]
        valid_introns=0
        introns=0
        tr_seq=str(seq_dict[tr])
        for i in range(df.shape[0]-1):
            introns=introns+1
            intron_start=start-df.iloc[i,1]+1-1# '-1' used because in "bedtools getfasta" locations start from '0'
            intron_end=start-df.iloc[i+1,2]-1-1# '-1' used because in "bedtools getfasta" locations start from '0'
            intron_acceptor=tr_seq[intron_start:intron_start+2]
            intron_donor=tr_seq[intron_end-1:intron_end+1]
            if (intron_acceptor=='CT' or intron_acceptor=='ct') and (intron_donor=='AC' or intron_donor=='ac'):
#            if (intron_acceptor=='AC' or intron_acceptor=='ac') and (intron_donor=='CT' or intron_donor=='ct'):
                valid_introns=valid_introns+1
        introns_info=introns_info.append(pd.DataFrame({'tr': tr, 'number_of_introns': introns,'number_of_valid_introns':valid_introns,'strand':df.iloc[0,3]}, index=[0]), ignore_index=True)
        introns_info = introns_info[['tr','number_of_introns','number_of_valid_introns','strand']]
    return(introns_info)

def paralle_chunk_input(chunk):
    info=pd.DataFrame([])
    for tr in range(len(chunk)):
        if len(chunk[tr][1])>1:
            intron_info=get_intron_validation(chunk[tr][1],chunk[tr][0])
            info=info.append(intron_info)
    return(info)
    
def ilen(it):
    """Yield number of items in generator."""
    return len(list(it))

        
ex_coords = nested_dict()
for line in open(input):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        tid=raw[-1].split('; ')[1].split()[1][1:-2]
        exnumber=raw[3]
        ex_coords[tid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
        
        
seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

items = list(ex_coords.items())
p = Pool() 
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, 400))
end = time.time()
print('total time (s)= ' + str(end-start)) 
introns_info=pd.concat(processed_values)

test=[x for x in introns_info.iterrows() if x[1][1]==x[1][2] ]
out=pd.DataFrame([])
for i in range(len(test)):
    if test[i][1][3]=="-":
        out=out.append(pd.DataFrame({'tr': test[i][1][0]}, index=[0]), ignore_index=True)


introns_info.to_csv(r'introns_info',index = None, header=False,sep="\t")
