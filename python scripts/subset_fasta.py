#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 15:42:22 2020

@author: beiki
"""

from Bio import SeqIO
import sys
import pandas as pd
import functools
import operator

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])
 
input_fasta=sys.argv[1]#combined_transcripts.fa
keep=sys.argv[2]# a file containing list of ids to keep, nc-transcripts
keep_df=pd.read_csv(keep,sep=' ')
l=keep_df.values.tolist()
keep_l=flatten_nested_list(l)

for rec in SeqIO.parse(input_fasta, "fasta"): 
    if rec.id in keep_l:
        SeqIO.write(rec, sys.stdout, "fasta")