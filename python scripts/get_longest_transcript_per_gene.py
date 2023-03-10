#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:32:22 2020

@author: beiki
"""

'''
module load python/3.6.3-u4oaxsb
module load py-pandas
module load py-biopython/1.70-py3-wos466g
'''

from Bio import SeqIO
import sys


lastGene = None
longest = (None, None)
res_dict={}
input_file=sys.argv[1]

for rec in SeqIO.parse(input_file, "fasta"):#combined_transcripts.fa
    a=rec.id.split("_")[0].split(".")
    gid=a[0] + '.' + a[1]
    l = len(rec)#same as len(rec.seq)
    if lastGene is not None:
        if gid == lastGene:
            if longest[0] < l:
                longest = (l, rec)
        else:
            res_dict[lastGene]=longest[1]
            lastGene = gid
            longest = (l, rec)
    else:
        lastGene = gid
        longest = (l, rec)
if gid not in res_dict:
    res_dict[gid]=rec


out_name='aaa'
outfile = open(out_name,"w")
for key,value in res_dict.items():
    value.id=key
    SeqIO.write(value, outfile, "fasta")
    
    
