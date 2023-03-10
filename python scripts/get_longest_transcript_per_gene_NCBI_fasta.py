#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 08:48:35 2020

@author: beiki
"""

'''
module load python/3.6.3-u4oaxsb
module load py-pandas
module load py-biopython/1.70-py3-wos466g
'''

from Bio import SeqIO
import sys
import re

lastGene = None
longest = (None, None)
res_dict={}
input_file=sys.argv[1]

for rec in SeqIO.parse(input_file, "fasta"):#'Bos_taurus.ARS-UCD1.2.ncrna.fa'
    if rec.description.find('LOC')!=-1:
        f=0
        s=rec.description.split(',')
        for i in s:
            f=f+1
            if i.find('LOC')!=-1:
                flag=f-1
        gid=re.split('\(|\)',s[flag])[1]
    else:
        gid=rec.id.split('.')[0]
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

out_name=input_file.split('.fa')[0] + '_primary_transcripts.fa'
outfile = open(out_name,"w")
for key,value in res_dict.items():
    value.id=key
    SeqIO.write(value, outfile, "fasta")
    
