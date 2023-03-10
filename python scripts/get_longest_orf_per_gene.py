#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:31:59 2020

@author: beiki
"""

from Bio import SeqIO
import sys

coding_file=sys.argv[1]#combined_coding_nonNMD_orfs.fa
noncoding_file=sys.argv[2]#combined_tnc_nonNMD_orfs.fa
lastGene = None
longest = (None, None)
res_dict={}
for rec in SeqIO.parse(coding_file, "fasta"): 
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

res_dict2={}
for rec in SeqIO.parse(noncoding_file, "fasta"): 
    a=rec.id.split("_")[0].split(".")
    gid=a[0] + '.' + a[1]
    if gid not in res_dict:
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
if gid not in res_dict2 and gid not in res_dict:
    res_dict2[gid]=rec

res_dict.update(res_dict2)

outfile = open('combined_genes_longest_orfs.fa',"w")
for key,value in res_dict.items():
    value.id=key
    SeqIO.write(value, outfile, "fasta")