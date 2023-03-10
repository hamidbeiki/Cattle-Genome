#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 12:06:06 2021

@author: beiki
"""

""" 
    NOTE: on nova use
    module load miniconda3/4.3.30-qdauveb
    source activate /home/beiki/.conda/envs/py-libs
"""
'''
module load python/3.6.3-u4oaxsb
module load py-biopython/1.70-py3-wos466g
'''

from Bio import SeqIO
import math
import sys

bunch_size=int(sys.argv[2])#number of sequences in each splitted fasta output file
out_name=sys.argv[3]
number_of_seqs=0
for rec in SeqIO.parse(sys.argv[1], "fasta"):# 3UTRs.fasta or protein_coding_lncRNAs_3UTRs.fasta
    number_of_seqs=number_of_seqs+1
lim=math.floor(number_of_seqs/float(bunch_size))
result_dict={}
bunch=[]
counter=0
flag=0
last_bunch=[]
for rec in SeqIO.parse(sys.argv[1], "fasta"):# 3UTRs.fasta or protein_coding_lncRNAs_3UTRs.fasta
    counter=counter+1
    flag=flag+1
    if flag<bunch_size:
        bunch.append(rec)
    elif flag>bunch_size:
        bunch=[]
        flag=1
        bunch.append(rec)
    elif flag==bunch_size:
        result_dict[rec.id]=bunch
    elif counter>(lim*bunch_size):
        last_bunch.append(rec)
result_dict[rec.id]=last_bunch
    

flag=1
for key,value in result_dict.items():
    out_name2=out_name + "_" + str(flag) + ".fasta"
    outfile = open(out_name2,"w")
    SeqIO.write(value, outfile, "fasta")
    flag=flag+1
    print("%s\t%s" %(out_name2,'file is written into the disc'))
        

