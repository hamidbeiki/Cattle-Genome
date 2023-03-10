#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 11:55:23 2020

@author: beiki
"""
from Bio import SeqIO
import sys

lastTranscript = None
longest = (None, None)
for rec in SeqIO.parse(sys.argv[1], "fasta"):#combined_coding_nonNMD_orfs.fa 
    transcript = rec.id.split("_")[0]
    l = len(rec)#same as len(rec.seq)
    if lastTranscript is not None:
        if transcript == lastTranscript:
            if longest[0] < l:
                longest = (l, rec)
        else:
            lastTranscript = transcript
            SeqIO.write(longest[1], sys.stdout, "fasta")
            longest = (l, rec)
    else:
        lastTranscript = transcript
        longest = (l, rec)
SeqIO.write(longest[1], sys.stdout, "fasta")