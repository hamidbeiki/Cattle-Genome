#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 16:19:16 2020

@author: beiki
"""

from Bio import SeqIO
import sys

input_file=sys.argv[1]#your dna file, Assuming your fasta file correspond exactly to the CDS of the protein sequences you want

for rec in SeqIO.parse(input_file, "fasta"): 
    rec.seq=rec.seq.translate(table=1)
    SeqIO.write(rec, sys.stdout, "fasta")