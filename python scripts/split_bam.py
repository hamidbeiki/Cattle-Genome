#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:14:17 2019

@author: beiki
"""

import pysam
import sys
infile = pysam.AlignmentFile(sys.argv[1])

chunk_size = 1000
outfile_pattern = "ptmp_segement%d.bam"

chunk = 0
reads_in_this_chunk = 0
old_name = None
outfile = pysam.AlignmentFile(outfile_pattern % chunk, "w", template = infile)

for read in infile.fetch(until_eof=True):

    if old_name != read.query_name and reads_in_this_chunk > chunk_size:
        reads_in_this_chunk = 0
        chunk += 1
        outfile.close()
        outfile = pysam.AlignmentFile(outfile_pattern % chunk, "w", template = infile)

    outfile.write(read)
    old_name = read.query_name
    reads_in_this_chunk += 1

outfile.close()
