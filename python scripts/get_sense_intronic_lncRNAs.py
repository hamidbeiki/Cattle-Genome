#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 09:38:55 2021

@author: beiki
"""

import pandas as pd

intersect=pd.read_csv('intersect',names=['chr','start','end','transcript_id','99','strand','tr'],sep='\t')
tr_ex_info=pd.read_csv('tr_ex_info',names=['transcript_id','number_of_exons'],sep='\t')
query_info=intersect.groupby(['transcript_id']).size().reset_index(name='count')
m=tr_ex_info.merge(query_info)
out=m.loc[m['number_of_exons']==m['count']]['transcript_id'].to_frame()
out.to_csv('sense_intronic_lncRNAs',index = None, header=False,sep="\t")