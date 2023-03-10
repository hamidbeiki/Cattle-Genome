#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:31:16 2021

@author: beiki
"""

'''
module load python/3.6.3-u4oaxsb
module load py-pandas
module load py-biopython/1.70-py3-wos466g
module load py-matplotlib/3.0.0-py3-xr6eijv
module load py-scipy/1.1.0-py3-pzig4lr
'''

from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
import pandas as pd
import matplotlib.pyplot as plt

out=sys.argv[1].split('.')[0]#5UTRs.fasta,3UTRs.fasta
res=[]
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    a=len(rec.seq)
    b=GC(rec.seq)
    res.append([rec.id , len(rec.seq) , GC(rec.seq)])
#    print("%s\t%i\t%i" % (rec.id, len(rec.seq), GC(rec.seq)))
    
df=pd.DataFrame(res)    
df=df.rename({0:'tr_id',1:'UTR_length',2:'GC_content'},axis=1)
df.to_csv(out+'_length_GCcontent',index = None, header=True,sep="\t")     

df.plot(kind="scatter", x="GC_content", y="UTR_length", alpha=0.05,s=2)
plt.savefig(out+'_length_GCcontent'+'.png')
plt.close()



"""
utr5=pd.read_csv('5UTRs_length_GCcontent',header=0,sep='\t')
utr5=utr5.rename({'UTR_length':'5UTR_length','GC_content':'5UTR_GC_content'},axis=1)
utr3=pd.read_csv('3UTRs_length_GCcontent',header=0,sep='\t')
utr3=utr3.rename({'UTR_length':'3UTR_length','GC_content':'3UTR_GC_content'},axis=1)
known_coding_nonNMD=pd.read_csv('coding_nonNMD_transcripts_with_known_protein',names=['tr_id'],sep='\t')
known_tnc_nonNMD=pd.read_csv('tnc_nonNMD_transcripts_with_known_protein',names=['tr_id'],sep='\t')
tr_biotypes=pd.read_csv('combined_transcript_biotypes',header=0,sep='\t')
tr_biotypes=tr_biotypes.rename({'ID':'tr_id'},axis=1)
protein_coding=tr_biotypes.loc[tr_biotypes['biotype']=='protein_coding']['tr_id'].to_frame()

m=utr5.merge(protein_coding)
m.plot(kind="scatter", x="5UTR_GC_content", y="5UTR_length", alpha=0.05,s=2)
plt.savefig('protein_coding_tr_5UTR_length_GCcontent'+'.png')
plt.close()

m=utr3.merge(protein_coding)
m.plot(kind="scatter", x="3UTR_GC_content", y="3UTR_length", alpha=0.05,s=2)
plt.savefig('protein_coding_tr_3UTR_length_GCcontent'+'.png')
plt.close()

m=utr5['tr_id'].to_frame().merge(known_coding_nonNMD,how='outer', left_on=['tr_id'],right_on=['tr_id'],indicator=True)
m=m.rename({'_merge':'homology'},axis=1)
m['homology']=m['homology'].str.replace('left_only','0')#0:unkown protein
m['homology']=m['homology'].str.replace('right_only','1')#1:know protein
m['homology']=m['homology'].str.replace('both','1')#1:know protein
m2=utr5.merge(m)
m2['homology']=m2['homology'].astype(int)
m3=m2.merge(utr3)
m3.plot(kind="scatter", x="5UTR_GC_content", y="5UTR_length", alpha=0.4,
    figsize=(10,7),
    c="homology", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,s=2)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.savefig('protein-coding_transcripts-known.vs.unkown_protein_5UTR.png')
plt.close()

m3.plot(kind="scatter", x="3UTR_GC_content", y="3UTR_length", alpha=0.4,
    figsize=(10,7),
    c="homology", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,s=2)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.savefig('protein-coding_transcripts-known.vs.unkown_protein_3UTR.png')
plt.close()

from scipy import stats
import numpy as np

s1=m3.loc[m3['homology']==1]['5UTR_GC_content']
x=np.array(s1) 
s2=m3.loc[m3['homology']==0]['5UTR_GC_content']
y=np.array(s2)
#stats.ttest_ind(x,y)[1]  
#stats.wilcoxon(x,y)
stats.mannwhitneyu(x,y)#Mann-Whitney test that is equivalent to Wilcoxon Rank Sum test with uequal sample sizes: https://stats.stackexchange.com/questions/368881/wilcoxon-rank-sum-test-unequal-sample-sizes

groups = m2.groupby("homology")
for name, group in groups:
    plt.plot(group["X Value"], group["Y Value"], marker="o", linestyle="", label=name)

df.plot(kind="scatter", x="GC_content", y="UTR_length", alpha=0.4,
    s=df["population"]/100, label="population", figsize=(10,7),
    c="median_house_value", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False)
plt.legend()
"""
