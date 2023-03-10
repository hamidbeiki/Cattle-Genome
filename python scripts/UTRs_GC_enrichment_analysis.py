#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:33:38 2021

@author: beiki
"""
'''
module load python/3.6.3-u4oaxsb
module load py-pandas
module load py-biopython/1.70-py3-wos466g
module load py-matplotlib/3.0.0-py3-xr6eijv
module load py-scipy/1.1.0-py3-pzig4lr
'''

import numpy as np
from scipy import stats
from scipy.stats.stats import tiecorrect, rankdata
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Result(object):
    pass

def mannwhitneyu(x, y, use_continuity=True, n_rep=100):
    """
    https://gist.github.com/josef-pkt/3866149
    Computes the Mann-Whitney rank test on samples x and y.
    Parameters
    ----------
    x, y : array_like
        Array of samples, should be one-dimensional.
    use_continuity : bool, optional
        Whether a continuity correction (1/2.) should be taken into
        account. Default is True.
    n_rep : int
        number of replication in permutation test
    Returns
    -------
    u : float
        The Mann-Whitney statistics.
    prob : float
        Two-sided p-value assuming a asymptotic normal distribution.
    res : result instance
        results from permutation test, contains all permutations
    Notes
    -----
    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.
    This test corrects for ties and by default uses a continuity correction.
    The reported p-value is for a two-sided hypothesis, in contrast to the
    original test in scipy.stats
    """
    n_rep = max([n_rep, 1])  #need at least one iteration through loop
    x = np.asarray(x)
    y = np.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x,y)))
    #T = np.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    T = tiecorrect(ranked)
    if T == 0:
        raise ValueError('All numbers are identical in amannwhitneyu')
    sd = np.sqrt(T*n1*n2*(n1+n2+1)/12.0)
    res = Result()
    res.ranked = ranked.copy()
    res.n1 = n1
    res.n2 = n2
    res.T = T
    res.sd = sd
    res.all_results = np.zeros((n_rep, 3))
    for ii in range(n_rep):
        rankx = ranked[0:n1]       # get the x-ranks
        #ranky = ranked[n1:]        # the rest are y-ranks
        u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx,axis=0)  # calc U for x
        u2 = n1*n2 - u1                            # remainder is U for y
        bigu = max(u1,u2)
        smallu = min(u1,u2)
        if use_continuity:
            # normal approximation for prob calc with continuity correction
            z = np.abs((bigu-0.5-n1*n2/2.0) / sd)
        else:
            z = np.abs((bigu-n1*n2/2.0) / sd)  # normal approximation for prob calc
        res.all_results[ii] = [u1, bigu, z]
        if ii == 0:
            smallu0 = smallu
            z0 = z
        np.random.shuffle(ranked)
    pval_perm = (res.all_results[:,-1] >= res.all_results[0, -1]).mean(0)
    #correct pvalue to two-sided compared to scipy.stats
    return smallu0, 2*stats.norm.sf(z0), pval_perm, res  #(1.0 - zprob(z))
utr5=pd.read_csv('5UTRs_length_GCcontent',header=0,sep='\t')
utr5=utr5.rename({'UTR_length':'5UTR_length','GC_content':'5UTR_GC_content'},axis=1)
utr3=pd.read_csv('3UTRs_length_GCcontent',header=0,sep='\t')
utr3=utr3.rename({'UTR_length':'3UTR_length','GC_content':'3UTR_GC_content'},axis=1)
known_coding_nonNMD=pd.read_csv('coding_nonNMD_transcripts_with_known_protein',names=['tr_id'],sep='\t')
known_tnc_nonNMD=pd.read_csv('tnc_nonNMD_transcripts_with_known_protein',names=['tr_id'],sep='\t')
tr_biotypes=pd.read_csv('combined_transcript_biotypes',header=0,sep='\t')
tr_biotypes=tr_biotypes.rename({'ID':'tr_id'},axis=1)
protein_coding=tr_biotypes.loc[tr_biotypes['biotype']=='protein_coding']['tr_id'].to_frame()
coding_nonNMD_cds=pd.read_csv('combined_coding_nonNMD_transcripts_longest_cds_length_GCcontent',header=0,sep='\t')
coding_nonNMD_cds=coding_nonNMD_cds.rename({'GC_content':'CDS_GC_content'},axis=1)
tr_expression=pd.read_csv('quantification/RSEM_FAANG_transcripts_MAX_fpkm', header=0, sep='\t')


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


s1=m3.loc[m3['homology']==1]['5UTR_GC_content']
x=np.array(s1) 
s2=m3.loc[m3['homology']==0]['5UTR_GC_content']
y=np.array(s2)
#stats.ttest_ind(x,y)[1]  
#stats.wilcoxon(x,y)
stats.mannwhitneyu(x,y)#Mann-Whitney test that is equivalent to Wilcoxon Rank Sum test with sample sizes are uequivalent: https://stats.stackexchange.com/questions/368881/wilcoxon-rank-sum-test-unequal-sample-sizes

su, pval, pval_perm, res =mannwhitneyu(x, y, use_continuity=True, n_rep=10000)

m4=m3.merge(coding_nonNMD_cds)

m4.plot(kind="scatter", x="CDS_GC_content", y="CDS_length", alpha=0.4,
    figsize=(10,7),
    c="homology", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,s=2)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.savefig('protein-coding_transcripts-known.vs.unkown_protein_CDS.png')
plt.close()

x=m4.loc[m4['homology']==1]['CDS_length']
y=m4.loc[m4['homology']==0]['CDS_length']
plt.hist(x, bins=50, alpha=0.5, label='Transcripts encoded homologous proteins')
plt.hist(y, bins=50, alpha=1, label='Transcripts encoded non-homologous proteins')
plt.xlabel('CDS length')
plt.ylabel('Number of transcripts')
plt.legend(loc='upper right')
plt.savefig('CDS_hist_of_protein_coding_nonNMD_transcripts.png')
plt.close()

x=m4.loc[m4['homology']==1]['CDS_GC_content']
y=m4.loc[m4['homology']==0]['CDS_GC_content']
plt.hist(x, bins=50, alpha=0.5, label='Transcripts encoded\nhomologous proteins')
plt.hist(y, bins=50, alpha=0.5, label='Transcripts encoded\nnon-homologous proteins')
plt.xlabel('CDS GC content')
plt.ylabel('Number of transcripts')
plt.legend(loc='upper left')
plt.savefig('CDS-GCcontent_hist_of_protein_coding_nonNMD_transcripts.png')
plt.close()

x=m4.loc[m4['homology']==1]['5UTR_GC_content']
y=m4.loc[m4['homology']==0]['5UTR_GC_content']
plt.hist(x, bins=50, alpha=0.5, label='Transcripts encoded\nhomologous proteins')
plt.hist(y, bins=50, alpha=0.5, label='Transcripts encoded\nnon-homologous proteins')
plt.xlabel('5UTR GC content')
plt.ylabel('Number of transcripts')
plt.legend(loc='upper left')
plt.savefig('UTR5-GCcontent_hist_of_protein_coding_nonNMD_transcripts.png')
plt.close()

x=m4.loc[m4['homology']==1]['3UTR_GC_content']
y=m4.loc[m4['homology']==0]['3UTR_GC_content']
plt.hist(x, bins=50, alpha=0.5, label='Transcripts encoded\nhomologous proteins')
plt.hist(y, bins=50, alpha=0.5, label='Transcripts encoded\nnon-homologous\nproteins')
plt.xlabel('3UTR GC content')
plt.ylabel('Number of transcripts')
plt.legend(loc='upper right')
plt.savefig('UTR3-GCcontent_hist_of_protein_coding_nonNMD_transcripts.png')
plt.close()


m5=m4.merge(tr_expression)
m5.to_csv('protein_coding_transcripts.info',index = None, header=True,sep="\t")

x=m5.loc[m5['homology']==1]['max_expression']
y=m5.loc[m5['homology']==0]['max_expression']
from scipy.stats import median_test
stat, p, med, tbl = median_test(x,y)#Test that two or more samples come from populations with the same median.

remove=pd.read_csv('remove_transcripts',names=['tr_id'])
known_trs=pd.read_csv('annotation/all-known-transcripts-main',names=['tr_id'])
novel_trs=pd.read_csv('annotation/all-novel-transcripts-main',names=['tr_id'])
info=pd.read_csv('protein_coding_transcripts.info',header=0,sep='\t')
known_trs=known_trs.merge(info)
known_trs=known_trs[~ known_trs.tr_id.isin(remove.tr_id)]
novel_trs=novel_trs.merge(info)
novel_trs=novel_trs[~ novel_trs.tr_id.isin(remove.tr_id)]
known_trs.to_csv('known_protein_coding_transcripts.info',index = None, header=True,sep="\t")
novel_trs.to_csv('novel_protein_coding_transcripts.info',index = None, header=True,sep="\t")

stat, p = stats.mannwhitneyu(known_trs['CDS_length'],novel_trs['CDS_length'],alternative='greater')
stat, p = stats.mannwhitneyu(novel_trs.loc[novel_trs['homology']==1]['CDS_length'],novel_trs.loc[novel_trs['homology']==0]['CDS_length'],alternative='greater')

stat, p = stats.mannwhitneyu(known_trs['CDS_GC_content'],novel_trs['CDS_GC_content'],alternative='greater')


known_trs['tr_length']=known_trs['5UTR_length']+known_trs['CDS_length']+known_trs['3UTR_length']
novel_trs['tr_length']=novel_trs['5UTR_length']+novel_trs['CDS_length']+novel_trs['3UTR_length']

a=known_trs['5UTR_length'].to_frame().rename({'5UTR_length':'length'},axis=1)
a['type1']='known'
a['type2']='5UTR'
b=known_trs['CDS_length'].to_frame().rename({'CDS_length':'length'},axis=1)
b['type1']='known'
b['type2']='CDS'
c=known_trs['3UTR_length'].to_frame().rename({'3UTR_length':'length'},axis=1)
c['type1']='known'
c['type2']='3UTR'
d=novel_trs['5UTR_length'].to_frame().rename({'5UTR_length':'length'},axis=1)
d['type1']='novel'
d['type2']='5UTR'
e=novel_trs['CDS_length'].to_frame().rename({'CDS_length':'length'},axis=1)
e['type1']='novel'
e['type2']='CDS'
f=novel_trs['3UTR_length'].to_frame().rename({'3UTR_length':'length'},axis=1)
f['type1']='novel'
f['type2']='3UTR'
df=pd.concat([a,b,c,d,e,f])
ax = sns.boxplot(x="type1", y="length", hue="type2", data=df, linewidth=2.5,showfliers=False)
plt.show()

a=known_trs['5UTR_length'].to_frame().rename({'5UTR_length':'length'},axis=1)
a['type1']='known'
a['type2']='5UTR'
b=known_trs['CDS_length'].to_frame().rename({'CDS_length':'length'},axis=1)
b['type1']='known'
b['type2']='CDS'
c=known_trs['3UTR_length'].to_frame().rename({'3UTR_length':'length'},axis=1)
c['type1']='known'
c['type2']='3UTR'
d=novel_trs.loc[novel_trs['homology']==1]['5UTR_length'].to_frame().rename({'5UTR_length':'length'},axis=1)
d['type1']='novel_homologous'
d['type2']='5UTR'
e=novel_trs.loc[novel_trs['homology']==1]['CDS_length'].to_frame().rename({'CDS_length':'length'},axis=1)
e['type1']='novel_homologous'
e['type2']='CDS'
f=novel_trs.loc[novel_trs['homology']==1]['3UTR_length'].to_frame().rename({'3UTR_length':'length'},axis=1)
f['type1']='novel_homologous'
f['type2']='3UTR'
g=novel_trs.loc[novel_trs['homology']==0]['5UTR_length'].to_frame().rename({'5UTR_length':'length'},axis=1)
g['type1']='novel_Nonhomologous'
g['type2']='5UTR'
h=novel_trs.loc[novel_trs['homology']==0]['CDS_length'].to_frame().rename({'CDS_length':'length'},axis=1)
h['type1']='novel_Nonhomologous'
h['type2']='CDS'
j=novel_trs.loc[novel_trs['homology']==0]['3UTR_length'].to_frame().rename({'3UTR_length':'length'},axis=1)
j['type1']='novel_Nonhomologous'
j['type2']='3UTR'
df=pd.concat([a,b,c,d,e,f,g,h,j])
ax = sns.boxplot(x="type1", y="length", hue="type2", data=df, linewidth=2.5,showfliers=False)
plt.savefig('comparision_of_known_novel_transcript_utr_cds_length.png')
plt.close()

a=known_trs['5UTR_GC_content'].to_frame().rename({'5UTR_GC_content':'GC_content'},axis=1)
a['type1']='known'
a['type2']='5UTR'
b=known_trs['CDS_GC_content'].to_frame().rename({'CDS_GC_content':'GC_content'},axis=1)
b['type1']='known'
b['type2']='CDS'
c=known_trs['3UTR_GC_content'].to_frame().rename({'3UTR_GC_content':'GC_content'},axis=1)
c['type1']='known'
c['type2']='3UTR'
d=novel_trs.loc[novel_trs['homology']==1]['5UTR_GC_content'].to_frame().rename({'5UTR_GC_content':'GC_content'},axis=1)
d['type1']='novel_homologous'
d['type2']='5UTR'
e=novel_trs.loc[novel_trs['homology']==1]['CDS_GC_content'].to_frame().rename({'CDS_GC_content':'GC_content'},axis=1)
e['type1']='novel_homologous'
e['type2']='CDS'
f=novel_trs.loc[novel_trs['homology']==1]['3UTR_GC_content'].to_frame().rename({'3UTR_GC_content':'GC_content'},axis=1)
f['type1']='novel_homologous'
f['type2']='3UTR'
g=novel_trs.loc[novel_trs['homology']==0]['5UTR_GC_content'].to_frame().rename({'5UTR_GC_content':'GC_content'},axis=1)
g['type1']='novel_Nonhomologous'
g['type2']='5UTR'
h=novel_trs.loc[novel_trs['homology']==0]['CDS_GC_content'].to_frame().rename({'CDS_GC_content':'GC_content'},axis=1)
h['type1']='novel_Nonhomologous'
h['type2']='CDS'
j=novel_trs.loc[novel_trs['homology']==0]['3UTR_GC_content'].to_frame().rename({'3UTR_GC_content':'GC_content'},axis=1)
j['type1']='novel_Nonhomologous'
j['type2']='3UTR'
df=pd.concat([a,b,c,d,e,f,g,h,j])
ax = sns.boxplot(x="type1", y="GC_content", hue="type2", data=df, linewidth=2.5,showfliers=False)
plt.savefig('comparision_of_known_novel_transcript_utr_cds_GC-contents.png')
plt.close()