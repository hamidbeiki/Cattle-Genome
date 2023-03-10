#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 13:35:39 2021

@author: beiki
"""

import pandas as pd
import re
import time
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib_venn import venn3,venn2, venn3_circles,venn2_circles,venn2_unweighted

def get_info(ex,biotype):
    ex = ex.replace(0, np.NaN)
    ex_num=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy
    row_medians = np.nanmedian(ex_num, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
    row_counts=np.count_nonzero(~np.isnan(ex_num),axis=1)#https://stackoverflow.com/questions/21778118/counting-the-number-of-non-nan-elements-in-a-numpy-ndarray-in-python
    df1=pd.DataFrame(row_counts).rename({0:'detected_tissues'},axis=1)
    df2=pd.DataFrame(row_medians).rename({0:'median_expression'},axis=1)
    df1[biotype]=ex.index
    df2[biotype]=ex.index
    m=df1.merge(df2)
    m=m[[biotype,'detected_tissues','median_expression']]
    return(m)
    
def get_venn_sections(sets):
    """
    Given a list of sets, return a new list of sets with all the possible
    mutually exclusive overlapping combinations of those sets.  Another way
    to think of this is the mutually exclusive sections of a venn diagram
    of the sets.  If the original list has N sets, the returned list will
    have (2**N)-1 sets.

    Parameters
    ----------
    sets : list of set

    Returns
    -------
    combinations : list of tuple
        tag : str
            Binary string representing which sets are included / excluded in
            the combination.
        set : set
            The set formed by the overlapping input sets.
    """
    num_combinations = 2 ** len(sets)
    bit_flags = [2 ** n for n in range(len(sets))]
    flags_zip_sets = [z for z in zip(bit_flags, sets)]
    combo_sets = []
    for bits in range(num_combinations - 1, 0, -1):
        include_sets = [s for flag, s in flags_zip_sets if bits & flag]
        exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
        combo = set.intersection(*include_sets)
        combo = set.difference(combo, *exclude_sets)
        tag = ''.join([str(int((bits & flag) > 0)) for flag in bit_flags])
        combo_sets.append((tag, combo))
    return combo_sets

def summariser(df,request):
    biotype_info={}
    if request=='transcript':
        info=df[['trA_biotype','trB_biotype']]
    elif request=='gene':
        info=df[['geneA_biotype','geneB_biotype']]
    for i in info.iterrows():
        key1=i[1][0] + '---' + i[1][1]
        key2=i[1][1] + '---' + i[1][0]
        if (key1 not in biotype_info) and (key2 not in biotype_info):
            key=i[1][0] + '---' + i[1][1]
            biotype_info[key]=1
        elif key1 in biotype_info:
            c=biotype_info[key1]
            n=c+1
            biotype_info[key1]=n
        elif key2 in biotype_info:
            c=biotype_info[key2]
            n=c+1
            biotype_info[key2]=n 
    return(biotype_info)
    
def get_divergent_biotypes(df,tissue,request):
    if request=='transcript':
        df_dict=summariser(df,'transcript')
        all_keys=['non_stop_decays---antisense_lncRNAs', 'non_stop_decays---NMDs', 'sncRNAs---non_stop_decays', 'protein_coding_transcripts---NMDs', 'non_stop_decays---non_stop_decays', 'non_stop_decays---protein_coding_transcripts', 'NMDs---NMDs', 'protein_coding_transcripts---sncRNAs', 'protein_coding_transcripts---protein_coding_transcripts', 'antisense_lncRNAs---sncRNAs', 'NMDs---antisense_lncRNAs', 'antisense_lncRNAs---antisense_lncRNAs', 'NMDs---sncRNAs', 'protein_coding_transcripts---antisense_lncRNAs']
        res1=[tissue]
        for key in all_keys:
            s=key.split('---')
            key2=s[1] + '---' + s[0]
            if key in df_dict:
                res1.append(df_dict[key])
            elif key2 in df_dict:
                res1.append(df_dict[key2])
            else:
                res1.append(0) 
    elif request=="gene":
        df_dict=summariser(df,'gene')
        all_keys=['protein-coding---protein-coding','protein-coding---non-coding','protein-coding---pseudogenes','non-coding---pseudogenes']
        res1=[tissue]
        for key in all_keys:
            s=key.split('---')
            key2=s[1] + '---' + s[0]
            if key in df_dict:
                res1.append(df_dict[key])
            elif key2 in df_dict:
                res1.append(df_dict[key2])
            else:
                res1.append(0) 
    return(res1)

def format_tissue_names(l):
    res=[]
    for tissue in l:
        f=tissue.replace('__',' ')
        f=f.lower()
        f=f.replace('_whole','')
        f=f.replace('_',' ')
        f=f.replace('-',' ')
        f=f.replace('mammarygland','mammary gland')
        if f.find('ipsilateral')>=0:
            f='uterine endometrium'
        elif f.find('follicular')>=0:
            f='follicular cells'
        elif f.find('frontal')>=0:
            f='brain frontal cortex' 
        elif f.find('ovarian')>=0:
            f='ovary' 
        elif f.find('longissimus')>=0:
            f='LD muscle'
        elif f.find('lobe')>=0:
            f='lung'
        elif f=='mammary gland ':
            f='virgin mammary gland'
        elif f.find('adipose')>=0:
            f='subcutaneous adipose'
        res.append(f)
    return(res)
                
tr_expressions=pd.read_csv("quantification/samples_RSEM_FAANG_transcripts_fpkm",header=0,sep="\t")
mir_targets=pd.read_csv('miRNA_target_prediction/miranda-PITA_mirna-transcripts3UTR-consensus_predicted_tragets',header=0,sep="\t")
mir_expressions=pd.read_csv('quantification/miRNA_normalised_counts-tissue_samples-duplicated_removed.txt',header=0,sep="\t")
tr_combined_expressions=pd.read_csv("quantification/RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
genes_combined_expressions=pd.read_csv("quantification/RSEM_FAANG_genes_modified_fpkm",header=0,sep="\t")
mir_combined_expressions=pd.read_csv('quantification/miRNA_normalised_counts-combined_samples-duplicates_removed.txt',header=0,sep='\t')
biotypes=pd.read_csv("final_transcript_biotypes",names=['transcript_id','biotype'],sep='\t')
TS_genes=pd.read_csv('quantification/tissue-specific_genes_info', header=0,sep="\t")
TS_transcripts=pd.read_csv('quantification/tissue-specific_transcripts_info',header=0,sep="\t")
genes_TSI=pd.read_csv("Combined_genes-tissue-specificity-scores.txt",header=0,sep="\t").rename({'genes':'gene_id'},axis=1)
tr_TSI=pd.read_csv("Combined_transcripts-tissue-specificity-scores.txt",header=0,sep="\t").reset_index().rename({'index':'transcript_id'},axis=1)
event_per_tr=pd.read_csv('AS-events/total_AS_event_per_transcript',header=0,sep='\t')
event_per_gene=pd.read_csv('AS-events/total_AS_event_per_gene',header=0,sep='\t')
ens_genes=pd.read_csv('annotation/combined_genes_ensembl_equivalent',names=['gene_id','ENS_id'],sep=' ')
trs_cds_length=pd.read_csv('coding_transcripts_longest_CDS_length',names=['transcript_id','CDS_length'],sep='\t')
trs_3utr_info=pd.read_csv('3UTRs_length_GCcontent',header=0,sep='\t').rename({'tr_id':'transcript_id'},axis=1) 

                  
tr_combined_expressions.columns=tr_combined_expressions.columns.str.replace('.', '-')
genes_combined_expressions.columns=genes_combined_expressions.columns.str.replace('.', '-')
tr_expressions.columns=tr_expressions.columns.str.replace('.', '-')
mir_combined_expressions.columns=mir_combined_expressions.columns.str.replace('.', '-')
mir_expressions.columns=mir_expressions.columns.str.replace('_sample', '')

a=set(tr_expressions.columns)
b=set(mir_expressions.columns)
shared_cols=list(a & b)
shared_cols.sort()

t=[re.sub(r'\_[1-9]', '', x) for x in shared_cols]
t=pd.DataFrame(t).rename({0:'tissue'},axis=1)
shared_tissue_sample_info=t.groupby(['tissue']).size().reset_index(name='number_of_samples')
useful_tissues=shared_tissue_sample_info.loc[shared_tissue_sample_info['number_of_samples']>=3]

useful_tissues=list(useful_tissues['tissue'])

useful_samples=[]
for tissue in useful_tissues:
    s=list(filter(lambda x: tissue in x, shared_cols))
    for i in s:
        useful_samples.append(i)

s=useful_samples.copy()
s=['transcript_id'] + s
tr_expressions=tr_expressions[s]
tr_expressions=tr_expressions.set_index('transcript_id')
s=useful_samples.copy()
s=['miRNA'] + s
mir_expressions=mir_expressions[s]
mir_expressions=mir_expressions.set_index('miRNA')
tr_combined_expressions=tr_combined_expressions.set_index('transcript_id')
mir_combined_expressions=mir_combined_expressions.set_index('miRNA')
genes_combined_expressions=genes_combined_expressions.set_index('gene_id')

transcripts_info=get_info(tr_combined_expressions,'transcript_id')
genes_info=get_info(genes_combined_expressions,'gene_id')
e=tr_combined_expressions[useful_tissues]
e['sum']=e.sum(axis=1)
e_t=e.loc[e['sum']>0]
del e_t['sum']
transcripts_info2=get_info(e_t,'transcript_id')
e_g=genes_combined_expressions[useful_tissues]
genes_info2=get_info(e_g,'gene_id')


interaction_corr_info={}
start_time = time.time()
for tissue in useful_tissues:
    info=pd.DataFrame(columns=['miRNA','transcript_id','correlations'])
    tr_ex=tr_combined_expressions[tissue].to_frame()
    tr_ex=tr_ex.loc[tr_ex[tissue]>0]
    trs=list(tr_ex.index)
    mir_ex=mir_combined_expressions[tissue].to_frame()
    mir_ex=mir_ex.loc[mir_ex[tissue]>0]
    mirs=list(mir_ex.index)
    a=mir_targets.loc[mir_targets.miRNA.isin(mirs)]
    intractions=a.loc[a.transcript_id.isin(trs)]
    s=list(filter(lambda x: tissue in x, shared_cols))
    df1=tr_expressions[s]
    df2=mir_expressions[s]
    tr_s_ex=df1.loc[df1.index.isin(intractions.transcript_id)]
    mir_s_ex=df2.loc[df2.index.isin(intractions.miRNA)]
    df=pd.concat([tr_s_ex,mir_s_ex])
    cor_mat=df.T.corr()#method='spearman'
    groups=intractions.groupby(['miRNA'])
    for g,v in groups:
        l=list(v['transcript_id'])
        c=cor_mat[g]
        d=cor_mat[g].to_frame()
        d=d.loc[d.index.isin(l)].rename({g:'correlations'},axis=1)
        d=d.loc[(d['correlations']>=0.6) | (d['correlations']<=-0.6)]#all_df created withoud this line
        d=d.reset_index().rename({'index':'transcript_id'},axis=1)
        d['miRNA']=g
        d['tissue']=tissue
        d=d[['tissue','miRNA','transcript_id','correlations']]
        info=pd.concat([info,d])
    interaction_corr_info[tissue]=info
    t=time.time() - start_time
    print([tissue,t])

flag=1
for key,value in interaction_corr_info.items():
    globals() ['df' + str(flag)]=value.copy()
    print(flag)
    flag=flag+1    
df=pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13])
df.to_csv('miRNA_target_prediction/miRNA_transcripts_targets_significant_correlations',index = None, header=True,sep="\t")
all_df=pd.read_csv('miRNA_target_prediction/miRNA_transcripts_targets_correlations',header=0,sep='\t')

df['id']=df['miRNA'] + '___' + df['transcript_id']
df['id2']=df['id'] + '___' + df['tissue']
intraction_replication=df.groupby(['id']).size().reset_index(name='replicated_tissues')
intraction_replication.groupby(['replicated_tissues']).size().reset_index(name='count')

pos=df.loc[df['correlations']>0]
neg=df.loc[df['correlations']<0]
pos_mirs=pos['miRNA'].drop_duplicates().to_frame()
neg_mirs=neg['miRNA'].drop_duplicates().to_frame()
pos_trs=pos['transcript_id'].drop_duplicates().to_frame()
neg_trs=neg['transcript_id'].drop_duplicates().to_frame()
i=pos['transcript_id'].str.split('.',2, expand=True)
pos_genes=(i[0] + '.' + i[1]).drop_duplicates().to_frame().rename({0:'gene_id'},axis=1)
i=neg['transcript_id'].str.split('.',2, expand=True)
neg_genes=(i[0] + '.' + i[1]).drop_duplicates().to_frame().rename({0:'gene_id'},axis=1)

plt.figure(figsize=(24,12),dpi=300)  
plt.subplot(2, 2, 1)# row 2, col 2 index 1
set1=set(pos['id'])
set2=set(neg['id'])
v=venn2([ set1, set2], ('Positive miRNA-transcript intractions', 'Negative miRNA-transcript intractions'))
c=venn2_circles([set1, set2], linestyle='dashed', color="grey")
plt.subplot(2, 2, 2)# row 2, col 2 index 2
set1=set(pos_mirs['miRNA'])
set2=set(neg_mirs['miRNA'])
v=venn2_unweighted([ set1, set2], ('miRNAs positively regulate their target transcripts', 'miRNAs negatively regulate their target transcripts'),set_colors=('purple', 'skyblue'))
plt.subplot(2, 2, 3)# row 2, col 2 index 3
set1=set(pos_trs['transcript_id'])
set2=set(neg_trs['transcript_id'])
v=venn2_unweighted([ set1, set2], ('Transcripts positively regulated by miRNAs', 'Transcripts negatively regulated by miRNAs'),set_colors=('c', 'pink'))
plt.subplot(2, 2, 4)# row 2, col 2 index 4
set1=set(pos_genes['gene_id'])
set2=set(neg_genes['gene_id'])
v=venn2_unweighted([ set1, set2], ('Genes positively regulated by miRNAs', 'Genes negatively regulated by miRNAs'),set_colors=('orange', 'greenyellow'))
plt.savefig('summary_of_positive_negative_miRNA-transcript_intractions.jpg',dpi=300)
plt.close()

sets=[set(pos_trs['transcript_id']),set(neg_trs['transcript_id'])]
ints=get_venn_sections(sets)
''' to get intersect classes '''
#for i in range(len(ints)):
#    print(i, ints[i][0])
#0 11
#1 01
#2 10
pos_only_trs=pd.DataFrame(list(ints[2][1])).rename({0:'transcript_id'},axis=1).merge(tr_TSI)
neg_only_trs=pd.DataFrame(list(ints[1][1])).rename({0:'transcript_id'},axis=1).merge(tr_TSI)
rest=tr_TSI.loc[(~tr_TSI.transcript_id.isin(pos_only_trs.transcript_id)) & (~tr_TSI.transcript_id.isin(neg_only_trs.transcript_id))]
stat, p1 = stats.mannwhitneyu(pos_only_trs['TSI'],rest['TSI'],alternative='greater')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_trs['TSI'],rest['TSI'],alternative='greater')#,alternative='greater'
pos_only_trs_AS=pos_only_trs.merge(event_per_tr)
neg_only_trs_AS=neg_only_trs.merge(event_per_tr)
rest=event_per_tr.loc[(~event_per_tr.transcript_id.isin(pos_only_trs.transcript_id)) & (~event_per_tr.transcript_id.isin(neg_only_trs.transcript_id))]
stat, p1 = stats.mannwhitneyu(pos_only_trs_AS['total_AS_event'],rest['total_AS_event'],alternative='less')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_trs_AS['total_AS_event'],rest['total_AS_event'],alternative='less')#,alternative='greater'
pos_only_trs_exp=pos_only_trs.merge(transcripts_info)
neg_only_trs_exp=neg_only_trs.merge(transcripts_info)
rest=transcripts_info2.loc[(~transcripts_info.transcript_id.isin(pos_only_trs.transcript_id)) & (~transcripts_info.transcript_id.isin(neg_only_trs.transcript_id))]
stat, p1 = stats.mannwhitneyu(pos_only_trs_exp['detected_tissues'],rest['detected_tissues'],alternative='less')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_trs_exp['detected_tissues'],rest['detected_tissues'],alternative='less')#,alternative='greater'
stat, p3 = stats.mannwhitneyu(pos_only_trs_exp['median_expression'],rest['median_expression'],alternative='less')#,alternative='greater'
stat, p4 = stats.mannwhitneyu(neg_only_trs_exp['median_expression'],rest['median_expression'],alternative='less')#,alternative='greater'


sets=[set(pos_genes['gene_id']),set(neg_genes['gene_id'])]
ints=get_venn_sections(sets)
pos_only_genes=pd.DataFrame(list(ints[2][1])).rename({0:'gene_id'},axis=1).merge(genes_TSI)
neg_only_genes=pd.DataFrame(list(ints[1][1])).rename({0:'gene_id'},axis=1).merge(genes_TSI)
rest=genes_TSI.loc[(~genes_TSI.gene_id.isin(pos_only_genes.gene_id)) & (~genes_TSI.gene_id.isin(neg_only_genes.gene_id))]
stat, p1 = stats.mannwhitneyu(pos_only_genes['TSI'],rest['TSI'],alternative='greater')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_genes['TSI'],rest['TSI'],alternative='greater')#,alternative='greater'
pos_only_genes_AS=pos_only_genes.merge(event_per_gene)
neg_only_genes_AS=neg_only_genes.merge(event_per_gene)
rest=event_per_gene.loc[(~event_per_gene.gene_id.isin(pos_only_genes.gene_id)) & (~event_per_gene.gene_id.isin(neg_only_genes.gene_id))]
stat, p1 = stats.mannwhitneyu(pos_only_genes_AS['total_AS_event'],rest['total_AS_event'],alternative='less')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_genes_AS['total_AS_event'],rest['total_AS_event'],alternative='less')#,alternative='greater'
pos_only_genes_exp=pos_only_genes.merge(genes_info)
neg_only_genes_exp=neg_only_genes.merge(genes_info)
rest=genes_info.loc[(~genes_info.gene_id.isin(pos_only_genes.gene_id)) & (~genes_info.gene_id.isin(neg_only_genes.gene_id))]
stat, p1 = stats.mannwhitneyu(pos_only_trs_exp['detected_tissues'],rest['detected_tissues'],alternative='less')#,alternative='greater'
stat, p2 = stats.mannwhitneyu(neg_only_trs_exp['detected_tissues'],rest['detected_tissues'],alternative='less')#,alternative='greater'
stat, p3 = stats.mannwhitneyu(pos_only_trs_exp['median_expression'],rest['median_expression'],alternative='less')#,alternative='greater'
stat, p4 = stats.mannwhitneyu(neg_only_trs_exp['median_expression'],rest['median_expression'],alternative='less')#,alternative='greater'

pos_only_genes_ens=pos_only_genes.merge(ens_genes)['ENS_id'].drop_duplicates().to_frame()
neg_only_genes_ens=neg_only_genes.merge(ens_genes)['ENS_id'].drop_duplicates().to_frame()


pos_only_genes_ens.to_csv('miRNA_target_prediction/genes_only_positively_regulated_by_miRNAs',index = None, header=False,sep="\t")
neg_only_genes_ens.to_csv('miRNA_target_prediction/genes_only_negatively_regulated_by_miRNAs',index = None, header=False,sep="\t")

plt.figure(figsize=(12,14),dpi=300)  
plt.subplot(2, 2, 1)# row 2, col 2 index 1
plt.hist(pos_only_trs_exp['detected_tissues'], bins=20,color='red',title='Transcripts just positively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of transcripts')
plt.title('Transcripts just positively \n regulated by miRNAs')
plt.subplot(2, 2, 2)
plt.hist(neg_only_trs_exp['detected_tissues'], bins=20,color='blue',label='Transcripts just negatively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of transcripts')
plt.title('Transcripts just negatively \n regulated by miRNAs')
plt.subplot(2, 2, 3)# row 2, col 2 index 1
plt.hist(pos_only_genes_exp['detected_tissues'], bins=20,color='darkred',label='Genes just positively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of genes')
plt.title('Genes just positively \n regulated by miRNAs')
plt.subplot(2, 2, 4)
plt.hist(neg_only_genes_exp['detected_tissues'], bins=20,color='darkblue',label='Genes just negatively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of genes')
plt.title('Genes just negatively \n regulated by miRNAs')
plt.savefig('genes_transcripts_either_positively_or_negatively_regulated_by_miRNAs_all_tissues.jpg',dpi=300)
plt.close()


pos_only_trs_exp2=pos_only_trs.merge(transcripts_info2)
neg_only_trs_exp2=neg_only_trs.merge(transcripts_info2)
pos_only_genes_exp2=pos_only_genes.merge(genes_info2)
neg_only_genes_exp2=neg_only_genes.merge(genes_info2)
plt.figure(figsize=(12,14),dpi=300)  
plt.subplot(2, 2, 1)# row 2, col 2 index 1
plt.hist(pos_only_trs_exp2['detected_tissues'], bins=20,color='red',title='Transcripts just positively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of transcripts')
plt.title('Transcripts just positively \n regulated by miRNAs')
plt.subplot(2, 2, 2)
plt.hist(neg_only_trs_exp2['detected_tissues'], bins=20,color='blue',label='Transcripts just negatively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of transcripts')
plt.title('Transcripts just negatively \n regulated by miRNAs')
plt.subplot(2, 2, 3)# row 2, col 2 index 1
plt.hist(pos_only_genes_exp2['detected_tissues'], bins=20,color='darkred',label='Genes just positively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of genes')
plt.title('Genes just positively \n regulated by miRNAs')
plt.subplot(2, 2, 4)
plt.hist(neg_only_genes_exp2['detected_tissues'], bins=20,color='darkblue',label='Genes just negatively \n regulated by miRNAs')#, weights=np.ones(len(x)) / len(x),alpha=0.5,label='+'
plt.xlabel('Number of detected tissues')
plt.ylabel('Number of genes')
plt.title('Genes just negatively \n regulated by miRNAs')
plt.savefig('genes_transcripts_either_positively_or_negatively_regulated_by_miRNAs_used_tissues.jpg',dpi=300)
plt.close()


pos_intraction_replication=pos.groupby(['id']).size().reset_index(name='POSreplicated_tissues')
neg_intraction_replication=neg.groupby(['id']).size().reset_index(name='NEGreplicated_tissues')

pos_intraction_replication['POSreplicated_tissues'].mean()
neg_intraction_replication['NEGreplicated_tissues'].mean()
stat, p = stats.mannwhitneyu(pos_intraction_replication['POSreplicated_tissues'],neg_intraction_replication['NEGreplicated_tissues'])#,alternative='greater'


all_df['id']=all_df['miRNA'] + '___' + all_df['transcript_id']
all_df['id2']=all_df['id'] + '___' + all_df['tissue']
detected_tissues=all_df.groupby(['id']).size().reset_index(name='detected_tissues')
intraction_detection=intraction_replication.merge(detected_tissues).rename({'replicated_tissues':'interacted_tissues'},axis=1)
intraction_detection.loc[intraction_detection['detected_tissues']>intraction_detection['interacted_tissues']].shape
intraction_detection.loc[intraction_detection['detected_tissues']==intraction_detection['interacted_tissues']].shape

interacted=df.copy()
not_interacted=all_df.loc[~all_df.id2.isin(interacted.id2)]

interacted['id'].drop_duplicates().shape[0]
not_interacted['id'].drop_duplicates().shape[0]
interacted['correlations'].median()
not_interacted['correlations'].median()

plt.figure(figsize=(18,6),dpi=300)  
plt.subplot(1, 2, 1)# row 1, col 2 index 1
a=pos['correlations'].to_frame().rename({'correlations':'info'},axis=1)
a['type1']='Significant positive miRNA-transcript \n correlation coeffiecients'
b=neg['correlations'].to_frame().rename({'correlations':'info'},axis=1)
b['type1']='Significant negative miRNA-transcript \n correlation coeffiecients'
c=not_interacted['correlations'].to_frame().rename({'correlations':'info'},axis=1)
c['type1']='Non-significant miRNA-transcript \n correlation coeffiecients'
df=pd.concat([a,b])
ax = sns.boxplot(x="type1", y="info", data=df, linewidth=0.6,showfliers=False)#, hue="type2"
plt.rc('xtick', labelsize=10)    # 8 fontsize of the tick labels
plt.rc('ytick', labelsize=10) 
plt.ylabel("Significant correlation coefficients \n between pairs of miRNA-transcript ",fontsize=10)
plt.subplot(1, 2, 2)# row 1, col 2 index 1
a=pos_intraction_replication['POSreplicated_tissues'].to_frame().rename({'POSreplicated_tissues':'info'},axis=1)
a['type1']='Positive miRNA-transcript \n correlation coeffiecients'
b=neg_intraction_replication['NEGreplicated_tissues'].to_frame().rename({'NEGreplicated_tissues':'info'},axis=1)
b['type1']='Negative miRNA-transcript \n correlation coeffiecients'
df=pd.concat([a,b])
ax = sns.boxplot(x="type1", y="info", data=df, linewidth=0.6,showfliers=False)#, hue="type2"
plt.rc('xtick', labelsize=10)    # 8 fontsize of the tick labels
plt.rc('ytick', labelsize=10) 
plt.ylabel("Number of tissues that each miRNA-transcript \n interaction has been detected",fontsize=10)
plt.savefig('positive_negative_miRNA-transcript_interactions.jpg',dpi=300)
plt.close()


'''
transcript biotype enrichment
'''
info=biotypes.groupby(['biotype']).size().reset_index(name='count')
interacted_trs=interacted['transcript_id'].drop_duplicates().to_frame()
interacted_trs=interacted_trs.merge(biotypes)
rest_trs=biotypes.loc[~biotypes.transcript_id.isin(interacted_trs.transcript_id)]
l=list(biotypes['biotype'].drop_duplicates())
enriched_biotypes=[]#enriched biotypes in transcripts interacted with miRNAs
for i in l:
    A=interacted_trs.loc[interacted_trs['biotype']==i].shape[0]
    C=interacted_trs.shape[0]-A
    B=info.loc[info['biotype']==i].iloc[0,1]-A
    D=biotypes.shape[0]
    obs = np.array([[A,B], [C,D]])
    oddsratio,pvalue = stats.fisher_exact(obs,alternative='greater')
    if pvalue<0.05:
        enriched_biotypes.append([i,pvalue])
#[['protein_coding_transcripts', 0.0], ['long_intergenic_lncRNAs', 0.0], ['long_intragenic_lncRNAs', 0.0], ['sense_intronic_lncRNAs', 4.6865777281592985e-163], ['antisense_lncRNAs', 0.0]] 
'''
itegration with divergent transcription:
'''
bidirection_single=[]
bidirection_both=[]
convergent_single=[]
convergent_both=[]
mir_divergent=[]
ppinfo=[]# comparision of coding trs targeted with miRNAs and those that were not (within coding-coding pairs)
for tissue in useful_tissues:
    bid=pd.read_csv(tissue + '_bidirection_promoter_events',header=0,sep='\t')
    del bid['trA_biotype']
    del bid['trB_biotype']
    bid=bid.rename({'trA':'transcript_id'},axis=1)
    m=bid.merge(biotypes)
    m=m.rename({'transcript_id':'trA','biotype':'trA_biotype','trB':'transcript_id'},axis=1)
    m=m.merge(biotypes)
    bid=m.rename({'transcript_id':'trB','biotype':'trB_biotype'},axis=1)
    conv=pd.read_csv(tissue + '_convergent_gene_events',header=0,sep='\t')
    del conv['trA_biotype']
    del conv['trB_biotype']
    conv=conv.rename({'trA':'transcript_id'},axis=1)
    m=conv.merge(biotypes)
    m=m.rename({'transcript_id':'trA','biotype':'trA_biotype','trB':'transcript_id'},axis=1)
    m=m.merge(biotypes)
    conv=m.rename({'transcript_id':'trB','biotype':'trB_biotype'},axis=1)
    trs=pd.read_csv('quantification/' + tissue + '_transcripts',names=['transcript_id'],sep='\t')
    a=bid['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b=bid['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    bid_trs=pd.concat([a,b]).drop_duplicates()
    a=conv['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b=conv['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    conv_trs=pd.concat([a,b]).drop_duplicates()
    t_intracted=interacted.loc[interacted['tissue']==tissue]['transcript_id'].drop_duplicates().to_frame()
    bid_interacted=t_intracted.loc[t_intracted.transcript_id.isin(bid_trs.transcript_id)]
    A=bid_interacted.shape[0]
    B=t_intracted.shape[0]-A
    C=bid_trs.shape[0]-A
    D=trs.shape[0]
    obs = np.array([[A,B], [C,D]])
    oddsratio,bid_pvalue = stats.fisher_exact(obs,alternative='greater')# the smallest pvalue that can be distinguished from zero in python is 7.1429666685167604e-293 :https://stackoverflow.com/questions/20530138/scipy-p-value-returns-0-0 
    conv_interacted=t_intracted.loc[t_intracted.transcript_id.isin(conv_trs.transcript_id)]
    A=conv_interacted.shape[0]
    B=t_intracted.shape[0]-A
    C=conv_trs.shape[0]-A
    D=trs.shape[0]
    obs = np.array([[A,B], [C,D]])
    oddsratio,conv_pvalue = stats.fisher_exact(obs,alternative='greater')# the smallest pvalue that can be distinguished from zero in python is 7.1429666685167604e-293 :https://stackoverflow.com/questions/20530138/scipy-p-value-returns-0-0 
    a=bid.loc[(bid.trA.isin(t_intracted.transcript_id)) & (~bid.trB.isin(t_intracted.transcript_id))]
    b=bid.loc[(~bid.trA.isin(t_intracted.transcript_id)) & (bid.trB.isin(t_intracted.transcript_id))]
    a2=a['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b2=b['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    bid_single_targeted_trs=pd.concat([a2,b2]).drop_duplicates().merge(biotypes)
    bid_single=pd.concat([a,b])
    bid_both=bid.loc[(bid.trA.isin(t_intracted.transcript_id)) & (bid.trB.isin(t_intracted.transcript_id))]
    a=bid_both['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b=bid_both['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    bid_both_targeted_trs=pd.concat([a,b]).drop_duplicates().merge(biotypes)
    a=conv.loc[(conv.trA.isin(t_intracted.transcript_id)) & (~conv.trB.isin(t_intracted.transcript_id))]
    b=conv.loc[(~conv.trA.isin(t_intracted.transcript_id)) & (conv.trB.isin(t_intracted.transcript_id))]
    a2=a['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b2=b['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    conv_single_targeted_trs=pd.concat([a2,b2]).drop_duplicates().merge(biotypes)
    conv_single=pd.concat([a,b])
    conv_both=conv.loc[(conv.trA.isin(t_intracted.transcript_id)) & (conv.trB.isin(t_intracted.transcript_id))]    
    a=conv_both['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b=conv_both['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    conv_both_targeted_trs=pd.concat([a,b]).drop_duplicates().merge(biotypes)
    bid_single_info=get_divergent_biotypes(bid_single,tissue,'transcript')
    bid_both_info=get_divergent_biotypes(bid_both,tissue,'transcript')
    conv_single_info=get_divergent_biotypes(conv_single,tissue,'transcript')
    conv_both_info=get_divergent_biotypes(conv_both,tissue,'transcript') 
    bid_single_p=(bid_single.shape[0]/bid.shape[0])*100
    bid_both_p=(bid_both.shape[0]/bid.shape[0])*100
    bid_p=bid_single_p+bid_both_p
    conv_single_p=(conv_single.shape[0]/conv.shape[0])*100
    conv_both_p=(conv_both.shape[0]/conv.shape[0])*100
    conv_p=conv_single_p+conv_both_p
    mir_divergent.append([tissue,bid_pvalue,conv_pvalue,bid_p,bid_single_p,bid_both_p,conv_p,conv_single_p,conv_both_p])
    bidirection_single.append(bid_single_info)
    bidirection_both.append(bid_both_info)
    convergent_single.append(conv_single_info)
    convergent_both.append(conv_both_info)    
    bid_single_pp=bid_single.loc[(bid_single['trA_biotype']=='protein_coding_transcripts') & (bid_single['trB_biotype']=='protein_coding_transcripts')]
    a=bid_single_pp.loc[(bid_single_pp.trA.isin(t_intracted.transcript_id)) & (~bid_single_pp.trB.isin(t_intracted.transcript_id))]
    b=bid_single_pp.loc[(~bid_single_pp.trA.isin(t_intracted.transcript_id)) & (bid_single_pp.trB.isin(t_intracted.transcript_id))]
    a2=a['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b2=b['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    bid_single_t=pd.concat([a2,b2]).drop_duplicates().merge(trs_cds_length).merge(trs_3utr_info)
    a2=a['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    b2=b['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    bid_single_r=pd.concat([a2,b2]).drop_duplicates().merge(trs_cds_length).merge(trs_3utr_info)
    conv_single_pp=conv_single.loc[(conv_single['trA_biotype']=='protein_coding_transcripts') & (conv_single['trB_biotype']=='protein_coding_transcripts')]
    a=conv_single_pp.loc[(conv_single_pp.trA.isin(t_intracted.transcript_id)) & (~conv_single_pp.trB.isin(t_intracted.transcript_id))]
    b=conv_single_pp.loc[(~conv_single_pp.trA.isin(t_intracted.transcript_id)) & (conv_single_pp.trB.isin(t_intracted.transcript_id))]
    a2=a['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    b2=b['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    conv_single_t=pd.concat([a2,b2]).drop_duplicates().merge(trs_cds_length).merge(trs_3utr_info)
    a2=a['trB'].to_frame().rename({'trB':'transcript_id'},axis=1)
    b2=b['trA'].to_frame().rename({'trA':'transcript_id'},axis=1)
    conv_single_r=pd.concat([a2,b2]).drop_duplicates().merge(trs_cds_length).merge(trs_3utr_info)
    stat, p1 = stats.mannwhitneyu(bid_single_t['CDS_length'],bid_single_r['CDS_length'],alternative='greater')#,alternative='greater'
    stat, p2 = stats.mannwhitneyu(bid_single_t['UTR_length'],bid_single_r['UTR_length'],alternative='greater')
    stat, p3 = stats.mannwhitneyu(bid_single_t['GC_content'],bid_single_r['GC_content'])
    stat, p4 = stats.mannwhitneyu(conv_single_t['CDS_length'],conv_single_r['CDS_length'],alternative='greater')#,alternative='greater'
    stat, p5 = stats.mannwhitneyu(conv_single_t['UTR_length'],conv_single_r['UTR_length'],alternative='greater')
    stat, p6 = stats.mannwhitneyu(conv_single_t['GC_content'],conv_single_r['GC_content'])
    ppinfo.append([tissue,p1,p2,p3,p4,p5,p6,bid_single_t['CDS_length'].median(),bid_single_t['UTR_length'].median(),bid_single_t['GC_content'].median(),bid_single_r['CDS_length'].median(),bid_single_r['UTR_length'].median(),bid_single_r['GC_content'].median(),conv_single_t['CDS_length'].median(),conv_single_t['UTR_length'].median(),conv_single_t['GC_content'].median(),conv_single_r['CDS_length'].median(),conv_single_r['UTR_length'].median(),conv_single_r['GC_content'].median()])
    print(tissue)    
    
ppinfo=pd.DataFrame(ppinfo).rename({0:'tissue',1:'bid_CDS_length_pvalue',2:'bid_UTR_length_pvalue',3:'bid_UTR_GC_pvalue',4:'conv_CDS_length_pvalue',5:'conv_UTR_length_pvalue',6:'conv_UTR_GC_pvalue',7:'targeted_bid_CDS_length',8:'targeted_bid_UTR_length',9:'targeted_bid_UTR_GC',10:'rest_bid_CDS_length',11:'rest_bid_UTR_length',12:'rest_bid_UTR_GC',13:'targeted_conv_CDS_length',14:'targeted_conv_UTR_length',15:'targeted_conv_UTR_GC',16:'rest_conv_CDS_length',17:'rest_conv_UTR_length',18:'rest_conv_UTR_GC'},axis=1)

mir_divergent=pd.DataFrame(mir_divergent).rename({0:'tissue',1:'miRNA_int_enrichment_bidirection',2:'miRNA_int_enrichment_convergent',3:'%_of_bidirection_targeted_by_miRNAs',4:'%_of_bidirection_ONLY_one_transcript_targeted_by_miRNAs',5:'%_of_bidirection_BOTH_transcripts_targeted_by_miRNAs',6:'%_of_convergent_targeted_by_miRNAs',7:'%_of_convergent_ONLY_one_transcript_targeted_by_miRNAs',8:'%_of_convergent_BOTH_transcripts_targeted_by_miRNAs'},axis=1)
bidirection_single=pd.DataFrame(bidirection_single).rename({0:'tissue',1:'non_stop_decays---antisense_lncRNAs', 2:'non_stop_decays---NMDs', 3:'sncRNAs---non_stop_decays', 4:'protein_coding_transcripts---NMDs', 5:'non_stop_decays---non_stop_decays', 6:'non_stop_decays---protein_coding_transcripts', 7:'NMDs---NMDs', 8:'protein_coding_transcripts---sncRNAs', 9:'protein_coding_transcripts---protein_coding_transcripts', 10:'antisense_lncRNAs---sncRNAs', 11:'NMDs---antisense_lncRNAs', 12:'antisense_lncRNAs---antisense_lncRNAs', 13:'NMDs---sncRNAs', 14:'protein_coding_transcripts---antisense_lncRNAs'},axis=1)
bidirection_both=pd.DataFrame(bidirection_both).rename({0:'tissue',1:'non_stop_decays---antisense_lncRNAs', 2:'non_stop_decays---NMDs', 3:'sncRNAs---non_stop_decays', 4:'protein_coding_transcripts---NMDs', 5:'non_stop_decays---non_stop_decays', 6:'non_stop_decays---protein_coding_transcripts', 7:'NMDs---NMDs', 8:'protein_coding_transcripts---sncRNAs', 9:'protein_coding_transcripts---protein_coding_transcripts', 10:'antisense_lncRNAs---sncRNAs', 11:'NMDs---antisense_lncRNAs', 12:'antisense_lncRNAs---antisense_lncRNAs', 13:'NMDs---sncRNAs', 14:'protein_coding_transcripts---antisense_lncRNAs'},axis=1)
convergent_single=pd.DataFrame(convergent_single).rename({0:'tissue',1:'non_stop_decays---antisense_lncRNAs', 2:'non_stop_decays---NMDs', 3:'sncRNAs---non_stop_decays', 4:'protein_coding_transcripts---NMDs', 5:'non_stop_decays---non_stop_decays', 6:'non_stop_decays---protein_coding_transcripts', 7:'NMDs---NMDs', 8:'protein_coding_transcripts---sncRNAs', 9:'protein_coding_transcripts---protein_coding_transcripts', 10:'antisense_lncRNAs---sncRNAs', 11:'NMDs---antisense_lncRNAs', 12:'antisense_lncRNAs---antisense_lncRNAs', 13:'NMDs---sncRNAs', 14:'protein_coding_transcripts---antisense_lncRNAs'},axis=1)
convergent_both=pd.DataFrame(convergent_both).rename({0:'tissue',1:'non_stop_decays---antisense_lncRNAs', 2:'non_stop_decays---NMDs', 3:'sncRNAs---non_stop_decays', 4:'protein_coding_transcripts---NMDs', 5:'non_stop_decays---non_stop_decays', 6:'non_stop_decays---protein_coding_transcripts', 7:'NMDs---NMDs', 8:'protein_coding_transcripts---sncRNAs', 9:'protein_coding_transcripts---protein_coding_transcripts', 10:'antisense_lncRNAs---sncRNAs', 11:'NMDs---antisense_lncRNAs', 12:'antisense_lncRNAs---antisense_lncRNAs', 13:'NMDs---sncRNAs', 14:'protein_coding_transcripts---antisense_lncRNAs'},axis=1)

bidirection_single['type']='bidirection_single'
bidirection_both['type']='bidirection_both'
convergent_single['type']='convergent_single'
convergent_both['type']='convergent_both'

divergent_miRNA_int=pd.concat([bidirection_single,bidirection_both,convergent_single,convergent_both])

mir_divergent.to_csv('miRNA_target_prediction/relation_between_divergent_transcription_and_miRNA_regulation',index = None, header=True,sep="\t")
divergent_miRNA_int.to_csv('miRNA_target_prediction/summary_of_divergent_transcription_intraction_with_miRNAs',index = None, header=True,sep="\t")
ppinfo.to_csv('miRNA_target_prediction/miRNA-targeted_vs_not-targeted_coding_tanscripts_in_coding_diveregent_pairs',index = None, header=True,sep="\t")


t=format_tissue_names(list(mir_divergent['tissue']))
bid_o=list(mir_divergent['%_of_bidirection_ONLY_one_transcript_targeted_by_miRNAs'])
bid_b=list(mir_divergent['%_of_bidirection_BOTH_transcripts_targeted_by_miRNAs'])
conv_o=list(mir_divergent['%_of_convergent_ONLY_one_transcript_targeted_by_miRNAs'])
conv_b=list(mir_divergent['%_of_convergent_BOTH_transcripts_targeted_by_miRNAs'])

N = 13
ind = np.arange(N) 
width = 0.1
bar1 = plt.bar(ind, bid_o, width, color = 'gold')
bar2 = plt.bar(ind+width, bid_b, width, color='red')
bar3 = plt.bar(ind+width*2, conv_o, width, color='cyan')
bar4 = plt.bar(ind+width*3, conv_b, width, color='navy')
plt.xlabel("Tissues")
plt.ylabel('Percentage of head2head/tail2tail organized \n transcripts tageted by miRNAs')
#plt.title("Title")
plt.xticks(ind+width,t)
plt.legend( (bar1, bar2, bar3,bar4), ('one of head2head transcripts targeted by miRNAs', 'both of head2head transcripts targeted by miRNAs', 'one of tail2tail transcripts targeted by miRNAs','both of tail2tail transcripts targeted by miRNAs') )
plt.show()
#relation_between_divergent_transcription_and_miRNA_regulatio.jpg


plt.figure(figsize=(18,18),dpi=300)  
plt.subplot(3, 2, 1)# row 3, col 2 figure 1
data=divergent_miRNA_int.loc[divergent_miRNA_int['type']=='bidirection_single']
a=data['non_stop_decays---antisense_lncRNAs'].to_frame().rename({'non_stop_decays---antisense_lncRNAs':'info'},axis=1)
a['type']=1
b=data['non_stop_decays---NMDs'].to_frame().rename({'non_stop_decays---NMDs':'info'},axis=1)
b['type']=2
c=data['sncRNAs---non_stop_decays'].to_frame().rename({'sncRNAs---non_stop_decays':'info'},axis=1)
c['type']=3
d=data['protein_coding_transcripts---NMDs'].to_frame().rename({'protein_coding_transcripts---NMDs':'info'},axis=1)
d['type']=4
e=data['non_stop_decays---non_stop_decays'].to_frame().rename({'non_stop_decays---non_stop_decays':'info'},axis=1)
e['type']=5
f=data['non_stop_decays---protein_coding_transcripts'].to_frame().rename({'non_stop_decays---protein_coding_transcripts':'info'},axis=1)
f['type']=6
g=data['NMDs---NMDs'].to_frame().rename({'NMDs---NMDs':'info'},axis=1)
g['type']=7
h=data['protein_coding_transcripts---sncRNAs'].to_frame().rename({'protein_coding_transcripts---sncRNAs':'info'},axis=1)
h['type']=8
j=data['protein_coding_transcripts---protein_coding_transcripts'].to_frame().rename({'protein_coding_transcripts---protein_coding_transcripts':'info'},axis=1)
j['type']=9
k=data['antisense_lncRNAs---sncRNAs'].to_frame().rename({'antisense_lncRNAs---sncRNAs':'info'},axis=1)
k['type']=10
l=data['NMDs---antisense_lncRNAs'].to_frame().rename({'NMDs---antisense_lncRNAs':'info'},axis=1)
l['type']=11
m=data['antisense_lncRNAs---antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs---antisense_lncRNAs':'info'},axis=1)
m['type']=12
n=data['NMDs---sncRNAs'].to_frame().rename({'NMDs---sncRNAs':'info'},axis=1)
n['type']=13
o=data['protein_coding_transcripts---antisense_lncRNAs'].to_frame().rename({'protein_coding_transcripts---antisense_lncRNAs':'info'},axis=1)
o['type']=14
df=pd.concat([a,b,c,d,e,f,g,h,j,k,l,m,n,o])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('one of head-to-head promoter transcripts targeted by miRNAs')
plt.subplot(3, 2, 2)
data=divergent_miRNA_int.loc[divergent_miRNA_int['type']=='bidirection_both']
a=data['non_stop_decays---antisense_lncRNAs'].to_frame().rename({'non_stop_decays---antisense_lncRNAs':'info'},axis=1)
a['type']=1
b=data['non_stop_decays---NMDs'].to_frame().rename({'non_stop_decays---NMDs':'info'},axis=1)
b['type']=2
c=data['sncRNAs---non_stop_decays'].to_frame().rename({'sncRNAs---non_stop_decays':'info'},axis=1)
c['type']=3
d=data['protein_coding_transcripts---NMDs'].to_frame().rename({'protein_coding_transcripts---NMDs':'info'},axis=1)
d['type']=4
e=data['non_stop_decays---non_stop_decays'].to_frame().rename({'non_stop_decays---non_stop_decays':'info'},axis=1)
e['type']=5
f=data['non_stop_decays---protein_coding_transcripts'].to_frame().rename({'non_stop_decays---protein_coding_transcripts':'info'},axis=1)
f['type']=6
g=data['NMDs---NMDs'].to_frame().rename({'NMDs---NMDs':'info'},axis=1)
g['type']=7
h=data['protein_coding_transcripts---sncRNAs'].to_frame().rename({'protein_coding_transcripts---sncRNAs':'info'},axis=1)
h['type']=8
j=data['protein_coding_transcripts---protein_coding_transcripts'].to_frame().rename({'protein_coding_transcripts---protein_coding_transcripts':'info'},axis=1)
j['type']=9
k=data['antisense_lncRNAs---sncRNAs'].to_frame().rename({'antisense_lncRNAs---sncRNAs':'info'},axis=1)
k['type']=10
l=data['NMDs---antisense_lncRNAs'].to_frame().rename({'NMDs---antisense_lncRNAs':'info'},axis=1)
l['type']=11
m=data['antisense_lncRNAs---antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs---antisense_lncRNAs':'info'},axis=1)
m['type']=12
n=data['NMDs---sncRNAs'].to_frame().rename({'NMDs---sncRNAs':'info'},axis=1)
n['type']=13
o=data['protein_coding_transcripts---antisense_lncRNAs'].to_frame().rename({'protein_coding_transcripts---antisense_lncRNAs':'info'},axis=1)
o['type']=14
df=pd.concat([a,b,c,d,e,f,g,h,j,k,l,m,n,o])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('both of head-to-head promoter transcripts targeted by miRNAs')
plt.subplot(3, 2, 3)# row 3, col 2 figure 1
data=divergent_miRNA_int.loc[divergent_miRNA_int['type']=='convergent_single']
a=data['non_stop_decays---antisense_lncRNAs'].to_frame().rename({'non_stop_decays---antisense_lncRNAs':'info'},axis=1)
a['type']=1
b=data['non_stop_decays---NMDs'].to_frame().rename({'non_stop_decays---NMDs':'info'},axis=1)
b['type']=2
c=data['sncRNAs---non_stop_decays'].to_frame().rename({'sncRNAs---non_stop_decays':'info'},axis=1)
c['type']=3
d=data['protein_coding_transcripts---NMDs'].to_frame().rename({'protein_coding_transcripts---NMDs':'info'},axis=1)
d['type']=4
e=data['non_stop_decays---non_stop_decays'].to_frame().rename({'non_stop_decays---non_stop_decays':'info'},axis=1)
e['type']=5
f=data['non_stop_decays---protein_coding_transcripts'].to_frame().rename({'non_stop_decays---protein_coding_transcripts':'info'},axis=1)
f['type']=6
g=data['NMDs---NMDs'].to_frame().rename({'NMDs---NMDs':'info'},axis=1)
g['type']=7
h=data['protein_coding_transcripts---sncRNAs'].to_frame().rename({'protein_coding_transcripts---sncRNAs':'info'},axis=1)
h['type']=8
j=data['protein_coding_transcripts---protein_coding_transcripts'].to_frame().rename({'protein_coding_transcripts---protein_coding_transcripts':'info'},axis=1)
j['type']=9
k=data['antisense_lncRNAs---sncRNAs'].to_frame().rename({'antisense_lncRNAs---sncRNAs':'info'},axis=1)
k['type']=10
l=data['NMDs---antisense_lncRNAs'].to_frame().rename({'NMDs---antisense_lncRNAs':'info'},axis=1)
l['type']=11
m=data['antisense_lncRNAs---antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs---antisense_lncRNAs':'info'},axis=1)
m['type']=12
n=data['NMDs---sncRNAs'].to_frame().rename({'NMDs---sncRNAs':'info'},axis=1)
n['type']=13
o=data['protein_coding_transcripts---antisense_lncRNAs'].to_frame().rename({'protein_coding_transcripts---antisense_lncRNAs':'info'},axis=1)
o['type']=14
df=pd.concat([a,b,c,d,e,f,g,h,j,k,l,m,n,o])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('one of tail-to-tail promoter transcripts targeted by miRNAs')
plt.subplot(3, 2, 4)# row 3, col 2 figure 1
data=divergent_miRNA_int.loc[divergent_miRNA_int['type']=='convergent_both']
a=data['non_stop_decays---antisense_lncRNAs'].to_frame().rename({'non_stop_decays---antisense_lncRNAs':'info'},axis=1)
a['type']=1
b=data['non_stop_decays---NMDs'].to_frame().rename({'non_stop_decays---NMDs':'info'},axis=1)
b['type']=2
c=data['sncRNAs---non_stop_decays'].to_frame().rename({'sncRNAs---non_stop_decays':'info'},axis=1)
c['type']=3
d=data['protein_coding_transcripts---NMDs'].to_frame().rename({'protein_coding_transcripts---NMDs':'info'},axis=1)
d['type']=4
e=data['non_stop_decays---non_stop_decays'].to_frame().rename({'non_stop_decays---non_stop_decays':'info'},axis=1)
e['type']=5
f=data['non_stop_decays---protein_coding_transcripts'].to_frame().rename({'non_stop_decays---protein_coding_transcripts':'info'},axis=1)
f['type']=6
g=data['NMDs---NMDs'].to_frame().rename({'NMDs---NMDs':'info'},axis=1)
g['type']=7
h=data['protein_coding_transcripts---sncRNAs'].to_frame().rename({'protein_coding_transcripts---sncRNAs':'info'},axis=1)
h['type']=8
j=data['protein_coding_transcripts---protein_coding_transcripts'].to_frame().rename({'protein_coding_transcripts---protein_coding_transcripts':'info'},axis=1)
j['type']=9
k=data['antisense_lncRNAs---sncRNAs'].to_frame().rename({'antisense_lncRNAs---sncRNAs':'info'},axis=1)
k['type']=10
l=data['NMDs---antisense_lncRNAs'].to_frame().rename({'NMDs---antisense_lncRNAs':'info'},axis=1)
l['type']=11
m=data['antisense_lncRNAs---antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs---antisense_lncRNAs':'info'},axis=1)
m['type']=12
n=data['NMDs---sncRNAs'].to_frame().rename({'NMDs---sncRNAs':'info'},axis=1)
n['type']=13
o=data['protein_coding_transcripts---antisense_lncRNAs'].to_frame().rename({'protein_coding_transcripts---antisense_lncRNAs':'info'},axis=1)
o['type']=14
df=pd.concat([a,b,c,d,e,f,g,h,j,k,l,m,n,o])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('both of tail-to-tail promoter transcripts targeted by miRNAs')
plt.savefig('summary_of_divergent_transcription_intraction_with_miRNAs.jpg',dpi=300)
plt.close()


plt.figure(figsize=(18,6),dpi=300)  
plt.subplot(1, 3, 1)# row 1, col 3 figure 1
a=ppinfo['targeted_bid_CDS_length'].to_frame().rename({'targeted_bid_CDS_length':'info'},axis=1)
a['type1']='head-to-head'
a['type2']='targeted'
b=ppinfo['rest_bid_CDS_length'].to_frame().rename({'rest_bid_CDS_length':'info'},axis=1)
b['type1']='head-to-head'
b['type2']='rest'
c=ppinfo['targeted_conv_CDS_length'].to_frame().rename({'targeted_conv_CDS_length':'info'},axis=1)
c['type1']='tail-to-tail'
c['type2']='targeted'
d=ppinfo['rest_conv_CDS_length'].to_frame().rename({'rest_conv_CDS_length':'info'},axis=1)
d['type1']='tail-to-tail'
d['type2']='rest'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type1", y="info", hue="type2",data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('CDS length')
plt.subplot(1, 3, 2)# row 1, col 3 figure 2
a=ppinfo['targeted_bid_UTR_length'].to_frame().rename({'targeted_bid_UTR_length':'info'},axis=1)
a['type1']='head-to-head'
a['type2']='targeted'
b=ppinfo['rest_bid_UTR_length'].to_frame().rename({'rest_bid_UTR_length':'info'},axis=1)
b['type1']='head-to-head'
b['type2']='rest'
c=ppinfo['targeted_conv_UTR_length'].to_frame().rename({'targeted_conv_UTR_length':'info'},axis=1)
c['type1']='tail-to-tail'
c['type2']='targeted'
d=ppinfo['rest_conv_UTR_length'].to_frame().rename({'rest_conv_UTR_length':'info'},axis=1)
d['type1']='tail-to-tail'
d['type2']='rest'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type1", y="info", hue="type2",data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('3UTR length')
plt.subplot(1, 3, 3)# row 1, col 3 figure 3
a=ppinfo['targeted_bid_UTR_GC'].to_frame().rename({'targeted_bid_UTR_GC':'info'},axis=1)
a['type1']='head-to-head'
a['type2']='targeted'
b=ppinfo['rest_bid_UTR_GC'].to_frame().rename({'rest_bid_UTR_GC':'info'},axis=1)
b['type1']='head-to-head'
b['type2']='rest'
c=ppinfo['targeted_conv_UTR_GC'].to_frame().rename({'targeted_conv_UTR_GC':'info'},axis=1)
c['type1']='tail-to-tail'
c['type2']='targeted'
d=ppinfo['rest_conv_UTR_GC'].to_frame().rename({'rest_conv_UTR_GC':'info'},axis=1)
d['type1']='tail-to-tail'
d['type2']='rest'
df=pd.concat([a,b,c,d])
ax = sns.boxplot(x="type1", y="info", hue="type2",data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.title('3UTR_GC')
plt.savefig('miRNA-targeted_vs_not-targeted_coding_tanscripts_in_coding_diveregent_pairs.jpg',dpi=300)
plt.close()











