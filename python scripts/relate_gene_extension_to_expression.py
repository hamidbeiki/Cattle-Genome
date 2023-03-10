#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 09:25:17 2021

@author: beiki
"""
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def compare_expressions(query,ref):
    results={}
    f=query.merge(faang_counts)
    r=query.merge(ref)
    for tissue in tissues:
        a=f.loc[f[tissue]>0]
        b=r.loc[r.refrence_gene.isin(a.refrence_gene)][tissue].to_frame()+1
        a=a[tissue].to_frame()+1
        a=a.rename({tissue:"a"},axis=1)
        b=b.rename({tissue:"b"},axis=1)
        c=pd.concat([a,b],axis=1)
        d=c.loc[c['a']>c['b']]
        dc=d.copy()
        dc['diff']=dc['a']-dc['b']
        dc['fc']=np.log2(dc['a']/dc['b'])#log2 fold change
        fc=2**dc['fc'].median()#everage fold change
        reads=dc['diff'].median()# average increase of rna-seq reads due to gene extension
        results[tissue]=[reads,fc]
    return(results)

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

def get_extension_reads(annotation,info):
    global faang_counts
    global tissue
    df=annotation.reset_index()
    df=info.merge(df)
    df_ref=faang_counts.loc[faang_counts.FAANG_gene.isin(df.FAANG_gene)]
    df=df.set_index('FAANG_gene')
    df_ref=df_ref.set_index('FAANG_gene')
    df=df.loc[:, df.columns.isin(df_ref.columns)]
    df_np=df.to_numpy()
    df_ref_np=df_ref.to_numpy()
    diff=df_ref_np - df_np
    df_np=df_np + 0.0001# to remove divide by zero
    fc=df_ref_np / df_np
    diff_df=pd.DataFrame(diff)
    diff_df.columns=tissues
    diff_df['FAANG_gene']=df_ref.index
    fc_df=pd.DataFrame(fc)
    fc_df.columns=tissues
    fc_df['FAANG_gene']=df_ref.index
    return(diff_df,df_ref.reset_index(),fc_df)

def summarizer(diff_df,expr,fc):
    global tissues
    res_dict={}
    res_dict2={}
    for tissue in tissues:
        ex=expr[['FAANG_gene',tissue]]
        ex=ex.loc[ex[tissue]>0]['FAANG_gene'].to_frame()
        df=diff_df[['FAANG_gene',tissue]]
        df=ex.merge(df)
        res_df=df.loc[df[tissue]>=0]
        res_dict[tissue]=res_df
        df=fc[['FAANG_gene',tissue]]
        df=ex.merge(df)
        res_df=df.loc[df[tissue]>=1]
        res_dict2[tissue]=res_df        
    return(res_dict,res_dict2)


def summarizer2(info1,info2):
    info=get_extension_reads(info1,info2)
    diff=info[0]
    expr=info[1]
    fc=info[2]
    info=summarizer(diff,expr,fc)
    res_dict=info[0]#reads mapped to the extended regions
    res_dict2=info[1]#fold-change resulted from extended regions
    return(res_dict,res_dict2)  
    
ens_same=pd.read_csv('genes_atleast_1_exact_same_transcript_with_Ensembl_gene',names=['FAANG_gene'],sep='\t')
ens_ex_5_df=pd.read_csv('Ensembl_genes_5end_extensions.txt',header=0,sep='\t')
ens_ex_5_df=ens_same.merge(ens_ex_5_df)
ens_ex_3_df=pd.read_csv('Ensembl_genes_3end_extensions.txt',header=0,sep='\t')
ens_ex_3_df=ens_same.merge(ens_ex_3_df)
ens_both_ex_df=pd.read_csv('Ensembl_genes_bothEnd_extensions.txt',header=0,sep='\t')

ncbi_same=pd.read_csv('genes_atleast_1_exact_same_transcript_with_NCBI_gene',names=['FAANG_gene'],sep='\t')
ncbi_ex_5_df=pd.read_csv('NCBI_genes_5end_extensions.txt',header=0,sep='\t')
ncbi_ex_5_df=ncbi_same.merge(ncbi_ex_5_df)
ncbi_ex_3_df=pd.read_csv('NCBI_genes_3end_extensions.txt',header=0,sep='\t')
ncbi_ex_3_df=ncbi_same.merge(ncbi_ex_3_df)
ncbi_both_ex_df=pd.read_csv('NCBI_genes_bothEnd_extensions.txt',header=0,sep='\t')

remove=pd.read_csv('remove_genes',header=0,names=['FAANG_gene'])
faang_counts=pd.read_csv('../quantification/RSEM_FAANG_genes_counts',header=0,sep="\t").rename({'gene_id':'FAANG_gene'},axis=1)
ens_counts=pd.read_csv('../quantification/RSEM_ENSEMBL_genes_counts',header=0,sep="\t").rename({'gene_id':'refrence_gene'},axis=1)
ncbi_counts=pd.read_csv('../quantification/RSEM_NCBI_genes_counts',header=0,sep="\t").rename({'gene_id':'gene_name'},axis=1)
ncbi_info=pd.read_csv('../quantification/RSEM-quantified_ncbi_gene_name_to_geneID',names=['refrence_gene','gene_name'],sep="\t")
ncbi_counts=ncbi_info.merge(ncbi_counts)
tissues=list(faang_counts.columns)[1:]

faang_counts=faang_counts.set_index('FAANG_gene')
faang_counts=faang_counts.loc[(faang_counts!=0).any(1)]#remove genes that couldn't be quantified
faang_counts=faang_counts.loc[~ faang_counts.index.isin(remove.FAANG_gene)]
faang_counts=faang_counts.reset_index()



    
for i in ['ensembl','ncbi']:
    if i=='ensembl':
        info=summarizer2(ens_counts,ens_ex_5_df)
        ens_5_dict=info[0]#reads mapped to the extended region of expressed estended genes
        ens_5_FC_dict=info[1]#fold-change resulted from extended regions
        info=summarizer2(ens_counts,ens_ex_3_df)
        ens_3_dict=info[0]
        ens_3_FC_dict=info[1]
        info=summarizer2(ens_counts,ens_both_ex_df)
        ens_both_dict=info[0]
        ens_both_FC_dict=info[1]
    else:
        info=summarizer2(ncbi_counts,ncbi_ex_5_df)
        ncbi_5_dict=info[0]#reads mapped to the extended region of expressed estended genes
        ncbi_5_FC_dict=info[1]#fold-change resulted from extended regions
        info=summarizer2(ncbi_counts,ncbi_ex_3_df)
        ncbi_3_dict=info[0]
        ncbi_3_FC_dict=info[1]
        info=summarizer2(ncbi_counts,ncbi_both_ex_df)
        ncbi_both_dict=info[0]
        ncbi_both_FC_dict=info[1]      

flag=1
for key,value in ncbi_3_FC_dict.items():     
    globals() ['a' + str(flag)]=value[key].to_frame().rename({key:"info"},axis=1)
    globals() ['a' + str(flag)]['type']=flag
    flag=flag+1
        
df=pd.concat([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47])
df.to_csv('results.txt',index = None, header=True,sep="\t")
#df=df.loc[df['info']>0]
df=df.loc[df['info']>1]
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
#plt.savefig('distribution_of_reads_mapped_to_extended_genes_3end.png')
#plt.close()

ncbi_3_res=compare_expressions(ncbi_ex_3_df,ncbi_counts)
ncbi_5_res=compare_expressions(ncbi_ex_5_df,ncbi_counts)
ncbi_both_res=compare_expressions(ncbi_both_ex_df,ncbi_counts)
ens_3_res=compare_expressions(ens_ex_3_df,ens_counts)
ens_5_res=compare_expressions(ens_ex_5_df,ens_counts)
ens_both_res=compare_expressions(ens_both_ex_df,ens_counts)

ncbi_3_res_df=pd.DataFrame.from_dict(ncbi_3_res, orient='index').reset_index().rename({"index":"tissue",0:'ncbi_3_reads',1:'ncbi_3_FC'},axis=1)
ncbi_5_res_df=pd.DataFrame.from_dict(ncbi_5_res, orient='index').reset_index().rename({"index":"tissue",0:'ncbi_5_reads',1:'ncbi_5_FC'},axis=1).iloc[0:,1:]
ncbi_both_res_df=pd.DataFrame.from_dict(ncbi_both_res, orient='index').reset_index().rename({"index":"tissue",0:'ncbi_both_reads',1:'ncbi_both_FC'},axis=1).iloc[0:,1:]

ens_3_res_df=pd.DataFrame.from_dict(ens_3_res, orient='index').reset_index().rename({"index":"tissue",0:'ens_3_reads',1:'ens_3_FC'},axis=1).iloc[0:,1:]
ens_5_res_df=pd.DataFrame.from_dict(ens_5_res, orient='index').reset_index().rename({"index":"tissue",0:'ens_5_reads',1:'ens_5_FC'},axis=1).iloc[0:,1:]
ens_both_res_df=pd.DataFrame.from_dict(ens_both_res, orient='index').reset_index().rename({"index":"tissue",0:'ens_both_reads',1:'ens_both_FC'},axis=1).iloc[0:,1:]

df=pd.concat([ncbi_3_res_df,ncbi_5_res_df,ncbi_both_res_df,ens_3_res_df,ens_5_res_df,ens_both_res_df],axis=1)
df.to_csv('effect_of_gene_extension_on_expression.txt',index = None, header=True,sep="\t")

stat, p=stats.ttest_ind(df['ncbi_3_reads'], df['ncbi_5_reads'], equal_var=True,alternative='greater')#https://stackoverflow.com/questions/25064506/what-scipy-statistical-test-do-i-use-to-compare-sample-means

df=df.reset_index()

df.plot(kind="scatter", x="tissue", y="ncbi_3_reads", alpha=1,
    figsize=(10,8),
    c="ncbi_5_reads", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,label='ncbi_3_FC', s=10)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
#plt.xticks(rotation=60)
plt.legend()
plt.savefig('test.png')
plt.close()

df.plot(kind="scatter", x="index", y="ens_both_reads", alpha=1,
    figsize=(10,8),
    c="ens_both_reads", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,label='ens_both_FC', s=10)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.legend()
plt.savefig('test.png')
plt.close()

ncbi_counts=ncbi_counts.drop(['gene_name'],axis=1)
ncbi_counts=ncbi_counts.set_index('refrence_gene')
ens_counts=ens_counts.set_index('refrence_gene')
faang_counts=faang_counts.set_index('FAANG_gene')

ens_ex_5_ex=ens_counts.loc[ens_counts.index.isin(ens_ex_5_df.refrence_gene)]
ens_ex_3_ex=ens_counts.loc[ens_counts.index.isin(ens_ex_3_df.refrence_gene)]
ens_both_ex=ens_counts.loc[ens_counts.index.isin(ens_both_ex_df.refrence_gene)]

ncbi_ex_5_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_5_df.refrence_gene)]
ncbi_ex_3_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_3_df.refrence_gene)]
ncbi_both_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_both_ex_df.refrence_gene)]

ens_ex_5_info=get_info(ens_ex_5_ex,'refrence_gene').rename({'detected_tissues':'ens_detected_tissues','median_expression':'ens_median_expression'},axis=1)
ens_ex_5_info=ens_ex_5_df.merge(ens_ex_5_info)
ens_ex_3_info=get_info(ens_ex_3_ex,'refrence_gene').rename({'detected_tissues':'ens_detected_tissues','median_expression':'ens_median_expression'},axis=1)
ens_ex_3_info=ens_ex_3_df.merge(ens_ex_3_info)
ens_both_info=get_info(ens_both_ex,'refrence_gene').replace(np.NaN,0).rename({'detected_tissues':'ens_detected_tissues','median_expression':'ens_median_expression'},axis=1)# waring is due to genes that were not expressed in any tissue
ens_both_info=ens_both_ex_df.merge(ens_both_info)

ncbi_ex_5_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_5_df.refrence_gene)]
ncbi_ex_3_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_3_df.refrence_gene)]
ncbi_both_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_both_ex_df.refrence_gene)]

ncbi_ex_5_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_5_df.refrence_gene)]
ncbi_ex_3_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_ex_3_df.refrence_gene)]
ncbi_both_ex=ncbi_counts.loc[ncbi_counts.index.isin(ncbi_both_ex_df.refrence_gene)]

ncbi_ex_5_info=get_info(ncbi_ex_5_ex,'refrence_gene').rename({'detected_tissues':'ncbi_detected_tissues','median_expression':'ncbi_median_expression'},axis=1)
ncbi_ex_5_info=ncbi_ex_5_df.merge(ncbi_ex_5_info)
ncbi_ex_3_info=get_info(ncbi_ex_3_ex,'refrence_gene').rename({'detected_tissues':'ncbi_detected_tissues','median_expression':'ncbi_median_expression'},axis=1)
ncbi_ex_3_info=ncbi_ex_3_df.merge(ncbi_ex_3_info)
ncbi_both_info=get_info(ncbi_both_ex,'refrence_gene').replace(np.NaN,0).rename({'detected_tissues':'ncbi_detected_tissues','median_expression':'ncbi_median_expression'},axis=1)# waring is due to genes that were not expressed in any tissue
ncbi_both_info=ncbi_both_ex_df.merge(ncbi_both_info)


faang_info=get_info(faang_counts,'FAANG_gene').rename({'detected_tissues':'faang_detected_tissues','median_expression':'faang_median_expression'},axis=1)


ens_ex_5_info=ens_ex_5_info.merge(faang_info)
ens_ex_3_info=ens_ex_3_info.merge(faang_info)
ens_both_info=ens_both_info.merge(faang_info)

ncbi_ex_5_info=ncbi_ex_5_info.merge(faang_info)
ncbi_ex_3_info=ncbi_ex_3_info.merge(faang_info)
ncbi_both_info=ncbi_both_info.merge(faang_info)

ens_ex_5_info.to_csv('Ensembl_genes_5end_extensions_info.txt',index = None, header=True,sep="\t")
ens_ex_3_info.to_csv('Ensembl_genes_3end_extensions_info.txt',index = None, header=True,sep="\t")
ens_both_info.to_csv('Ensembl_genes_bothEnd_extensions_info.txt',index = None, header=True,sep="\t")

ncbi_ex_5_info.to_csv('NCBI_genes_5end_extensions_info.txt',index = None, header=True,sep="\t")
ncbi_ex_3_info.to_csv('NCBI_genes_3end_extensions_info.txt',index = None, header=True,sep="\t")
ncbi_both_info.to_csv('NCBI_genes_bothEnd_extensions_info.txt',index = None, header=True,sep="\t")

ens_ex_5_info=ens_ex_5_info.loc[ens_ex_5_info['faang_median_expression']>ens_ex_5_info['ens_median_expression']]
ens_ex_3_info=ens_ex_3_info.loc[ens_ex_3_info['faang_median_expression']>ens_ex_3_info['ens_median_expression']]
ens_both_info=ens_both_info.loc[ens_both_info['faang_median_expression']>ens_both_info['ens_median_expression']]
ncbi_ex_5_info=ncbi_ex_5_info.loc[ncbi_ex_5_info['faang_median_expression']>ncbi_ex_5_info['ncbi_median_expression']]
ncbi_ex_3_info=ncbi_ex_3_info.loc[ncbi_ex_3_info['faang_median_expression']>ncbi_ex_3_info['ncbi_median_expression']]
ncbi_both_info=ncbi_both_info.loc[ncbi_both_info['faang_median_expression']>ncbi_both_info['ncbi_median_expression']]

ens_ex_5_info.plot(kind="scatter", x="ens_median_expression", y="faang_median_expression", alpha=1,
    figsize=(10,8),
    c="extension_length", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,label='ens_both_FC', s=10)#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.legend()

ll=plt.scatter(ens_ex_5_info["ens_median_expression"], ens_ex_5_info["faang_median_expression"], marker='o',s=4,alpha=0.3,color='blue')#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
l=plt.plot(ens_ex_5_info["ens_median_expression"],ens_ex_5_info["ens_median_expression"], color='red')#s=out['total'], cmap=plt.cm.get_cmap('cubehelix', 6) or cmap=plt.get_cmap("jet")
plt.xlabel("Number of detected tissues")
plt.ylabel("Median expression level (RPKM) in detected tissues")
plt.show()









"""
test=ens_counts.loc[(ens_counts!=0).any(1)]
zeros=ens_counts[~ ens_counts.index.isin(test.index)]
ens_not_quantified_5=zeros.loc[zeros.index.isin(ens_ex_5_df.refrence_gene)]
ens_not_quantified_3=zeros.loc[zeros.index.isin(ens_ex_3_df.refrence_gene)]
ens_not_quantified_both=zeros.loc[zeros.index.isin(ens_both_ex_df.refrence_gene)].reset_index()

test=ncbi_counts.loc[(ncbi_counts!=0).any(1)]
zeros=ncbi_counts[~ ncbi_counts.index.isin(test.index)]
ncbi_not_quantified_5=zeros.loc[zeros.index.isin(ncbi_ex_5_df.refrence_gene)]
ncbi_not_quantified_3=zeros.loc[zeros.index.isin(ncbi_ex_3_df.refrence_gene)]
ncbi_not_quantified_both=zeros.loc[zeros.index.isin(ncbi_both_ex_df.refrence_gene)].reset_index()
"""


"""
from pandas.plotting import table
fig, ax = plt.subplots(1, 1)
table(ax, np.round(df['ncbi_3_FC'].describe(), 2), loc="upper right", colWidths=[0.2, 0.2, 0.2])
ax.xaxis.tick_top()  # Display x-axis ticks on top.
df.plot(ax=ax,kind="scatter", x="tissue", y="ncbi_3_reads", alpha=1,
    figsize=(10,7),
    c="ncbi_5_reads", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,label='ncbi_3_FC', s=df["ncbi_3_FC"])#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
#plt.xticks(rotation=60)
plt.legend()
plt.savefig('test.png')
plt.close()

fig, ax = plt.subplots(1, 1)
ax.xaxis.tick_top()  # Display x-axis ticks on top.
df.plot(table=np.round(df['ncbi_3_FC'].describe(), 1),ax=ax,kind="scatter", x="tissue", y="ncbi_3_reads", alpha=1,
    figsize=(10,8),
    c="ncbi_5_reads", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,label='ncbi_3_FC', s=df["ncbi_3_FC"])#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
#plt.xticks(rotation=60)
plt.legend()
plt.savefig('test.png')
plt.close()
"""
