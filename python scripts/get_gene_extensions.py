#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 08:58:21 2021

@author: beiki
"""

import pandas as pd
import functools
import operator
import statistics

def flatten_nested_list(lis):#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return functools.reduce(operator.iconcat, lis, [])

def get_borders(df):
    strand=pd.DataFrame(df[['gene','strand']]).drop_duplicates()
    genes_start=df.groupby(['gene'], sort=False)['ex_s'].min()
    genes_end=df.groupby(['gene'], sort=False)['ex_e'].max()
    genes_borders=pd.concat([genes_start,genes_end],axis=1).reset_index()
    m=genes_borders.merge(strand)
    m=m.rename({'ex_s':'g_s','ex_e':'g_e'},axis=1)
    return(m)

def get_intervals(df):
    l=[]
    l1=list(df['ex_s'])
    l2=list(df['ex_e'])
    for i in range(len(l1)):
        l.append([l1[i],l2[i]])
    return(l)

def merge_intervals(arr):
    # Sorting based on the increasing order  
    # of the start intervals
    arr.sort(key = lambda x: x[0]) 
    # array to hold the merged intervals 
    m = []
    s = -10000
    max = -100000
    for i in range(len(arr)):
        a = arr[i]
        if a[0] > max:
            if i != 0:
                m.append([s,max])
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]
    #'max' value gives the last point of
    # that particular interval
    # 's' gives the starting point of that interval
    # 'm' array contains the list of all merged intervals
    if max != -100000 and [s, max] not in m:
        m.append([s, max])
    return(m)
    
def polish_dfs5(df2,df3,strand):
    if strand=="+":
        if df3.shape[0]==0:
            return(df2['length'].sum())
        elif df3.shape[0]==df2.shape[0]:
            res=(abs(df3['ex_s']-df3['ref_g_s']+1)).to_list()[0]
            return(res)
        else:
            a=df2[~ df2.ex_s.isin(df3.ex_s)]
            if a.shape[0]>0:
                arr=get_intervals(df3)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r1=abs((d.iloc[0,0]-df3['ref_g_s']+1).to_list()[0])
                arr=get_intervals(a)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()            
                return(r1+r2)
            else:
                arr=get_intervals(df2)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r2)
    elif strand=="-":
        if df3.shape[0]==0:
            return(df2['length'].sum())
        elif df3.shape[0]==df2.shape[0]:
            res=abs(df3['ex_e']-df3['ref_g_e']+1).to_list()[0]
            return(res)
        else:
            a=df2[~ df2.ex_s.isin(df3.ex_s)]
            if a.shape[0]>0:
                arr=get_intervals(df3)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r1=(d.iloc[0,1]-df3['ref_g_e']+1).to_list()[0]
                arr=get_intervals(a)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r1+r2)
            else:
                arr=get_intervals(df2)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r2)                

def get_extensions5(query):
    out=[]#extension length
    out2=[]#refrence gene that was extended by FAANG gene
    out3=[]#FAANG gene
    df=query.copy()
    genes=list(df['ref_gene'].drop_duplicates())
    genes2=list(df['gene'].drop_duplicates())
    df['length']=abs(df['ex_e']-df['ex_s'])+1#exon length
    flag=0
    for gene in genes:
        gene2=genes2[flag]
        flag=flag+1
        info=df.loc[df['ref_gene']==gene]
        strand=info.iloc[0,5]
        if strand=="+":
            df2=info.loc[(info['ex_s'] < info['ref_g_s'])]
            df3=df2.loc[(df2['ex_e'] > df2['ref_g_s'])]
            out.append(polish_dfs5(df2,df3,strand))
            out2.append(gene)
            out3.append(gene2)
        elif strand=="-":
            df2=info.loc[(info['ex_e'] > info['ref_g_e'])]
            df3=df2.loc[(df2['ex_s'] < df2['ref_g_e'])]
            out.append(polish_dfs5(df2,df3,strand))
            out2.append(gene)
            out3.append(gene2)
    return(out,out2,out3)

def polish_dfs3(df2,df3,strand):
    if strand=="+":
        if df3.shape[0]==0:
            return(df2['length'].sum())
        elif df3.shape[0]==df2.shape[0]:
            res=(abs(df3['ex_e']-df3['ref_g_e']+1)).to_list()[0]
            return(res)
        else:
            a=df2[~ df2.ex_s.isin(df3.ex_s)]
            if a.shape[0]>0:
                arr=get_intervals(df3)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r1=(d.iloc[0,1]-df3['ref_g_e']+1).to_list()[0]
                arr=get_intervals(a)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r1+r2)
            else:
                arr=get_intervals(df2)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r2)                                 
    elif strand=="-":
        if df3.shape[0]==0:
            return(df2['length'].sum())
        elif df3.shape[0]==df2.shape[0]:
            res=abs(df3['ex_s']-df3['ref_g_s']+1).to_list()[0]
            return(res)
        else:
            a=df2[~ df2.ex_s.isin(df3.ex_s)]
            if a.shape[0]>0:
                arr=get_intervals(df3)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r1=abs((d.iloc[0,0]-df3['ref_g_s']+1).to_list()[0])
                arr=get_intervals(a)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()            
                return(r1+r2)
            else:
                arr=get_intervals(df2)
                l=merge_intervals(arr)
                d=pd.DataFrame(l)
                r2=(d[1]-d[0]+1).sum()
                return(r2)                 
            
def get_extensions3(query):
    out=[]#extension length
    out2=[]#refrence gene that was extended by FAANG gene
    out3=[]#FAANG gene
    df=query.copy()
    genes=list(df['ref_gene'].drop_duplicates())
    genes2=list(df['gene'].drop_duplicates())
    df['length']=abs(df['ex_e']-df['ex_s'])+1#exon length
    flag=0
    for gene in genes:
        gene2=genes2[flag]
        flag=flag+1
        info=df.loc[df['ref_gene']==gene]
        strand=info.iloc[0,5]
        if strand=="+":
            df2=info.loc[(info['ex_e'] > info['ref_g_e'])]
            df3=df2.loc[(df2['ex_s'] < df2['ref_g_e'])]
            out.append(polish_dfs3(df2,df3,strand))
            out2.append(gene)
            out3.append(gene2)
        elif strand=="-":
            df2=info.loc[(info['ex_s'] < info['ref_g_s'])]
            df3=df2.loc[(df2['ex_e'] > df2['ref_g_s'])]
            out.append(polish_dfs3(df2,df3,strand))
            out2.append(gene)
            out3.append(gene2)
    return(out,out2,out3)
            
def get_extensions(query,ex_type):#query could be ens_5,ens_3,...
    out=[]#extension length
    out2=[]#refrence gene that was extended by FAANG gene
    out3=[]#FAANG gene
    m=faang.merge(query)
    p=m.loc[m['strand']=="+"]
    n=m.loc[m['strand']=="-"]
    if ex_type=="5_end":
        out.append(get_extensions5(p)[0])
        out.append(get_extensions5(n)[0])
        out2.append(get_extensions5(p)[1])
        out2.append(get_extensions5(n)[1])
        out3.append(get_extensions5(p)[2])
        out3.append(get_extensions5(n)[2])
    elif ex_type=="3_end":
        out.append(get_extensions3(p)[0])
        out.append(get_extensions3(n)[0])
        out2.append(get_extensions5(p)[1])
        out2.append(get_extensions5(n)[1])
        out3.append(get_extensions5(p)[2])
        out3.append(get_extensions5(n)[2])
    re=flatten_nested_list(out)
    d=pd.DataFrame(re).rename({0:"extension_length"},axis=1)
    re2=flatten_nested_list(out2)
    d2=pd.DataFrame(re2).rename({0:"refrence_gene"},axis=1)
    re3=flatten_nested_list(out3)
    d3=pd.DataFrame(re3).rename({0:"FAANG_gene"},axis=1)
    out_df=pd.concat([d3,d2,d],axis=1)
    return(re,out_df)
    
ncbi=pd.read_csv('ncbi-genes.bed',names=['chr','ex_s','ex_e','gene','num','strand','g'],sep="\t")
ensembl=pd.read_csv('ensembl-genes.bed',names=['chr','ex_s','ex_e','gene','num','strand','g'],sep="\t")
remove=pd.read_csv('remove_genes',names=['gene'])
faang=pd.read_csv('combined-genes.bed',names=['chr','ex_s','ex_e','gene','num','strand','g'],sep="\t")
faang=faang[~faang.gene.isin(remove.gene)]
m_ensmbl=pd.read_csv('genes_covered_multiple_ensembl',names=['gene'],sep="\t")
m_ncbi=pd.read_csv('genes_covered_multiple_ncbi',names=['gene'],sep="\t")
faang=faang[~faang.gene.isin(m_ensmbl.gene)]
faang=faang[~faang.gene.isin(m_ncbi.gene)]
ens_eq=pd.read_csv('combined_genes_ensembl_equivalent',names=['gene','ens_g'],sep=" ")
ncbi_eq=pd.read_csv('combined_genes_ncbi_equivalent',names=['gene','ncbi_g'],sep=" ")
ens_eq=ens_eq[~ens_eq.gene.isin(remove.gene)]
ncbi_eq=ncbi_eq[~ncbi_eq.gene.isin(remove.gene)]

rampage=pd.read_csv('gene_transcript_validation/genes_supported_by_CAGE',names=['FAANG_gene'])
rampage=rampage.loc[~ rampage.FAANG_gene.isin(remove.gene)]
wtts=pd.read_csv('gene_transcript_validation/genes_supported_by_WTTS',names=['FAANG_gene'])
wtts=wtts.loc[~ wtts.FAANG_gene.isin(remove.gene)]

faang_b=get_borders(faang)
ens_b=get_borders(ensembl).rename({'gene':'ens_g','g_s':'ens_g_s','g_e':'ens_g_e','strand':'ens_strand'},axis=1)
ncbi_b=get_borders(ncbi).rename({'gene':'ncbi_g','g_s':'ncbi_g_s','g_e':'ncbi_g_e','strand':'ncbi_strand'},axis=1)

ens_m=ens_b.merge(ens_eq)
ens_m=ens_m.merge(faang_b)
ncbi_m=ncbi_b.merge(ncbi_eq)
ncbi_m=ncbi_m.merge(faang_b)

a=ens_m.loc[(ens_m['ens_strand']=="+") & (ens_m['g_s'] < ens_m['ens_g_s']) & (ens_m['g_e'] > ens_m['ens_g_e'])]
b=ens_m.loc[(ens_m['ens_strand']=="-") & (ens_m['g_s'] < ens_m['ens_g_s']) & (ens_m['g_e'] > ens_m['ens_g_e'])]
ens_both=pd.concat([a,b])
a=ens_m.loc[(ens_m['ens_strand']=="+") & (ens_m['g_s'] < ens_m['ens_g_s'])]
b=ens_m.loc[(ens_m['ens_strand']=="-") & (ens_m['g_e'] > ens_m['ens_g_e'])]
c=pd.concat([a,b])
ens_5=c[~ c.gene.isin(ens_both.gene)]
a=ens_m.loc[(ens_m['ens_strand']=="+") & (ens_m['g_e'] > ens_m['ens_g_e'])]
b=ens_m.loc[(ens_m['ens_strand']=="-") & (ens_m['g_s'] < ens_m['ens_g_s'])]
c=pd.concat([a,b])
ens_3=c[~ c.gene.isin(ens_both.gene)]

a=ncbi_m.loc[(ncbi_m['ncbi_strand']=="+") & (ncbi_m['g_s'] < ncbi_m['ncbi_g_s']) & (ncbi_m['g_e'] > ncbi_m['ncbi_g_e'])]
b=ncbi_m.loc[(ncbi_m['ncbi_strand']=="-") & (ncbi_m['g_s'] < ncbi_m['ncbi_g_s']) & (ncbi_m['g_e'] > ncbi_m['ncbi_g_e'])]
ncbi_both=pd.concat([a,b])
a=ncbi_m.loc[(ncbi_m['ncbi_strand']=="+") & (ncbi_m['g_s'] < ncbi_m['ncbi_g_s'])]
b=ncbi_m.loc[(ncbi_m['ncbi_strand']=="-") & (ncbi_m['g_e'] > ncbi_m['ncbi_g_e'])]
c=pd.concat([a,b])
ncbi_5=c[~ c.gene.isin(ncbi_both.gene)]
a=ncbi_m.loc[(ncbi_m['ncbi_strand']=="+") & (ncbi_m['g_e'] > ncbi_m['ncbi_g_e'])]
b=ncbi_m.loc[(ncbi_m['ncbi_strand']=="-") & (ncbi_m['g_s'] < ncbi_m['ncbi_g_s'])]
c=pd.concat([a,b])
ncbi_3=c[~ c.gene.isin(ncbi_both.gene)]


ensembl=ensembl.rename({'gene':'ref_gene'},axis=1)
ens_5=ens_5.rename({'ens_g':'ref_gene','ens_g_s':'ref_g_s','ens_g_e':'ref_g_e'},axis=1)
ens_3=ens_3.rename({'ens_g':'ref_gene','ens_g_s':'ref_g_s','ens_g_e':'ref_g_e'},axis=1)  
ens_both=ens_both.rename({'ens_g':'ref_gene','ens_g_s':'ref_g_s','ens_g_e':'ref_g_e'},axis=1)  
res=get_extensions(ens_5,"5_end")
ens_ex_5=res[0]
ens_ex_5_df=res[1]
res=get_extensions(ens_3,"3_end")
ens_ex_3=res[0]
ens_ex_3_df=res[1]
res1=get_extensions(ens_both,"5_end")
res2=get_extensions(ens_both,"3_end")
d1=res1[1].rename({"extension_length":"5_extension_length"},axis=1)
d2=res2[1].rename({"extension_length":"3_extension_length"},axis=1)
m=d1.merge(d2)
ens_both_ex_df=m
ens_both_ex_5=res1[0]
ens_both_ex_3=res2[0]

ncbi=ncbi.rename({'gene':'ref_gene'},axis=1)
ncbi_5=ncbi_5.rename({'ncbi_g':'ref_gene','ncbi_g_s':'ref_g_s','ncbi_g_e':'ref_g_e'},axis=1)
ncbi_3=ncbi_3.rename({'ncbi_g':'ref_gene','ncbi_g_s':'ref_g_s','ncbi_g_e':'ref_g_e'},axis=1)
ncbi_both=ncbi_both.rename({'ncbi_g':'ref_gene','ncbi_g_s':'ref_g_s','ncbi_g_e':'ref_g_e'},axis=1)
res=get_extensions(ncbi_5,"5_end")
ncbi_ex_5=res[0]
ncbi_ex_5_df=res[1]
res=get_extensions(ncbi_3,"3_end")
ncbi_ex_3=res[0]
ncbi_ex_3_df=res[1]
res1=get_extensions(ncbi_both,"5_end")
res2=get_extensions(ncbi_both,"3_end")
d1=res1[1].rename({"extension_length":"5_extension_length"},axis=1)
d2=res2[1].rename({"extension_length":"3_extension_length"},axis=1)
m=d1.merge(d2)
ncbi_both_ex_df=m
ncbi_both_ex_5=res1[0]
ncbi_both_ex_3=res2[0]

statistics.median(ens_ex_5)
statistics.median(ens_ex_3)
statistics.median(ens_both_ex_5)
statistics.median(ens_both_ex_3)

statistics.median(ncbi_ex_5)
statistics.median(ncbi_ex_3)
statistics.median(ncbi_both_ex_5)
statistics.median(ncbi_both_ex_3) 

ens_ex_5_df.to_csv('Ensembl_genes_5end_extensions.txt',index = None, header=True,sep="\t")
ens_ex_3_df.to_csv('Ensembl_genes_3end_extensions.txt',index = None, header=True,sep="\t")
ens_both_ex_df.to_csv('Ensembl_genes_bothEnd_extensions.txt',index = None, header=True,sep="\t")

ncbi_ex_5_df.to_csv('NCBI_genes_5end_extensions.txt',index = None, header=True,sep="\t")
ncbi_ex_3_df.to_csv('NCBI_genes_3end_extensions.txt',index = None, header=True,sep="\t")
ncbi_both_ex_df.to_csv('NCBI_genes_bothEnd_extensions.txt',index = None, header=True,sep="\t")

""""
out=[]
query=n
df=query.copy()
genes=list(df['ref_gene'].drop_duplicates())
df['length']=abs(df['ex_e']-df['ex_s'])+1#exon length
for gene in genes:
    info=df.loc[df['ref_gene']==gene]
    strand=info.iloc[0,5]
    if strand=="+":
        df2=info.loc[(info['ex_s'] < info['ref_g_s'])]
        df3=df2.loc[(df2['ex_e'] > df2['ref_g_s'])]
        a=df2[~ df2.ex_s.isin(df3.ex_s)]
        res=polish_dfs5(df2,df3,strand)
        out.append(res)       
    elif strand=="-":
        df2=info.loc[(info['ex_e'] > info['ref_g_e'])]
        df3=df2.loc[(df2['ex_s'] < df2['ref_g_e'])]
        res=polish_dfs5(df2,df3,strand)
        if res==0:
            print(gene)
        out.append(res)
            
    
info=df.loc[df['ref_gene']=='ENSBTAG00000014945']#ENSBTAG00000053072 ENSBTAG00000049054 ENSBTAG00000014945 
"""
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




































    