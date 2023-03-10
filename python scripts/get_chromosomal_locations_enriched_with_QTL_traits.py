#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 13:57:37 2021

@author: beiki
"""
import pandas as pd
import pybedtools
from scipy import stats

""" requiered modules:
    module load py-pandas/0.21.1-py2-326uzkn
    module load py-scipy/1.0.0-py2-flqcuxg
    module load bedtools2
    module load py-pybedtools/0.7.10-py2-a7gbt3r
"""

def get_windows_enriched_with_qtls(df,chrom):
    global qtl_num
    res=[]    
    for index, row in df.iterrows(): 
        window=chrom + '_' + row['window']
        trait=row['trait']
        CA=row['windowCount']
        AB=row['AllQTLcount']
        A=row['windowQTLcount']
        B=AB-A
        C=CA-A
        D=qtl_num-A-B-C
        if A>=3:
            oddsratio,pvalue = stats.fisher_exact([[A, B], [C, D]],alternative='greater')
            if pvalue<0.0001:
                res.append([window,trait,A,pvalue])
    if len(res)>0:
        res=pd.DataFrame(res).rename({0:'chr_window',1:'trait',2:'traits_Qtls',3:'QTL_enrich_pvalue'},axis=1)
    else:
        res=pd.DataFrame()
    return(res)

def reformat_chr_window_names(df):
    res=[]
    for row in df.iterrows():
        b=int(row[1][0].split('_window_')[1])+2000000
        b=str(b)[::-1].replace('000000','bM')[::-1]
        b=str(b)[::-1].replace('000','bK')[::-1]
        a='chr' + row[1][0].replace('_window_',':')
        a=a[::-1].replace('000000','bM')[::-1]
        a=a[::-1].replace('000','bK')[::-1]
        name=a+'-'+str(b)
        res.append([name])
    df=pd.DataFrame(res).rename({0:'chr_window'},axis=1)
    return(df)
    
qtl_bed=pd.read_csv('CattleQTLdb_April2020.bed',sep="\t")
chr_info=pd.read_csv('bt_ref_ARS-UCD1.2.RepeadMasked-chromosome_info',sep="\t",names=['chr','chr_length'])
chr_info=chr_info[~chr_info["chr"].str.contains('ref|MT')]#remove unplaced contigs and MT chromosome

qtl_bed_df=qtl_bed.copy()
qtl_bed_df['trait'] = qtl_bed_df['Body_weight_(mature)_QTL_(65327).1'].str.split('\_QTL\_\(', 1).str[0]
qtl_num=qtl_bed_df.shape[0]
general_info=qtl_bed_df.groupby(['trait']).size().reset_index(name='AllQTLcount').rename({3:'window'},axis=1)
qtl_bed=pybedtools.BedTool.from_dataframe(qtl_bed)

window_size=2000000 ## see section "window size selection down below", this size had the highest number of regions enriched with event genes
regions={}
for i in range(chr_info.shape[0]):
     regions[chr_info.iloc[i,0]]=[]
upper_lim=chr_info['chr_length'].max()+1000000#chr length upper limit
for i in range(chr_info.shape[0]):
    chr_len=chr_info.iloc[i,1]
    for j in range(0,upper_lim,window_size):
        flag=j+window_size-1
        if j<=chr_len:
            regions[chr_info.iloc[i,0]].append([chr_info.iloc[i,0],j,flag,'window_' + str(j),99,'.',chr_info.iloc[i,0] + '_' + 'window_' + str(j)])
        elif flag>chr_len and abs(flag-chr_len)<window_size:
            regions[chr_info.iloc[i,0]].append([chr_info.iloc[i,0],j,chr_len,'window_' + str(j),99,'.','window_' + str(j)])
        else:
            break           
for i in range(chr_info.shape[0]):
    my_list=regions[chr_info.iloc[i,0]]
    regions[chr_info.iloc[i,0]]=pd.DataFrame(my_list).rename({0:'chr',1:'window_start',2:'window_end',3:'window_name',4:'info',5:'strand',6:'window_info'},axis=1)
    
windows_enriched_with_qtls={}
windows_qtl_details={}
for i in range(chr_info.shape[0]):
    chrom=chr_info.iloc[i,0]
    df=regions[chrom]
    df_bed=pybedtools.BedTool.from_dataframe(df)
    intersections=df_bed.intersect(qtl_bed,wa=True,wb=True).to_dataframe()
    intersections['trait'] = intersections[10].str.split('\_QTL\_\(', 1).str[0]
    chrom_info=intersections.groupby([3,'trait']).size().reset_index(name='windowQTLcount').rename({3:'window'},axis=1)
    window_info=chrom_info.groupby(['window'])['windowQTLcount'].sum().reset_index(name='windowCount')
    info=chrom_info.merge(general_info)
    info=info.merge(window_info)
    enr=get_windows_enriched_with_qtls(info,chrom)
    if enr.shape[0]>0:
        windows_enriched_with_qtls[chrom]=enr
    df=intersections[[6,10,'trait']].rename({6:'chr_window',10:'QTL'},axis=1)
    windows_qtl_details[chrom]=df

flag=1
for key,value in windows_enriched_with_qtls.items():
    globals() ['df' + str(flag)]=value.copy()
    print(flag)
    flag=flag+1    
windows_enriched_with_qtls=pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,df23,df24,df25,df26,df27,df28,df29,df30])
windows_enriched_with_qtls=windows_enriched_with_qtls.merge(general_info)
df=windows_enriched_with_qtls['chr_window'].to_frame()
f=reformat_chr_window_names(df)
del windows_enriched_with_qtls['chr_window']
windows_enriched_with_qtls['chr_window']=f
windows_enriched_with_qtls=windows_enriched_with_qtls[['chr_window','trait','traits_Qtls','QTL_enrich_pvalue','AllQTLcount']]

flag=1
for key,value in windows_qtl_details.items():
    globals() ['df' + str(flag)]=value.copy()
    print(flag)
    flag=flag+1    
windows_qtl_details=pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,df23,df24,df25,df26,df27,df28,df29,df30])
df=windows_qtl_details['chr_window'].to_frame()
f=reformat_chr_window_names(df)
del windows_qtl_details['chr_window']
windows_qtl_details['chr_window']=f
windows_qtl_details=windows_qtl_details[['chr_window','trait','QTL']]


windows_enriched_with_qtls.to_csv('chromosomal_locations_enriched_with_QTL_traits',index = None, header=True,sep="\t")
windows_qtl_details.to_csv('QTLs_located_in_chromosomal_windows',index = None, header=True,sep="\t")