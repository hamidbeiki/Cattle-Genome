#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 08:29:39 2021

@author: beiki
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:33:56 2020

@author: beiki
"""

## This program works with python3
import pandas as pd
import sys
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

#intervals = [[3,9],[2,6],[8,10],[15,18]]
#intervals = [[1,2],[4,5],[7,8],[3,4],[7,10],[1,2],[4,5],[6,10]]

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

def get_geneLength(gene_structure):
    df=pd.DataFrame(gene_structure)
    df['length']=df[1]-df[0]
    return(df['length'].sum())


def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def paralle_length_calculator(chunk):
    gene_length=pd.DataFrame([])
    for i in range(len(chunk)):
        gene=chunk[i][0]
        gene_structure=merge_intervals(chunk[i][1][0])
        length=get_geneLength(gene_structure)
        gene_length=gene_length.append(pd.DataFrame({'gene':gene,'exonic_length':length},index=[0]), ignore_index=True)
    return(gene_length)        

def parallle_format_df(chunk):
    """Yield successive n-sized chunks from grouped_df."""
    df=pd.DataFrame([])
    for i in range(len(chunk)):
        df=df.append(chunk[i][1])
    a=df.iloc[:,1:3].apply(list,axis=1)
    df['list']=a
    grouped=df.groupby('gene').agg({'list':lambda x: list(x)})
    return(grouped)

def get_list(d):
    a=d.values# instead of "pd.DataFrame.to_numpy(d)" in python 2.7
    l=a.tolist()
    return(l)    

def get_tr_lenght(l):
    length=0
    for i in l:
        ex_le=i[1]-i[0]+1
        length=length+ex_le
    return(length)

''' run single CPU '''
df1=pd.read_csv(sys.argv[1],sep=' ',names=["gene_id","ex_start","ex_end","strand"])#awk '$3=="exon"{print $10, $4,$5, $7}' combined.gff | sed 's/\"//g;s/\;//g' | sort | uniq > out1
df2=pd.read_csv(sys.argv[2],sep=' ',names=["tr_id","ex_start","ex_end","strand"])# awk '$3=="exon"{print $12, $4,$5,$7}' combined.gff | sed 's/\"//g;s/\;//g' | sort | uniq > out2
tr_biotypes=pd.read_csv("final_transcript_biotypes",names=['transcript_id','tr_biotype'],sep='\t')
gene_biotypes=pd.read_csv('final_gene_biotypes',names=['gene_id','gene_biotype'],sep='\t')
## Step1 get transcript length
groups=df1.groupby('gene_id') 
gene_info={}
for gene, group in groups:
    d=group[['ex_start', 'ex_end']]
    l=get_list(d)
    res=merge_intervals(l)
    gene_info[gene]=get_tr_lenght(res) 
gene_length=pd.DataFrame.from_dict(gene_info, orient='index').reset_index().rename({'index':'gene_id',0:'gene_exonic_length'},axis=1)#gene_length=pd.read_csv('combined_exonic_gene_length',header=0,sep='\t')
gene_length.to_csv('combined_exonic_gene_length',index = None, header=True,sep="\t")
## Step2 get transcript length
groups=df2.groupby('tr_id') 
tr_info={}
for tr, group in groups:
    d=group[['ex_start', 'ex_end']]
    l=get_list(d)
    res=merge_intervals(l)
    tr_info[tr]=get_tr_lenght(res) 

tr_length=pd.DataFrame.from_dict(tr_info, orient='index').reset_index().rename({'index':'tr_id',0:'tr_exonic_length'},axis=1)#tr_length=pd.read_csv('combined_exonic_transcript_length',header=0,sep='\t')
tr_length.to_csv('combined_exonic_transcript_length',index = None, header=True,sep="\t")


pos=df1.loc[df1['strand']=="+"]
neg=df1.loc[df1['strand']=="-"]
pos_start=pos.groupby(['gene_id'], sort=False)['ex_start'].min().reset_index().rename({'ex_start':'gene_start'},axis=1)
neg_start=pos.groupby(['gene_id'], sort=False)['ex_end'].max().reset_index().rename({'ex_end':'gene_start'},axis=1)
genes_start=pd.concat([pos_start,neg_start])

pos=df2.loc[df2['strand']=="+"]
neg=df2.loc[df2['strand']=="-"]
pos_start=pos.groupby(['tr_id'], sort=False)['ex_start'].min().reset_index().rename({'ex_start':'tr_start'},axis=1)
pos_start['ex_start']=pos_start['tr_start']
pos_start=pos_start.merge(df2,on=['tr_id','ex_start'])
neg_start=neg.groupby(['tr_id'], sort=False)['ex_end'].max().reset_index().rename({'ex_end':'tr_start'},axis=1)
neg_start['ex_end']=neg_start['tr_start']
neg_start=neg_start.merge(df2,on=['tr_id','ex_end'])
trs_start=pd.concat([pos_start,neg_start])


groups=df1.groupby('gene_id') 
gene_info={}
for gene, group in groups:
    d=group[['ex_start', 'ex_end']]
    l=get_list(d)
    res=merge_intervals(l)
    gene_info[gene]=res

d=df1[['gene_id','strand']].drop_duplicates()
d.set_index("gene_id", drop=True, inplace=True)
strand_dict=d.to_dict(orient="index")

d=gene_length.copy()
d.set_index("gene_id", drop=True, inplace=True)
gene_length_dict=d.to_dict(orient="index")


l1=list(trs_start['tr_id'])
l1_2=list(trs_start['tr_start'])
l1_3=list(trs_start['ex_start'])
l1_4=list(trs_start['ex_end'])

info={}#position transcript start site on gene exonic length
for i in range(len(l1)):
    tr=l1[i]
    s=l1_2[i]
    a=l1_3[i]
    b=l1_4[i]
    gene="PB" + "." + tr.split('.')[1]
    strand=strand_dict[gene]['strand']
    l=gene_info[gene]
    l2=[]
    if strand=="+":
        for ex in l:
            ex_s,ex_e=ex
            if ex_s<s and ex_e<s:
                l2.append(ex)
            elif ex_s<s and ex_e>s:
                l2.append([ex_s,s])
            elif ex_s==s and ex_e>s:
                l2.append([a,b])
    elif strand=="-":
        for ex in l:
            ex_s,ex_e=ex
            if ex_s>s and ex_e>s:
                l2.append(ex)
            elif ex_s>s and ex_e<s:
                l2.append([ex_s,s])
            elif ex_e==s and ex_s<s:
                l2.append([a,b])
    if len(l2)>0:
        res=merge_intervals(l2)
        length=get_tr_lenght(res)
        if length !=gene_length_dict[gene]['gene_exonic_length']:
            info[tr]=length
        else:
            info[tr]=0
    else:
        info[tr]=0             
l2=[]   
for ex in l:
    ex_s,ex_e=ex
    if ex_s<s and ex_e<s:
        l2.append(ex)
    elif ex_s<s and ex_e>s:
        l2.append([ex_s,s])
    elif ex_s==s and ex_e>s:
        l2.append([a,b])
                
positions1=pd.DataFrame.from_dict(info,orient='index').reset_index().rename({'index':"tr_id",0:"strart_pos_on_exonic_gene_length"},axis=1)#transcript start site on gene exonic length

pos_end=pos.groupby(['tr_id'], sort=False)['ex_end'].max().reset_index().rename({'ex_end':'tr_end'},axis=1)
pos_end['ex_end']=pos_end['tr_end']
pos_end=pos_end.merge(df2,on=['tr_id','ex_end'])
neg_end=neg.groupby(['tr_id'], sort=False)['ex_start'].min().reset_index().rename({'ex_start':'tr_end'},axis=1)
neg_end['ex_start']=neg_end['tr_end']
neg_end=neg_end.merge(df2,on=['tr_id','ex_start'])
trs_end=pd.concat([pos_end,neg_end])

l1=list(trs_end['tr_id'])
l1_2=list(trs_end['tr_end'])
l1_3=list(trs_end['ex_start'])
l1_4=list(trs_end['ex_end'])

info2={}#position transcript terminal site on gene exonic length
for i in range(len(l1)):
    tr=l1[i]
    s=l1_2[i]
    a=l1_3[i]
    b=l1_4[i]
    gene="PB" + "." + tr.split('.')[1]
    strand=strand_dict[gene]['strand']
    l=gene_info[gene]
    l2=[]
    if strand=="+":
        for ex in l:
            ex_s,ex_e=ex
            if ex_s<s and ex_e<s:
                l2.append(ex)
            elif ex_s<s and ex_e>s:
                l2.append([ex_s,s])
            elif ex_s==s and ex_e>s:
                l2.append([a,b])
            else:
                l2.append([ex_s,ex_e])
    elif strand=="-":
        for ex in l:
            ex_s,ex_e=ex
            if ex_s>s and ex_e>s:
                l2.append(ex)
            elif ex_s>s and ex_e<s:
                l2.append([ex_s,s])
            elif ex_e==s and ex_s<s:
                l2.append([a,b])
            else:
                l2.append([ex_s,ex_e])
    tr_s=info[tr]
    res=merge_intervals(l2)
    length=get_tr_lenght(res)
    if length==tr_s:
        info2[tr]=gene_length_dict[gene]['gene_exonic_length']
    else:
        info2[tr]=length

l2=[]
for ex in l:
    ex_s,ex_e=ex
    if ex_s<s and ex_e<s:
        l2.append(ex)
    elif ex_s<s and ex_e>s:
        l2.append([ex_s,s])
    elif ex_s==s and ex_e>s:
        l2.append([a,b])
    else:
        l2.append([ex_s,ex_e])          
   
positions2=pd.DataFrame.from_dict(info2,orient='index').reset_index().rename({'index':"tr_id",0:"end_pos_on_exonic_gene_length"},axis=1)#transcript start site on gene exonic length

positions=positions1.merge(positions2)

g=positions['tr_id'].str.split('.',n = 2, expand = True)[[0,1]]
genes=(g[0] + "." + g[1]).to_frame().rename({0:'gene_id'},axis=1)
positions['gene_id']=genes

gene_length=gene_length.rename({'gene':'gene_id'},axis=1)
m=positions.merge(gene_length).rename({'tr_id':'transcript_id'},axis=1)
m['start%']=(m['strart_pos_on_exonic_gene_length']/m['gene_exonic_length'])*100
m['end%']=(m['end_pos_on_exonic_gene_length']/m['gene_exonic_length'])*100
m['gene_length%']=m['end%'] - m['start%']# % of gene length that covered by each transcript
m=m.merge(tr_biotypes)
m=m.merge(gene_biotypes)
m.to_csv('details_of_transcript_position_on_genes',index = None, header=True,sep="\t")

df=m.loc[(m['gene_biotype']=='protein-coding')]
protein_coding=df.loc[df['tr_biotype']=='protein_coding_transcripts']
NMDs=df.loc[df['tr_biotype']=='NMDs']
non_stop_decays=df.loc[df['tr_biotype']=='non_stop_decays']
sncRNAs=df.loc[df['tr_biotype']=='sncRNAs']
d=pd.concat([protein_coding,NMDs,non_stop_decays,sncRNAs])
long_intragenic_lncRNAs=df.loc[~ df.transcript_id.isin(d.transcript_id)]

protein_coding['start%'].median()
NMDs['start%'].median()
non_stop_decays['start%'].median()
sncRNAs['start%'].median()
long_intragenic_lncRNAs['start%'].median()

a=protein_coding['gene_length%'].to_frame().rename({'gene_length%':'percentage'},axis=1)
a['type1']='protein_coding'
a['type2']='gene_length%'
b=protein_coding['start%'].to_frame().rename({'start%':'percentage'},axis=1)
b['type1']='protein_coding'
b['type2']='start%'
c=protein_coding['end%'].to_frame().rename({'end%':'percentage'},axis=1)
c['type1']='protein_coding'
c['type2']='end%'
d=NMDs['gene_length%'].to_frame().rename({'gene_length%':'percentage'},axis=1)
d['type1']='NMDs'
d['type2']='gene_length%'
e=NMDs['start%'].to_frame().rename({'start%':'percentage'},axis=1)
e['type1']='NMDs'
e['type2']='start%'
f=NMDs['end%'].to_frame().rename({'end%':'percentage'},axis=1)
f['type1']='NMDs'
f['type2']='end%'
g=non_stop_decays['gene_length%'].to_frame().rename({'gene_length%':'percentage'},axis=1)
g['type1']='non_stop_decays'
g['type2']='gene_length%'
h=non_stop_decays['start%'].to_frame().rename({'start%':'percentage'},axis=1)
h['type1']='non_stop_decays'
h['type2']='start%'
j=non_stop_decays['end%'].to_frame().rename({'end%':'percentage'},axis=1)
j['type1']='non_stop_decays'
j['type2']='end%'
k=sncRNAs['gene_length%'].to_frame().rename({'gene_length%':'percentage'},axis=1)
k['type1']='sncRNAs'
k['type2']='gene_length%'
l=sncRNAs['start%'].to_frame().rename({'start%':'percentage'},axis=1)
l['type1']='sncRNAs'
l['type2']='start%'
m=sncRNAs['end%'].to_frame().rename({'end%':'percentage'},axis=1)
m['type1']='sncRNAs'
m['type2']='end%'
n=long_intragenic_lncRNAs['gene_length%'].to_frame().rename({'gene_length%':'percentage'},axis=1)
n['type1']='long_intragenic_lncRNAs'
n['type2']='gene_length%'
o=long_intragenic_lncRNAs['start%'].to_frame().rename({'start%':'percentage'},axis=1)
o['type1']='long_intragenic_lncRNAs'
o['type2']='start%'
p=long_intragenic_lncRNAs['end%'].to_frame().rename({'end%':'percentage'},axis=1)
p['type1']='long_intragenic_lncRNAs'
p['type2']='end%'
df2=pd.concat([a,b,c,d,e,f,g,h,j,k,l,m,n,o,p])
ax = sns.boxplot(x="type1", y="percentage", hue="type2", data=df2, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()


''' run parallel 
## Step1 get gene length
df=pd.read_csv(sys.argv[1],sep=' ',names=["gene","ex_start","ex_end"])#awk '$3=="exon"{print $10, $4,$5, $7}' combined.gff | sed 's/\"//g;s/\;//g' | sort | uniq > out1
df2=pd.read_csv(sys.argv[2],sep=' ',names=["tr_id","ex_start","ex_end"])# awk '$3=="exon"{print $12, $4,$5}' combined.gff | sed 's/\"//g;s/\;//g' | sort | uniq > out2
grouped_df=df.groupby('gene') 
p = Pool() 
dflist = []
for group in grouped_df:
    dflist.append(group)
processed_values= p.map(parallle_format_df, get_chunks(dflist, 100))
grouped=pd.concat(processed_values)
## Step2 get exonic gene length
p = Pool() 
items = list(grouped.iterrows())
processed_values= p.map(paralle_length_calculator, get_chunks(items, 100))
gene_length=pd.concat(processed_values)#gene_length=pd.read_csv('combined_exonic_gene_length',header=0,sep='\t')
'''