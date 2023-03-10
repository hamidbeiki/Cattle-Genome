#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:59:54 2021

@author: beiki
"""

from scipy import stats
import pandas as pd
import numpy as np
from scipy.stats import median_test
import matplotlib.pyplot as plt
import seaborn as sns

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
    
expressions=pd.read_csv("RSEM_FAANG_transcripts_modified_fpkm",header=0,sep="\t")
df=expressions.iloc[:,1:]
df=df.loc[(df!=0).any(1)]
ids=expressions['transcript_id']
expressions=df
expressions['transcript_id']=ids
tissues=list(expressions.columns)[0:-1]
biotypes=pd.read_csv("final_transcript_biotypes",names=['transcript_id','biotype'],sep='\t')
ncbi_known=pd.read_csv("gene_transcript_validation/all_validated_transcripts_by_ncbi",names=['transcript_id'],sep='\t')
ens_known=pd.read_csv("gene_transcript_validation/all_validated_transcripts_by_ensembl",names=['transcript_id'],sep='\t')
m=ncbi_known.merge(ens_known,how='outer',indicator=True)
m.loc[m['_merge']=='right_only'].shape[0]
m.loc[m['_merge']=='left_only'].shape[0]
#know_trs=m['transcript_id'].to_frame()
know_trs=m
know_trs=know_trs.merge(biotypes)
novel_trs=biotypes[~ biotypes.transcript_id.isin(know_trs.transcript_id)]
exon_info=pd.read_csv('transcripts_exons_info',names=['transcript_id','number_of_exons'],sep='\t')
know_trs=know_trs.merge(exon_info)
novel_trs=novel_trs.merge(exon_info)

tissues_info=[]
for tissue in tissues:
    df=expressions[['transcript_id',tissue]]
    df=df.merge(biotypes)
    df=df.loc[(df[tissue]>0)]
    k=df.merge(know_trs)
    n=df.merge(novel_trs)
    tissues_info.append([tissue,k.shape[0],n.shape[0]])

tissues_info=pd.DataFrame(tissues_info).rename({0:'Tissue',1:'known',2:'novel'},axis=1)
a=tissues_info['known'].to_frame().rename({'known':'info'},axis=1)
a['type']='known'
b=tissues_info['novel'].to_frame().rename({'novel':'info'},axis=1)
b['type']='novel'
df=pd.concat([a,b])
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_known-novel_transcript_number.png',dpi=300)
plt.close()  



results=[]
results2=[]
results3=[]
results4=[]
results5=[]
types=list(biotypes['biotype'].drop_duplicates())
test_set=['antisense_lncRNAs','NMDs']
for tissue in tissues:
    df=expressions[['transcript_id',tissue]]
    df=df.merge(biotypes)
    df=df.loc[df[tissue]>0]
    known=df.merge(know_trs)
    novel=df[~ df.transcript_id.isin(m.transcript_id)]
    protein_coding=df.loc[df['biotype']=="protein_coding_transcripts"]
    protein_coding['log']=np.log2(protein_coding[tissue])
    protein_coding_transcripts=protein_coding
    nc=df[~ df.transcript_id.isin(protein_coding.transcript_id)]
    nc['log']=np.log2(nc[tissue])
    sncRNAs=df.loc[df['biotype']=="sncRNAs"]
    sncRNAs['log']=np.log2(sncRNAs[tissue])
    NMDs=df.loc[df['biotype']=="NMDs"]
    NMDs['log']=np.log2(NMDs[tissue])
    non_stop_decays=df.loc[df['biotype']=="non_stop_decays"]
    non_stop_decays['log']=np.log2(non_stop_decays[tissue])
    long_intergenic_lncRNAs=df.loc[df['biotype']=="long_intergenic_lncRNAs"]
    long_intergenic_lncRNAs['log']=np.log2(long_intergenic_lncRNAs[tissue])
    long_intragenic_lncRNAs=df.loc[df['biotype']=="long_intragenic_lncRNAs"]
    long_intragenic_lncRNAs['log']=np.log2(long_intragenic_lncRNAs[tissue])
    sense_intronic_lncRNAs=df.loc[df['biotype']=="sense_intronic_lncRNAs"]
    sense_intronic_lncRNAs['log']=np.log2(sense_intronic_lncRNAs[tissue])
    antisense_lncRNAs=df.loc[df['biotype']=="antisense_lncRNAs"]
    antisense_lncRNAs['log']=np.log2(antisense_lncRNAs[tissue])
    stat, p, med, tbl = median_test(protein_coding[tissue],nc[tissue])
    results.append([tissue,protein_coding[tissue].median(),nc[tissue].median(),p])
    results2.append([tissue,protein_coding[tissue].median(),sncRNAs[tissue].median(),NMDs[tissue].median(),non_stop_decays[tissue].median(),long_intergenic_lncRNAs[tissue].median(),long_intragenic_lncRNAs[tissue].median(),sense_intronic_lncRNAs[tissue].median(),antisense_lncRNAs[tissue].median()])
    results4.append([tissue,protein_coding[tissue].shape[0],sncRNAs[tissue].shape[0],NMDs[tissue].shape[0],non_stop_decays[tissue].shape[0],long_intergenic_lncRNAs[tissue].shape[0],long_intragenic_lncRNAs[tissue].shape[0],sense_intronic_lncRNAs[tissue].shape[0],antisense_lncRNAs[tissue].shape[0]])   
    results5.append([tissue,known.shape[0],novel.shape[0]])
    stat, p, med, tbl = median_test(known[tissue],novel[tissue])
    results3.append([tissue,known[tissue].median(),novel[tissue].median(),p])
    print(tissue)
df=pd.DataFrame(results)
stat, p, med, tbl = median_test(df[1],df[2])
df2=pd.DataFrame(results2)
df2=df2.rename({0:'tissue',1:'protein_coding',2:'sncRNAs',3:'NMDs',4:'non_stop_decays',5:'long_intergenic_lncRNAs',6:'long_intragenic_lncRNAs',7:'sense_intronic_lncRNAs',8:'antisense_lncRNAs'},axis=1)
df2.to_csv("results.txt",index = False, header=True,sep="\t") 

a=df2['protein_coding'].to_frame().rename({'protein_coding':"expr"},axis=1)
a['type']='coding'
b=df2['sncRNAs'].to_frame().rename({'sncRNAs':"expr"},axis=1)
b['type']='sncRNAs'
c=df2['NMDs'].to_frame().rename({'NMDs':"expr"},axis=1)
c['type']='NMDs'
d=df2['non_stop_decays'].to_frame().rename({'non_stop_decays':"expr"},axis=1)
d['type']='non_stop_decays'
e=df2['long_intergenic_lncRNAs'].to_frame().rename({'long_intergenic_lncRNAs':"expr"},axis=1)
e['type']='long_intergenic_lncRNAs'
f=df2['long_intragenic_lncRNAs'].to_frame().rename({'long_intragenic_lncRNAs':"expr"},axis=1)
f['type']='long_intragenic_lncRNAs'
g=df2['sense_intronic_lncRNAs'].to_frame().rename({'sense_intronic_lncRNAs':"expr"},axis=1)
g['type']='sense_intronic_lncRNAs'
h=df2['antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs':"expr"},axis=1)
h['type']='antisense_lncRNAs'
df=pd.concat([a,b,c,d,e,f,g,h])
ax = sns.boxplot(x="type", y="expr", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()

df4=pd.DataFrame(results4)
df4=df4.rename({0:'tissue',1:'protein_coding',2:'sncRNAs',3:'NMDs',4:'non_stop_decays',5:'long_intergenic_lncRNAs',6:'long_intragenic_lncRNAs',7:'sense_intronic_lncRNAs',8:'antisense_lncRNAs'},axis=1)
df4['sum']=df4.loc[:,'protein_coding':'antisense_lncRNAs'].sum(axis=1)
df4_p=(df4.loc[:,'protein_coding':'antisense_lncRNAs'].div(df4["sum"], axis=0))*100

a=df4_p['protein_coding'].to_frame().rename({'protein_coding':"p"},axis=1)
a['type']='coding'
b=df4_p['sncRNAs'].to_frame().rename({'sncRNAs':"p"},axis=1)
b['type']='sncRNAs'
c=df4_p['NMDs'].to_frame().rename({'NMDs':"p"},axis=1)
c['type']='NMDs'
d=df4_p['non_stop_decays'].to_frame().rename({'non_stop_decays':"p"},axis=1)
d['type']='non_stop_decays'
e=df4_p['long_intergenic_lncRNAs'].to_frame().rename({'long_intergenic_lncRNAs':"p"},axis=1)
e['type']='long_intergenic_lncRNAs'
f=df4_p['long_intragenic_lncRNAs'].to_frame().rename({'long_intragenic_lncRNAs':"p"},axis=1)
f['type']='long_intragenic_lncRNAs'
g=df4_p['sense_intronic_lncRNAs'].to_frame().rename({'sense_intronic_lncRNAs':"p"},axis=1)
g['type']='sense_intronic_lncRNAs'
h=df4_p['antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs':"p"},axis=1)
h['type']='antisense_lncRNAs'
df=pd.concat([a,b,c,d,e,f,g,h])
ax = sns.boxplot(x="type", y="p", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
#plt.rcParams["figure.figsize"] = (36,12)
#plt.savefig('transcript_biotypes_percentage.jpg',dpi=300)
#plt.close()

a=df4['protein_coding'].to_frame().rename({'protein_coding':"p"},axis=1)
a['type']='coding'
b=df4['sncRNAs'].to_frame().rename({'sncRNAs':"p"},axis=1)
b['type']='sncRNAs'
c=df4['NMDs'].to_frame().rename({'NMDs':"p"},axis=1)
c['type']='NMDs'
d=df4['non_stop_decays'].to_frame().rename({'non_stop_decays':"p"},axis=1)
d['type']='non_stop_decays'
e=df4['long_intergenic_lncRNAs'].to_frame().rename({'long_intergenic_lncRNAs':"p"},axis=1)
e['type']='long_intergenic_lncRNAs'
f=df4['long_intragenic_lncRNAs'].to_frame().rename({'long_intragenic_lncRNAs':"p"},axis=1)
f['type']='long_intragenic_lncRNAs'
g=df4['sense_intronic_lncRNAs'].to_frame().rename({'sense_intronic_lncRNAs':"p"},axis=1)
g['type']='sense_intronic_lncRNAs'
h=df4['antisense_lncRNAs'].to_frame().rename({'antisense_lncRNAs':"p"},axis=1)
h['type']='antisense_lncRNAs'
df=pd.concat([a,b,c,d,e,f,g,h])
ax = sns.boxplot(x="type", y="p", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.rcParams["figure.figsize"] = (36,12)
plt.savefig('transcript_biotypes_number.jpg',dpi=300)
plt.close()

results=[]
for tissue in tissues:
    df=expressions[['transcript_id',tissue]]
    df=df.loc[df[tissue]>0]
    df=df.merge(biotypes)
    results.append([tissue,df.shape[0]])
results_df=pd.DataFrame(results).rename({0:'tissue',1:'number_transcripts'},axis=1)
df=results_df['number_transcripts'].to_frame().rename({'number_transcripts':'info'},axis=1)
df['type']='transcript'
ax = sns.boxplot(x="type", y="info", data=df, linewidth=0.6,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('tissue_comparision_transcript_number.png',dpi=300)
plt.close()


coding=biotypes.loc[biotypes['biotype']=="protein_coding_transcripts"]['transcript_id'].to_frame()
sncRNAs=biotypes.loc[biotypes['biotype']=="sncRNAs"]['transcript_id'].to_frame()
NMDs=biotypes.loc[biotypes['biotype']=="NMDs"]['transcript_id'].to_frame()
non_stop_decays=biotypes.loc[biotypes['biotype']=="non_stop_decays"]['transcript_id'].to_frame()
long_intergenic_lncRNAs=biotypes.loc[biotypes['biotype']=="long_intergenic_lncRNAs"]['transcript_id'].to_frame()
long_intragenic_lncRNAs=biotypes.loc[biotypes['biotype']=="long_intragenic_lncRNAs"]['transcript_id'].to_frame()
sense_intronic_lncRNAs=biotypes.loc[biotypes['biotype']=="sense_intronic_lncRNAs"]['transcript_id'].to_frame()
antisense_lncRNAs=biotypes.loc[biotypes['biotype']=="antisense_lncRNAs"]['transcript_id'].to_frame()

expressions=expressions.set_index('transcript_id')
coding_ex=expressions.loc[expressions.index.isin(coding.transcript_id)]
sncRNAs_ex=expressions.loc[expressions.index.isin(sncRNAs.transcript_id)]
NMDs_ex=expressions.loc[expressions.index.isin(NMDs.transcript_id)]
non_stop_decays_ex=expressions.loc[expressions.index.isin(non_stop_decays.transcript_id)]
long_intergenic_lncRNAs_ex=expressions.loc[expressions.index.isin(long_intergenic_lncRNAs.transcript_id)]
long_intragenic_lncRNAs_ex=expressions.loc[expressions.index.isin(long_intragenic_lncRNAs.transcript_id)]
sense_intronic_lncRNAs_ex=expressions.loc[expressions.index.isin(sense_intronic_lncRNAs.transcript_id)]
antisense_lncRNAs_ex=expressions.loc[expressions.index.isin(antisense_lncRNAs.transcript_id)]

coding_info=get_info(coding_ex,'coding')
sncRNAs_info=get_info(sncRNAs_ex,'sncRNAs')
NMDs_info=get_info(NMDs_ex,'NMDs')
non_stop_decays_info=get_info(non_stop_decays_ex,'non_stop_decays')
long_intergenic_lncRNAs_info=get_info(long_intergenic_lncRNAs_ex,'long_intergenic_lncRNAs')
long_intragenic_lncRNAs_info=get_info(long_intragenic_lncRNAs_ex,'long_intragenic_lncRNAs')
sense_intronic_lncRNAs_info=get_info(long_intragenic_lncRNAs_ex,'sense_intronic_lncRNAs')
antisense_lncRNAs_info=get_info(antisense_lncRNAs_ex,'antisense_lncRNAs')

coding_info['log']=np.log2(coding_info['median_expression'])
sncRNAs_info['log']=np.log2(sncRNAs_info['median_expression'])
NMDs_info['log']=np.log2(NMDs_info['median_expression'])
non_stop_decays_info['log']=np.log2(non_stop_decays_info['median_expression'])
long_intergenic_lncRNAs_info['log']=np.log2(long_intergenic_lncRNAs_info['median_expression'])
long_intragenic_lncRNAs_info['log']=np.log2(long_intragenic_lncRNAs_info['median_expression'])
sense_intronic_lncRNAs_info['log']=np.log2(sense_intronic_lncRNAs_info['median_expression'])
antisense_lncRNAs_info['log']=np.log2(antisense_lncRNAs_info['median_expression'])

a=coding_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='coding'
a['type2']='detected_tissues'
b=coding_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='coding'
b['type2']='median_expression'
c=sncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='sncRNAs'
c['type2']='detected_tissues'
d=sncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='sncRNAs'
d['type2']='median_expression'
e=NMDs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
e['type1']='NMDs'
e['type2']='detected_tissues'
f=NMDs_info['log'].to_frame().rename({'log':"info"},axis=1)
f['type1']='NMDs'
f['type2']='median_expression'
g=non_stop_decays_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
g['type1']='non_stop_decays'
g['type2']='detected_tissues'
h=non_stop_decays_info['log'].to_frame().rename({'log':"info"},axis=1)
h['type1']='non_stop_decays'
h['type2']='median_expression'
i=long_intergenic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
i['type1']='long_intergenic_lncRNAs'
i['type2']='detected_tissues'
j=long_intergenic_lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
j['type1']='long_intergenic_lncRNAs'
j['type2']='median_expression'
k=long_intragenic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
k['type1']='long_intragenic_lncRNAs'
k['type2']='detected_tissues'
l=long_intragenic_lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
l['type1']='long_intragenic_lncRNAs'
l['type2']='median_expression'
m=sense_intronic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
m['type1']='sense_intronic_lncRNAs'
m['type2']='detected_tissues'
n=sense_intronic_lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
n['type1']='sense_intronic_lncRNAs'
n['type2']='median_expression'
o=antisense_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
o['type1']='antisense_lncRNAs'
o['type2']='detected_tissues'
p=antisense_lncRNAs_info['log'].to_frame().rename({'log':"info"},axis=1)
p['type1']='antisense_lncRNAs'
p['type2']='median_expression'
#df=pd.concat([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p])
#ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
df=pd.concat([a,c,e,g,i,k,m,o])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)
plt.show()
df2=pd.concat([b,d,f,h,j,l,n,p])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df2, linewidth=2.5,showfliers=True,fliersize=0.1)
plt.savefig('test.png',dpi=300)
plt.close() 

results=[]
types=list(df2.columns)[1:]
for i in types:
    r=[]
    for j in types:
        stat, p1 = stats.ttest_ind(df2[i],df2[j],equal_var=True,alternative='greater')
        stat, p2 = stats.ttest_ind(df2[i],df2[j],equal_var=True,alternative='less')
        if p1<=0.05:
            r.append(1)
        elif p2<=0.05:
            r.append(-1)
        elif p1>0.05 and p2>0.05:
            r.append(0)
    results.append(r)
df=pd.DataFrame(results)
df.columns=types
df.index=types            
df.to_csv("results.txt",index = True, header=True,sep="\t") 

ex=expressions
df=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy
detections=(pd.DataFrame(np.count_nonzero(df, axis=1))-1)

ids=pd.DataFrame(expressions['transcript_id'])
ids.index = np.arange(0, ids.shape[0])
df=pd.concat([ids,detections],ignore_index=True,axis=1).rename({0:'transcript_id',1:"number_of_expressed_tissue"},axis=1)
m=df.merge(biotypes)
m.loc[m['biotype']=='sncRNAs']
m.loc[m['biotype']=='sncRNAs']['number_of_expressed_tissue'].median()
rest=df[~ df.transcript_id.isin(m.transcript_id)]

df=pd.DataFrame(results3)
df=df.rename({0:'tissue',1:'novel',2:'known',3:'p-value'},axis=1)
df['-log2(p_value)']=(-np.log2(df["p-value"]+1e-306)) 
df=df.reset_index()
ll = plt.scatter(df["known"], df["novel"], marker='o', c=df["-log2(p_value)"],cmap=plt.get_cmap("jet"))#, s=df["-log2(p_value)"]
plt.colorbar(label="-log2(p-value) for the comparision of expression values" "\n" "between known and novel transcripts in each tissue" "\n")
plt.plot(df["known"], df["known"],linestyle = 'dashed',color = 'grey')
plt.xlabel("median expression of known transcripts (RPKM)")
plt.ylabel("median expression of novel transcripts (RPKM)")
#plt.show()
plt.savefig('known_novel_transcript_expression_comparision.png')
plt.close()

df.plot(kind="scatter", x="known", y="novel", alpha=0.4,
    figsize=(10,7),
    c="index", cmap=plt.get_cmap("jet"), colorbar=True,
    sharex=False,s=df["-log2(p_value)"])#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"
plt.show()

known=m.merge(biotypes)
novel=df[~ df.transcript_id.isin(m.transcript_id)]

expressions=expressions.set_index('transcript_id')
remove=pd.read_csv('remove_transcripts',names=['transcript_id'])
known=pd.read_csv('all-known-transcripts-main',names=['transcript_id'])
expressions=expressions.loc[~ expressions.index.isin(remove.transcript_id)]

known_ex=expressions.loc[expressions.index.isin(known.transcript_id)]
novel_ex=expressions.loc[~ expressions.index.isin(known.transcript_id)]

known_info=get_info(known_ex,'transcript_id')
novel_info=get_info(novel_ex,'transcript_id')

known_info['log']=np.log2(known_info['median_expression'])
novel_info['log']=np.log2(novel_info['median_expression'])

a=known_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='known'
a['type2']='detected_tissues'
b=known_info['log'].to_frame().rename({'log':"info"},axis=1)
b['type1']='known'
b['type2']='median_expression'
c=novel_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='novel'
c['type2']='detected_tissues'
d=novel_info['log'].to_frame().rename({'log':"info"},axis=1)
d['type1']='novel'
d['type2']='median_expression'
df=pd.concat([a,c])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()
df=pd.concat([b,d])
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=True,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.savefig('test.png',dpi=300)
plt.close()    
    
stat, p = stats.mannwhitneyu(known_info['median_expression'],novel_info['median_expression'],alternative='greater')
stat, p = stats.mannwhitneyu(known_info['detected_tissues'],novel_info['detected_tissues'],alternative='greater')
    

#ex = ex.replace(0, np.NaN)
#df=pd.DataFrame.to_numpy(ex)#convert dataframe to numpy
#df_row_medians = np.nanmedian(df, axis=1)# use use numpy.nanmedian() to get medain of no-zero values per row https://stackoverflow.com/questions/33217636/mean-calculation-in-pandas-excluding-zeros
#df_row_means = np.nanmean(df, axis=1)
""" ALTERNATIVE PLOT
coding=biotypes.loc[biotypes['biotype']=="protein_coding_transcripts"]['transcript_id'].to_frame()
sncRNAs=biotypes.loc[biotypes['biotype']=="sncRNAs"]['transcript_id'].to_frame()
NMDs=biotypes.loc[biotypes['biotype']=="NMDs"]['transcript_id'].to_frame()
non_stop_decays=biotypes.loc[biotypes['biotype']=="non_stop_decays"]['transcript_id'].to_frame()
long_intergenic_lncRNAs=biotypes.loc[biotypes['biotype']=="long_intergenic_lncRNAs"]['transcript_id'].to_frame()
long_intragenic_lncRNAs=biotypes.loc[biotypes['biotype']=="long_intragenic_lncRNAs"]['transcript_id'].to_frame()
sense_intronic_lncRNAs=biotypes.loc[biotypes['biotype']=="sense_intronic_lncRNAs"]['transcript_id'].to_frame()
antisense_lncRNAs=biotypes.loc[biotypes['biotype']=="antisense_lncRNAs"]['transcript_id'].to_frame()

expressions=expressions.set_index('transcript_id')
coding_ex=expressions.loc[expressions.index.isin(coding.transcript_id)]
sncRNAs_ex=expressions.loc[expressions.index.isin(sncRNAs.transcript_id)]
NMDs_ex=expressions.loc[expressions.index.isin(NMDs.transcript_id)]
non_stop_decays_ex=expressions.loc[expressions.index.isin(non_stop_decays.transcript_id)]
long_intergenic_lncRNAs_ex=expressions.loc[expressions.index.isin(long_intergenic_lncRNAs.transcript_id)]
long_intragenic_lncRNAs_ex=expressions.loc[expressions.index.isin(long_intragenic_lncRNAs.transcript_id)]
sense_intronic_lncRNAs_ex=expressions.loc[expressions.index.isin(sense_intronic_lncRNAs.transcript_id)]
antisense_lncRNAs_ex=expressions.loc[expressions.index.isin(antisense_lncRNAs.transcript_id)]

coding_info=get_info(coding_ex,'coding')
sncRNAs_info=get_info(sncRNAs_ex,'sncRNAs')
NMDs_info=get_info(NMDs_ex,'NMDs')
non_stop_decays_info=get_info(non_stop_decays_ex,'non_stop_decays')
long_intergenic_lncRNAs_info=get_info(long_intergenic_lncRNAs_ex,'long_intergenic_lncRNAs')
long_intragenic_lncRNAs_info=get_info(long_intragenic_lncRNAs_ex,'long_intragenic_lncRNAs')
sense_intronic_lncRNAs_info=get_info(long_intragenic_lncRNAs_ex,'sense_intronic_lncRNAs')
antisense_lncRNAs_info=get_info(antisense_lncRNAs_ex,'antisense_lncRNAs')

a=coding_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
a['type1']='coding'
a['type2']='detected_tissues'
b=coding_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
b['type1']='coding'
b['type2']='median_expression'
c=sncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
c['type1']='sncRNAs'
c['type2']='detected_tissues'
d=sncRNAs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
d['type1']='sncRNAs'
d['type2']='median_expression'
e=NMDs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
e['type1']='NMDs'
e['type2']='detected_tissues'
f=NMDs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
f['type1']='NMDs'
f['type2']='median_expression'
g=non_stop_decays_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
g['type1']='non_stop_decays'
g['type2']='detected_tissues'
h=non_stop_decays_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
h['type1']='non_stop_decays'
h['type2']='median_expression'
i=long_intergenic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
i['type1']='long_intergenic_lncRNAs'
i['type2']='detected_tissues'
j=long_intergenic_lncRNAs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
j['type1']='long_intergenic_lncRNAs'
j['type2']='median_expression'
k=long_intragenic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
k['type1']='long_intragenic_lncRNAs'
k['type2']='detected_tissues'
l=long_intragenic_lncRNAs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
l['type1']='long_intragenic_lncRNAs'
l['type2']='median_expression'
m=sense_intronic_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
m['type1']='sense_intronic_lncRNAs'
m['type2']='detected_tissues'
n=sense_intronic_lncRNAs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
n['type1']='sense_intronic_lncRNAs'
n['type2']='median_expression'
o=antisense_lncRNAs_info['detected_tissues'].to_frame().rename({'detected_tissues':"info"},axis=1)
o['type1']='antisense_lncRNAs'
o['type2']='detected_tissues'
p=antisense_lncRNAs_info['median_expression'].to_frame().rename({'median_expression':"info"},axis=1)
p['type1']='antisense_lncRNAs'
p['type2']='median_expression'
df=pd.concat([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p])
#df['info'][ df['info'] > 100] =100
ax = sns.boxplot(x="type1", y="info", hue="type2", data=df, linewidth=2.5,showfliers=False,fliersize=0.1)#showfliers=True,fliersize=0.5
plt.show()

"""
'''
    flag=0
    for test in test_set:
        test=globals()[test]
        name=test_set[flag]
        flag=flag+1
        rest=[x for x in types if x !=name and x!='protein_coding_transcripts']
        flag2=0
        for target in rest:
            name2=rest[flag2]
            flag2=flag2+1
            target=globals()[target]
            stat, p, med, tbl = median_test(test[tissue],target[tissue])
            results2.append([tissue,name,name2,test[tissue].median(),target[tissue].median(),p])
    
'''



'''
types=list(biotypes['biotype'].drop_duplicates())
test_set=['antisense_lncRNAs','non_stop_decays']
protein_coding_transcripts=protein_coding
results=[]
flag=0
for test in test_set:
    test=globals()[test]
    name=test_set[flag]
    flag=flag+1
    rest=[x for x in types if x !=name]
    flag2=0
    for target in rest:
        name2=rest[flag2]
        flag2=flag2+1
        target=globals()[target]
#        stat, p, med, tbl = median_test(test['max'],target['max'])
        stat, p = stats.mannwhitneyu(test['mean'],target['mean'],alternative='greater')
        results.append([name,name2,p,test['mean'].mean(),target['mean'].mean()])
df=pd.DataFrame(results)
df=df.rename({0:"biotype_A",1:"biotype_B",2:"p-value",3:"median_expr_A",4:"median_expr_B"},axis=1)
df['-log2(p_value)']=(-np.log2(df["p-value"]+1e-163))      
df.biotype_A = pd.Categorical(df.biotype_A)
df['code'] = df.biotype_A.cat.codes

df1=df.loc[df['biotype_A']=="antisense_lncRNAs"]
df2=df.loc[df['biotype_A']=="non_stop_decays"]

colors = ['b', 'c', 'y', 'm', 'r']
ll = plt.scatter(df1["median_expr_A"], df1["median_expr_B"], s=df1["-log2(p_value)"], marker='o', color=colors[0])
l = plt.scatter(df2["median_expr_A"], df2["median_expr_B"], s=df2["-log2(p_value)"], marker='o', color=colors[1])
gll = plt.scatter([], [], s=3, marker='o', color='#555555')
gl = plt.scatter([], [], s=30, marker='o', color='#555555')
ga = plt.scatter([],[], s=300, marker='o', color='#555555')
                 
plt.legend((gll,gl,ga),
       ('3', '30', '300'),
       scatterpoints=1,
       loc='upper center',
       ncol=1,
       fontsize=12)                 



fig=df.plot(kind="scatter", x="median_expr_A", y="median_expr_B", alpha=0.4,
    figsize=(10,7),
    c="code", cmap=plt.get_cmap("jet"), colorbar=False,
    sharex=False,s=df["-log2(p_value)"])#to change dot size based on 3UTR_GC_content add: s=m3["3UTR_GC_content"],label="3UTR_GC_content"

'''







'''


plt.legend(fig,
       ('5', '50', '500'),
       scatterpoints=1,
       loc='lower left',
       ncol=1,
       fontsize=8)

plt.show()

plt.savefig('test.png')
plt.close()
'''
"""
antisense_lncRNAs['max'].median()
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],sncRNAs['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],NMDs['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],non_stop_decays['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],long_intergenic_lncRNAs['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],long_intragenic_lncRNAs['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],sense_intronic_lncRNAs['max'])
stat, p, med, tbl = median_test(antisense_lncRNAs['max'],sense_intronic_lncRNAs['max'])

non_stop_decays['max'].median()
stat, p, med, tbl = median_test(non_stop_decays['max'],sncRNAs['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],NMDs['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],non_stop_decays['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],long_intergenic_lncRNAs['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],long_intragenic_lncRNAs['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],sense_intronic_lncRNAs['max'])
stat, p, med, tbl = median_test(non_stop_decays['max'],sense_intronic_lncRNAs['max'])
"""
