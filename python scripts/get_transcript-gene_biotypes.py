#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:25:58 2019

@author: beiki
"""

""" this program has been written on python 3.7.5, except parallel processing, all other parts work on python 2.7.5
On Ceres: please issue "module load miniconda, before start python: 'source activate /KEEP/cattle_genome_assemblies/bioconda' 
On Nove issue "module load python/3.6.3-u4oaxsb" and "module load py-biopython/1.70-py3-wos466g" and "module load py-pandas/0.23.4-py3-sx6iffy"
""" 

from Bio import SeqIO
import re 
import pandas as pd
from multiprocessing import Pool
import time
import sys
import collections


class coding_evaluator:
    def is_lncrna(self,startposition,stopposition):
        if len(startposition)==0:
            self.lncrna=True
        else:
            self.lncrna=False
        return(self.lncrna)

def get_orfs(record,sequence,startposition,stopposition):
    min_ORFtranscript_coverage=0.35 #based on Ensembl pipeline, please see: Functional impacts of non-coding RNA processing on enhancer activity and target gene expression.pdf
    position=pd.DataFrame([])
    flag=pd.DataFrame([])
    position2=pd.DataFrame([])
    biotype=pd.DataFrame([])
    maxwithstop=pd.DataFrame([])
    maxwithstop=maxwithstop.append(pd.DataFrame({'length':0},index=[0]), ignore_index=True)
    maxwithNOstop=pd.DataFrame([])
    maxwithNOstop=maxwithNOstop.append(pd.DataFrame({'length':0},index=[0]), ignore_index=True)
    tr_length=len(record.seq)
    tr_end=tr_length-1
    for start in startposition:
        flag=[x for x in stopposition if x>start and (x-start)%3==0]
        if len(flag)>0 and maxwithstop.loc[maxwithstop['length']<(flag[0]-start)].shape[0]>0:
            maxwithstop=maxwithstop.append(pd.DataFrame({'length':flag[0]-start},index=[0]), ignore_index=True)
            tr_coverage=(flag[0]-start+1)/tr_length# ratio of transcript covered by cds region
            position=position.append(pd.DataFrame({'ID': record.id, 'start': start, 'stop': flag[0], 'length': flag[0]-start, 'tr_coverage': tr_coverage, 'orf_id':  str(record.id) + str("_") + str(start) + str("_") + str(flag[0]), 'protein': str(record.seq[start:flag[0]].translate(table=1)), 'CDS_sequence' : str(record.seq[start:flag[0]])}, index=[0]), ignore_index=True)
        elif len(flag)==0 and maxwithstop.iloc[0,0]==0 and position2.shape[0]==0:
            end=get_closest_dividable(start,tr_end,3)
            maxwithNOstop=maxwithNOstop.append(pd.DataFrame({'length':end-start},index=[0]), ignore_index=True)
            tr_coverage=(end-start+1)/tr_length
            position2=position2.append(pd.DataFrame({'ID': record.id, 'start': start, 'stop': end, 'length': maxwithNOstop.iloc[-1,0],'tr_coverage': tr_coverage, 'orf_id':  str(record.id) + str("_") + str(start) + str("_") + str(end), 'protein': str(record.seq[start:end].translate(table=1)), 'CDS_sequence' : str(record.seq[start:end])}, index=[0]), ignore_index=True)
    if position.shape[0]>0:
        position=position[['ID','start','stop','length','tr_coverage','orf_id','protein', 'CDS_sequence']]
    elif position2.shape[0]>0:
        position2=position2[['ID','start','stop','length','tr_coverage','orf_id','protein', 'CDS_sequence']]
    if position.shape[0]>0 and position.loc[position['tr_coverage']>=min_ORFtranscript_coverage].shape[0]>0:
        position=position.loc[position['tr_coverage']>=min_ORFtranscript_coverage].sort_values(by=['length'],ascending=False)
        biotype=biotype.append(pd.DataFrame({'type': 'coding'}, index=[0]), ignore_index=True)
        return(position.iloc[0:3,0:8],biotype)
    elif position.shape[0]>0 and position.loc[position['tr_coverage']>=min_ORFtranscript_coverage].shape[0]==0:
        position=position.sort_values(by=['length'],ascending=False)
        biotype=biotype.append(pd.DataFrame({'type': 'lowCoverage_coding'}, index=[0]), ignore_index=True) 
        return(position.iloc[0:3,0:8],biotype)
    elif position.shape[0]==0 and position2.loc[position2['tr_coverage']>=min_ORFtranscript_coverage].shape[0]>0:
        biotype=biotype.append(pd.DataFrame({'type': 'coding_with_no_stop'}, index=[0]), ignore_index=True)
        return(position2.iloc[0:3,0:8],biotype)
    elif position.shape[0]==0 and position2.loc[position2['tr_coverage']>=min_ORFtranscript_coverage].shape[0]==0:
        biotype=biotype.append(pd.DataFrame({'type': 'low_coverage_coding_with_no_stop'}, index=[0]), ignore_index=True)
        return(position2.iloc[0:3,0:8],biotype)

def get_closest_dividable(start,end,target):
    for i in reversed(range(end+1)):
        if (i-start)%target==0:
            answer=i
            break
    return(answer)
            
def get_biotypes(record,minimum_number_of_aminoAcids,start_codon_type):
    lncRNAs=pd.DataFrame([])
    best_orfs_cds=pd.DataFrame([])
    best_orfs_short_cds=pd.DataFrame([])
    no_stop_best_orfs=pd.DataFrame([])
    no_stop_best_orfs_short_cds=pd.DataFrame([])
    best_orfs_lowCoverage_cds=pd.DataFrame([])
    best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
    short_lowCoverage_cds=pd.DataFrame([])
    no_stop_short_lowCoverage_cds=pd.DataFrame([])
    sequence=str(record.seq)
    if start_codon_type=="canonical":
        startposition =[m.start() for m in re.finditer('(ATG)', sequence,  re.IGNORECASE)]#without Kozak motif https://www.cell.com/trends/biochemical-sciences/pdf/S0968-0004(19)30146-X.pdf
#        startposition =[m.start() for m in re.finditer('(GCCACCATGG)|(GCCGCCATGG)', sequence,  re.IGNORECASE)]#with Kozak motif
#        startposition= [x+6 for x in startposition]#with Kozak motif
    elif start_codon_type=="noncanonical":
        startposition =[m.start() for m in re.finditer('(ATG)|(CTG)', sequence,  re.IGNORECASE)]# near-cognate start codon list has been gotten from "http://genesdev.cshlp.org/content/31/17/1717.full.pdf" . also, it is expected that lncRNAs predominantly use near-cognate start codons such as CUG compared to protein coding mRNAs (Coding functions of “noncoding” RNAs)  
#        startposition =[m.start() for m in re.finditer('(GCCACCATGG)|(GCCGCCATGG)|(GCCACCCTGG)|(GCCGCCCTGG)', sequence,  re.IGNORECASE)]
#        startposition= [x+6 for x in startposition]
    stopposition=[m.start() for m in re.finditer('(TAA)|(TAG)|(TGA)', sequence,  re.IGNORECASE)]
    seq=coding_evaluator()
    lncRNA_info=seq.is_lncrna(startposition,stopposition)
    if lncRNA_info==False:
        cds_info=get_orfs(record,sequence,startposition,stopposition)
    if (lncRNA_info==True):
        lncRNAs=lncRNAs.append(pd.DataFrame({'transcript': record.id}, index=[0]), ignore_index=True)
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='coding_with_no_stop') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]>0:
        no_stop_best_orfs=pd.concat([no_stop_best_orfs,cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)]])
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='coding_with_no_stop') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]==0:
        no_stop_best_orfs_short_cds=pd.concat([no_stop_best_orfs_short_cds,cds_info[0]])
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='coding') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]>0:
        best_orfs_cds=pd.concat([best_orfs_cds,cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)]])
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='coding') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]==0:
        best_orfs_short_cds=pd.concat([best_orfs_short_cds,cds_info[0]])
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='lowCoverage_coding') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]>0:
        best_orfs_lowCoverage_cds=pd.concat([best_orfs_lowCoverage_cds,cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)]])
    elif (lncRNA_info==False and cds_info[1].iloc[0,0]=='low_coverage_coding_with_no_stop') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]>0:
        best_orfs_no_stop_lowCoverage_cds=pd.concat([best_orfs_no_stop_lowCoverage_cds,cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)]])
    elif lncRNA_info==False and cds_info[1].iloc[0,0]=='lowCoverage_coding' and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]==0:
        short_lowCoverage_cds=pd.concat([short_lowCoverage_cds,cds_info[0]])
    elif lncRNA_info==False and cds_info[1].iloc[0,0]=='low_coverage_coding_with_no_stop' and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*3)].shape[0]==0:
        no_stop_short_lowCoverage_cds=pd.concat([no_stop_short_lowCoverage_cds,cds_info[0]])      
    return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
 

def sort_df(df):
    df=df.sort_values(by=['length'],ascending=False)
    if df.shape[0]<=3:
        return(df)
    elif df.shape[0]>3:
        return(df.head(n=3))

def run_parallel(record):
    lncRNAs=pd.DataFrame([])
    best_orfs_cds=pd.DataFrame([])
    best_orfs_short_cds=pd.DataFrame([])
    no_stop_best_orfs=pd.DataFrame([])
    no_stop_best_orfs_short_cds=pd.DataFrame([])
    best_orfs_lowCoverage_cds=pd.DataFrame([])
    best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
    short_lowCoverage_cds=pd.DataFrame([])
    no_stop_short_lowCoverage_cds=pd.DataFrame([])
    minimum_number_of_aminoAcids=44 ##The smallest human protein is 44 amino acids but it could be an abortive translation from the 5' UTR of another mRNA. The smallest functional polypeptide is glutathione with only three amino acids. https://www.science20.com/princerain/blog/whats_biggest_and_whats_smallest AND https://www.ncbi.nlm.nih.gov/pmc/articles/PMC528919/; ALSO Proteins generally contain from 50 to 1000 amino acid residues (AAs) per polypeptide chain https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3864261/; ALSO To form even the simplest transmembrane α-helix (TMH) structure, 30 amino acids are needed "Peptides encoded by noncoding genes: challenges needed" 
    results=get_biotypes(record,minimum_number_of_aminoAcids,"canonical")
    lncRNAs=pd.concat([lncRNAs,results[0]])
    no_stop_best_orfs=pd.concat([no_stop_best_orfs,results[1]])
    best_orfs_cds=pd.concat([best_orfs_cds,results[2]])
    best_orfs_short_cds=pd.concat([best_orfs_short_cds,results[3]])
    no_stop_best_orfs_short_cds=pd.concat([no_stop_best_orfs_short_cds,results[4]])
    best_orfs_lowCoverage_cds=pd.concat([best_orfs_lowCoverage_cds,results[5]])
    best_orfs_no_stop_lowCoverage_cds=pd.concat([ best_orfs_no_stop_lowCoverage_cds,results[6]])
    short_lowCoverage_cds=pd.concat([ short_lowCoverage_cds,results[7]])
    no_stop_short_lowCoverage_cds=pd.concat([no_stop_short_lowCoverage_cds,results[8]])
    if lncRNAs.shape[0]>0:
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=lncRNAs.append(pd.DataFrame({'transcript': record.id}, index=[0]), ignore_index=True)
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]>0:
        best_orfs_cds=sort_df(best_orfs_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=sort_df(best_orfs_lowCoverage_cds)
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)        
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]>0:
        best_orfs_short_cds=sort_df(best_orfs_short_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=sort_df(short_lowCoverage_cds) 
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)                 
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=sort_df(no_stop_best_orfs)
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=sort_df(best_orfs_no_stop_lowCoverage_cds)
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])        
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs_short_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=sort_df(no_stop_best_orfs_short_cds)
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])        
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs_short_cds.shape[0]==0 and no_stop_short_lowCoverage_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=sort_df(no_stop_short_lowCoverage_cds)       
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         

def run_parallel2(record):
    lncRNAs=pd.DataFrame([])
    best_orfs_cds=pd.DataFrame([])
    best_orfs_short_cds=pd.DataFrame([])
    no_stop_best_orfs=pd.DataFrame([])
    no_stop_best_orfs_short_cds=pd.DataFrame([])
    best_orfs_lowCoverage_cds=pd.DataFrame([])
    best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
    short_lowCoverage_cds=pd.DataFrame([])
    no_stop_short_lowCoverage_cds=pd.DataFrame([])
    minimum_number_of_aminoAcids=100 ## all transcripts passes to this function already have canonical orfs<44 aa, this function specidicly set to explore lncRNAs wich mostle have non_cononical start sites with orfs<100 aa "Coding functions of “noncoding” RNAs" and "Pervasive functional translation of noncanonical human open reading frames"
    results=get_biotypes(record,minimum_number_of_aminoAcids,"noncanonical")
    lncRNAs=pd.concat([lncRNAs,results[0]])
    no_stop_best_orfs=pd.concat([no_stop_best_orfs,results[1]])
    best_orfs_cds=pd.concat([best_orfs_cds,results[2]])
    best_orfs_short_cds=pd.concat([best_orfs_short_cds,results[3]])
    no_stop_best_orfs_short_cds=pd.concat([no_stop_best_orfs_short_cds,results[4]])
    best_orfs_lowCoverage_cds=pd.concat([best_orfs_lowCoverage_cds,results[5]])
    best_orfs_no_stop_lowCoverage_cds=pd.concat([ best_orfs_no_stop_lowCoverage_cds,results[6]])
    short_lowCoverage_cds=pd.concat([ short_lowCoverage_cds,results[7]])
    no_stop_short_lowCoverage_cds=pd.concat([no_stop_short_lowCoverage_cds,results[8]])
    if lncRNAs.shape[0]>0:
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=lncRNAs.append(pd.DataFrame({'transcript': record.id}, index=[0]), ignore_index=True)
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]>0:
        best_orfs_cds=sort_df(best_orfs_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=sort_df(best_orfs_lowCoverage_cds)
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)        
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]>0:
        best_orfs_short_cds=sort_df(best_orfs_short_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=sort_df(short_lowCoverage_cds) 
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)                 
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]>0:
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=sort_df(no_stop_best_orfs)
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=sort_df(best_orfs_no_stop_lowCoverage_cds)
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])        
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs_short_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=sort_df(no_stop_best_orfs_short_cds)
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=pd.DataFrame([])        
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)
    elif best_orfs_cds.shape[0]==0 and best_orfs_lowCoverage_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and short_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs.shape[0]==0 and best_orfs_no_stop_lowCoverage_cds.shape[0]==0 and no_stop_best_orfs_short_cds.shape[0]==0 and no_stop_short_lowCoverage_cds.shape[0]>0:     
        best_orfs_short_cds=pd.DataFrame([])
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        no_stop_best_orfs_short_cds=pd.DataFrame([])
        best_orfs_lowCoverage_cds=pd.DataFrame([])
        best_orfs_no_stop_lowCoverage_cds=pd.DataFrame([])
        short_lowCoverage_cds=pd.DataFrame([])
        no_stop_short_lowCoverage_cds=sort_df(no_stop_short_lowCoverage_cds)       
        return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds,no_stop_best_orfs_short_cds,best_orfs_lowCoverage_cds,best_orfs_no_stop_lowCoverage_cds,short_lowCoverage_cds,no_stop_short_lowCoverage_cds)         

def nested_dict():
    return collections.defaultdict(nested_dict)        

def get_genes(transcripts):
    global tr_to_gene
    genes=transcripts.merge(tr_to_gene).iloc[0:,1].to_frame()
    genes.sort_values("gene_id", inplace = True)
    genes=genes.drop_duplicates(subset ="gene_id")
    return(genes)

def get_last_splice_junction(input,tr):
    """ no need to set minimum_exon_length as it has been applied in previous steps"""
    introns_info=pd.DataFrame([])
    last_splice_junction=pd.DataFrame([])
    df=pd.DataFrame([v.split(';') for k, v in input.items()],columns=['chr',
                         'ex_start', 'ex_end', 'strand'])# on python 2.7, like on Ceres: change "items" to "iteritems"
    df['ex_start']=df['ex_start'].astype(int)
    df['ex_end']=df['ex_end'].astype(int)
    df['ex_length']=df['ex_end'] - df['ex_start'] +1
    if df.shape[0]>1 and df.iloc[0,3]=="+":
        df=df.sort_values(by='ex_start', ascending=True)
        exonic_length=0
        for i in range(df.shape[0]-1):
            exonic_length=exonic_length + df.iloc[i,2] - df.iloc[i,1] + 1
            introns_info=introns_info.append(pd.DataFrame({'intron_end':exonic_length+1}, index=[0]), ignore_index=True)
    elif df.shape[0]>1 and df.iloc[0,3]=="-":
        df=df.sort_values(by='ex_start', ascending=False)
        exonic_length=0
        for i in range(df.shape[0]-1):
            exonic_length=exonic_length + df.iloc[i,2] - df.iloc[i,1] + 1
            introns_info=introns_info.append(pd.DataFrame({'intron_end':exonic_length+1}, index=[0]), ignore_index=True)
        introns_info = introns_info[['intron_end']]
    info=introns_info.iloc[-1,0]
    last_splice_junction=last_splice_junction.append(pd.DataFrame({'ID': tr,'last_splice_junction':info}, index=[0]), ignore_index=True)
    return(last_splice_junction)

def paralle_chunk_input(chunk):
    last_splice_junction_info=pd.DataFrame([])
    for tr in range(len(chunk)):
        if len(chunk[tr][1])>1:
            info=get_last_splice_junction(chunk[tr][1],chunk[tr][0])#RNAbased method
            last_splice_junction_info=last_splice_junction_info.append(info)
    return(last_splice_junction_info)

def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]    

def get_NMDs(input_file,last_splice_junction_info):
    min_downstream_distanse=50 # based on my pig paper
    spliced_best_orfs_cds=input_file.merge(last_splice_junction_info)
    NMD_orfs=spliced_best_orfs_cds[(spliced_best_orfs_cds['stop'] - spliced_best_orfs_cds['last_splice_junction'])<-min_downstream_distanse]
    df2 = pd.DataFrame(NMD_orfs['ID'].value_counts().values, index=NMD_orfs['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)### numner of NMD orfs detected per transcript
    df1 = pd.DataFrame(input_file['ID'].value_counts().values, index=input_file['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)### numner of orfs detected per transcript
    NMD_transcripts=df2.merge(df1).iloc[0:,0].to_frame()
    all_orf_ids=input_file.iloc[0:,5].to_frame()
    NMD_orf_ids=NMD_orfs.iloc[0:,5].to_frame()
    nonNMDs_orf_ids=all_orf_ids[~all_orf_ids.orf_id.isin(NMD_orf_ids.orf_id)]
    nonNMDs_orf=input_file.merge(nonNMDs_orf_ids)
    return(NMD_transcripts,nonNMDs_orf)
    
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

def get_unique_genes_across_dfs(df_target,df_rest):
    unique_genes=df_target.merge(df_rest,how="outer",indicator=True).loc[lambda x : x['_merge']=='left_only'].iloc[0:,0].to_frame()
    return(unique_genes)
    
    
''' STEP1 find biotypes based on ORFs '''    
start = time.time() 
"""In Python 2.x and 3.0, 3.1 and 3.2, multiprocessing.Pool() objects are not context managers. You cannot use them in a with statement.
    Only in Python 3.3 and up can you use them as such. """
input_fasta=sys.argv[1]# your transcripts fasta file, like input_fasta='Cerebral_Cortex_final.collapsed_transcripts.fa' 
input_gff=sys.argv[2]# your gff file, like input_gff='Cerebral_Cortex_final.collapsed.gff'
tissue=input_gff.split("_final")[0]
with open(input_fasta) as handle, Pool() as pool:
    results = pool.map(
        run_parallel,
        (record for record in SeqIO.parse(handle, 'fasta')),
        chunksize=30,
    )
end = time.time()
print('total time (s)= ' + str(end-start)) 
for i in range(0,9):
    if i==0:
        outputs = [result[i] for result in results]
        lncRNAs = pd.concat(outputs)
    elif i==1:
        outputs = [result[i] for result in results]
        no_stop_best_orfs = pd.concat(outputs)
    elif i==2:
        outputs = [result[i] for result in results]
        best_orfs_cds = pd.concat(outputs)
    elif i==3:
        outputs = [result[i] for result in results]
        best_orfs_short_cds = pd.concat(outputs)
    elif i==4:
        outputs = [result[i] for result in results]
        no_stop_best_orfs_short_cds = pd.concat(outputs)
    elif i==5:
        outputs = [result[i] for result in results]
        best_orfs_lowCoverage_cds = pd.concat(outputs)
    elif i==6:
        outputs = [result[i] for result in results]
        best_orfs_no_stop_lowCoverage_cds = pd.concat(outputs)
    elif i==7:
        outputs = [result[i] for result in results]
        short_lowCoverage_cds = pd.concat(outputs)
    elif i==8:
        outputs = [result[i] for result in results]
        no_stop_short_lowCoverage_cds=pd.concat(outputs)

lncRNAs=lncRNAs.rename({'transcript':'ID'},axis=1)
short_lowCoverage_cds=short_lowCoverage_cds.rename({'transcript':'ID'},axis=1)
""" codidng transcripts """        
#best_orfs_cds.iloc[0:,0].nunique()
""" questionable coding transcripts? """
#best_orfs_lowCoverage_cds.iloc[0:,0].nunique()
""" non_stop_decay transcripts """
#no_stop_best_orfs.iloc[0:,0].nunique()
#best_orfs_no_stop_lowCoverage_cds.iloc[0:,0].nunique()
#no_stop_best_orfs_short_cds.iloc[0:,0].nunique()
#no_stop_short_lowCoverage_cds.iloc[0:,0].nunique()
""" non-codidng transcripts """
#lncRNAs.shape[0]
#best_orfs_short_cds.iloc[0:,0].nunique()
#short_lowCoverage_cds.iloc[0:,0].nunique()

''' STEP2 detect coding transcripts, their related NMDs and ncRNAs (non_coding_candidates)'''
ex_coords = nested_dict()
for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'exon':
        tid=raw[-1].split('; ')[1].split()[1][1:-2]
        exnumber=raw[3]
        ex_coords[tid][exnumber] = "{0};{1};{2};{3}".format(raw[0], raw[3], raw[4], raw[6])
            
items = list(ex_coords.items())
p = Pool() 
start = time.time()
# out=get_chunks(items, 70) breaks items into 320 chunks (ilen(out)) and 70 genesets in each chunk ##
processed_values= p.map(paralle_chunk_input, get_chunks(items, 400))
end = time.time()
print('total time (s)= ' + str(end-start)) 
last_splice_junction_info=pd.concat(processed_values)
coding_candidates=["best_orfs_cds", "best_orfs_lowCoverage_cds"]
NMD_transcripts=pd.DataFrame([])
nonNMD_orfs=pd.DataFrame([])
for input_file in coding_candidates:
    info=get_NMDs(globals()[input_file],last_splice_junction_info)
    NMD_transcripts=pd.concat([NMD_transcripts, info[0]])
    nonNMD_orfs=pd.concat([nonNMD_orfs, info[1]]) ## this includes both spliced and unspliced cds

tnc_orfs=pd.DataFrame([])
nc_info=pd.concat([best_orfs_short_cds,short_lowCoverage_cds])
candidates=nc_info.loc[nc_info['length']<9*3].iloc[0:,0].to_frame().drop_duplicates(subset="ID") ## shortest non-coding gene (NGC) peptite length is 9 aa: Peptides encoded by noncoding genes- challenges and perspectives (Nature)
candidates_info=nc_info.loc[nc_info['length']<9*3]
df1 = pd.DataFrame(nc_info['ID'].value_counts().values, index=nc_info['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)
df2 = pd.DataFrame(candidates_info['ID'].value_counts().values, index=candidates_info['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)
lncRNAs_2=df1.merge(df2).iloc[0:,0].to_frame().drop_duplicates(subset="ID")
tnc_candidates=nc_info[~nc_info.orf_id.isin(candidates_info.orf_id)]## translated non-coding (TNC) trnascript
lncRNA_transcripts=pd.concat([lncRNAs,lncRNAs_2])
info=get_NMDs(tnc_candidates,last_splice_junction_info)
NMD_transcripts=pd.concat([NMD_transcripts, info[0]])
tnc_orfs=pd.concat([tnc_orfs, info[1]]) ## this includes both spliced and unspliced transcripts
tnc_transcripts=tnc_orfs.iloc[0:,0].to_frame().drop_duplicates(subset="ID")
coding_transcripts=nonNMD_orfs.iloc[0:,0].to_frame().drop_duplicates(subset="ID")
non_coding_tr_with_cds=tnc_transcripts
non_coding_tr_without_cds=lncRNA_transcripts
non_coding_candidates=pd.concat([non_coding_tr_with_cds,non_coding_tr_without_cds])
#non_stop_decay_transcripts=pd.DataFrame([])
#if no_stop_best_orfs_short_cds.shape[0]>0:
#    non_stop_decay_transcripts=pd.concat([non_stop_decay_transcripts,no_stop_best_orfs_short_cds.iloc[0:,0].to_frame()])
#if no_stop_short_lowCoverage_cds.shape[0]>0:
#    non_stop_decay_transcripts=pd.concat([non_stop_decay_transcripts,no_stop_short_lowCoverage_cds.iloc[0:,0].to_frame()])
#if best_orfs_no_stop_lowCoverage_cds.shape[0]>0:
#    non_stop_decay_transcripts=pd.concat([non_stop_decay_transcripts,best_orfs_no_stop_lowCoverage_cds.iloc[0:,0].to_frame()])
#if no_stop_best_orfs.shape[0]>0:
#    non_stop_decay_transcripts=pd.concat([non_stop_decay_transcripts,no_stop_best_orfs.iloc[0:,0].to_frame()])
    
#non_stop_decay_transcripts=non_stop_decay_transcripts.drop_duplicates(subset="ID")    
non_stop_decay_transcripts=pd.concat([no_stop_best_orfs_short_cds.iloc[0:,0],no_stop_short_lowCoverage_cds.iloc[0:,0],best_orfs_no_stop_lowCoverage_cds.iloc[0:,0],no_stop_best_orfs.iloc[0:,0]]).to_frame().drop_duplicates(subset="ID")
''' total transcripts '''
#coding_transcripts.shape[0]+non_coding_tr_with_cds.shape[0]+non_coding_tr_without_cds.shape[0]+non_stop_decay_transcripts.shape[0]+NMD_transcripts.shape[0]

""" STEP3 adjust ncRNA ORF based on non-cononocal start sited and detect NMD transcripts among them """
with open(input_fasta) as handle, Pool() as pool:
    results = pool.map(
        run_parallel2,
        (record for record in SeqIO.parse(handle, 'fasta') if non_coding_candidates.loc[non_coding_candidates["ID"]==record.id].shape[0]>0),
        chunksize=30,
    )
    
for i in range(0,9):
    if i==0:
        outputs = [result[i] for result in results]
        lncRNAs = pd.concat(outputs)
    elif i==1:
        outputs = [result[i] for result in results]
        no_stop_best_orfs2 = pd.concat(outputs)
    elif i==2:
        outputs = [result[i] for result in results]
        best_orfs_cds = pd.concat(outputs)
    elif i==3:
        outputs = [result[i] for result in results]
        best_orfs_short_cds = pd.concat(outputs)
    elif i==4:
        outputs = [result[i] for result in results]
        no_stop_best_orfs_short_cds2 = pd.concat(outputs)
    elif i==5:
        outputs = [result[i] for result in results]
        best_orfs_lowCoverage_cds = pd.concat(outputs)
    elif i==6:
        outputs = [result[i] for result in results]
        best_orfs_no_stop_lowCoverage_cds2 = pd.concat(outputs)
    elif i==7:
        outputs = [result[i] for result in results]
        short_lowCoverage_cds = pd.concat(outputs)
    elif i==8:
        outputs = [result[i] for result in results]
        no_stop_short_lowCoverage_cds2=pd.concat(outputs)

lncRNAs=lncRNAs.rename({'transcript':'ID'},axis=1)
short_lowCoverage_cds=short_lowCoverage_cds.rename({'transcript':'ID'},axis=1)
candidates=["no_stop_best_orfs_short_cds2","best_orfs_cds","best_orfs_short_cds","no_stop_best_orfs_short_cds2",
            "best_orfs_lowCoverage_cds","best_orfs_no_stop_lowCoverage_cds2","short_lowCoverage_cds","no_stop_short_lowCoverage_cds2"]
NMD_transcripts2=pd.DataFrame([])
nonNMD_orfs2=pd.DataFrame([])
for input_file in candidates:
    if globals()[input_file].shape[0]>0:
        info=get_NMDs(globals()[input_file],last_splice_junction_info)
        NMD_transcripts2=pd.concat([NMD_transcripts2, info[0]])
        nonNMD_orfs2=pd.concat([nonNMD_orfs2, info[1]]) ## this includes both spliced and unspliced cds

tnc_orfs=pd.DataFrame([])
nc_info=pd.concat([best_orfs_short_cds,short_lowCoverage_cds,best_orfs_cds,best_orfs_lowCoverage_cds])
nc_info=nc_info[~nc_info.ID.isin(NMD_transcripts2.ID)]
candidates=nc_info.loc[nc_info['length']<9*3].iloc[0:,0].to_frame().drop_duplicates(subset="ID") ## shortest non-coding gene (NGC) peptite length is 9 aa: Peptides encoded by noncoding genes- challenges and perspectives (Nature)
candidates_info=nc_info.loc[nc_info['length']<9*3]
df1 = pd.DataFrame(nc_info['ID'].value_counts().values, index=nc_info['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)
df2 = pd.DataFrame(candidates_info['ID'].value_counts().values, index=candidates_info['ID'].value_counts().index, columns=['Count']).reset_index().rename({'index':'ID'},axis=1)
lncRNAs_2=df1.merge(df2).iloc[0:,0].to_frame().drop_duplicates(subset="ID")
lncRNA_transcripts=pd.concat([lncRNAs,lncRNAs_2])
tnc_orfs=nc_info[~nc_info.orf_id.isin(candidates_info.orf_id)]## translated non-coding (TNC) trnascript
tnc_transcripts=tnc_orfs.iloc[0:,0].to_frame().drop_duplicates(subset="ID")
non_coding_tr_with_cds=tnc_transcripts
non_coding_tr_without_cds=lncRNA_transcripts
non_coding_transcripts=pd.concat([non_coding_tr_with_cds,non_coding_tr_without_cds])
non_stop_decay_transcripts2=pd.DataFrame([])

if no_stop_best_orfs_short_cds2.shape[0]>0:
    non_stop_decay_transcripts2=pd.concat([non_stop_decay_transcripts2,no_stop_best_orfs_short_cds2.iloc[0:,0].to_frame()])
if no_stop_short_lowCoverage_cds2.shape[0]>0:
    non_stop_decay_transcripts2=pd.concat([non_stop_decay_transcripts2,no_stop_short_lowCoverage_cds2.iloc[0:,0].to_frame()])
if best_orfs_no_stop_lowCoverage_cds2.shape[0]>0:
    non_stop_decay_transcripts2=pd.concat([non_stop_decay_transcripts2,best_orfs_no_stop_lowCoverage_cds2.iloc[0:,0].to_frame()])
if no_stop_best_orfs2.shape[0]>0:
    non_stop_decay_transcripts2=pd.concat([non_stop_decay_transcripts2,no_stop_best_orfs2.iloc[0:,0].to_frame()])
    
non_stop_decay_transcripts=pd.concat([non_stop_decay_transcripts,non_stop_decay_transcripts2]).drop_duplicates(subset="ID")

NMD_transcripts=pd.concat([NMD_transcripts,NMD_transcripts2])
nonNMD_orfs2=pd.concat([nonNMD_orfs,nonNMD_orfs2])

""" total transcriptps check """
## NMD_transcripts.shape[0]+non_stop_decay_transcripts.shape[0]+non_coding_transcripts.shape[0]+coding_transcripts.shape[0]

""" STEP4 get genes biotypes """
tr_to_gene=pd.DataFrame([])   
for line in open(input_gff):
    raw=line.strip().split("\t")
    if raw[2] == 'transcript':
        tid=raw[-1].split('; ')[1].split()[1][1:-2]
        gid=raw[-1].split('; ')[0].split()[1][1:-1]
        tr_to_gene=tr_to_gene.append(pd.DataFrame({'ID': tid, 'gene_id': gid}, index=[0]), ignore_index=True)
        
coding_genes=set(get_genes(coding_transcripts).iloc[0:,0])
genes_with_non_coding_tr_with_cds=set(get_genes(non_coding_tr_with_cds).iloc[0:,0])
genes_with_non_coding_tr_without_cds=set(get_genes(non_coding_tr_without_cds).iloc[0:,0])
genes_with_non_stop_decay_transcripts=set(get_genes(non_stop_decay_transcripts).iloc[0:,0])
genes_with_NMD_transcripts=set(get_genes(NMD_transcripts).iloc[0:,0])
sets=[coding_genes,genes_with_non_coding_tr_with_cds,genes_with_non_coding_tr_without_cds,genes_with_non_stop_decay_transcripts,genes_with_NMD_transcripts]
genes_intersetcs=get_venn_sections(sets)
''' to get intersect classes '''
#for i in range(len(genes_intersetcs)):
#    print(i, genes_intersetcs[i][0])
genes_with_just_coding_tr=pd.DataFrame(list(genes_intersetcs[30][1])).rename({0:'gene_ID'},axis=1)
genes_with_just_non_coding_with_cds_tr=pd.DataFrame(list(genes_intersetcs[29][1])).rename({0:'gene_ID'},axis=1)
genes_with_just_non_coding_without_cds_tr=pd.DataFrame(list(genes_intersetcs[27][1])).rename({0:'gene_ID'},axis=1)
genes_with_just_non_stop_decay_tr=pd.DataFrame(list(genes_intersetcs[23][1])).rename({0:'gene_ID'},axis=1)
genes_with_just_NMD_tr=pd.DataFrame(list(genes_intersetcs[15][1])).rename({0:'gene_ID'},axis=1)
genes_with_all_type_tr=pd.DataFrame(list(genes_intersetcs[0][1])).rename({0:'gene_ID'},axis=1)
genes_with_NMD_or_non_stop_decay_tr=pd.concat([get_genes(NMD_transcripts).iloc[0:,0],get_genes(non_coding_tr_without_cds).iloc[0:,0]]).to_frame().drop_duplicates(subset="gene_id")

genes_with_noncoding_tr=set(get_genes(pd.concat([non_coding_tr_with_cds,non_coding_tr_without_cds,non_stop_decay_transcripts,NMD_transcripts])).iloc[0:,0])
sets=[coding_genes,genes_with_noncoding_tr]
genes_intersetcs=get_venn_sections(sets)
genes_with_coding_noncoding_tr=pd.DataFrame(list(genes_intersetcs[0][1])).rename({0:'gene_ID'},axis=1)
genes_with_just_noncoding_tr=pd.DataFrame(list(genes_intersetcs[1][1])).rename({0:'gene_ID'},axis=1)


''' STEP4 export results '''
coding_transcripts['biotype']='protein_coding'
non_coding_tr_with_cds['biotype']='non_coding_with_cds'
non_coding_tr_without_cds['biotype']='non_coding_without_cds'
non_stop_decay_transcripts['biotype']='non_stop_decay'
NMD_transcripts['biotype']='nonsense_mediated_decay'
transcript_biotypes=pd.concat([coding_transcripts,non_coding_tr_with_cds,non_coding_tr_without_cds,non_stop_decay_transcripts,NMD_transcripts])

genes_with_just_coding_tr['biotype']="genes_with_just_coding_tr"
genes_with_just_noncoding_tr['biotype']="genes_with_just_noncoding_tr"
genes_with_coding_noncoding_tr['biotype']="genes_with_coding_noncoding_tr"
genes_biotypes=pd.concat([genes_with_just_coding_tr,genes_with_just_noncoding_tr,genes_with_coding_noncoding_tr])

transcript_biotypes.to_csv(tissue + '_transcript_biotypes',index = None, header=True,sep="\t")
genes_biotypes.to_csv(tissue + '_genes_biotypes',index = None, header=True,sep="\t")
nonNMD_orfs.to_csv(tissue + '_coding_nonNMD_orfs',index = None, header=True,sep="\t")
tnc_orfs.to_csv(tissue + '_tnc_nonNMD_orfs',index = None, header=True,sep="\t")
genes_with_NMD_or_non_stop_decay_tr.to_csv(tissue + '_genes_with_NMD_or_non_stop_decay_tr',index = None, header=True,sep="\t")
genes_with_coding_noncoding_tr.to_csv(tissue + '_genes_with_coding_noncoding_tr',index = None, header=True,sep="\t")

#near_cognates=[]
#for i in range(tnc_orfs.shape[0]):
#    if tnc_orfs.iloc[i,7][0:3]!="ATG":#check the start codon
#        near_cognates.append(nonNMD_orfs.iloc[i,5])
############################## END ##############################
""" to run on single CPU 
handle=open("Cerebral_Cortex_final.collapsed_transcripts.fa","r")
lncRNAs=pd.DataFrame([])
best_orfs_cds=pd.DataFrame([])
best_orfs_short_cds=pd.DataFrame([])
no_stop_best_orfs=pd.DataFrame([])
processed=pd.DataFrame([])
minimum_number_of_aminoAcids=44
for record in SeqIO.parse(handle,"fasta"):
    results=get_biotypes(record,minimum_number_of_aminoAcids)
    lncRNAs=pd.concat([lncRNAs,results[0]])
    no_stop_best_orfs=pd.concat([no_stop_best_orfs,results[1]])
    best_orfs_cds=pd.concat([best_orfs_cds,results[2]])
    best_orfs_short_cds=pd.concat([best_orfs_short_cds,results[3]])
    if lncRNAs.shape[0]==3:
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
        lncRNAs=lncRNAs.append(pd.DataFrame({'transcript': record.id}, index=[0]), ignore_index=True)
    elif best_orfs_cds.shape[0]>0:
        best_orfs_cds=sort_df(best_orfs_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
    elif lncRNAs.shape[0]<3 and best_orfs_cds.shape[0]==0 and best_orfs_short_cds.shape[0]>0:
        best_orfs_short_cds=sort_df(best_orfs_short_cds)
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        no_stop_best_orfs=pd.DataFrame([])
    elif lncRNAs.shape[0]<3 and best_orfs_cds.shape[0]==0 and best_orfs_short_cds.shape[0]==0 and no_stop_best_orfs.shape[0]>0:
        no_stop_best_orfs=sort_df(no_stop_best_orfs)
        lncRNAs=pd.DataFrame([])
        best_orfs_cds=pd.DataFrame([])
        best_orfs_short_cds=pd.DataFrame([])
"""
        

################################################
"""slow program, just for edujaction purposes
class slow_coding_evaluator:
    def is_lncrna(self,startposition,stopposition):
        if len(startposition)==0:
            self.lncrna=True
        else:
            self.lncrna=False
        return(self.lncrna)
        
    def get_orfs(self,ID,tr_end,startposition,stopposition):
        self.position=pd.DataFrame([])
        self.flag=pd.DataFrame([])
        self.starts_with_no_stop_position=pd.DataFrame([])
        self.type=pd.DataFrame([])
        for start in startposition:
            for stop in stopposition:
              if (stop>start) and ((stop-start) % 3==0):
                  self.position=self.position.append(pd.DataFrame({'ID': ID, 'start': start, 'stop': stop, 'length': stop-start}, index=[0]), ignore_index=True)
                  break
              elif (stop<start):
                  self.flag=self.flag.append(pd.DataFrame({'check': 'yes'}, index=[0]), ignore_index=True)
                  self.starts_with_no_stop_position=self.starts_with_no_stop_position.append(pd.DataFrame({'ID': ID, 'start': start, 'stop': tr_end, 'length': tr_end-start}, index=[0]), ignore_index=True)
              elif (stop>start) and not ((stop-start) % 3==0):
                  self.flag=self.flag.append(pd.DataFrame({'check': 'yes'}, index=[0]), ignore_index=True)
                  self.starts_with_no_stop_position=self.starts_with_no_stop_position.append(pd.DataFrame({'ID': ID, 'start': start, 'stop': tr_end, 'length': tr_end-start}, index=[0]), ignore_index=True)
        if self.position.shape[0]<2:
            self.position=self.position.sort_values(by=['length'],ascending=False)
            self.type=self.type.append(pd.DataFrame({'type': 'coding'}, index=[0]), ignore_index=True)
            return(self.position,self.type)
        elif self.position.shape[0]>2:
            self.position=self.position.sort_values(by=['length'],ascending=False)
            self.type=self.type.append(pd.DataFrame({'type': 'coding'}, index=[0]), ignore_index=True)
            return(self.position.head(n=3),self.type)
        elif (self.position.shape[0]==0) and (self.flag.loc[self.flag['check']=='yes'].shape[0]>0):
            self.starts_with_no_stop_position=self.starts_with_no_stop_position.sort_values(by=['length'],ascending=False)
            self.type=self.type.append(pd.DataFrame({'type': 'coding_with_no_stop'}, index=[0]), ignore_index=True)
            if self.starts_with_no_stop_position.shape[0]<2:
                return(self.starts_with_no_stop_position,self.type)
            elif self.starts_with_no_stop_position.shape[0]>2:
                return(self.starts_with_no_stop_position.head(n=3),self.type)


def run_singleFrame_at_time(record):
    lncRNAs=pd.DataFrame([])
    best_orfs_cds=pd.DataFrame([])
    best_orfs_short_cds=pd.DataFrame([])
    no_stop_best_orfs=pd.DataFrame([])
    minimum_number_of_aminoAcids=44
    first=str(record.seq)[0:]
    startposition =[m.start() for m in re.finditer('ATG', first)]
    stopposition=[m.start() for m in re.finditer('(TAA)|(TAG)|(TGA)', first)]
    seq=coding_evaluator()
    cds_info=seq.get_orfs(record,len(record.seq)-1,startposition,stopposition)
    if (seq.is_lncrna(startposition,stopposition)==True):
        lncRNAs=lncRNAs.append(pd.DataFrame({'transcript': record.id}, index=[0]), ignore_index=True)
    elif (cds_info[1].iloc[0,0]=='coding_with_no_stop'):
        no_stop_best_orfs=pd.concat([no_stop_best_orfs,cds_info[0]])
    elif (cds_info[1].iloc[0,0]=='coding') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*4)].shape[0]>0:
        best_orfs_cds=pd.concat([best_orfs_cds,cds_info[0]])
    elif (cds_info[1].iloc[0,0]=='coding') and cds_info[0].loc[cds_info[0]['length']>=(minimum_number_of_aminoAcids*4)].shape[0]==0:
        best_orfs_short_cds=pd.concat([best_orfs_short_cds,cds_info[0]])
    return(lncRNAs,no_stop_best_orfs,best_orfs_cds,best_orfs_short_cds)
"""
############################        
"""               
df=evaluate_length(startposition,stopposition)
df.loc[df.groupby('start')['stop'].idxmin(), :].reset_index()
"""