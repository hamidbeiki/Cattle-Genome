#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:56:33 2021

@author: beiki
"""

from matplotlib_venn import venn3,venn2, venn3_circles,venn2_circles
import pandas as pd
import matplotlib.pyplot as plt

protein_coding_genes=pd.read_csv('/Users/beiki/Desktop/test/protein-coding_genes',sep='\t',names=['gene_id'])
protein_coding_genes_with_protein_homology=pd.read_csv('/Users/beiki/Desktop/test/protein-coding_genes_with_protein_homology',sep='\t',names=['gene_id'])
protein_coding_genes_with_ncRNA_homology=pd.read_csv('/Users/beiki/Desktop/test/protein-coding_genes_with_ncRNA_homology',sep='\t',names=['gene_id'])
non_coding_genes=pd.read_csv('/Users/beiki/Desktop/test/non-coding_genes',sep='\t',names=['gene_id'])
non_coding_genes_with_peptide_homology=pd.read_csv('/Users/beiki/Desktop/test/non-coding_genes_with_peptide_homology',sep='\t',names=['gene_id'])
non_coding_genes_with_ncRNA_homology=pd.read_csv('/Users/beiki/Desktop/test/non-coding_genes_with_ncRNA_homology',sep='\t',names=['gene_id'])

set1=set(non_coding_genes['gene_id'])
set2=set(non_coding_genes_with_peptide_homology['gene_id'])
set3=set(non_coding_genes_with_ncRNA_homology['gene_id'])

v=venn3([set1, set2, set3])#('Set1', 'Set2', 'Set3')
c=venn3_circles([set1, set2, set3], linestyle='dashed', color="grey")
plt.show()


'''
v=venn2([ set2, set3], ('Set2', 'Set3'))
c=venn2_circles([set2, set3], linestyle='dashed', color="grey")
plt.show()
'''