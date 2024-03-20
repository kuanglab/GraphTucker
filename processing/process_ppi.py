import os
import csv
import gzip
# import torch
# import anndata

import numpy as np
import pandas as pd
# import scanpy as sc
# import tensorly as tl
# import sklearn.neighbors

# from scipy.io import mmread
from scipy.sparse import coo_matrix
from scipy.io import savemat

PPI_data_path = '../data/BRCA1/BIOGRID-ORGANISM-Homo_sapiens-4.4.231.tab3.txt'  # path to downloaded PPI file
ppi_name = 'HSA_ppi'  # this is the name of the PPI it will be saved as, much shorter than the original name
save_path = '../data/BRCA1/' # change this to where you want PPI.mat file saved
ppi_save_path = save_path + ppi_name + '.mat'
print('PPI will be saved to: ', ppi_save_path)

# load list of gene names present in Visium dataset
features_path = '../data/BRCA1/filtered_feature_bc_matrix/features.tsv.gz'

gene_names = np.char.lower(np.array([row[1] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter="\t")]))

# Load PPI graph
biogrid_PPI = pd.read_csv(PPI_data_path, delimiter='\t')

# Construct PPI graph
gene_names_A = biogrid_PPI['Official Symbol Interactor A'].str.lower().to_numpy().astype('<U16')
gene_names_B = biogrid_PPI['Official Symbol Interactor B'].str.lower().to_numpy().astype('<U16')
matched_indices = np.where(np.in1d(gene_names_A, gene_names) & np.in1d(gene_names_B, gene_names))[0]
gene_names_A = gene_names_A[matched_indices]; gene_names_B = gene_names_B[matched_indices]
gene_names_PPI = np.unique(np.union1d(gene_names_A, gene_names_B))
row = np.array([np.where(gene_names_PPI == gene_name)[0][0] for gene_name in gene_names_A])
col = np.array([np.where(gene_names_PPI == gene_name)[0][0] for gene_name in gene_names_B])
data = np.ones(len(row) + len(col)); n_g = len(gene_names_PPI)
A_g = coo_matrix((data, (np.concatenate([row, col]),
                         np.concatenate([col, row]))),
                  shape=(n_g, n_g)).toarray()
A_g[np.where(A_g > 0)] = 1  # convert weighted adjacent matrix to adjacent matrix
np.fill_diagonal(A_g, 0)  # remove self connections

# the variable names HSA_UNIQ_BIOGRID_GENE and HSA_BIOGRID can be changed,
# but make sure they are matched with loaded variable names in data_prep.m script (or update those variables to whatever)
# you name them here
vardict = {'HSA_UNIQ_BIOGRID_GENE': list(gene_names_PPI), 'HSA_BIOGRID': np.array(A_g)}
savemat(ppi_save_path, vardict)