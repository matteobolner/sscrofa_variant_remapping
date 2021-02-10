import pandas as pd
import os
import requests
import json

#import custom rest functions
from ensembl_rest_client import EnsemblRestClient

#set the snakemake input and output variables
starting_file = snakemake.input[0]
vars_coords_v10_file = snakemake.output.vars_file
complete_coords_file = snakemake.output.complete_coords
genes_info = snakemake.output.genes_info

starting_df = pd.read_csv(starting_file, sep='|')
starting_df = starting_df.drop_duplicates()
ensembl_ids = starting_df['ensembl_id'].unique().tolist()

genes_info_df = EnsemblRestClient().get_genes_info(ensembl_ids, server_version='89')
genes_info_df = genes_info_df.dropna(subset=['id'])
strand_df = genes_info_df[['id','strand']]
starting_df = starting_df.merge(strand_df, left_on='ensembl_id', right_on='id')
starting_df = starting_df.drop(columns = 'id')
starting_df = starting_df[['gene_name', 'ensembl_id', 'chr', 'strand', 'gene_start', 'gene_end', 'var_pos']]

forward_genes = starting_df[starting_df['strand']== 1]
reverse_genes = starting_df[starting_df['strand']== -1]
reverse_genes['var_pos_absolute'] = reverse_genes['gene_end'] - reverse_genes['var_pos'] + 1
forward_genes['var_pos_absolute'] = forward_genes['gene_start'] + forward_genes['var_pos'] - 1
df = pd.concat([forward_genes, reverse_genes]).sort_index()

#format the variant coordinates in a ncbi-remap readable format
vars_coords = starting_df['chr'].astype(str) + ":" + df['var_pos_absolute'].astype(str)+ "-" + df['var_pos_absolute'].astype(str)

#save the files
vars_coords.to_csv(vars_coords_v10_file, header = False, index = False)
df.to_csv(complete_coords_file, index=False, sep = ',')
genes_info_df.to_csv(genes_info, index=False, sep = ',')
