import pandas as pd
import requests
import sys
import json
import os
#os.chdir('/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/scripts')
from ensembl_rest_client import EnsemblRestClient

remapped_ids= snakemake.input[0]
merged_file = snakemake.input[1]
#remapped_ids = '/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/data/starting_files/remapped_ids.csv'
#merged_file = '/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/data/complete_coords/merged_coords.csv'

remapped_df = pd.read_csv(remapped_ids)#, index_col=0)
merged_df = pd.read_csv(merged_file)
ensembl_ids_v11 = remapped_df['updated_id']
ensembl_ids_v11 = ensembl_ids_v11.tolist()

#for id in ensembl_ids_v11:
#    if id == 'NOT_FOUND':
#        ensembl_ids_v11.remove(id)
#    elif id == 'PROBLEMATIC':
#        ensembl_ids_v11.remove(id)
#remapped_df = remapped_df[remapped_df['updated_id'] != 'NOT_FOUND']
#remapped_df = remapped_df[remapped_df['updated_id'] != 'PROBLEMATIC']

genes_info_v11_df = EnsemblRestClient().get_genes_info(ensembl_ids_v11)
merged_df = merged_df.merge(remapped_df, left_on='ensembl_id_10', right_on='id_10_2')
final_merged_df = merged_df.merge(genes_info_v11_df, left_on='updated_id', right_on='id')
final_merged_df = final_merged_df.drop(columns=['id_10_2', 'id_11_1', 'id', 'assembly_name', 'biotype', 'description'])
final_merged_df.columns=['feat_name', 'gene_name_10', 'ensembl_id_10', 'chr_10_ncbi', 'chr_11_ncbi', 'strand_10', 'gene_start_10', 'gene_end_10', 'rel_var_pos_10', 'var_pos_10', 'var_pos_11', 'ensembl_id_11', 'gene_name_11', 'strand_11', 'chr_11_ensembl', 'gene_start_11', 'gene_end_11']
final_merged_df = final_merged_df[final_merged_df['var_pos_11'] != 'NOMAP']

forward_genes_11 = final_merged_df[final_merged_df['strand_11']== 1]
reverse_genes_11 = final_merged_df[final_merged_df['strand_11']== -1]
forward_genes_11['rel_var_pos_11'] = forward_genes_11['var_pos_11'].astype(int) - forward_genes_11['gene_start_11'].astype(int) + 1
reverse_genes_11['rel_var_pos_11'] = reverse_genes_11['gene_end_11'].astype(int) - reverse_genes_11['var_pos_11'].astype(int) + 1
df = pd.concat([forward_genes_11, reverse_genes_11]).sort_index()

#some variants were remapped outside the gene coordinates, leading to negative relative variant positions
df = df[df['rel_var_pos_11'] > 0] #remove rows with negative variant relative positions
df = df[['feat_name', 'gene_name_10', 'gene_name_11', 'ensembl_id_10', 'ensembl_id_11', 'chr_10_ncbi', 'chr_11_ncbi', 'chr_11_ensembl', 'strand_10', 'strand_11', 'gene_start_10', 'gene_end_10', 'rel_var_pos_10', 'var_pos_10', 'gene_start_11', 'gene_end_11', 'rel_var_pos_11', 'var_pos_11']]

genes_info_v11_df.to_csv(snakemake.output[0])
df.to_csv(snakemake.output[1], index=False)
