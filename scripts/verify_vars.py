import pandas as pd
import json
import sys
import numpy as np
import os

#os.chdir("/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/scripts")
from ensembl_rest_client import EnsemblRestClient
from difflib import SequenceMatcher
#from rest_functions import get_gene_seqs_v10

input_df = snakemake.input[0]
#input_df = '/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/data/complete_coords/merged_updated_coords_v11.csv'
genes_info_v10 = snakemake.input[1]
#genes_info_v10 = '/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/data/genes_info/genes_info_v10.csv'
genes_info_v11 = snakemake.input[2]
#genes_info_v11 = '/home/pelmo/data_and_pipelines/sscrofa_variant_remapping/data/genes_info/genes_info_v11.csv'

input_df = pd.read_csv(input_df)
genes_df_v10 = pd.read_csv(genes_info_v10)
genes_df_v11 = pd.read_csv(genes_info_v11)

gene_list_v10 = genes_df_v10['id'].unique().tolist()
gene_list_v11 = genes_df_v11['id'].unique().tolist()
len(gene_list_v10)
seq_df_v10 = EnsemblRestClient().get_gene_seqs(gene_list_v10, server_version='89')
seq_df_v11 = EnsemblRestClient().get_gene_seqs(gene_list_v11)

temp_df = input_df.merge(seq_df_v10, left_on='ensembl_id_10', right_on='id')
temp_df.rename(columns={'seq':'seq_v10'}, inplace=True)
temp_df = temp_df.drop(columns=['id'])
df_with_seqs = temp_df.merge(seq_df_v11, left_on='ensembl_id_11', right_on='id')
df_with_seqs.rename(columns={'seq':'seq_v11'}, inplace=True)
df_with_seqs = df_with_seqs.drop(columns=['id'])

df_with_seqs['slice_coords_v10_start'] = (df_with_seqs['rel_var_pos_10'] -1 -4)
df_with_seqs['slice_coords_v10_end'] = (df_with_seqs['rel_var_pos_10'] -1 +5)
df_with_seqs['slice_coords_v11_start'] = (df_with_seqs['rel_var_pos_11'] -1 -4)
df_with_seqs['slice_coords_v11_end'] = (df_with_seqs['rel_var_pos_11'] -1 +5)

df_with_seqs['var_seq_v10'] = df_with_seqs.apply(lambda x: x['seq_v10'][x['slice_coords_v10_start']:x['slice_coords_v10_end']], axis=1)
df_with_seqs['var_base_v10'] = df_with_seqs.apply(lambda x: x['seq_v10'][x['rel_var_pos_10']-1], axis=1)
df_with_seqs['var_seq_v11'] = df_with_seqs.apply(lambda x: x['seq_v11'][x['slice_coords_v11_start']:x['slice_coords_v11_end']], axis=1)
df_with_seqs['len_seq_v11'] = df_with_seqs['seq_v11'].str.len()
df_with_seqs = df_with_seqs[df_with_seqs['rel_var_pos_11'] <= df_with_seqs['len_seq_v11'] ]
df_with_seqs = df_with_seqs[df_with_seqs['var_seq_v10'].str.len() != 0].reset_index()
df_with_seqs = df_with_seqs[df_with_seqs['var_seq_v11'].str.len() != 0].reset_index()
df_with_seqs['var_base_v11'] = df_with_seqs.apply(lambda x: x['seq_v11'][x['rel_var_pos_11']-1], axis=1)


df_with_seqs = df_with_seqs.drop(columns = ['seq_v10', 'seq_v11', 'slice_coords_v10_start', 'slice_coords_v10_end', 'slice_coords_v11_start', 'slice_coords_v11_end', 'level_0', 'index', 'feat_name', 'len_seq_v11'])
df_with_seqs['identical_base'] = df_with_seqs['var_base_v10'] == df_with_seqs['var_base_v11']

def sim_metric(df, col1, col2):
    return SequenceMatcher(None, df[col1], df[col2]).ratio()

df_with_seqs['seq_similarity'] = df_with_seqs.apply(sim_metric, args = ('var_seq_v10', 'var_seq_v11'), axis=1)
unsure_remappings = df_with_seqs[df_with_seqs['seq_similarity'] <= 0.5].reset_index()
df_with_seqs = df_with_seqs[df_with_seqs['seq_similarity'] >= 0.5].reset_index()
df_with_seqs.to_csv(snakemake.output[0], index=False)
unsure_remappings.to_csv(snakemake.output[1], index=False)
