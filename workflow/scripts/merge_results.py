import pandas as pd
import numpy as np

#define snakemake variables
report_file = snakemake.input[0]
complete_file = snakemake.input[1]
merged_file = snakemake.output[0]
unmapped_vars_file = snakemake.output[1]

#read the files
report_df = pd.read_csv(report_file, sep='\t')
complete_df = pd.read_csv(complete_file)

#trim the report file to remove useless columns
report_df_trimmed = report_df[['#feat_name', 'source_id', 'mapped_id', 'source_start', 'mapped_start']]
report_df_trimmed.columns = ['feat_name', 'chr_v10_ncbi', 'chr_v11_ncbi', 'var_pos_10', 'var_pos_11']
report_df_trimmed = report_df_trimmed.drop_duplicates(subset=['feat_name']) #removes multiple remappings of the same variant; recheck this if problems with the remapping arise
report_df_trimmed = report_df_trimmed.reset_index().drop(columns = 'index')

#its important to use the merging_column since if you use the variant position without the chromosome you find duplicate positions
report_df_trimmed['merging_column'] = report_df_trimmed['chr_v10_ncbi'].astype(str) + report_df_trimmed['var_pos_10'].astype(str)
complete_df['merging_column'] = complete_df['chr'].astype(str) + complete_df['var_pos_absolute'].astype(str)

#merge the dataframes and rename the columns
final_df = complete_df.merge(report_df_trimmed, left_on='merging_column', right_on='merging_column')
final_df = final_df.drop(columns=['var_pos_absolute', 'merging_column'])
final_df.columns = ['gene_name_10','ensembl_id_10', 'chr_10_ensembl', 'strand', 'gene_start_10', 'gene_end_10', 'relative_var_pos_10', 'feat_name', 'chr_v10_ncbi', 'chr_v11_ncbi', 'var_pos_10', 'var_pos_11']
final_df = final_df[['feat_name', 'gene_name_10','ensembl_id_10', 'chr_10_ensembl', 'chr_v10_ncbi', 'chr_v11_ncbi', 'strand', 'gene_start_10', 'gene_end_10', 'relative_var_pos_10', 'var_pos_10', 'var_pos_11']]
final_df['chr_v10_ncbi'].equals(final_df['chr_10_ensembl'])
final_df = final_df.drop(columns=['chr_10_ensembl'])

#select the variant positions which ncbi was not able to remap
unmapped_df = final_df[final_df['var_pos_11'].isna()]
unmapped_df['var_pos_11'] = unmapped_df['var_pos_11'].fillna('NOMAP')

#fill the NaN values with NOMAP to better mark the unmapped vars
final_df['var_pos_11'] = final_df['var_pos_11'].fillna('NOMAP')
final_df['chr_v11_ncbi'] = final_df['chr_v11_ncbi'].fillna('NOMAP')
final_df['var_pos_11'] = final_df['var_pos_11'].astype(str).str.split('.', expand = True)[0]


#save the files
unmapped_df.to_csv(unmapped_vars_file, index = False, na_rep = 'NOMAP')
final_df.to_csv(merged_file, index = False)
