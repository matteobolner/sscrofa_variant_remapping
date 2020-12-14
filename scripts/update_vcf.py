import pandas as pd
input_vcf = snakemake.input[0]
remapped_df = snakemake.input[1]

#remapped_df = "/home/pelmo/sscrofa_variant_remapping/data/final_files/remapped_and_verified_vars.csv"
#input_vcf = "/home/pelmo/sscrofa_variant_remapping/data/vcf_files/whole_POP_all.vcf"
df = pd.read_csv(input_vcf, sep = "\t", comment = '#', header = None)
remapped_df = pd.read_csv(remapped_df)
remapped_df = remapped_df[['gene_name_10', 'ensembl_id_10', 'ensembl_id_11', 'chr_11_ensembl', 'strand_11', 'gene_start_11', 'gene_end_11', 'rel_var_pos_10', 'rel_var_pos_11', 'var_pos_11']].reset_index()
remapped_df['id_pos'] = remapped_df['ensembl_id_10'] + "_" + remapped_df['rel_var_pos_10'].astype(str)
df['id_pos'] = df[0].str.split('|').str[1] + "_" + df[1].astype(str)
merged = df.merge(remapped_df, left_on='id_pos', right_on='id_pos')
merged['header'] = merged['gene_name_10'] +"|"+ merged['ensembl_id_11'] +"|"+ merged['chr_11_ensembl'].astype(str) +"|"+ merged['gene_start_11'].astype(str) +"|"+ merged['gene_end_11'].astype(str)
merged = merged.drop(columns = [0,'gene_name_10', 'ensembl_id_10', 'ensembl_id_11', 'chr_11_ensembl', 'strand_11', 'gene_start_11', 'gene_end_11', 'rel_var_pos_10', 'var_pos_11', 'id_pos'])
cols_to_move = ['header', 'rel_var_pos_11']
merged = merged[cols_to_move + [col for col in merged.columns if col not in cols_to_move]]
merged = merged.drop(columns = 1)

merged.to_csv(snakemake.output[0], index =False, header = False,  sep = "\t")
