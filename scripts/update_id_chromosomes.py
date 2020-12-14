import pandas as pd
import requests, sys

def update(remapped_file, updated_df):
    #remapped_df = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/remapped_vars/ENSSSCG00000016940.remapped', sep = ',')
    updated_ids = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/remapped_ids.csv', sep = ',')
    temp_df = updated_ids.loc[updated_ids['id_10_2'] == remapped_df['ensembl_id_10'][0]]
    updated_id = temp_df['updated_id'].iloc[0]
    remapped_df['ensembl_id_11'] = updated_id
    if updated_id != 'NOT_FOUND':
        server = "https://rest.ensembl.org"
        ext = "/lookup/id/"+ updated_id + "?"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        decoded = r.json()
        chr_11 = decoded['seq_region_name']
        remapped_df['chromosome_11'] = chr_11
        remapped_df
    else:
        remapped_df['chromosome_11'] = 'NOT_FOUND'
    return(updated_df)
remapped_df
