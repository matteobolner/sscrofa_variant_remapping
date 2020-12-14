import pandas as pd
import requests
import sys
import json


                                        #################################################
                                        ###REST FUNCTIONS FOR ENSEMBL 89 - Sscrofa10.2###
                                        #################################################

server_v10 = "https://may2017.rest.ensembl.org" #ensembl rest interface, release 89

def lookup_single_id_v10(ensembl_id):
    #get a dictionary containing information about an ensembl id
    ext = "/lookup/id/" + ensembl_id + "?"
    r = requests.get(server_v10 + ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()
    return(decoded)

def lookup_multiple_ids_v10(ensembl_ids):
    #get a dictionary containing information about a list of ensembl ids
    id_list = json.dumps(ensembl_ids)
    ext = "/lookup/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.post(server_v10 + ext, headers=headers, data='{ "ids" : ' + id_list + ' }')
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return(decoded)

def get_genes_info_v10(ensembl_ids): #TODO : remove empty lines
    #get a dataframe containing selected informations from a list of ensembl ids
    #requires lookup_multiple_ids function to work
    genes_info_dict = lookup_multiple_ids_v10(ensembl_ids)
    genes_info_df = pd.DataFrame.from_dict(genes_info_dict).transpose().reset_index()
    genes_info_df = genes_info_df.drop(columns = ['source', 'object_type', 'logic_name', 'version', 'species', 'db_type'])
    genes_info_df = genes_info_df[['display_name', 'id', 'assembly_name', 'biotype', 'strand', 'seq_region_name', 'start', 'end', 'description' ]]
    return(genes_info_df)


def get_gene_seqs_v10(ensembl_ids):
    chunks = [ensembl_ids[x:x+50] for x in range(0, len(ensembl_ids), 50)]
    seqs_df = pd.DataFrame([], columns = ['id', 'seq'])

    for chunk in chunks:
        id_list = json.dumps(chunk)
        ext = "/sequence/id"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        r = requests.post(server_v10 + ext, headers=headers, data='{ "ids" : ' + id_list + ' }')
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        temp_df = pd.DataFrame.from_dict(decoded)
        temp_df = temp_df[['id','seq']]
        seqs_df = pd.concat([seqs_df, temp_df ], ignore_index=True)

    return(seqs_df)



                                        #################################################
                                        ###REST FUNCTIONS FOR ENSEMBL 102 - Sscrofa11.1##
                                        #################################################

server_v11 = "https://rest.ensembl.org" #ensembl rest interface, release 102

def lookup_single_id_v11(ensembl_id):
    #get a dictionary containing information about an ensembl id
    ext = "/lookup/id/" + ensembl_id + "?"
    r = requests.get(server_v11 + ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()
    return(decoded)

def lookup_multiple_ids_v11(ensembl_ids):
    #get a dictionary containing information about a list of ensembl ids
    id_list = json.dumps(ensembl_ids)
    ext = "/lookup/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.post(server_v11 + ext, headers=headers, data='{ "ids" : ' + id_list + ' }')
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return(decoded)

def get_genes_info_v11(ensembl_ids): #TODO : remove empty lines
    #get a dataframe containing selected informations from a list of ensembl ids
    #requires lookup_multiple_ids function to work
    genes_info_dict = lookup_multiple_ids_v11(ensembl_ids)
    genes_info_df = pd.DataFrame.from_dict(genes_info_dict).transpose().reset_index()
    genes_info_df = genes_info_df.drop(columns = ['source', 'object_type', 'logic_name', 'version', 'species', 'db_type'])
    genes_info_df = genes_info_df[['display_name', 'id', 'assembly_name', 'biotype', 'strand', 'seq_region_name', 'start', 'end', 'description' ]]
    return(genes_info_df)

def get_gene_seqs_v11(ensembl_ids):
    chunks = [ensembl_ids[x:x+50] for x in range(0, len(ensembl_ids), 50)]
    seqs_df = pd.DataFrame([], columns = ['id', 'seq'])

    for chunk in chunks:
        id_list = json.dumps(chunk)
        ext = "/sequence/id"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        r = requests.post(server_v11 + ext, headers=headers, data='{ "ids" : ' + id_list + ' }')
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        temp_df = pd.DataFrame.from_dict(decoded)
        temp_df = temp_df[['id','seq']]
        seqs_df = pd.concat([seqs_df, temp_df ], ignore_index=True)

    return(seqs_df)
