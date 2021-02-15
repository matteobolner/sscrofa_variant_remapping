#!/usr/bin/env python
import sys
import json
import time
import requests
import pandas as pd


#Define the class EnsemblRestClient
class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=15):
        #the two versions of Ensembl REST are those corresponding with Sscrofa11.1 (102) and Sscrofa10.2 (89)
        self.server = {'102': 'http://rest.ensembl.org', '89' : 'https://may2017.rest.ensembl.org'}
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.species = 'sus scrofa'
    #Define a generic function to retrieve data from the server based on the endpoint and server version specified
    def perform_rest_action(self, endpoint, server_version='102', hdrs=None, params=None, data=None):
        server = self.server[server_version]
        endpoint_url = server+endpoint

        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        #if params:
        #    endpoint += '?' + params

        # check if rate limit is needed

        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            if data:
                request = requests.post(endpoint_url, params = params, headers=hdrs, data = data )
            else:
                request = requests.get(endpoint_url, params = params, headers=hdrs)

            if request:
                output_data = request.json()
            self.req_count += 1

        except requests.exceptions.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return output_data


    #ENDPOINT-SPECIFIC FUNCTIONS:


    #Retrieve Ensembl stable ID from the input symbol (gene name or stable ID)
    def get_cross_ref_ids(self, symbols, server_version='102'):
        dict = {}
        for symbol in symbols:
            genes = self.perform_rest_action(
                endpoint='/xrefs/symbol/{0}/{1}'.format(self.species, symbol),
                params={'object_type': 'gene'},
                server_version = server_version
                )
            if genes:
                stable_id = genes[0]['id']
                dict[symbol] = stable_id
        return dict

    #Retrieve genes information (name,version,coordinates,strand etc.)
    def get_genes_info(self, ids, server_version='102'):
        #convert IDs to list and convert to json format
        #if IDs already in list just convert to json format
        try:
            ids = json.dumps(ids.tolist())
        except:
            ids = json.dumps(ids)
        genes = self.perform_rest_action(
            endpoint = '/lookup/id',
            data = '{ "ids" : ' + ids + ' }',
            server_version = server_version
        )
        genes_info_df = pd.DataFrame.from_dict(genes).transpose().reset_index()
        genes_info_df = genes_info_df.drop(columns = ['source', 'object_type', 'logic_name', 'version', 'species', 'db_type'])
        genes_info_df = genes_info_df[['display_name', 'id', 'assembly_name', 'biotype', 'strand', 'seq_region_name', 'start', 'end', 'description' ]]
        return(genes_info_df)

    #Retrieve the sequence of the input Ensembl Stable IDs
    def get_gene_seqs(self, ids, server_version = '102'):
        ensembl_ids = ids
        chunks = [ensembl_ids[x:x+50] for x in range(0, len(ensembl_ids), 50)]
        seqs_df = pd.DataFrame([], columns = ['id', 'seq'])
        for chunk in chunks:
            id_list = json.dumps(chunk)
            seqs = self.perform_rest_action(
            endpoint = "/sequence/id",
            data = '{ "ids" : ' + id_list + ' }',
            server_version = server_version)
            temp_df = pd.DataFrame.from_dict(seqs)
            temp_df = temp_df[['id','seq']]
            seqs_df = pd.concat([seqs_df, temp_df ], ignore_index=True)
        return(seqs_df)


    #Retrieve all variants overlapping the input gene
    #Old function, useful but not necessary for the project
    def get_variants(self, symbol, server_version ='102'):
        genes = self.perform_rest_action(
            endpoint='/xrefs/symbol/{0}/{1}'.format(self.species, symbol),
            params={'object_type': 'gene'},
            server_version = server_version
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'},
                server_version = server_version
            )
            return variants
        return None


#ids =  ['ENSSSCG00000000038', 'ENSSSCG00000002385']

'''
#CALL THE CLASS FUNCTIONS
def get_ids(species, symbols, server_version='102'):
    client = EnsemblRestClient()
    stable_ids = client.get_cross_ref_ids(species, symbols, server_version=server_version)
    return(stable_ids)

def get_variants(species, symbol, server_version='102'):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol, server_version='102')
    if variants:
        for v in variants:
            print('{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))
'''
