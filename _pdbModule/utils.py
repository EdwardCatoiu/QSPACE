import os
import os.path as op
import json
from tqdm import tqdm_notebook as tqdm
import ast
import pandas as pd
import numpy as np
import urllib
# import urllib2
import requests
import copy

import sys
# sys.path.append('../')
# import qspace_directories as qspaceDirs

#import ssbio Functions
# sys.path.append(qspaceDirs.ssbioDir)
with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import logging

import ssbio
from ssbio.protein.structure.structprop import StructProp
from ssbio.databases.pdb import best_structures as ssbio_bs
import ssbio.protein.sequence.utils.alignment as ssbioaln

#import qspace functions
# sys.path.append(qspaceDirs.FunctionsDir)
import PDB_search
import API_BLAST_pdb as blast_api
import download_structures
import handle_bioassemblies
import handle_pdbs
import pseudo_structures as pseudo

#

from bioservices import UniProt
bs_unip = UniProt()

from slugify import Slugify
custom_slugify = Slugify(safe_chars='-_.')

sys.path.append('/usr/local/lib/python3.7/dist-packages/')
import pdbecif
from pdbecif.mmcif_io import CifFileReader

import Bio
from Bio import SeqIO


def textService_getPdbInfo(list_of_PDB_entries):
#     print '\n------------------------------------------------'
#     print "Running textService PDB API to get structure info..."
    appender = []
    for pdb_id in tqdm(list_of_PDB_entries):

        url_encoded = PDB_search.encode_query(pdb_id)
        url = 'https://data.rcsb.org/graphql?query={}'.format(url_encoded)

        req = requests.get(url)
        if req.status_code != 200:
            print (pdb_id, req.status_code)
        response = req.text
        response = response.replace('null','"null"')
        response_dict = ast.literal_eval(response)['data']['entries']

        for data in (response_dict):
            pdb_id = data['rcsb_id']

            polymer_entities = data['polymer_entities']

            for entity in polymer_entities:
                info = [pdb_id] + PDB_search.parse_pdb_search(entity)
                appender.append(info)

    cols = ['pdb_entry', 'entity_id','asym_ids','auth_asym_ids','databaseName','databaseId',
            'seq','polymer_entity_seq_len','entity_macro_type','formula_weight'] 
    dfpdb_search = pd.DataFrame.from_records(appender, columns = cols)
    return dfpdb_search
   


