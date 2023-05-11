import urllib
import numpy as np
import pandas as pd



def get_info_from_best_structures(uniprot_id, best_structures,    appender = []):

    rank = 0
    pdblist = []
    for bs_info in best_structures:

        pdb_id = bs_info["pdb_id"]

        if pdb_id not in pdblist:
            pdblist += [pdb_id]
            rank +=1


        resolution = bs_info["resolution"]
        experimental_method = bs_info["experimental_method"]

        start = bs_info["start"]
        end = bs_info["end"]
        unp_start = bs_info["unp_start"]
        unp_end = bs_info["unp_end"]
        coverage = bs_info["coverage"]

        chain_id = str(bs_info['chain_id'])

#         print rank, pdb_id ,  chain_id,resolution, '\t', experimental_method,uniprot_id, start,end, unp_start,unp_end,coverage
        appender.append([uniprot_id, pdb_id, chain_id, rank,coverage, resolution, start,end, unp_start,unp_end, experimental_method])
    return appender



def encode_query(string_of_pdbs):

    query_dict = """{
      entries(entry_ids: ["PDBLIST"]) {
        rcsb_id
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entity_id
            asym_ids
            auth_asym_ids
            reference_sequence_identifiers {
              database_name
              database_accession
            }
          }
          entity_poly {
            pdbx_seq_one_letter_code_can
            rcsb_sample_sequence_length
            type
          }
          rcsb_polymer_entity {
            formula_weight
          }
        }
      }
    }
    """
    query_dict = query_dict.replace('PDBLIST', string_of_pdbs)
    url_encoded = urllib.parse.quote(query_dict, safe = '"(entry_ids)"')
    return url_encoded

def parse_pdb_search(data):
    seq = data['entity_poly']['pdbx_seq_one_letter_code_can']
    polymer_entity_seq_len = data['entity_poly']['rcsb_sample_sequence_length']
    entity_macro_type = data['entity_poly']['type']
    
    formula_weight = data['rcsb_polymer_entity']['formula_weight']
    
    asym_ids = data['rcsb_polymer_entity_container_identifiers']['asym_ids']
    auth_asym_ids = data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']
    entity_id = data['rcsb_polymer_entity_container_identifiers']['entity_id']
    reference = data['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers']
    
    if reference == 'null':
        databaseId = entity_macro_type
        databaseName = np.nan
    else:
        databaseId = reference[0]['database_accession']
        databaseName = reference[0]['database_name']
    
    return [entity_id,asym_ids,auth_asym_ids,databaseName,databaseId,seq,polymer_entity_seq_len,entity_macro_type,formula_weight]


def encode_query_structure_info(string_of_pdbs):

    query_dict = """{
    entries(entry_ids: ["PDBLIST"])
    {
    rcsb_id
    struct {
    title
    }
    rcsb_entry_info {
    structure_determination_methodology
    resolution_combined
    molecular_weight
    deposited_polymer_monomer_count
    deposited_atom_count
    }
    exptl {
    method
    }
    rcsb_accession_info {
    initial_release_date
    }
    struct_keywords {
    pdbx_keywords
    }
    audit_author {
    name
    }
    }
    }
    """
    query_dict = query_dict.replace('PDBLIST', string_of_pdbs)
    url_encoded = urllib.parse.quote(query_dict, safe = '"(entry_ids)"')
    return url_encoded