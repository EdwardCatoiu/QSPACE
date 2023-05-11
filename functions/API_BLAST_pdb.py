import json
import urllib
from tqdm import tqdm_notebook as tqdm
# def encode_seq_request(seq, return_type = 'polymer_entity'):
#     query_dict = {"query": {"type": "terminal",
#                      "service": "sequence",   
#                      "parameters": {  
#                          "evalue_cutoff": 0.0001,  
#                          "identity_cutoff": 0.20,
#                          "target": "pdb_protein_sequence", 
#                          "value": seq 
#                      }  
#                     }, 
#            "request_options": {   
#                "scoring_strategy": "sequence",
#                "pager": {
#                    "start": 0,
#                    "rows": 200
            
#                },
#            },  
#            "return_type": return_type
#           }
#     query_dict = json.dumps(query_dict)
#     url_encoded = urllib.quote(query_dict)
#     url = 'https://search.rcsb.org/rcsbsearch/v2/query?json={}'.format(url_encoded)
# #     print url
#     return url


def encode_seq_request_new(seq, return_type = 'polymer_entity'):
    query_dict= {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.0001,
                "identity_cutoff": 0.20,
                "sequence_type": "protein",
                "value": seq,        
            }
        },
        "return_type": return_type,
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 250
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "sequence"
        }
    }
    query_dict = json.dumps(query_dict)
    url_encoded = urllib.parse.quote(query_dict)
    url = 'https://search.rcsb.org/rcsbsearch/v2/query?json={}'.format(url_encoded)
#     print url
    return url

# def get_results_dict(response,seq):
#     query_result = {}
#     for result in response['result_set']:
#         pdb_id = result['identifier']

#         service_info = result['services'][0]['nodes'][0]
#         original_score = service_info['original_score']
#         norm_score = service_info['norm_score']
#         match = service_info['match_context'][0]
#         seq_ident = float(match['alignment_length']) - match['mismatches']
#         percent_seq_ident = seq_ident / len(seq)
#         gaps = match['gaps_opened']
#         percent_gaps = gaps / len(seq)
#         evalue = match['evalue']
#         query_result.update( {pdb_id : {"original_score":original_score,
#                                        "hit_score":norm_score,
#                                        "hit_num_ident":seq_ident,
#                                        "hit_percent_ident":percent_seq_ident,
#                                        "hit_num_gaps":gaps,
#                                        "hit_percent_gaps":percent_gaps,
#                                        "hit_evalue":evalue,
#                                        }
#                              })
#     return query_result


def get_results_dict_new(response):
    query_result = {}
    for result in response['result_set']:
        pdb_id, enitity_id = result['identifier'].split('_')
        score = result['score']

        if pdb_id not in query_result:
            query_result.update( {pdb_id : {enitity_id : {"score":score}}})
        elif enitity_id not in query_result[pdb_id]:
            query_result[pdb_id].update( {enitity_id : {"score":score}})
    return query_result





def interpret_api_result(api_blast_result):
    pdb_list = {}
    for gene, info in tqdm(api_blast_result.items()):
        for pdb_entity in info.keys():
            [pdb_id , entity_blast]= pdb_entity.split('_')
            pdb_id = pdb_id.lower()
            if pdb_id not in pdb_list:
                pdb_list.update({pdb_id: [gene]})
            elif gene not in pdb_list[pdb_id]:
                pdb_list[pdb_id] += [gene]
    return pdb_list



