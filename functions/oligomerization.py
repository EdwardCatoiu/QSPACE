import pandas as pd
import numpy as np
from tqdm import tqdm_notebook as tqdm
import ast
import copy
import json
import os
import os.path as op
import logging
log = logging.getLogger(__name__)

import sys
with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)

swiss_metric_thresholds = {"QMN4": -4.0, "QSPRD": 0.5, "GMQE" : 0.5}

def get_dfhetero(parent,real_structure_quality_pdb,oa_dict_pdb,o_dict_pdb,o_dict_swiss, blast_cutoff = 80):
    columns = ['gene', 'gene_stoichiometry','manual_stoich','manual', 'Hetero_kmer_PDB','Homo_kmer_PDB','Homo_kmer_SWISS']
    dfoa = pd.DataFrame(columns = columns)
    for cplx, slist in oa_dict_pdb.items():
        gene = list(parent[cplx].keys())[0]
        dfoa.loc[cplx , 'gene'] = gene
        dfoa.loc[cplx , 'gene_stoichiometry'] = str({gene : 1})

        quality = {}
        for s in slist:
            quality.update({s: real_structure_quality_pdb[s]})
        dfoa.loc[cplx , 'Hetero_kmer_PDB'] = str(quality)

        if cplx in o_dict_pdb:
            dfoa.loc[cplx , 'Homo_kmer_PDB'] = str(o_dict_pdb[cplx])
        if cplx in o_dict_swiss:
            dfoa.loc[cplx , 'Homo_kmer_SWISS'] = str(o_dict_swiss[cplx])


    #qcqa        
    for cplx in tqdm(dfoa.index.tolist()):

        gene = dfoa.loc[cplx, 'gene']
        oa_structures = copy.deepcopy(ast.literal_eval(dfoa.loc[cplx,'Hetero_kmer_PDB']))
        del_list = []
        for structure, gs in oa_structures.items():
            gene_blast = np.mean(np.array(list(gs[gene].values()))[:,1])
            if gene_blast < blast_cutoff:
                del_list += [structure]
        #     print structure ,gs[gene]
        #     break
    #     print len(del_list), len(oa_structures), cplx

        for structure in del_list:
            del oa_structures[structure]

        if len(oa_structures) >0:
            dfoa.loc[cplx,'Hetero_kmer_PDB'] = str(oa_structures)
        else:
            dfoa.loc[cplx,'Hetero_kmer_PDB'] = np.nan

    keep = []
    for index, row in tqdm(dfoa.iterrows()):
        if len(row.dropna().drop('gene')) > 0:
            keep += [index]
    dfoa = dfoa[dfoa.index.isin(keep)]
#     dfoa=dfoa[dfoa.oa_dict.isna() == False]
    return dfoa



def auto_qcqa_swiss(parent,
                    o_dict_swiss,
                    oa_dict_pdb,
                    o_dict_pdb,
                    df_swiss,
                    real_structure_stoich_swiss, 
                   **kwargs):
    
    o_dict_keys_swiss = set(o_dict_swiss)
    oa_dict_keys_pdb = set(oa_dict_pdb)
    o_dict_keys_pdb = set(o_dict_pdb)
    swiss_only = o_dict_keys_swiss-o_dict_keys_pdb-oa_dict_keys_pdb

    swiss_only_dict = {}
    for cplx in swiss_only:
        swiss_only_dict.update({cplx : o_dict_swiss[cplx]})



    swiss_results = {}
    failed_swiss = {}
    for k , swiss_model_list in swiss_only_dict.items():
        for swiss_id in swiss_model_list:
    #         df_swiss.loc[swiss_id]
            
        
            gs_swiss = df_swiss.loc[swiss_id,'gene_stoichiometry']
            if type(gs_swiss) == str:
                gs_swiss = ast.literal_eval(gs_swiss)
            is_monomer = np.sum(list(gs_swiss.values())) == 1
            is_good_model, metrics = check_swiss_model(df_swiss,
                                                       swiss_id = swiss_id, 
                                                       total_quality_cutoff = 70,
                                                       is_monomer= is_monomer)
            if is_good_model:
                kmer = real_structure_stoich_swiss[swiss_id]
                kmer = np.sum(list(kmer.values()))

                if k not in swiss_results:
                    swiss_results.update({k : {kmer : [swiss_id]}})
                elif kmer not in swiss_results[k]:
                    swiss_results[k].update({kmer : [swiss_id]})
                elif swiss_id not in swiss_results[k][kmer]:
                    swiss_results[k][kmer] += [swiss_id]
            else:
                kmer = real_structure_stoich_swiss[swiss_id]
                kmer = np.sum(list(kmer.values()))

                if k not in failed_swiss:
                    failed_swiss.update({k : {kmer : [swiss_id]}})
                elif kmer not in failed_swiss[k]:
                    failed_swiss[k].update({kmer : [swiss_id]})
                elif swiss_id not in failed_swiss[k][kmer]:
                    failed_swiss[k][kmer] += [swiss_id]
                



    df_swiss_results = pd.DataFrame(columns = ["total_quality","QSPRD","QMN4","GMQE"])
    for cplx, results in swiss_results.items():
        if len(results) > 1 :
            print (cplx)
            print (results)
            break
        for kmer, swiss_id_list in results.items():
            if len(swiss_id_list) > 1:
                print (cplx, kmer, swiss_id_list)


        for kmer, swiss_id_list in results.items():
            swiss_id = swiss_id_list[0]
            df_swiss_results.loc[cplx, 'total_quality'] = df_swiss.loc[swiss_id, 'total_quality']
            df_swiss_results.loc[cplx, kmer] = swiss_id
            is_good_model, metrics = check_swiss_model(df_swiss,
                                                       swiss_id = swiss_id,
                                                       total_quality_cutoff = 70,
                                                      )
            for metric, value in metrics.items():
                df_swiss_results.loc[cplx, metric] = value
    print ("Confirmed {} of {} Genes as Homo-Oligomers".format(len(df_swiss_results), len(swiss_only)))
    
    ordered_colums = []
    for c in df_swiss_results.columns:
        try:
            int(c)
            ordered_colums +=[ c]
        except ValueError:
            continue
            
    ordered_colums.sort()
    
    updated_gene_stoich = {}
    for index, row in df_swiss_results.iterrows():
        row = row.drop(["total_quality","QSPRD","QMN4","GMQE"]).dropna()
        if len(row) != 1:
            log.warn(index,row)
            raise ValueError ('Multiple Good Swiss Models, cannot find the right stoichiometry')
        gene_stoich = row.keys()[0]
        updated_gene_stoich.update({index : {list(parent[index].keys())[0] : gene_stoich}})
    df_swiss_results['gene_stoich'] = pd.Series(updated_gene_stoich)
    df_swiss_results = df_swiss_results[['gene_stoich'] + ordered_colums + ["total_quality","QSPRD",'QMN4','GMQE']]

    return df_swiss_results,swiss_only_dict,failed_swiss
    
    
def auto_qcqa_swiss_and_pdb(parent,
                            o_dict_swiss,       
                            oa_dict_pdb,
                            o_dict_pdb,
                            df_swiss,
                            df_pdb,
                            dfbest_trim,
                            real_structure_stoich_swiss,
                           **kwargs):
    
    keys  = set(o_dict_swiss.keys()).intersection(set(o_dict_pdb.keys()))-set(oa_dict_pdb.keys())

    o_dict  = {}
    for cplx in keys:
        o_dict.update({cplx : o_dict_swiss[cplx] +  o_dict_pdb[cplx]})



    swiss_results = {}
    print ('PDB homo-k-mer results')

    for k , swiss_model_list in tqdm(o_dict.items()):
        for swiss_id in swiss_model_list:
    # #         df_swiss.loc[swiss_id]
            try:
                is_good_model, metrics = check_swiss_model(df_swiss,
                                                           swiss_id = swiss_id,
                                                           total_quality_cutoff = 70
                                                          )
            except KeyError:
                try:
                    kmer = df_pdb.loc[swiss_id,'gene_stoichiometry']
                except KeyError:
                    log.warn("The current Manual Curation contains structure {} that is NOT in your current version of the SWISS-Model Repository. Please re-do manual curation for enzyme {}".format(swiss_id, k))
                    continue
                if type(kmer) == str:
                    kmer = ast.literal_eval(kmer)
                kmer = np.sum(list(kmer.values()))

                if k not in swiss_results:
                    swiss_results.update({k : {kmer : [swiss_id]}})
                elif kmer not in swiss_results[k]:
                    swiss_results[k].update({kmer : [swiss_id]})
                elif swiss_id not in swiss_results[k][kmer]:
                    swiss_results[k][kmer] += [swiss_id]

                continue

            if is_good_model:
                kmer = real_structure_stoich_swiss[swiss_id]
                kmer = np.sum(list(kmer.values()))

                if k not in swiss_results:
                    swiss_results.update({k : {kmer : [swiss_id]}})
                elif kmer not in swiss_results[k]:
                    swiss_results[k].update({kmer : [swiss_id]})
                elif swiss_id not in swiss_results[k][kmer]:
                    swiss_results[k][kmer] += [swiss_id]



    total_quality = {}           
    df_swiss_results = pd.DataFrame(columns = ['gene_stoich','manual_stoich','total_quality','calculated','manual'])
    print ('SWISS homo-k-mer results')
    for cplx, results in tqdm(swiss_results.items()):
  

        for kmer, swiss_id_list in results.items():
            if cplx not in total_quality:
                total_quality.update({cplx : dfbest_trim[dfbest_trim.index.isin(swiss_id_list)].total_quality.to_dict()})
            else:
                total_quality[cplx].update(dfbest_trim[dfbest_trim.index.isin(swiss_id_list)].total_quality.to_dict())

            df_swiss_results.loc[cplx, kmer] = str(swiss_id_list)


    df_swiss_results['total_quality'] = pd.Series(total_quality)
    print ("Found SWISS & PDBD Homo-Oligomerization evidence for {} Genes".format(len(df_swiss_results)))

    for index, row in df_swiss_results.iterrows():
        row = row.drop('total_quality').dropna()
        if len(row) ==1:
            df_swiss_results.loc[index, 'calculated'] = row.keys()[0]
            df_swiss_results.loc[index, 'manual'] = row.keys()[0]

    print ("\t{} of {} Oligomers need Manual Curation".format(len(df_swiss_results[df_swiss_results.manual.isna() == True]), len(keys)))
    print ("\t{} of {} Oligomers AUTO Confirmed".format(len(df_swiss_results[df_swiss_results.manual.isna() == False]), len(keys)))

    # df_swiss_results.to_csv('df_339_swiss_and_pdb.csv')
    return df_swiss_results[df_swiss_results.manual.isna()],df_swiss_results[df_swiss_results.manual.isna() == False]


def check_swiss_model(dfpseudo_best_pdbswiss,
                      swiss_model_metrics=False,
                      swiss_id = None, 
                      total_quality_cutoff = 70,
                      is_monomer = False):
    
    if not swiss_model_metrics:
        src = op.join(qspaceDirs['DataOutput_dir'], '001B-swiss_model_metrics.json') 
        with open(src,'rb') as f:
            swiss_model_metrics = json.load(f)
    
    is_good_model = True
    
    manual_swiss_ok = ['P0CF26_41_143_1czb_&_NOBFS_0',
                      'P0CF80_3_275_5hoo_&_NOBFS_0',
                      'P0CF82_3_275_5hoo_&_NOBFS_0',
                      'P0CF60_3_194_5hoo_&_NOBFS_0',
                      'P76508_49_172_4mtm_&_NOBFS_0',
                      'P0CF81_3_275_5hoo_&_NOBFS_0',
                      'P33227_19_102_2xgf_&_NOBFS_0',
                      'P0CF89_57_380_4fcy_&_NOBFS_0',
                      'P75980_32_137_3kdr_&_NOBFS_0',
                      'P75684_3_130_5fhk_&_NOBFS_0',
                      'P0CF79_3_275_5hoo_&_NOBFS_0',
        
                      ]
    
    if swiss_id in manual_swiss_ok:
        return is_good_model,swiss_model_metrics[ swiss_id.split('_&_')[0]]
        
        
    if dfpseudo_best_pdbswiss.loc[swiss_id].total_quality < total_quality_cutoff:
        is_good_model = False
        return is_good_model,swiss_model_metrics[ swiss_id.split('_&_')[0]]
    
    if swiss_id is None:
        raise KeyError ("no swiss id given")
    swiss_id = swiss_id.split('_&_')[0]
    for metric in [u'GMQE', u'QMN4', u'QSPRD']:
        try:
            metric_value = swiss_model_metrics[swiss_id][metric]
        except KeyError:
            continue
        
        if metric_value < swiss_metric_thresholds[metric]:
            if metric == 'QSPRD' and is_monomer:
                continue
            
            is_good_model = False
            break
            
    return is_good_model, swiss_model_metrics[swiss_id]