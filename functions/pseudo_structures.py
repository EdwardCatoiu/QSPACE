import pandas as pd
from tqdm import tqdm_notebook as tqdm
import copy
import operator
import ast
import numpy as np
from bfs_algo import bfs_algo

import logging
log = logging.getLogger(__name__)


def remake_chains_from_quality(gene_quality_dict):
    quality_input_chains = {}
    
    
    for structure, gene_info in tqdm(gene_quality_dict.items()):

        chain_list = []
        mapping_dict = {}
        for gene, chain_info in gene_info.items():
            if gene in ['GMQE','QMN4','QSPRD']:
                continue
            
            gene_chain_map = chain_info.keys()
            for chain in chain_info.keys():
                chain_list += [chain]
            chain_list = list(set(chain_list))
            mapping_dict.update({gene : gene_chain_map})        
        quality_input_chains.update({structure : {'mapping': mapping_dict}})
        quality_input_chains[structure].update({'chains'  : chain_list})
    return quality_input_chains

def get_quals(pseudo_stoich_chains,input_quality):
    pseudo_quality = {}

    for k , v in tqdm(pseudo_stoich_chains.items()):
        s_id = k.split('_&_')[0]
        squal = {}
        for g, chain_list in v.items():
            gqual = input_quality[s_id][g]
            chain_quailty = {}
            for c in chain_list:            
                chain_qual = gqual[c]
                chain_quailty.update({c : chain_qual})
            squal.update({g : chain_quailty})
        pseudo_quality.update({k : squal})
    return pseudo_quality



########################PSEUDO STRUCTURE CODE BELOW#####################################################



def add_zero_value_to_subcomplexes(child, parent):
    """
    for the algorithm to work, a child class must be assigned a 
    key:value of {'gene' : 0} if that gene is missing
    """
    
    new_child = copy.deepcopy(child)
    for key in parent:
        if key not in child.keys():
            new_child.update({key : 0})
    return new_child


def split_shared_bfs_mappings____pseudo_mathches(structure, Input_stoich_chains):

    """
    Takes the prefiltered PDB:chain mapping and determines which chains are shared among only one gene and which parts of the 
    genes need to be run with BFS matching
    
    Inputs 
        -->Input_stoich_chains <dict> : {u'4kqh': {u'chains': [u'A'], u'mapping': {u'b1991': [u'A']}}}
                                    this dictionary is directly from nathan and blasts gene to chains in pdb

        --> structure <str> : pdb ID or bio assembly ID
    
    
    Returns 
    
        --> shared_stoich <dict> : {gene : # of times mapped to structure} (trivial)
        --> to_bfs_dict <dict> : {chain : [list of genes blasted to chain per structure]} (non-trivial)
    """

    def prefilter_pdb____pseudo_mathches(structure, Input_stoich_chains=Input_stoich_chains):
        """

        Inputs 
        --> see parent function


        Returns  
        -->   chain_gene_list_dict  <dict> : {chain_ID : [list of blasted genes]}

        """
        chain_gene_list_dict = {}
        for chain in Input_stoich_chains[structure]['chains']:
            chain_gene_list_dict.update({chain: []})
            chain_count = 0.
            for gene_map, chain_list in Input_stoich_chains[structure]['mapping'].items():
                if chain in chain_list:
                    chain_gene_list_dict[chain] += [gene_map]
        #print chain_gene_list_dict
        return chain_gene_list_dict

    
    shared_stoich = {}
    to_bfs_dict = {}
    filtered = prefilter_pdb____pseudo_mathches(structure,Input_stoich_chains)
    
    used_gene_chains = {}
    for chain, gene_list in filtered.items():
        if len(gene_list) ==1:
            
            if gene_list[0] in shared_stoich:
                shared_stoich[gene_list[0]] += 1
                
            else:
                shared_stoich.update({gene_list[0] : 1})
                
                
            if gene_list[0] in used_gene_chains:
                used_gene_chains[gene_list[0]] += [chain]
            else:
                used_gene_chains.update({gene_list[0] : [chain]})
              
        else:
            to_bfs_dict.update({chain : gene_list})
    return shared_stoich, to_bfs_dict,used_gene_chains


def create_parents_and_children____pseudo_mathches(structure, to_bfs_dict,Input_stoich_chains,used_gene_chains):
    """
    Splits the to_bfs_chain mapping dict into the appropriate children dictionaries
    using the chains in Input_biostoich_chains[structure]. 
    
    Returns
    parent <dict> : {chain : chain_stoich}
    child_list <list[<dict>]> : {all the children and their chain mappings}
    childnum_gene <dict> : {child # : gene_ID}
    
    
    
    """
    parent = {}
    for chain in to_bfs_dict:
        parent.update({chain : 1})
    
    count = 0
    child_list = []
    childnum_gene = {}
    checked_children = []
    for chain, gene_list in to_bfs_dict.items():
        for g in gene_list:
            if g in checked_children:
                continue
            mapping = Input_stoich_chains[structure]['mapping'][g]
            if g in used_gene_chains:
                use_gene_mapping = used_gene_chains[g]
                mapping = list(set(mapping) - set(use_gene_mapping))
            temp = {}
            for c in mapping:
                temp.update({c : 1.})
            temp.update({'child' : count})
        #     temp.update({'gene' : g})
            temp = add_zero_value_to_subcomplexes(child = temp, parent = parent)
            childnum_gene.update({ count: g})
            count +=1
            child_list.append(temp)
            checked_children.append(g)
    return parent, child_list, childnum_gene

def parse_BFS_answer____pseudo_mathches(structure, s,childnum_gene,child_list, Input_stoich_chains,used_gene_chains):
    answer_count = 0.
    bfs_answer = {}
    bfs_answer_chains = {}
    for answer in s:
        gene_stoich = {}
#         print '\n '
        answer = answer.split('_')
    #     print answer
        answer_chains = {}
        for child_num, child_count in enumerate(answer):
            child_count = int(child_count)
            
            if int(child_count) == 0:
                continue
            
            g = childnum_gene[child_num]
#             print g
            
            for child in child_list:
                if child['child'] == child_num:
                    chains = []
                    for chain, value in child.items():
                        if chain == 'child':
                            continue
                        if value == 0:
                            continue
                        chains += [chain]
                    if g in answer_chains:
                        answer_chains[g] += chains
                    else:
                        answer_chains.update({g : chains})
                        
                        
                        
            used_stoich = 0 
            if g in used_gene_chains:
                used_stoich = len(used_gene_chains[g])
            
            
            
            if g in gene_stoich:
                gene_stoich[g] += child_count * (len(Input_stoich_chains[structure]['mapping'][g])-used_stoich)
            else:
                gene_stoich.update({ g : child_count * (len(Input_stoich_chains[structure]['mapping'][g])-used_stoich)})
        bfs_answer.update({structure + '_&_MATCH_%i' %answer_count : gene_stoich})
        bfs_answer_chains.update({structure + '_&_MATCH_%i' %answer_count : answer_chains})
        answer_count += 1
    return bfs_answer,bfs_answer_chains


def add_shared_stoich_to_BFS(structure, bfs_answer,bfs_answer_chains, shared_stoich, used_gene_chains):
    if len(bfs_answer) == 0:
        bfs_answer = {structure +'_&_NOBFS_0' : shared_stoich}
        bfs_answer_chains = {structure +'_&_NOBFS_0' : used_gene_chains}
        
        return bfs_answer,bfs_answer_chains
    for match, bfs_genes in bfs_answer.items():
        bfs_chains = bfs_answer_chains[match]
        
        for g, gs in shared_stoich.items():
            if g in bfs_genes:
                bfs_genes[g] += gs
            else:
                bfs_genes.update({g:gs})
                
        for g, g_chains in used_gene_chains.items():
            if g in bfs_chains:
                bfs_chains[g] += g_chains
            else:
                bfs_chains.update({g:g_chains})
                 
        
        #bfs_genes.update(shared_stoich)
    return bfs_answer, bfs_answer_chains

def find_pseudo_matches(Input_stoich_chains):

    final_answer  ={}
    final_answer_chains = {}
    for structure in tqdm(Input_stoich_chains):
        #if structure != '5ng5':
            #continue

        #split the pdb into shared and non-shared gene stoichs
        shared_stoich, to_bfs_dict,used_gene_chains = split_shared_bfs_mappings____pseudo_mathches(structure,Input_stoich_chains)

        #create parent, childlist for non-shared region (input to BFS)
        parent, child_list, childnum_gene= create_parents_and_children____pseudo_mathches(structure, to_bfs_dict,Input_stoich_chains,used_gene_chains)

        #run the BFS on non-shared inputs
        a,s, i ,o = bfs_algo(parent=parent, children_list=child_list, answers=[], incomplete_matches= [],over_matches = [])

        #parse the BFS answer, indexing from above step
        bfs_answer,bfs_answer_chains =  parse_BFS_answer____pseudo_mathches(structure, s,childnum_gene,child_list, Input_stoich_chains,used_gene_chains)

        #add the shared stoich to easch match
        bfs_answer,bfs_answer_chains  = add_shared_stoich_to_BFS(structure, bfs_answer,bfs_answer_chains, shared_stoich, used_gene_chains)

        #add this to the final output
        final_answer.update(bfs_answer)
        final_answer_chains.update(bfs_answer_chains)
    return final_answer,final_answer_chains


############################ Finding and Ranking Structures code Below ################
def make_dataframe_pseudo_structures(input_stoich, input_quality,df_gene_length):

    df = pd.DataFrame(columns= ["stype", 'gene_stoichiometry', 'pdb_quality','identical_structures','AA_length'])
    try:
        print ("Making Dataframe...........\n-------------------")
        errors_from_gempro = {}
        gene_stoich={}
        gene_quality = {}
        stype_dict = {}
        # input_stoich = pseudo_stoich_filtered
        for pdb in tqdm(input_stoich.keys()):

            try:
                gene_stoich.update({pdb : input_stoich[pdb]})
                gene_quality.update({pdb : input_quality[pdb]})
                s = pdb.split('_&_')[0]
                if 'ECOLI' in s or s[0] =='E' or 'clean' in s:
                    stype_dict.update({pdb : 'ITASSER'})
                elif 'assembly' in s or len(s) == 4 or 'bio' in s:
                    stype_dict.update({pdb : 'PDB'})
                elif 'AF' in s and 'model_v2' in s:
                    stype_dict.update({pdb : 'ALPHAFOLD'})
                else:
                    stype_dict.update({pdb : 'SWISS'})

            except KeyError:
                log.warn("Could not find Structure Quality : {}".format(pdb))
#                 errors_from_gempro.update({pdb : 'missing quality'})
                continue

        df['gene_stoichiometry'] = pd.Series(gene_stoich)
        df['pdb_quality'] = pd.Series(gene_quality)
        df['stype'] = pd.Series(stype_dict)

        print ('Checking Dataframe...........\n-------------------')
        for pdb in tqdm(set(input_quality.keys()) - set(df.index.tolist())):
            errors_from_gempro.update({pdb : 'missing stoich'})

            
#         for gene_stoich
            
            
        #assign identical structures
        print ('Finding Unique Gene Stoich...........\n-------------------')
        unique_gs_list = []
        for gs in tqdm(df.get("gene_stoichiometry").values.tolist()):
            if gs not in unique_gs_list:
                unique_gs_list.append(gs)
            
        print ('Assigning Unique Gene Stoich as String...........\n-------------------')
        gene_stoichiometry_string = {}
        for gs in tqdm(unique_gs_list):
            dfg = df[df.gene_stoichiometry == gs]
            for index in dfg.index:
                gene_stoichiometry_string.update({index : str(gs)})
        df['gene_stoichiometry'] = pd.Series(gene_stoichiometry_string)

            
        
        print ('Assigning Identical Structures & AA Enzyme Length...........\n-------------------')
        gene_length_structures ={}
        
        
        rseq_len = df_gene_length["UniProtSeqLen"].to_dict()
        cseq_len = df_gene_length["AlleleomeSeqLen"].to_dict()
        
        
        ident_gene_stoich = {}
        for structure, gene_stoich in tqdm(df['gene_stoichiometry'].to_dict().items()):
            if gene_stoich in ident_gene_stoich:
                ident_gene_stoich[gene_stoich] += [structure]
            else:
                ident_gene_stoich.update({gene_stoich : [structure]})
                
                
        identical_structure_dict ={}

        for unique_gene_stoich , identical_structures in tqdm(ident_gene_stoich.items()):
        
#         #old
# #         for unique_gene_stoich in tqdm(unique_stoich_list):
# #             temp_df = df[df.gene_stoichiometry == unique_gene_stoich]
# #             identical_structures = temp_df.index.tolist()
            
                
            gs = ast.literal_eval(unique_gene_stoich)
            total_length = 0
            for g , g_stoich in gs.items():
                try:
                    seqlen = np.array([rseq_len[g],cseq_len[g]])
                    seqlen = seqlen[~np.isnan(seqlen)]

                   
                    length = seqlen.max()
                    total_length += length * g_stoich
                except KeyError:
                    log.warn("{} not in MG1655 genome".format(g))

            for i in identical_structures:
                gene_length_structures.update({i  : total_length})
                identical_structure_dict.update({i : str(identical_structures)})

        df['AA_length'] = pd.Series(gene_length_structures)
        df['identical_structures'] = pd.Series(identical_structure_dict)

        return df,errors_from_gempro
    except KeyboardInterrupt:
        df['AA_length'] = pd.Series(gene_length_structures)
        df['identical_structures'] = pd.Series(identical_structure_dict)
        return df, errors_from_gempro
    
def assign_structure_quality(df,df_gene_length,quality_key = 'pdb_quality'):
    # figure out the global quality and false/true percentage for each structure
    
    total_quality_dict = {}
    percent_true_dict = {}
    len_rseq = df_gene_length.UniProtSeqLen.to_dict()
    len_cseq = df_gene_length.AlleleomeSeqLen.to_dict()
    
    df = copy.deepcopy(df)
    print ("Finding Total Match Quality...........\n-------------------")
    for index, row in tqdm(df.iterrows()):
        try:

#             if index != 'P00561_1_461_3c1n_&_NOBFS_0':
#                 continue
#             print index
#             gs = ast.literal_eval(row.gene_stoichiometry)
            
            quality = row.get(quality_key)
            if type(quality) == str:
                quality = ast.literal_eval(row.get(quality_key))
            
            quality = filter_empty_chains(quality)
            gs = find_gs_from_pdb_quality(quality)
            
            
            
            
            percent_length_per_gene = {}

            for g , s in gs.items():
                
                seqlen = np.array([len_rseq[g],len_cseq[g]])
                seqlen = seqlen[~np.isnan(seqlen)]
                
                aa_length = seqlen.max()

                gl = aa_length
                percent_gl = float(gl * s) / float(row.AA_length) 
                percent_length_per_gene.update({g : percent_gl})

            total_quality = 0
            total_falsecount=0
            for g, info in quality.items():
        #         print '\n', g , info
                avg_quality = []
                false_count = []
                for chain, chain_quality in info.items():
        #             print chain_quality
                    avg_quality += [chain_quality[1]]
                    false_count += [chain_quality[0]]

                avg_quality = np.mean(avg_quality)
                false_count= np.mean(false_count)
        #         print g, avg_quality,'\n'
                total_quality += avg_quality * percent_length_per_gene[g]
                total_falsecount += false_count * percent_length_per_gene[g]
#             df.at[index, 'total_quality'] = total_quality
            total_quality_dict.update({index : total_quality})
#             df.at[index, 'percent_true'] = total_falsecount
            percent_true_dict.update({index : total_falsecount})

        except (KeyError, ValueError,AttributeError, IndexError) as e:
            log.warn("Cannot Assign Structure Quality : {}".format(index))
            continue
    
    df['total_quality'] = pd.Series(total_quality_dict)
    df['percent_true'] = pd.Series(percent_true_dict)

    return df

def filter_empty_chains(pdb_quality):
    
    if type(pdb_quality) == str:
        pdb_quality = ast.literal_eval(pdb_quality)
    
    new_gene_quality = {}
    for gene, gene_quality in pdb_quality.items():
        new_chain_info = {}
        for chain, chain_info in gene_quality.items():
            if chain in ['-','_','X'] and chain_info[1] < 5:
                continue
            new_chain_info.update({chain : chain_info})
        if new_chain_info == {}:
            continue
        new_gene_quality.update({gene : new_chain_info})
        
    return new_gene_quality

def find_gs_from_pdb_quality(pdb_quality):
    if type(pdb_quality) == str:
        pdb_quality = ast.literal_eval(pdb_quality)
    
    new_gs = {}
    for gene, gene_info in pdb_quality.items():
        new_gs.update({gene : len(gene_info)})
    return new_gs

def assign_quality_ranks(df):
    print ("Ranking Matches...........\n-------------------")

    ident_gene_stoich = {}
    for structure, gene_stoich in tqdm(df['gene_stoichiometry'].to_dict().items()):
        if type(gene_stoich) == dict:
            gene_stoich  = str(gene_stoich)
        
        if gene_stoich in ident_gene_stoich:
            ident_gene_stoich[gene_stoich] += [structure]
        else:
            ident_gene_stoich.update({gene_stoich : [structure]})
           
    pdb_quality_dict_all = df['total_quality'].to_dict()
    percent_true_dict_all = df['percent_true'].to_dict()

    
    
    for gs, i_list in tqdm(ident_gene_stoich.items()):
#         i_list = ast.literal_eval(i_list)
#         temp_df = df[df.index.isin(i_list)]
        pdb_quality_dict = {}
        percent_true = {}
        for structure in i_list:
            pdb_quality_dict.update({structure :pdb_quality_dict_all[structure]})
            percent_true.update({structure :percent_true_dict_all[structure]})
        
        
        
#     for list_string in tqdm(df.identical_structures.unique().tolist()):
#         temp_df = df[df.index.isin(i_list)]

#         pdb_quality_dict = temp_df['total_quality'].to_dict()
        sorted_x = sorted(pdb_quality_dict.items(), key=operator.itemgetter(1), reverse=True)
        ranked_x = rank_pdb_structures(sorted_x)
        for pdb, ranking in ranked_x.items():
            df.at[pdb, 'quality_rank'] = ranking[0]

        sorted_x = sorted(percent_true.items(), key=operator.itemgetter(1), reverse=True)
        ranked_x = rank_pdb_structures(sorted_x)
        for pdb, ranking in ranked_x.items():
            df.at[pdb, 'percent_true_rank'] = ranking[0]
    return df

def rank_pdb_structures(sorted_x):
    """
    ranks a sorted dictionary of {pdb_id : pdb_quality}
    
    returns {pdb_id: pdb_rank}
    
    is used for quality and false count rankings
    
    """
    
    rank = 1
    rank_dict = {}


    for i, quality in enumerate(sorted_x):
        #print i, quality[0], quality[1]
        if i == 0:
            pass
        else:
            if quality[1] == sorted_x[i-1][1]:
                rank = rank
            else:
                rank += 1

#         print i, rank, quality[0], quality[1]
        rank_dict.update({quality[0] : (rank, quality[1])})
    return rank_dict


def find_best_resolution(ilist_min,structure_info,stype):
    dfresolution  = pd.DataFrame()
         
    if stype == 'PDB':
        for structId in ilist_min:
            pdbId = structId.split('-assembly')[0].upper()
            dfresolution.loc[structId,'pdb_entry'] = pdbId
            resolution = structure_info[structure_info.pdb_id == pdbId].resolution.values[0]
            dfresolution.loc[structId,'resolution'] = resolution

        dfresolution = dfresolution.sort_values(by ='resolution', ascending=True)
        return dfresolution.first_valid_index()
    
    elif stype =='SWISS':
        index_keep={}

        for structId in ilist_min:
#             print structId     
            uniprotId, from_, to_, template = structId.split('_&')[0].split('_')
            
            dfs_uni = structure_info[structure_info.UniProtKB_ac == uniprotId]
            dfs_uni_from = dfs_uni[dfs_uni.get('from') == int(from_)]
            dfs_uni_from_to = dfs_uni_from[dfs_uni_from.get('to') == int(to_)]
            for index, templateId in dfs_uni_from_to.get('template').to_dict().items():
                if template in templateId:
                    index_keep.update({index : structId})
                    
        dfs = structure_info[structure_info.index.isin(index_keep.keys())]
        dfs=dfs.sort_values(by = 'qmean_norm',ascending=False)
        best_index = dfs.first_valid_index()
        best_structure = index_keep[best_index]
        
#         print (index_keep.keys(), best_structure)
        return best_structure

            
#     for index
#             resolution = structure_info[structure_info.pdb_id == pdbId].resolution.values[0]
#             dfresolution.loc[structId,'resolution'] = resolution
            
            
    

def find_best_structures(df4, query_type =['PDB','ITASSER', 'SWISS','ALPHAFOLD'], structure_info = pd.DataFrame() ):

    best_structures = {}
    quality_rank_df4  = df4['quality_rank'].to_dict()

    
    print ("Finding Best Matches...........\n-------------------")

    
    
    for stype in query_type:

        temp_stype = df4[df4.stype == stype]
        temp_stype_index_list = temp_stype.index.tolist()
        best_list = []
        ilist_total = temp_stype.identical_structures.unique()
        print ("{} Unique gene stoichiometries from {} {} structures...".format(len(ilist_total),len(temp_stype),stype))
        for ilist in tqdm(ilist_total):
            ilist = ast.literal_eval(ilist)
            if stype == 'PDB':
                ilist_quality_dict = {k : quality_rank_df4[k] for k in ilist if k[0:3] != 'AF-' and 'ECOLI' not in k and 'clean' not in k and len(k.split('_&')[0]) in [18,19]}
    #         elif stype == 'ITASSER':
    #             ilist_quality_dict = {k : quality_rank_df4[k] for k in ilist if k[0:3] != 'AF-' and 'ECOLI' in k or 'clean'  in k }


            else:
                ilist_quality_dict = {k : quality_rank_df4[k] for k in set(ilist).intersection(set(temp_stype_index_list))}

            ilist_min = []
            for k,v in ilist_quality_dict.items():
#                 print (v,np.min(list(ilist_quality_dict.values())))
                if v == np.min(list(ilist_quality_dict.values())):
                    ilist_min += [k]
#                     print (k)
#                     break
                    
         
                
                
            if len(ilist_min) > 1:
                ilist_min = [find_best_resolution(ilist_min,structure_info,stype)]

         
            best_list += ilist_min

        best_structures.update({stype: best_list})
    return best_structures
