import pandas as pd
import numpy as np
from tqdm import tqdm_notebook as tqdm
import ast
import copy

def get_data(df, parent_input, stype = ['PDB'], quality_cutoff = 0,passed_index = [], quality_source = "pdb_quality"):

    dfpass = df[df.stype.isin(stype)]
    dfpass = dfpass[dfpass.total_quality >= quality_cutoff]
    keep_index = dfpass.index.tolist() + passed_index
    df=df[df.index.isin(keep_index)]
    
    if len(df) == 0:
        return pd.DataFrame(), {},{},{} ,  {},  {}, {}, {}
    real_structure_stoich,real_structure_quality =create_pseudo_structures_dicts(df,quality_source = quality_source)
    output_dict= match_pdbs_to_me_model_stoich(parent_stoich = parent_input, input_stoich = real_structure_stoich)
    to_bfs,  full_1_1, o_dict, oa_dict= filter_output_dict_for_children(output_dict= output_dict, parent_stoich=parent_input)

    
    
    return df, real_structure_stoich,real_structure_quality,output_dict ,  to_bfs,  full_1_1, o_dict, oa_dict


def create_pseudo_structures_dicts(df, quality_source = "pdb_quality"):
    structure_quality = df.get(quality_source).to_dict()
    for k, v in structure_quality.items():
        if type(v) == str:
            v = ast.literal_eval(v)
        
        structure_quality.update({k : v})


    real_structure_quality = {}

    for k, gene_info in structure_quality.items():
        new_gene_info = {}
        for gene , chain_info in gene_info.items():
            new_chain_info = {}
            for chain, chain_quality in chain_info.items():
                if chain in ['_','-','',' '] and chain_quality[1] < 3:
                    continue
                new_chain_info.update({chain : chain_quality})
            if new_chain_info != {}:
                new_gene_info.update({gene : new_chain_info})

        if new_gene_info != {}:
            real_structure_quality.update({k : new_gene_info})

    real_structure_stoich = {}
    for sid, stoich in tqdm(real_structure_quality.items()):
        gene_newinfo  = {}
        for gene, geneinfo in stoich.items():
            gene_newinfo.update({gene : len(geneinfo)})
        real_structure_stoich.update({sid : gene_newinfo})   
    return real_structure_stoich,real_structure_quality 


def match_pdbs_to_me_model_stoich(parent_stoich , input_stoich ):
    output_dict ={}
    found_a_child = []
    for enzyme in tqdm(parent_stoich):
    #     if enzyme != 'RpoS_mono':
    #         continue
        parent = parent_stoich[enzyme]
        all_me_genes = set(parent.keys()) 
        enzymes_children_dict = {}
        for pdb in input_stoich:
            child = input_stoich[pdb] #use dict because it is way faster than df.loc[pdb, "gene_stoichiometry"]
            all_pdb_genes = set(child.keys())

            if len(all_me_genes.intersection(all_pdb_genes)) == 0:
                continue
            else:
    #             this will pring the enzyme, its me model stoich, and then all associated children this loop finds
    #             if enzyme not in found_a_child:
    #                print enzyme, parent
                found_a_child.append(enzyme)
    #             print '\t -->', pdb,  INPUT_GEM_STOICH[pdb]
                child = add_zero_value_to_subcomplexes(parent=parent, child=child)
                enzymes_children_dict.update({pdb:child})
        output_dict.update({enzyme:enzymes_children_dict})
    return output_dict

def add_zero_value_to_subcomplexes(child, parent):
    """
    for the algorithm to work, a child class must be assigned a 
    key:value of {'gene' : 0} if that gene is missing
    """
    
    new_child = copy.deepcopy(child)
    for key in parent:
        if key not in list(child.keys()):
            new_child.update({key : 0})
    return new_child

def filter_output_dict_for_children(output_dict, parent_stoich):

    def obtain_sets(input_dict):
        output_set = []
        for k,v in input_dict.items():
            if len(v)  == 0:
                continue
            output_set.append(k)
        return set(output_set)
    
    
    
    bfs_children_dict = {}
    full_1_1_match_dict = {}
    o_dict = {}
    oa_dict = {}
    
    
    

    for enzyme, pdb_dict in tqdm(output_dict.items()):

        parent = parent_stoich[enzyme]
        enzyme_df, oa_list, o_list,one_to_one_list =create_dataframe_per_enzyme(enzyme=enzyme, parent_stoich=parent_stoich, output_dict=output_dict)

        children = {}
        child_count = 0
        for pdb, pdb_stoich in pdb_dict.items():

            if pdb not in enzyme_df.index:
                continue
            else:
                pdb_stoich = add_zero_value_to_subcomplexes(child=pdb_stoich, parent=parent)
                pdb_stoich.update({'child' : child_count})
                children.update({pdb : pdb_stoich})
                child_count +=1
        if children != {}:
            bfs_children_dict.update({enzyme : children})
        full_1_1_match_dict.update({enzyme : one_to_one_list})
        o_dict.update({enzyme : o_list})
        oa_dict.update({enzyme : oa_list})
    
    
    f_set = obtain_sets(full_1_1_match_dict)
    o_set = obtain_sets(o_dict)
    oa_set = obtain_sets(oa_dict)
    
    new_full_1_1_match_dict = {}
    for k in f_set:
        new_full_1_1_match_dict.update({k : full_1_1_match_dict[k]}) 
        
    new_o_dict = {}
    for k in o_set:
        new_o_dict.update({k : o_dict[k]}) 
        
    new_oa_dict = {}
    for k in oa_set:
        new_oa_dict.update({k : oa_dict[k]}) 
    
    
    
    return bfs_children_dict,  new_full_1_1_match_dict, new_o_dict, new_oa_dict


def create_dataframe_per_enzyme(enzyme,  output_dict,parent_stoich):
    
    #this df will filter out overmatched children

    v = output_dict[enzyme] #the children matches to the enzyme
    enzyme_df = pd.DataFrame(columns = list(parent_stoich[enzyme].keys()))


    for b_id,s in parent_stoich[enzyme].items():
        enzyme_df.at[enzyme, b_id] = s
    oa_list = []
    
    
    for k, dict_k in v.items():
        
        if 'child' in list(dict_k.keys()):
            dict_k_genes = len(dict_k)- 1
        else:
            dict_k_genes = len(dict_k)
        
        if dict_k_genes > len(enzyme_df.columns):
            oa_list.append(k)
#             print k, dict_k
            continue

        for b ,  s in dict_k.items():
            if b == 'child':
                continue
    #         if b not in ribo_df.columns:
    #             print k, b, s
    #             continue
            enzyme_df.at[k,b] = s
    enzyme_df = enzyme_df.fillna(value = 0 )
    
    o_list = []
    one_to_one_list = []
    is_over = []
    for index in enzyme_df.index:
        

        if index == enzyme:
            continue
            
            
        if  list(enzyme_df.loc[index].values) == list(enzyme_df.loc[enzyme].values):        
#             print '1 -->' , index
            one_to_one_list.append(index)
            continue    
            
        found_overset = True
        for c in enzyme_df.columns:
            if enzyme_df.at[index, c] <= enzyme_df.at[enzyme,c]:
                found_overset= False
                continue
            is_over += [index]
#             print (index)
                
#             print c, enzyme_df.at[index, c]  , '>>' , enzyme_df.at[enzyme,c], '\t', index
        if found_overset:
           o_list.append(index)
            #these are the indexes that have are overmatches
            

    enzyme_df=enzyme_df[enzyme_df.index.isin(o_list + is_over) == False] 
#     enzyme_df=enzyme_df.drop(is_over) 
    enzyme_df = enzyme_df.drop(one_to_one_list)
    
    return     enzyme_df, oa_list, o_list,one_to_one_list


