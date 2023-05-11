import pandas as pd
import numpy as np
import ast
import os
import sys
import logging
log = logging.getLogger(__name__)


#from structure_filtering import *
#from figure_functions import *

import datetime

#import ssbio
import copy
import operator
from copy import deepcopy
import json
import cloudpickle as cPickle
from tqdm import tqdm_notebook as tqdm
#import seaborn as sns

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

def create_dataframe_per_enzyme(enzyme,  output_dict,parent_stoich):
    
    #this df will filter out overmatched children

    v = output_dict[enzyme] #the children matches to the enzyme
    enzyme_df = pd.DataFrame(columns = parent_stoich[enzyme].keys())


    for b_id,s in parent_stoich[enzyme].items():
        enzyme_df.at[enzyme, b_id] = s
    oa_list = []
    
    
    for k, dict_k in v.items():
        
        if 'child' in dict_k.keys():
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
            
            
        if  list(enzyme_df.loc[index].get_values()) == list(enzyme_df.loc[enzyme].get_values()):        
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

def obtain_sets(input_dict):
    output_set = []
    for k,v in input_dict.items():
        if len(v)  == 0:
            continue
        output_set.append(k)
    return set(output_set)

def filter_output_dict_for_children(output_dict, parent_stoich):

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
#         if children != {}:
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


###################################################################################
# The BFS STUFF



#take level above it and make the next one
def buildLevel(queue):
    levelHash = {}
    level = []
    for node in tqdm(queue):
        #create level by adding children to each queued item
        for child in children:
            #checks go here
            leaf = add(node, child)
            leafHash = buildHash(leaf)
            #check if existis else add
            # print leafHash
            # print levelHash
            if leafHash not in levelHash:
                #add leafHash & exec code
                levelHash[leafHash] = 1
                if alive(leaf):
                    level.append(leaf)
                    alive_leafs.append(leafHash)

    # print level;
    return level;

def buildHash(node):
    hasher = [0]*len(children)
    for child in node['children']:
        hasher[child['child']] += 1
    # print ''.join(map(lambda x: str(x),children))
    # print(hasher)
    return ''.join(map(lambda x: str(x),hasher))

def main():
    global parent
    sum = {}
    for key in parent:
        sum[key] = 0
    root = {
        'sum':sum,
        'children': []
        }
    queue = [root]
    while len(queue) > 0:
        #process queue
        # print(queue)
        queue = buildLevel(queue)

def add(a, b):
    # print a
    # print b
    # a = a['sum']
    sum = {}
    for key in parent:
        sum[key] = a['sum'][key] + b[key]
    return {
            'sum': sum,
            'children': a['children'] + [b]
            }

def alive(node):
    # under = True
    match = True

    if buildHash(node) in seen:
        return False
    else:
        seen[buildHash(node)] = 1
    # compare to Parent
    for key in parent:
        if node['sum'][key] == parent[key]:
            # under = under and True
            match = match and True
        elif node['sum'][key] < parent[key]:
            # under = under and True
            match = match and False
        elif node ['sum'][key] > parent[key]:
            return False

    if match:
        answers.append(node['children'])
        return False
    return True
############################################################################################
def convert_answer_to_pdb_dict(answers_string , children_dict):

    subset_pdb_answers ={}
    answer_count = 0
    for ans in answers_a:
        pdb_ans = {}

        for child_num, child_count in enumerate(ans.split('.')):
            if int(child_count) == 0:
                continue
            else:
                pdb = get_pdb_from_childnum(children_dict, child_num)
                pdb_ans.update({pdb : int(child_count)})
        subset_pdb_answers.update({'subset_match_%i' % answer_count : pdb_ans})
        answer_count += 1
    return subset_pdb_answers

def get_pdb_from_childnum(children_dict, child_num):
    for pdb, genes in children_dict.items():
        if genes['child']== child_num:
            break
    return pdb

def convert_answer_to_pdb_dict_incompletes(answers_string ,children_dict ):

    subset_pdb_answers ={}
    answer_count = 0
    for ans in answers_string:
        pdb_ans = {}

        for child_num, child_count in enumerate(ans):
            if int(child_count) == 0:
                continue
            else:
                pdb = get_pdb_from_childnum(children_dict, child_num)
                pdb_ans.update({pdb : int(child_count)})
        subset_pdb_answers.update({'incomplete_match_%i' % answer_count : pdb_ans})
        answer_count += 1
    return subset_pdb_answers



def run_bfs(bfs_children_dict,parent_stoich,  test_pyr  = False, test_ribo = False,    bfs_answers = {}):
    
#     bfs_answers = {}
    alive_dict = {}
    global parent
    global children
    global seen
    global alive_leafs
    global answers
    global answers_a
    
    try:
        for enzyme, children_dict in tqdm(bfs_children_dict.items()):
            if enzyme in bfs_answers:
                continue
            if children_dict == {}:
                continue
                
                
            if not test_ribo and 'ribosome' in enzyme:
#                 print enzyme, '<<<<skipped'
                continue
            
            

            if not test_pyr and enzyme ==  'PYRUVATEDEH-CPLX':
#                 print enzyme, '<<<<skipped'
                continue
            
            alive_leafs = []

#             print enzyme
            parent = parent_stoich[enzyme]
            children = children_dict.values()
        #     for key in children_dict:
        #         children.append(children_dict[key])

            answers = []

            seen={}

            main()

            answers_a = []

            for ans in answers:
                hasher = [0]*len(children)
                for child in ans:
                    hasher[child['child']] += 1
                answers_a.append('.'.join(map(lambda x: str(x),hasher)))

            answer_dict = convert_answer_to_pdb_dict(answers_string = answers_a, children_dict = children_dict)
            if answer_dict != {}:
                bfs_answers.update({enzyme : answer_dict})
            else:
                alive_dict.update({enzyme : alive_leafs})
        return bfs_answers, alive_dict
    except KeyboardInterrupt:
        return bfs_answers, alive_dict




def filter_incomplete_matches(i, e):
    """
    filters the solution space of i for the "mostcomplete" incomplete match
    """
    
    def build_I_df_from_string(i,e):
        """
        this just builds a pandas dataframe from the string solutions for incomplete matches
        it allows the algo to parse the df efficiently
        """
        matches = i
        if len(matches) > 700:
#             print '\n' , e
            print ('Before Filtering Incomplete :  %i ' %len(matches))
        df = pd.DataFrame(columns=range(len(matches[0])))

        for index in matches:
            for x, num in enumerate(list(index)):
                df.at[index,x  ] = int(num)
            non_zero = list(df.loc[index].get_values()[0:len(matches[0])].nonzero()[0])
            non_zero = [str(s) for s in non_zero]
        #     print non_zero
            s = '-'
            s=s.join(non_zero)
            df.loc[index, 'non_zero'] = s
            df.loc[index, 'total_children'] = np.sum([int(s) for s in list(index)])
        return df

    def iteration_step(df,input_queue):
        #this will filter for
        # 1) the best per each unique MI
        # 2) all the children of the MI that are not up to par
        
        """
        the iterative step in reducing the solution space of incomplete matches
        in each step, there is a solution space and a queue
        
        this step will check the solution space for each Mutual Information in the queue
        it will 1) pick only the best solution that shares the same MI of interest
        it will 2) find all other solutions with MI that are subsets of the MI of interests 
        and then eliminate those that do not provide any additional information
        """

        previous_queue = []
        queue = copy.deepcopy(input_queue)
        while queue!= previous_queue:
            previous_queue = queue
            for a in queue:
                a_split = a.split('-')

                children_of_a = []

                for b in queue:
                    b_split = b.split('-')
                    if a == b:
                        continue
                    if set(a_split).intersection(set(b_split)) == set(b_split):
                        #b is a child of a
                        children_of_a.append(b)
                        continue
                    if set(a_split).intersection(set(b_split)) == set(a_split):
                        children_of_a.append(a)
                        # a is a child of b

                queue = set(queue) - set(children_of_a)


        children = set(input_queue) - set(queue)


        #create the dictionary to check parents
        parent_checks = {}
        for c in children:
            for p in queue:
                if set(c.split('-')).intersection(set(p.split('-'))) == set(c.split('-')):# in p:
                    if p not in parent_checks:
                        parent_checks.update({p : [c]})
                    else:
                        parent_checks[p]+= [c]

        v_l = []
        for k,v in parent_checks.items():
            v_l += v



        #this should technically pick the best indices for each MI parsed
        for mi in df.non_zero.unique().tolist():
            t = df[df.non_zero == mi]

            max_pbd_index = t[t.total_children == max(t.total_children.tolist())].index.tolist()[0]

            #delete each row in the df unless the child stoich is greater than the max pdbchild stoich for that column
            for index, row in t.iterrows():
                if index == max_pbd_index:
                    continue
                else:
                    keep = False
                    for column in df.columns:
                        try:
                            int(column)    
                        except:
                            continue
                        if t.loc[index, column] > t.loc[max_pbd_index, column]:
                            keep = True
                            break
                if not keep:
                    df = df.drop([index])
    #     print '\t Queue Solution (MI) : %i' %len(df)
        for mi, mi_subset in parent_checks.items():

            t = df[df.non_zero == mi]
            if len(t) == 0:
                continue
            max_pbd_index = t[t.total_children == max(t.total_children.tolist())].index.tolist()[0]

            #delete each row in the df unless the child stoich is greater than the max pdbchild stoich for that column

            df_left = df[df.non_zero.isin(mi_subset)]
            for index, row in df_left.iterrows():
                if index == max_pbd_index:
                    continue
                else:
                    keep = False
                    for column in df_left.columns:
                        try:
                            int(column)    
                        except:
                            continue
                        if df_left.loc[index, column] > t.loc[max_pbd_index, column]:
                            keep = True
                            break
                if not keep:
                    df = df.drop([index])
    #     print '\t Queue Solution (PC) : %i' %len(df)

        return df, parent_checks.keys()
    
    
    def reduce_solution_space(df):
        """
        each queue is built up of the most complete mutual information among the solution space
        after filtering each through the iterative step (above), the set of most relevant mutual information
        is removed from the queue and the best MIs are found from the remaining MIs avaialble.
        in principle, if the length of children is the maximum, it will use all MIs of that same length, perform the 
        iterative step, and then remove those MI for the list. the next iteration will used MIs of length (Max - 1)        
        """
        
        
        queue_total = set(df.non_zero.unique().tolist())
        previous_queue = []
        printing = False
        if len(df) > 700:
            printing = True
        while  queue_total!= previous_queue:

            df, parent_checks = iteration_step(df, input_queue=queue_total)

            previous_queue = copy.deepcopy(queue_total)
            queue_total -= set(parent_checks)
        
        if printing:
            print ('After Filtering Incomplet :  %i ' %len(df))

        return df.index.tolist()
    
    df = build_I_df_from_string(i, e)
    i = reduce_solution_space(df)
    return i

def remove_1_1_matches_from_incomplete_matches(alive_dict , full_dict ):
    new_alive_dict = {}
    for k,v in alive_dict.items():
        if k not in full_dict:
            new_alive_dict.update({k : v})
    return new_alive_dict

def find_unique_genes(df):
    glist = []
    for gene_stoichiometry in df.gene_stoichiometry.unique():
        gene_stoichiometry = ast.literal_eval(gene_stoichiometry)
        glist += gene_stoichiometry.keys()
    glist = set(glist)
#     print len(glist)
    return glist



