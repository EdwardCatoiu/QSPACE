# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def run_004A_get_matches(dfbest,
                         proteinTargetsNew,
                         quality_cutoff  = {'PDB': 70,'SWISS': 70,'ITASSER': 70,'ALPHAFOLD': 70, 'ALPHAFOLD_MULTIMER' : 0}, #multimer results already checked,
                         remove_bad_structures = True,
                        ):
    
    keep_index = []
    for stype, cutoff in quality_cutoff.items():
        dfb = dfbest[dfbest.stype == stype]
        dfb = dfb[dfb.total_quality >= cutoff]
        keep_index += dfb.index.tolist()
    if remove_bad_structures:
        print ("Removing {} Pseudo-Structures that do not meet the quality thresholds...".format(len(dfbest[dfbest.index.isin(keep_index) == False])))

        dfbest = dfbest[dfbest.index.isin(keep_index) == True]

    parent_gene_stoich = proteinTargetsNew.gene_stoich.to_dict()
    for k, v in parent_gene_stoich.items():
        if type(v) == str:
            v = ast.literal_eval(v)
        parent_gene_stoich.update({k : v})

    results = prepBFS.get_data(df=dfbest, 
                               parent_input=parent_gene_stoich,
                               stype = list(quality_cutoff.keys()),
                               quality_cutoff = 0,  #already removed from dfbest structures above^
                               passed_index = [], 
                               quality_source = "pdb_quality")

    [df, real_structure_stoich,real_structure_quality,output_dict ,  to_bfs,  full_1_1, o_dict, oa_dict] = results



    to_bfs_filtered = copy.deepcopy(to_bfs)
    c = 0
    for k , slist in full_1_1.items():

        is_good_match = False
        for structure in slist:
            #test to see if the full match structure quality is better than the 80% cutoff

            stype = dfbest.loc[structure,"stype"]
            if dfbest.loc[structure,"total_quality"] > quality_cutoff[stype]:
                is_good_match = True
                break
        if is_good_match and k in to_bfs_filtered:
            c +=1

            del to_bfs_filtered[k]

#     print "{} protein complexes removed from BFS".format(c)


    data = {'ProteinTarget' : parent_gene_stoich ,
            'FullMatch': full_1_1,
            'BFSAvailable': to_bfs,
            'BFSRequired': to_bfs_filtered}
    return data


def renumber_children(childrenDict):
    for i, structure in enumerate(list(childrenDict.keys())):
        childrenDict[structure].update({"child" : i})
    return childrenDict
        
def simplify_children(cplx, childrenDict,dfbest):
    simplified = {}
    for structure, gene_stoich in childrenDict.items():
        for gene, stoich in gene_stoich.items():
            if gene =='child':
                continue
            if stoich == 0 :
                continue
            if structure not in simplified:
                simplified.update({structure : {gene: stoich}})
            elif gene not in simplified[structure]:
                simplified[structure].update({gene: stoich})
                       
    dfsimplified = pd.DataFrame()
    dfsimplified['gene_stoich']= pd.Series(simplified)
    for index in dfsimplified.index:
        dfsimplified.loc[index,'total_quality'] = dfbest.loc[index,'total_quality']

    for index, row in dfsimplified.iterrows():
        dfsimplified.loc[index,'gene_stoich'] = str(dfsimplified.loc[index,'gene_stoich'] )

    indexKeep = []
    for gs in dfsimplified.gene_stoich.unique():
        dfgs = dfsimplified[dfsimplified.gene_stoich == gs]
        dfgs = dfgs.sort_values(by ='total_quality', ascending = False)
        indexKeep += [dfgs.first_valid_index()]

    dfsimplified = dfsimplified[dfsimplified.index.isin(indexKeep)]

    new_childrenDict = {}
    for index in dfsimplified.index.tolist():
        new_childrenDict.update({index : childrenDict[index]})

    return new_childrenDict


manual_bfs_answers = {}
cplx = 'FLAGELLAR-MOTOR-COMPLEX'
new_match = {'7nvg-assembly1.cif_&_MATCH_0': 1,
             'AF-P15070-F1-model_v2_&_NOBFS_0': 102,
             "AF-P0ABZ1-F1-model_v2_&_NOBFS_0" : 34,
             "AF-P06974-F1-model_v2_&_NOBFS_0" : 34,
             'P0AF06_101_274_2zvz_&_NOBFS_0' : 1,
             'P09348_124_253_5zfp_&_NOBFS_0' : 1
}
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})

cplx = 'CPLX0-7452'
new_match = {"P75937_2_402_6jzt_&_NOBFS_0" : 6,
             "AF-P06974-F1-model_v2_&_NOBFS_0":34,
             "AF-P0ABZ1-F1-model_v2_&_NOBFS_0" : 26,
             "7cg0-assembly1.cif_&_NOBFS_0": 1,
             'P75938_1_247_6jzr_&_NOBFS_0': 1,
             "P0ABX5_1_260_6jzr_&_NOBFS_0": 1,
             '7cbl-assembly1.cif_&_NOBFS_0' : 1,
             'AF-P33235-F1-model_v2_&_NOBFS_0': 1,
             'AF-P29744-F1-model_v2_&_NOBFS_0': 1,
             'AF-P0AF06-F1-model_v2_&_NOBFS_0':1,
             'AF-P09348-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P04949-F1-model_v2_&_NOBFS_0':1,
             'P24216_43_411_5h5v_&_NOBFS_0' : 1,
             'AF-P25798-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P15070-F1-model_v2_&_NOBFS_0' : 1             
            }
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})


cplx = 'PYRUVATEDEH-CPLX'
new_match = {'P0AFG8_57_887_1l8a_&_NOBFS_0':12,
             '4jdr-assembly1.cif_&_MATCH_1' : 6, 
             'ODP2_ECOLI_model1_clean_residues_removed_&_NOBFS_0':24}
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})

cplx = 'CPLX0-3803'
new_match = {"3sxu-assembly1.cif_&_NOBFS_0" : 4, 
             "P0A988_1_366_1mmi_&_NOBFS_0" : 2,
             'CPLX0-2361_AlphaMulti' : 3,
             "P28631_1_329_1a5t_&_NOBFS_0" :1,
             "1jqj-assembly3.cif_&_NOBFS_0":1,
             "AF-P06710-F1-model_v2_&_NOBFS_0" : 8 
            }
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})


cplx = 'ribosome'
new_match = {'3j9y-assembly1.cif_&_MATCH_1' : 1,
            '2mlx-assembly1.cif_&_NOBFS_0' : 1,
            'AF-P0AG67-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P68191-F1-model_v2_&_NOBFS_0' : 1,
             '7n2c-assembly1.cif_&_NOBFS_0' : 1,
             'AF-P0A7V0-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P0AG59-F1-model_v2_&_NOBFS_0' : 1
            }
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})


cplx = 'CPLX0-3382'
new_match = {'6hcg-assembly1.cif_&_MATCH_0' : 1,
             'P45759_106_484_4kss_&_NOBFS_0' : 1,
             'P36678_7_149_3jc9_&_NOBFS_0' : 6,
             'AF-P45763-F1-model_v2_&_NOBFS_0' : 6,
             'AF-P41441-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P41442-F1-model_v2_&_NOBFS_0' : 1,
             'P41443_30_169_2knq_&_NOBFS_0' : 1,
             'AF-P45760-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P45761-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P45762-F1-model_v2_&_NOBFS_0' : 1,
             'AF-P25960-F1-model_v2_&_NOBFS_0' : 1,
            }
manual_bfs_answers.update({cplx : {'subset_match_0': new_match}})

def run_004A_MultiStructureMatching(to_bfs_filtered,
                                    proteinTargetsNew,
                                    dfbest,
                                    test_pyr  = True,
                                    test_ribo = True,
                                    bfs_answers = {},
                                   ):
    parent_gene_stoich = proteinTargetsNew.gene_stoich.to_dict()
    for k, v in parent_gene_stoich.items():
        if type(v) == str:
            v = ast.literal_eval(v)
        parent_gene_stoich.update({k : v})
        
        
    for cplx, childrenDict in tqdm(to_bfs_filtered.items()):
        simplified = simplify_children(cplx, childrenDict,dfbest)
        simplified = renumber_children(simplified)
        to_bfs_filtered.update({cplx : simplified})
        
   
    if len(bfs_answers) > 0:
#         print "Found manually added multi-structure matching results for {} protein complexes...".format(len(bfs_answers) )
        
        del_list = []
        for cplx in bfs_answers:
            if cplx not in parent_gene_stoich:
                del_list +=[cplx]
        for cplx in del_list:
            del bfs_answers[cplx]
#         print "Removed {} protein complexes that fall outside of the current query...".format(len(del_list) )

#         print "Speeding up multi-structure matching by using manually found best matches for {} protein complexes...\n>{}".format(len(bfs_answers),bfs_answers.keys() )
        
#     print "Running BFS multi-structure matching..."
    bfs_answers, alive_dict = bfs_utils.run_bfs(bfs_children_dict=to_bfs_filtered,
                                        parent_stoich = parent_gene_stoich,
                                        bfs_answers=bfs_answers ,
                                        test_pyr  = test_pyr,
                                        test_ribo = test_ribo,
                                       )
#     print "Found Multi-structure matches for {} proteins".format(len(bfs_answers) )
    return bfs_answers

def find_min_num(e,match_dict):
    
    temp = {}
    min_num = 100
    for m , det in match_dict[e].items():
        if np.sum(list(det.values())) <min_num:
            min_num = np.sum(list(det.values())) 
            
    for m , det in match_dict[e].items():
        if np.sum(list(det.values())) > min_num + 2:
            continue
        temp.update({m : det})
    return temp

   
def run_004A_GatherMatches(bfs_answers,
                        full_1_1,
                        dfparent,
                        dfbest,
                       ):
    
    column_list= ['cplx','match_number','cplx_gene_stoich','cplx_structure_stoich','num_pdbs','match_quality','type_of_match']
    
    ####################################################
#     print '---------------------------------------------------'
#     print "Gathering multi-structure matches..."
    checked = []
    match_dict = bfs_answers


    dfmatch = pd.DataFrame(columns= column_list)


    for cplx in tqdm(dfparent.index.tolist()):
        try:
            match_info = match_dict[cplx]
        except KeyError:
            continue
        if cplx in checked:
            continue
#         print len(match_info), cplx

        #shorten the BFS answer matches with this
        if len(match_info) >400:
            match_info = find_min_num(cplx,match_dict)
#             print cplx, len(match_info)

        for match, match_structure_info in tqdm(match_info.items()):
            
            index = "{}_&_{}".format(cplx, match)
            total_AA_length = 0
            dfmatch.at[index,"cplx"] = cplx
            dfmatch.at[index,"match_number"] = match
            dfmatch.at[index,"type_of_match"] = 'Subset'
            dfmatch.at[index,"cplx_gene_stoich"] = str(dfparent.loc[cplx,'gene_stoich'])

            num_pdbs = np.sum(list(match_structure_info.values()))
            dfmatch.at[index,"num_pdbs"] = num_pdbs
            dfmatch.at[index,"cplx_structure_stoich"] = str(match_structure_info )


            dfs_info = dfbest[dfbest.index.isin(match_structure_info)]
    #         if len(dfs_info) != len(match_structure_info):


            type_of_structures = dfs_info['stype'].to_dict()
            quality = dfs_info['total_quality'].to_dict()
            aa_covered = dfs_info['AA_length'].to_dict()
            gene_stoichiometry = dfs_info['gene_stoichiometry'].to_dict()
            dftemp = pd.DataFrame()
            dftemp['type_of_structures'] = pd.Series(type_of_structures)
            dftemp['quality'] = pd.Series(quality)   
            dftemp['aa_covered'] = pd.Series(aa_covered)
            dftemp['structure_stoich'] = pd.Series(match_structure_info)
            dftemp['gene_stoichiometry'] = pd.Series(gene_stoichiometry)
            total_AA = dftemp['aa_covered'] * dftemp['structure_stoich']
            total_AA= total_AA.sum()
            dfmatch.loc[index, 'total_AA'] = total_AA
            for i in dftemp.index:
                dftemp.at[i, 'total_AA'] = total_AA

            dftemp['ratio'] = dftemp['aa_covered'] * dftemp['structure_stoich']  / dftemp['total_AA'] 

            match_quality =  dftemp['ratio'] * dftemp['quality'] 

            for stype in dftemp.type_of_structures.unique():
                dfmatch.at[index,str(stype)] = dftemp[dftemp.type_of_structures == stype].ratio.sum()

            dfmatch.at[index,"match_quality"] = match_quality.sum()
    #         index +=1
        checked+=[cplx]
        
        
   
    ############################################################
#     print '---------------------------------------------------'
#     print "Gathering single-structure matches..."
    checked= []
    match_dict = full_1_1

    for cplx in tqdm(dfparent.index.tolist()):
        try:
            fullmatch_list = match_dict[cplx]
        except KeyError:
            continue

        for i, fullmatch in enumerate(fullmatch_list):
            total_AA_length = 0
            index = "{}_&_full_match_{}".format(cplx, i)
            dfmatch.at[index,"cplx"] = cplx
            dfmatch.at[index,"match_number"] = 'full_match_{}'.format(i)
            dfmatch.at[index,"type_of_match"] = 'Full'
            dfmatch.at[index,"cplx_gene_stoich"] = str(dfparent.loc[cplx,'gene_stoich'])
            dfmatch.at[index,"num_pdbs"] = 1

            match_structure_info = {fullmatch : 1}
            dfmatch.at[index,"cplx_structure_stoich"] = str(match_structure_info )


            dfs_info = dfbest[dfbest.index.isin(match_structure_info)]
    #         if len(dfs_info) != len(match_structure_info):


            type_of_structures = dfs_info['stype'].to_dict()
            quality = dfs_info['total_quality'].to_dict()
            aa_covered = dfs_info['AA_length'].to_dict()
            gene_stoichiometry = dfs_info['gene_stoichiometry'].to_dict()
            dftemp = pd.DataFrame()
            dftemp['type_of_structures'] = pd.Series(type_of_structures)
            dftemp['quality'] = pd.Series(quality)   
            dftemp['aa_covered'] = pd.Series(aa_covered)
            dftemp['structure_stoich'] = pd.Series(match_structure_info)
            dftemp['gene_stoichiometry'] = pd.Series(gene_stoichiometry)
            total_AA = dftemp['aa_covered'] * dftemp['structure_stoich']
            total_AA= total_AA.sum()
            dfmatch.loc[index, 'total_AA'] = total_AA
            for i in dftemp.index:
                dftemp.at[i, 'total_AA'] = total_AA

            dftemp['ratio'] = dftemp['aa_covered'] * dftemp['structure_stoich']  / dftemp['total_AA'] 

            match_quality =  dftemp['ratio'] * dftemp['quality'] 

            for stype in dftemp.type_of_structures.unique():
                dfmatch.at[index,stype] = dftemp[dftemp.type_of_structures == stype].ratio.sum()

            dfmatch.at[index,"match_quality"] = match_quality.sum()
        checked+=[cplx]

#     print "Determining number of alphafold multimers in each match..."
    for index, row in tqdm(dfmatch.iterrows()):
        numAlphaMultimers = row.cplx_structure_stoich.count( 'AlphaMulti')
        dfmatch.loc[index,'numAlphaMultimers'] = numAlphaMultimers
    return dfmatch


def run_004A_rankMatches(dfmatch):
#     print '-------------------------------------------'
#     print "Ranking matches by quality, number of structures used, and number of AF multimer models used..."

    quality_rank_dict = {}
    num_pdb_rank_dict = {}
    numAlpha_rank_dict= {}
    for c in tqdm(dfmatch.cplx.unique()):
        dftemp = dfmatch[dfmatch.cplx == c]
        i = 1
        dftemp = dftemp.sort_values(by = ['match_quality'], ascending = [False])
        for value in dftemp.get("match_quality").unique():
            dfnum = dftemp[dftemp.get("match_quality") == value]
            for index, row in dfnum.iterrows():
    #             dfmatch.loc[index, 'quality_rank'] = i
                quality_rank_dict.update({index : i})
            i+=len(dfnum)




        i = 1
        dftemp = dftemp.sort_values(by = ['num_pdbs'], ascending = [True])
        for value in dftemp.get("num_pdbs").unique():
            dfnum = dftemp[dftemp.get("num_pdbs") == value]
            for index, row in dfnum.iterrows():
    #             dfmatch.loc[index, 'numpdbs_rank'] = i
                num_pdb_rank_dict.update({index : i})

            i+=len(dfnum)
            
        i = 1
        dftemp = dftemp.sort_values(by = ['numAlphaMultimers'], ascending = [False])
        for value in dftemp.get("numAlphaMultimers").unique():
            dfnum = dftemp[dftemp.get("numAlphaMultimers") == value]
            for index, row in dfnum.iterrows():
    #             dfmatch.loc[index, 'numpdbs_rank'] = i
                numAlpha_rank_dict.update({index : i})

            i+=len(dfnum)
            
            
    dfmatch[ 'quality_rank'] = pd.Series(quality_rank_dict)
    dfmatch['numpdbs_rank'] = pd.Series(num_pdb_rank_dict)
    dfmatch['numAlphaMulti_rank'] = pd.Series(numAlpha_rank_dict)
    def simplify_match(match_dict):
        new_dict = {}
        for k,v in match_dict.items():
            new_k = k.split('_&_')[0].split('_merged')[0]
            if new_k in new_dict:
                new_dict[new_k] += v
            else:
                new_dict.update({new_k : v})
        return new_dict

    def simplify_match_extra(match_dict):
        new_dict = {}
        for k,v in match_dict.items():
            new_k = k.split('_&_')[0].split('_merged')[0].split('_bio')[0]
            if new_k in new_dict:
                new_dict[new_k] += v
            else:
                new_dict.update({new_k : v})
  
        return new_dict


    match_simple_dict ={}
    match_simple_extra_dict = {}
    for index, row in tqdm(dfmatch.iterrows()):
        match = ast.literal_eval(row.cplx_structure_stoich)

        match_simple = simplify_match(match)
        match_simple_dict.update({index : match_simple})
    #     dfmatch.loc[index, 'simple_match'] = str(match_simple)
        match_simple = simplify_match_extra(match)
        match_simple_extra_dict.update({index : match_simple})

    #     dfmatch.loc[index, 'simple_match_extra'] = str(match_simple)

    dfmatch['simple_match'] = pd.Series(match_simple_dict)
    dfmatch['simple_match_extra'] = pd.Series(match_simple_extra_dict)
    return dfmatch

def sort_alpha(dfc):
    dfc= dfc.sort_values(by = ['numAlphaMulti_rank','quality_rank','numpdbs_rank'], ascending=['False','False','False'])
    return dfc

def sort_quality(dfc):
    dfc= dfc.sort_values(by = ['quality_rank','numAlphaMulti_rank','numpdbs_rank'], ascending=['False','False','False'])
    return dfc
def sort_few(dfc):
    dfc= dfc.sort_values(by = ['numpdbs_rank', 'numAlphaMulti_rank','quality_rank',], ascending=['False','False','False'])
    return dfc


def run_004A_chooseBestMatch(dfmatch_ranked):
#     print '-------------------------------------------'
#     print "Choosing Best Matches..."

    use_cplx = {}
    for cplx in tqdm(dfmatch_ranked.cplx.unique()):

        dfc = dfmatch_ranked[dfmatch_ranked.cplx == cplx]
        if [0] !=  dfc.numAlphaMultimers.unique().tolist():
            dfc = sort_alpha(dfc)
            use = dfc.first_valid_index()
            use_cplx.update({cplx :use})
#             if cplx == 'CPLX0-3802':
#                     print '1'
            ##### Use the match with the most alphafold multimers and highest quality
            continue

        else:

            #get the highest quality matches first
            dfc_qual = sort_quality(dfc)
            quality_match = dfc_qual.first_valid_index()
            quality_quality = dfc_qual.match_quality.values[0]
            quality_num_pdbs = dfc_qual.num_pdbs.values[0]

            #get the fewest structures matches first
            dfc_few = sort_few(dfc)
            few_match = dfc_few.first_valid_index()

            if few_match == quality_match:
                use = few_match
                use_cplx.update({cplx :use})
#                 if cplx == 'CPLX0-3802':
#                     print '2'
                #if the fewest structures is the same as the best quality, we've found out match
                continue

            index_keep_dfc = []
            for index, row in dfc_few.iterrows():
                few_quality = row.match_quality
                few_num_pdbs = row.num_pdbs
                if few_num_pdbs-quality_num_pdbs >=0:
                    #if the number of structures is the same as the best quality, but the quality is worse, DONT use
                    continue

                if quality_quality - few_quality>=70:
                    #if the quality loss is more than 70%, do not use
                    continue

                index_keep_dfc += [index]
                dfc_few.loc[index,'num_diff'] = few_num_pdbs-quality_num_pdbs
                dfc_few.loc[index,'qual_diff'] = few_quality-quality_quality

            dfc_few = dfc_few[dfc_few.index.isin(index_keep_dfc)]
            if len(dfc_few) == 0:
                use = quality_match
                use_cplx.update({cplx :use})
#                 if cplx == 'CPLX0-3802':
#                     print '3'
                # if no matches with fewer structures meet criteria, use the best quality match
                #### USE BEST MATCH###################################
                continue

            index_keep_dfc = []
            for num_diff in dfc_few.num_diff.unique():
                t = dfc_few[dfc_few.num_diff == num_diff]
                t = t.sort_values(by = 'match_quality', ascending = False)
    #             print t.match_quality.values
                index_keep_dfc += [t.first_valid_index()]
                #if multiple matches have same number of stuctures loss, only use the one with the best quality

            dfc_few = dfc_few[dfc_few.index.isin(index_keep_dfc)]
    #         print '\n'


            diff_max_loss_quality = dfc_few.loc[dfc_few.first_valid_index(), 'qual_diff']
            diff_max_loss_num = dfc_few.loc[dfc_few.first_valid_index(), 'num_pdbs']
            if diff_max_loss_quality > -25 or diff_max_loss_num == 1:
                use = dfc_few.first_valid_index()
                use_cplx.update({cplx :use})

#                 if cplx == 'CPLX0-3802':
#                     print '4','\t', diff_max_loss_quality,diff_max_loss_num
#                 print '\t', "{:.2f}\t{}\t{}".format(diff_max_loss_quality, diff_max_loss_num, cplx,use)
                #USE FEW MATCH###################################
                continue



            index_remove_dfc = []
            for index, row in dfc_few.iterrows():
                not_needed = '<>'
                if row.qual_diff < diff_max_loss_quality:
                    not_needed = 'X'
                    index_remove_dfc += [index]

    #         for index, row in dfc_few.iterrows():
    #         for index, row in dfc_few.iterrows():


    #             x = few_quality-quality_quality
    #             y = few_num_pdbs-quality_num_pdbs
                if not_needed =='X':
                    continue
    #             print "{} \t {:.1f}\t{:.1f}\t{:.1f} \t\t {}\t{}\t{} \t {}".format(not_needed,quality_quality , row.match_quality, row.qual_diff, quality_num_pdbs, row.num_pdbs, row.num_diff, cplx)

            dfc_few = dfc_few[dfc_few.index.isin(index_remove_dfc) == False]

            use = dfc_few.first_valid_index()
            use_cplx.update({cplx :use})
#             if cplx == 'CPLX0-3802':
#                     print '5'
    dfqual = pd.DataFrame()
    for cplx, match in tqdm(use_cplx.items()):
        dfqual.loc[cplx, 'final_match_id'] = match
        dfqual.loc[cplx, 'final_qual'] = dfmatch_ranked.loc[match,'match_quality']
        dfqual.loc[cplx, 'final_match'] = dfmatch_ranked.loc[match,'cplx_structure_stoich']
        dfqual.loc[cplx, 'num_pdbs'] = dfmatch_ranked.loc[match,'num_pdbs']
    #     dfqual.loc[cplx, 'simple_match'] = dfmatch.loc[match,'simple_match']
        for stype in ['PDB','ALPHAFOLD','ALPHAFOLD_MULTIMER','SWISS','ITASSER']:
            dfqual.loc[cplx, stype] = dfmatch_ranked.loc[match,stype]
    outfile = op.join(qspaceDirs['DataOutput_dir'], '004A-BFS_best_match_for_protein_complexes.csv')
    dfqual.to_csv(outfile)
#     print "Saving....\n\t> {}".format(outfile)
    return dfqual
