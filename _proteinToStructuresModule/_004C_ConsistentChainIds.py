# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def get_len_chains_w_same_letter(list_of_chain_ids):
    first_letter = {}
    for chain in list_of_chain_ids:
        if chain[0] not in first_letter:
            first_letter.update({chain[0] : 1})
        else:
            first_letter[chain[0]] += 1
    
    return first_letter

def get_locationOfAllStructures(dfalldata):   
    structure_dict = {}

    for structure in dfalldata.structureId.unique():
        if 'assembly' in structure:
            stype = 'PDB'
        elif 'AlphaMulti' in structure:
            stype = 'ALPHAFOLD_MULTIMER'
        elif 'AF-' in structure and 'model' in structure:
            stype= 'ALPHAFOLD'
        elif 'ECOLI' in structure and 'clean_residues' in structure:
            stype= 'ITASSER'
        else:
            stype= 'SWISS'

        structure_dict.update({structure : stype})

    dfstructure_file_locations = pd.DataFrame()
    dfstructure_file_locations['stype'] = pd.Series(structure_dict)
    for index, row in dfstructure_file_locations.iterrows():
        sfile = find_sfile(index, row.stype)
        dfstructure_file_locations.loc[index, 'sfile'] =  sfile

    sfileDict = dfstructure_file_locations.sfile.to_dict()
    return sfileDict

def UseReference_chains(dfallchains):
    
    print ("Finding exact matches to structure file chains...\n-------------------------------------")
    use_reference_dict = {}
    for index, row in dfallchains.iterrows():
        sfileChainsIds = row.sfileChains
        if type(sfileChainsIds) == str:
            sfileChainsIds = ast.literal_eval(sfileChainsIds)

        chainType_dict = {}
        for chainType, chainIds in row.drop(['sfileChains','use_reference','mismatch_reference','chain_translation']).items():
            if type(chainIds) == str:
                chainType_dict.update({chainType : ast.literal_eval(chainIds)})
            else:
                chainType_dict.update({chainType : chainIds})

        for chainType, chainIds in chainType_dict.items():
            if chainIds == sfileChainsIds:
                use_reference_dict.update({index : chainType})
                break
    dfallchains['use_reference'] = pd.Series(use_reference_dict)
    return dfallchains



def find_mismatches(dfallchains,PDBchainSeqsAPI, dfalldata):
    print ("Finding chain mismatches...\n-------------------------------------")
    
    ChainInDataFrame_to_sfileChain = {
#         '2scu-assembly1': {'E-0': 'B-2', 'A-0': 'A', 'B-0': B'', 'D-0': 'A-2'},
#                                       "6ag8-assembly1" : {"A-0":"C-2", "C-0":"C"},
                                     }
    
    fill_in_use_reference = {
#         '2scu-assembly1' : 'structChainId',
#                              "6ag8-assembly1" : 'structChainId'
                            }
    err= []
    for index, row in dfallchains[dfallchains.use_reference.isna()].iterrows():
        if index in fill_in_use_reference:
            continue
        sfileChainsIds = row.sfileChains
        if type(sfileChainsIds) == str:
            sfileChainsIds = ast.literal_eval(sfileChainsIds)

        sfileChainsIds = row.sfileChains
        if type(sfileChainsIds) == str:
            sfileChainsIds = ast.literal_eval(sfileChainsIds)

        chainType_dict = {}
        for chainType, chainIds in row.drop(['sfileChains','use_reference','mismatch_reference','chain_translation']).items():
            if type(chainIds) == str:
                chainType_dict.update({chainType : ast.literal_eval(chainIds)})
            else:
                chainType_dict.update({chainType : chainIds})

        compact_sfileChains = get_len_chains_w_same_letter(row.sfileChains)

        keepChainType = []
        for chainType, chainIds in chainType_dict.items():
            compact_dfChains = get_len_chains_w_same_letter(chainIds)
            if compact_dfChains == compact_sfileChains:
                keepChainType += [chainType]
        if len(keepChainType) > 0:

            use_reference_id = keepChainType[0]
            log.info("Using {} chains {} | Mismatched w/same chainId[0]".format(use_reference_id,index))
            fill_in_use_reference.update({index : use_reference_id})
            chainIds = chainType_dict[use_reference_id] 

            
            
            
            for i, chainId_in_df in enumerate(chainIds):
                if index not in ChainInDataFrame_to_sfileChain:
                    ChainInDataFrame_to_sfileChain.update({index : {chainId_in_df : list(sfileChainsIds)[i]} })
                else:
                    ChainInDataFrame_to_sfileChain[index].update( {chainId_in_df : list(sfileChainsIds)[i]})
        else:
            found_index = False
            if "assembly" in index:
                pdb_entry = index.split('-assembly')[0]

                dfpdb_seq = PDBchainSeqsAPI[PDBchainSeqsAPI.pdb_entry ==str(pdb_entry).upper()]
                if len(dfpdb_seq) == 1 and len(chainType_dict['structChainId']) == len(sfileChainsIds): 
                    found_index = True
                    log.info("Using {} chains {} | Homomeric --> mismatch not important".format('structChainId',index))
                    fill_in_use_reference.update({index : 'structChainId'})

                    #homomeric and same number of mapped chains in DF and sfile
#                     for i , chainId_in_df in enumerate(chainType_dict['structChainId']):
    #                     print (index, i)
#                         if index not in ChainInDataFrame_to_sfileChain:
#                             ChainInDataFrame_to_sfileChain.update({index : {chainId_in_df : list(sfileChainsIds)[i]} })
#                         else:
#                             ChainInDataFrame_to_sfileChain[index].update( {chainId_in_df : list(sfileChainsIds)[i]})
                    continue
            if found_index:
                continue

            found_subset = False
            for chainType, chainIds in chainType_dict.items():
                if type(chainIds) == str:
                    chainIds = ast.literal_eval(chainIds)
                if chainIds.intersection(sfileChainsIds) == chainIds:
                    #the structure file has addtional chains (i.e. antibody, etc. just ignore it)
                    found_subset = True
                    found_index = True
                    fill_in_use_reference.update({index : chainType})
                    log.info(" {} addtional chains (i.e. antibody) | use subset {}".format(index, chainType))

    #                 print ('found index : {}'.format(index))
#                     for chainId_in_df in chainIds:
#                         if index not in ChainInDataFrame_to_sfileChain:
#                             ChainInDataFrame_to_sfileChain.update({index : {chainId_in_df : chainId_in_df} })
#                         else:
#                             ChainInDataFrame_to_sfileChain[index].update( {chainId_in_df : chainId_in_df})
                    break
            if found_subset:
                continue


            log.warn('>> chains mismatched, using structChainId : {}'.format(index))
            err += [index]
            fill_in_use_reference.update({index : 'structChainId'})

    dfallchains['mismatch_reference'] = pd.Series(fill_in_use_reference)    
    return dfallchains#, ChainInDataFrame_to_sfileChain



from Bio.PDB.Polypeptide import three_to_one

def map_dfChain_to_sfileChain(dfallchains,dfallchains_w_seq, dfQuatProteome):
    print ("Resolving chains in structures vs chains in QSPACE dataframe...\n-------------------------------------")

    translation_glob = {}
    for index, row in tqdm (dfallchains.iterrows()):
        dfs = dfQuatProteome[dfQuatProteome.structureId == index]
        dfs = dfs[dfs.structNum.isna() == False]
        dfs_allchains_w_seq = dfallchains_w_seq[dfallchains_w_seq.structureId == index]
        dfs_allchains_w_seq = dfs_allchains_w_seq.sort_values(by = 'chainId')
        translation = {}
        for cplx in dfs.cplx.unique():
            dfc = dfs[dfs.cplx == cplx]

            if 'mismatch_reference' in row.dropna().keys():
                use_reference = row.get('mismatch_reference')

            elif 'use_reference' in row.dropna().keys():
                use_reference = row.get('use_reference')

            else:
                log.warn('Unknown Reference Chains to use : {}'.format(index))
                continue

            reference_chains = row.get(use_reference)
            if type(reference_chains) == str:
                reference_chains= ast.literal_eval(reference_chains)
            used_chains = []
            reference_chains = list(reference_chains)
            reference_chains.sort()

            for ref_chain in reference_chains:
                dfc_chain= dfc[dfc.get(use_reference) == ref_chain]
                dfc_chain_AAseq = "".join(dfc_chain.structAA1.tolist())

                same_seq = dfs_allchains_w_seq[dfs_allchains_w_seq.chainAAseq == dfc_chain_AAseq]
                if len(same_seq) == 0:
                    log.warn("No Structure AA seq matches QSPACE proteome {}-{}".format(ref_chain, index))
                    continue
                for same_chain in same_seq.chainId.tolist():
                    if same_chain in used_chains:
                        continue

                    translation.update({ref_chain : same_chain})
                    used_chains += [same_chain]
                    break


        translation_glob.update({index : translation})
    
    return translation_glob

def run_004C_get_all_chains(dfalldata,PDBchainSeqsAPI, force_rerun = False):
    print ("Getting chains in structures vs chains in QSPACE dataframe...\n-------------------------------------")

    outfile = op.join(qspaceDirs['DataOutput_dir'], '004C-structureFile_chains_mapped.csv')
    outfile2 = op.join(qspaceDirs['DataOutput_dir'], '004C-structureFile_chain_seqs.csv')
    
    if op.exists(outfile) and op.exists(outfile2) and not force_rerun:
        dfallchains = pd.read_csv(outfile,index_col= 0)
        dfallchains_w_seq = pd.read_csv(outfile2,index_col= 0)
    
    else:
        
        dfallchains = pd.DataFrame(columns = ['sfileChains', 'structChain', 'structChainId',
                                              'structChainId_mod','use_reference', 'mismatch_reference','chain_translation']
                                  )

        
        sfileChains = {}
        structChains = {}
        structChainIds = {}
        structChainId_mods = {}

        sfileDict = get_locationOfAllStructures(dfalldata)
        index = 0 
        
        glob_chain_aaseq = {}
        glob_chain_chain_id = {}
        glob_chain_struct_id = {}
        glob_chain_aaseq_len = {}
        for structureId, sfile in tqdm(sfileDict.items()):

            sfile = StructProp(ident=structureId, structure_path=sfile , file_type = op.basename(sfile).split(".")[-1])
            sfile = sfile.parse_structure()
            real_structure_chains = []
            # structure_chains = sfile.first_model.child_dict

            for sChainId, sChain in sfile.first_model.child_dict.items():
                realchainlen = 0 
                chain_aa_seq = ''
                for residueId, residue in sChain.child_dict.items():
                    
                    try:
                        chain_aa_seq += three_to_one(residue.resname)
                    except KeyError:
#                         log.info('{}-{}-{}-{} Not an AA'.format(structureId, sChainId, residueId[1], residue.resname))
                        continue
                    if residue.resname in amino_acid_ids_three:
                        realchainlen +=1
                if realchainlen > 5:
                    real_structure_chains += [sChainId]
                    
                    
                glob_chain_aaseq.update({index : chain_aa_seq})
                glob_chain_chain_id.update({index : sChainId})
                glob_chain_struct_id.update({index : structureId})
                glob_chain_aaseq_len.update({index : len(chain_aa_seq)})
                
                index +=1
            real_structure_chains = set(real_structure_chains)
            
            

            sfileChains.update({structureId : real_structure_chains})

            dfs = dfalldata[dfalldata.structureId == structureId]

            structChains.update({structureId : set(dfs.structChain.unique().tolist())})
            structChainIds.update({structureId : set(dfs.structChainId.unique().tolist())})
            structChainId_mods.update({structureId : set(dfs.structChainId_mod.unique().tolist())})
        #     break

        dfallchains['sfileChains'] = pd.Series(sfileChains)
        dfallchains['structChain'] = pd.Series(structChains)
        dfallchains['structChainId'] = pd.Series(structChainIds)
        dfallchains['structChainId_mod'] = pd.Series(structChainId_mods)
        dfallchains.to_csv(outfile)
        log.info("Saving dfChain-sfileChain map ...\n\t{}".format(outfile))

        
        dfallchains_w_seq = pd.DataFrame()
        dfallchains_w_seq['structureId'] = pd.Series(glob_chain_struct_id)
        dfallchains_w_seq['chainId'] = pd.Series(glob_chain_chain_id)
        dfallchains_w_seq['chainAAseq'] = pd.Series(glob_chain_aaseq)
        dfallchains_w_seq['polymer_entity_seq_len'] = pd.Series(glob_chain_aaseq_len)
         #fixing the chain mismatches
        dfallchains_w_seq.to_csv(outfile2)
        log.info("Saving sfileChain AA seqs ...\n\t{}".format(outfile2))
    
   
    dfallchains = UseReference_chains(dfallchains) 
    dfallchains  = find_mismatches(dfallchains,PDBchainSeqsAPI,dfalldata=dfalldata)
    
    translation_dict = map_dfChain_to_sfileChain(dfallchains=dfallchains,
                                                 dfallchains_w_seq=dfallchains_w_seq,
                                                 dfQuatProteome=dfalldata)
    
    dfallchains['chain_translation'] = pd.Series(translation_dict)
    dfallchains.to_csv(outfile)
    log.info("Saving dfChain-sfileChain map ...\n\t{}".format(outfile))

    return dfallchains, dfallchains_w_seq,translation_dict


def run_004C_global_chain_consistency(dfalldata, dfallchains):
    
    
#     dfalldata = dfalldata[dfalldata.structNum.isna() == False]
    global_sfile_chains = {}
    print ("\nMapping structureFileChainIds to QSPACE proteome.....\n---------------------------------")

    for index in tqdm(dfallchains.index.tolist()):
        row = dfallchains.loc[index]
        dfs = dfalldata[dfalldata.structureId == index]
        
        if 'use_reference' in row.dropna().keys():        
            sfileChainStandard = row.get('use_reference')
            global_sfile_chains.update(dfs.get(sfileChainStandard).to_dict())
        
        elif 'mismatch_reference' in row.dropna().keys(): 
            sfileChainStandard = row.get('mismatch_reference')
            
            df_to_struct_chain_map = row.get('chain_translation')
            if type(df_to_struct_chain_map) == str:
                df_to_struct_chain_map = ast.literal(df_to_struct_chain_map)
            
            for df_chain, sfileChainId in df_to_struct_chain_map.items():
                dfs_c = dfs[dfs.get(sfileChainStandard) == df_chain]
                for key, chain in dfs_c.get(sfileChainStandard).to_dict().items():
                    global_sfile_chains.update({key : sfileChainId})
    
    dfalldata['sfileChainId'] = pd.Series(global_sfile_chains)
    dfalldata_struct = dfalldata[dfalldata.structNum.isna() == False]
    resNums = dfalldata_struct.structNum.to_dict()
    sfileChains  = dfalldata_struct.sfileChainId.to_dict()
    
    sfileChain_residue = {}
    print ("Getting {structureFileChainId}_{Residue#}.....\n---------------------------------")
    for index, resNum in tqdm(resNums.items()):
        sfileChain_residue.update({index : "{}_{}".format(sfileChains[index], int(resNum))})
    dfalldata['sfileChain_Residue'] = pd.Series(sfileChain_residue)
    
    outfile = op.join(qspaceDirs['DataOutput_dir'] , '004C-alldata_skeleton.csv')
    dfalldata.to_csv(outfile)
    log.info("Saving QSPACE proteome with chains fixed ...\n\t{}".format(outfile))
    
    return dfalldata
            
            