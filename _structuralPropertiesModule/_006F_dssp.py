# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def get_dsspIndex_to_StructIndex(dssp_df,df_structure):

    dsspIndex_to_StructIndex = {}
    
    unique_chains= {}
    for chain in df_structure.sfileChainId.unique():
        chain1letter = chain[0]
        if chain1letter not in unique_chains:
            unique_chains.update({chain1letter : [chain]})
        else:
            unique_chains[chain1letter] += [chain]
    
    print (unique_chains)
    for dssp_chain, df_struct_chain_list in unique_chains.items():
        
    
        i = 0
        dfstruct_c = df_structure[df_structure.sfileChainId.isin(df_struct_chain_list)]
        dssp_df_c = dssp_df[dssp_df.dssp_chain == dssp_chain]
#         dssp_df_c = dssp_df_c[dssp_df_c.dssp_chain == dssp_chain]
        print ( dssp_chain, len(dfstruct_c),len(dssp_df_c))
        structure_index_list = dfstruct_c.index.tolist()
        
        for index , row in dssp_df_c.iterrows():
            dssp_aa = row.dssp_aa
            if dssp_aa == "X":
                continue
            if i >= len(structure_index_list):
                break
            structIndex = structure_index_list[i]
            structAA1 = dfstruct_c.loc[structIndex,'structAA1']
            while structAA1 != dssp_aa:
                i+=1
                if i >= len(structure_index_list):
                    break
                structIndex = structure_index_list[i]
                structAA1 = dfstruct_c.loc[structIndex,'structAA1']
                print (i,end =' ')
                
            dsspIndex_to_StructIndex.update({index:structIndex})
            i +=1
            
            
    return dsspIndex_to_StructIndex

def get_dssp_df(sfile,outdir, force_rerun = True):
        
    structureId = op.basename(sfile).split('.')[0]
    outfile = op.join(outdir, "{}_dssp.df".format(structureId))
    
    if op.exists(outfile) and not force_rerun:
        dssp_df = pd.read_csv(outfile, index_col = 0)
        return dssp_df,outfile
    
    if '.cif' in sfile:
        p = MMCIFParser()
    elif '.pdb' in sfile:
        p = PDBParser()
    structure = p.get_structure(structureId, sfile)
    model = structure[0]
    dssp = DSSP(model, sfile)

    appender = []
    for a_key in dssp.keys():
        a_value = dssp[a_key]
        info = [a_key[0]] + list(a_value) 
        appender.append(info)
    cols = ['dssp_chain','dssp_index','dssp_aa','dssp_ss','dssp_exposure_rsa','dssp_phi','dssp_psi','dssp_NH_O_1_relidx',
            'dssp_NH_O_1_energy','dssp_O_NH_1_relidx','dssp_O_NH_1_energy','dssp_NH_O_2_relidx',
            'dssp_NH_O_2_energy','dssp_O_NH_2_relidx','dssp_O_NH_2_energy']
    dssp_df = pd.DataFrame.from_records(appender, columns= cols )
    dssp_df.to_csv(outfile)
    
    return dssp_df,outfile

def parse_DSSP(dssp_file, dfstructure):
    dssp_properties = ['dssp_chain','dssp_index','dssp_aa','dssp_ss','dssp_exposure_rsa','dssp_phi',
                       'dssp_psi','dssp_NH_O_1_relidx','dssp_NH_O_1_energy','dssp_O_NH_1_relidx',
                       'dssp_O_NH_1_energy','dssp_NH_O_2_relidx','dssp_NH_O_2_energy','dssp_O_NH_2_relidx',
                       'dssp_O_NH_2_energy']
 
        
    dssp_df = pd.read_csv(dssp_file,index_col=0)
    
    dsspIndex_to_StructIndex =get_dsspIndex_to_StructIndex(dssp_df,dfstructure)
    
    
    dssp_dict = {}
    for prop in dssp_properties:
        prop_data =  dssp_df[dssp_df.index.isin(dsspIndex_to_StructIndex)][prop]
        
        
        temp = {}
#         print (prop)
        if prop  =='dssp_index' or 'relidx' in prop:
            prop_index_info = map(lambda index, residue_property : {int(dsspIndex_to_StructIndex[index]) :int(residue_property)}, prop_data.index, prop_data.values)
        else:
            prop_index_info = map(lambda index, residue_property : {int(dsspIndex_to_StructIndex[index]) :residue_property}, prop_data.index, prop_data.values)
        for index in prop_index_info:
            temp.update(index)
        dssp_dict.update({prop : temp})

    return dssp_dict 

def run_006F_runDSSP(dfalldata,
                     force_rerun = True,
                     dsspResultsFolder = qspaceDirs['dsspResultsDir'],
                    ):

    print ("Getting structure file locations...\n------------------------------------")
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
        
        
    print ("Running DSSP ...\n------------------------------------")
    err = 0
    err_list=[]
    for structure_id in tqdm(dfalldata.structureId.unique()):
        file_location = dfstructure_file_locations.loc[structure_id,'sfile']
#         print (file_location)
        try:
            dssp_df, dssp_outfile = get_dssp_df(sfile=file_location,outdir=dsspResultsFolder, force_rerun = force_rerun)
          
        except Exception:
            err += 1
#             print ('error', err, structure_id,'\t', file_location)
            log.warn("Error with DSSP : {}".format(file_location))
            err_list +=[structure_id]
      
          
#     dssp_properties = ['dssp_chain','dssp_index','dssp_aa','dssp_ss','dssp_exposure_rsa','dssp_phi','dssp_psi','dssp_NH_O_1_relidx',
#             'dssp_NH_O_1_energy','dssp_O_NH_1_relidx','dssp_O_NH_1_energy','dssp_NH_O_2_relidx',
#             'dssp_NH_O_2_energy','dssp_O_NH_2_relidx','dssp_O_NH_2_energy']

    dssp_global_properties_dict = {}
#     for k in dssp_properties:
#         dssp_global_properties_dict.update({k : {}})
        
    dssp_errorlist = []

    print ("Mapping DSSP to Global Indicies ...\n------------------------------------")
    cplx_list = dfalldata.cplx.unique()
    for cplx in tqdm(cplx_list):
        df_cplx = dfalldata[dfalldata.cplx ==cplx]

        for structure_id in df_cplx.structureId.unique():
            if structure_id in err_list:
                continue
            df_cplx_struct = df_cplx[df_cplx.structureId ==structure_id]

            df_structure = df_cplx_struct[df_cplx_struct.structNum.isna() == False]
            df_structure = df_structure[df_structure.structNum != np.inf]
            df_structure = df_structure[df_structure.sfileChainId.isna() == False]
            df_structure = df_structure[df_structure.sfileChain_Residue.isna() == False]
            
            if len(df_structure) == 0:
                log.warn("No viable protein structure coordinates in QSPACE : {}".format(structure_id))
            file_location = dfstructure_file_locations.loc[structure_id,'sfile']
            basename = op.basename(file_location).split('.')[0]

            dssp_file =op.join(dsspResultsFolder, '{}_dssp.df'.format(basename))
            if not op.exists(dssp_file):
                log.warn("No DSSP file : {}".format(structure_id))
                dssp_errorlist.append(structure_id)
                continue
            print (structure_id)
#             print (file_location)
#             print (dssp_file,'\n')
            
            info = parse_DSSP(dssp_file, df_structure)
            

            for dssp_prop, dssp_dict in info.items():
                if dssp_prop not in dssp_global_properties_dict:
                    dssp_global_properties_dict.update({dssp_prop : dssp_dict})
                else:
                    dssp_global_properties_dict[dssp_prop].update(dssp_dict)

    outfile = op.join(qspaceDirs['DataOutput_dir'],'006F-StructuralProteome_DSSP.json')    
    with open(outfile, 'w') as f:
        json.dump(dssp_global_properties_dict,f )
    log.info("Saving QSPACE Proteome - DSSP ...\n\t{}".format(outfile))
    return dssp_global_properties_dict            
        