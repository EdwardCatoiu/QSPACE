# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)



def get_df_coords(sfile, aa3 = amino_acid_ids_three ):
    dfp_interface = pd.DataFrame()
    s = StructureIO(sfile).structure[0]

    atom_ca_coord_dict = {}
    x_dict = {}
    y_dict = {}
    z_dict = {}
    chain_dict = {}
    residue_num_dict = {}
    for chain in s.get_chains():
        for residue in chain.get_residues():
            if residue.resname not in aa3:
                continue
            try:
                [x,y,z] = atom_ca_coord = residue['CA'].coord
            except KeyError:
                continue
            residue_id = '{}_{}'.format(chain.id,residue.id[1])
            atom_ca_coord_dict.update({residue_id : atom_ca_coord})
            x_dict.update({residue_id : x})
            y_dict.update({residue_id : y})
            z_dict.update({residue_id : z})
            chain_dict.update({residue_id : chain.id})
            residue_num_dict.update({residue_id : residue.id[1]})
    dfp_interface['chainId'] = pd.Series(chain_dict)
    dfp_interface['resNum'] = pd.Series(residue_num_dict)
    dfp_interface['atom_ca_coord'] = pd.Series(atom_ca_coord_dict)
    dfp_interface=dfp_interface.sort_values(by=['chainId','resNum'])
    return dfp_interface

def return_interfaces(dfp, threshold_list=[10,7.5,5,3]):

    interface_dict = {}

    all_interface_residues= dfp.index.tolist()
    for interface_threshold in (threshold_list):
        interface_dict.update({interface_threshold : {}})
        checked_chains = []
        checked_crs = []
        for chain in dfp.get('chainId').unique():
            dfchain = dfp[dfp.get('chainId') == chain]
            dfchain = dfchain[dfchain.index.isin(checked_crs) == False]
            #this is the non-redundant step of the process
            #only keep indicies that are within the threshold angstroms at each step
            dfchain = dfchain[dfchain.index.isin(all_interface_residues) == True]
#             print "Only considering {} Residues for {}A in chain {}".format(len(set(all_interface_residues)),interface_threshold,chain)


            dfchain_others = dfp[dfp.get('chainId') != chain]
            dfchain_others = dfchain_others[dfchain_others.get('chainId').isin(checked_chains) == False]
            #this is the non-redundant step of the process
            #only keep indicies that are within the threshold angstroms at each step
            dfchain_others = dfchain_others[dfchain_others.index.isin(all_interface_residues) == True]

            for index, coordinate in (dfchain['atom_ca_coord'].to_dict().items()):
                if '__' in index:
                    continue
                [c1,r1] = index.split('_')
                

                for other_index, other_coordinate in dfchain_others['atom_ca_coord'].to_dict().items():
                    if '__' in other_index:
                        continue
                    
                    [c2,r2] = other_index.split('_')

                    if protein_geometry.distance_two_points(coordinate, other_coordinate) < interface_threshold:
                        checked_crs+= ["{}_{}".format(c1, r1), "{}_{}".format(c2,r2)]
                        if interface_threshold not in interface_dict:
    #                         interface_dict.update({interface_threshold : ["{}_{}".format(c1, r1), "{}_{}".format(c2,r2)]})
                            interface_dict.update({interface_threshold :{"{}_{}".format(c1, r1) : 1}})
                            interface_dict.update({interface_threshold :{"{}_{}".format(c2, r2) : 1}})
                        else:
                            interface_dict[interface_threshold].update({"{}_{}".format(c1, r1) : 1})
                            interface_dict[interface_threshold].update({"{}_{}".format(c2, r2) : 1})
    #                         += ["{}_{}".format(c1, r1), "{}_{}".format(c2,r2)]


            checked_chains.append(chain)
        all_interface_residues = list(set(interface_dict[interface_threshold].keys()))
    return interface_dict


def run_006E_geometricInterfaces(dfalldata, outfolder= qspaceDirs['proteinInterfaceDir'],force_rerun= False):
#     print ("Getting structure file locations...\n------------------------------------")
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

        
    print ("Getting geometric protein-protein interfaces ...\n------------------------------------")
    for p in tqdm(dfalldata.structureId.unique()):
        sfile = dfstructure_file_locations.loc[p, 'sfile']
        outfile = op.join(outfolder, '{}_interfaces.csv'.format(p))

        if op.exists(outfile) and not force_rerun:
            continue

        dfp = get_df_coords(sfile, aa3 = amino_acid_ids_three )
        print (len(dfp), p)
        threshold_list = [10,7.5,5,3]
        interface_dict = return_interfaces(dfp,threshold_list = threshold_list)
        for threshold in threshold_list:
            interface = interface_dict[threshold]
            dfp['Interface_{}_Angstroms'.format(threshold)] = pd.Series(interface)
        dfp.to_csv(outfile) 
        
def run_006E_mapGeoInterfaceToProteome(dfalldata,interfaceFolder = qspaceDirs['proteinInterfaceDir']):
    global_interfaces = {3 :{},5:{},7.5:{},10:{}}
    print ("Mapping geometric protein-protein interfaces to proteome ...\n------------------------------------")

    for f in tqdm(os.listdir(interfaceFolder)):
        structureId = f.split('_interfaces')[0]
        if structureId not in dfalldata.structureId.unique():
            continue
        
        infile = op.join(interfaceFolder, f)
        dfinterface = pd.read_csv(infile, index_col=0)

        dfs = dfalldata[dfalldata.structureId == structureId]
        global_index_dict_structure = dfs.sfileChain_Residue.to_dict()

        for threshold in [3,5,7.5,10]:
            columnKey = "Interface_{}_Angstroms".format(threshold)
            interface_crs = dfinterface[columnKey].dropna().index.tolist()

            for index, cr in global_index_dict_structure.items():
                if cr not in interface_crs:
                    continue

                if index not in global_interfaces[threshold]:
                    global_interfaces[threshold].update({index : 1})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006E-StructuralProteome-interfaces_geometry.json')    
    with open(outfile, 'w') as f:
        json.dump(global_interfaces,f )
    log.info("Saving QSPACE Proteome - protein interfaces (geometric)...\n\t{}".format(outfile))
    return global_interfaces

def generate_ScanNet_command(dfalldata,
                             path_to_ScanNet,
                             commandFile ,  
                             force_rerun = False,
                             scannetResultsFolder = qspaceDirs['scannetResultsDir'],
                             outPDBfolder = qspaceDirs['scannetStructuresDir'],
                             
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
    
    for index, row in dfstructure_file_locations.iterrows():
        structureFileId = op.basename(row.sfile).split('.')[0]
        dfstructure_file_locations.loc[index, 'structureFileId'] = structureFileId
    
    print ("Writing ScanNet Command...\n------------------------------------")

    cmd = ''
    c = 0
    for sid in tqdm(dfalldata.structureId.unique()):
        sfile = dfstructure_file_locations.loc[sid, 'sfile']

        basename = op.basename(sfile).split('.')[0]

        results_file = "{}_single_ScanNet_interface_noMSA/predictions_{}.csv".format(basename,basename)

        outfileResults = op.join(scannetResultsFolder, results_file)
        if op.exists(outfileResults) and not force_rerun:
            continue
        print ("{}; ".format(sid),)
        s = StructProp(ident=sid, structure_path=sfile, file_type = op.basename(sfile).split('.')[-1])
        parsed = s.parse_structure()
        #     break
        outfilePDB = op.join(outPDBfolder, sid + '.pdb') 
        outfilePDB = write_pdb(s = s, p= parsed, outfile = outfilePDB)

        if c %100 == 0:
            cmd +='\n\n'

        cmd += '{}python predict_bindingsites.py {} --noMSA ; '.format(path_to_ScanNet, outfilePDB)
        c +=1

    with open(commandFile, 'w') as f2:
        f2.write("#!/bin/bash\n")
        f2.write(cmd)
    return commandFile, dfstructure_file_locations

def run_006E_readScannetResults(dfalldata,
                                dfstructure_file_locations,
                                scannetFolder = qspaceDirs['scannetResultsDir'],
                               ):
    scannetResults = {}

    for folder in tqdm(os.listdir(scannetFolder)):
        structure_id = folder.split('_single_ScanNet')[0]
        if structure_id not in dfalldata.structureId.unique():
            continue
            
        
        print (folder)

        prediction_file = op.join(scannetFolder,folder,"predictions_{}.csv".format(structure_id))
        dfscannet = pd.read_csv(prediction_file)

        scannetChain = dfscannet['Chain'].to_dict()
        scannetResidue = dfscannet.get('Residue Index').to_dict()
        scannetBindingSite = dfscannet.get('Binding site probability').to_dict()


        scannetChainResidueBinding = {}
        for index, chain in scannetChain.items():
            residue = scannetResidue[index]
            binding = scannetBindingSite[index]
            scannetChainResidueBinding.update({"{}_{}".format(chain, residue) : binding})

        scannetResults.update({structure_id : scannetChainResidueBinding})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006E-ScanNetResults.json')    
    with open(outfile, 'w') as f:
        json.dump(scannetResults,f )
    log.info("Saving ScanNet Results ... \n\t{}".format(outfile))
    
    t = 0 
    checked = []
    global_scannet = {}
    for sfileId, scannetBinding in tqdm(scannetResults.items()):
 
    # sfileId = '1-PFK_unrelaxed_rank_1_model_3'
    # scannetBinding = scannetResults[sfileId]
        if sfileId in dfalldata.structureId.unique():
            structureId = sfileId
        else:
            try:
                structureId = dfstructure_file_locations[dfstructure_file_locations.structureFileId == sfileId].index[0]
            except IndexError:
                continue
        dfs = dfalldata[dfalldata.structureId == structureId]

        t += len(dfs[dfs.sfileChain_Residue.isin(scannetBinding)])
        print (t, '\t\t',len(dfs[dfs.sfileChain_Residue.isin(scannetBinding)]),'\t',len(dfs), '\t',structureId)
        sfileChain_Residue_globIndex = dfs[dfs.sfileChain_Residue.isin(scannetBinding)].sfileChain_Residue.to_dict()
        for index, sfileChain_Residue in sfileChain_Residue_globIndex.items():
            global_scannet.update({index : scannetBinding[sfileChain_Residue]})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006E-StructuralProteome-interfaces_ScanNet.json')    
    with open(outfile, 'w') as f:
        json.dump(global_scannet,f )
    
    log.info("Saving QSPACE Proteome - protein interfaces (ScanNet)...\n\t{}".format(outfile))

    return scannetResults,global_scannet

    