# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def run_005B_writeOPMpdbs(potential_membrane_structures,
                          outfolder = qspaceDirs['opmStructuresToSendDir'],
                         force_rerun= True,
                         ):
    
    
    ###Get locations of original files
    membrane_structure_dict = {}
    for structure in potential_membrane_structures:
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

        membrane_structure_dict.update({structure : stype})
    
    dfstructure_file_locations = pd.DataFrame()
    dfstructure_file_locations['stype'] = pd.Series(membrane_structure_dict)
    for index, row in dfstructure_file_locations.iterrows():
        sfile = find_sfile(index, row.stype)
        dfstructure_file_locations.loc[index, 'sfile'] =  sfile
    
    structreOPMChains = {}
    infile_OPMchains = op.join(qspaceDirs['DataOutput_dir'], "005B-structChain_opmChain.json")
    if op.exists(infile_OPMchains):
        with open(infile_OPMchains ,'r' ) as f:
            structreOPMChains = json.load( f)
            
            
            
#     print '-------------------------------------------\nWriting .pdb files for OPM in ...\n\t> {}'.format(outfolder)
    for s_id  in tqdm(dfstructure_file_locations.index.tolist()):
        row = dfstructure_file_locations.loc[s_id]
        sfile = row.sfile

        outfile = op.join(outfolder, "{}.pdb".format(s_id))
        if op.exists(outfile) and not force_rerun:
            dfstructure_file_locations.loc[s_id, 'to_opm_file'] = outfile
            if s_id in structreOPMChains:
                dfstructure_file_locations.loc[s_id, 'origChain_opmChain'] = str(structreOPMChains[s_id])
            else:
                dfstructure_file_locations.loc[s_id, 'origChain_opmChain'] =  str('same_as_pdb')
#             dfstructure_file_locations.loc[s_id, 'origChain_opmChain'] = str(structreOPMChains[s_id])
            continue

        if ".pdb" in sfile:
            shutil.copy(src = sfile ,dst = outfile)
            dfstructure_file_locations.loc[s_id, 'to_opm_file'] = outfile
            dfstructure_file_locations.loc[s_id, 'origChain_opmChain'] = str('same_as_pdb')

            continue
#         print (sfile)
        s = StructProp(ident=s_id, structure_path=sfile, file_type = op.basename(sfile).split('.')[-1])
        p = s.parse_structure()

        with open(outfile ,'wb') as f:
            f.write('MODEL \n'.encode())

            i = 1

            used_chains = {}
            for chain in p.first_model.get_chains():
                #we need to select a new chaing for sending PDBs to OPM
                for new_chain in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ123456789':
                    if new_chain not in used_chains.values():
                        break

    #             print '\n\n' , chain.id
                for residue in chain.get_residues():
    #                 print residue.id[1],
                    for atom in residue.get_atoms():
            #             print chain.id, residue.id, atom.id
                        newline = pdbLine(new_chain, residue, atom, i)
                        f.write(newline.encode())

                        i +=1

                used_chains.update({chain.id : new_chain})
    #             break


                f.write('TER \n'.encode())
            f.write('END \n'.encode())
        structreOPMChains.update({s_id:used_chains})   
        dfstructure_file_locations.loc[s_id, 'to_opm_file'] = outfile
        dfstructure_file_locations.loc[s_id, 'origChain_opmChain'] = str(used_chains)
    
    
    
    with open(infile_OPMchains ,'w' ) as f:
        json.dump(structreOPMChains,f)
    log.info("Saving OPM chain translation ....\n\t{}".format(infile_OPMchains))
    outfile = op.join(qspaceDirs['DataOutput_dir'], "005B-opmStructureLocations.csv")
#     print "Saving....\n\t> {}".format(outfile)
    dfstructure_file_locations.to_csv(outfile)
    
    return dfstructure_file_locations,structreOPMChains
