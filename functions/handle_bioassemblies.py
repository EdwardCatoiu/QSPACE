import os
import os.path as op
import ast
import pandas as pd
from tqdm import tqdm_notebook as tqdm
import json
import sys

sys.path.append('/usr/local/lib/python3.7/dist-packages/')
import pdbecif
from pdbecif.mmcif_io import CifFileReader



ignore_list_bioassembly_mapping = ['_chem_comp',
 '_audit_author',
 '_struct_mon_prot_cis',
 '_diffrn',
 '_pdbx_struct_assembly_gen',
 '_atom_sites',
 '_atom_type',
 '_reflns',

 '_entity_poly',
 '_pdbx_audit_revision_details',
 '_struct_ref_seq',
 '_citation_author',
 '_pdbx_struct_oper_list',
 '_pdbx_refine_tls',
 '_citation',
 '_struct_sheet',
 '_symmetry',
 '_pdbx_validate_rmsd_bond',
 '_pdbx_audit_revision_category',
 '_pdbx_struct_special_symmetry',
 '_pdbx_struct_mod_residue',
 '_software',
 '_pdbx_audit_revision_history',
 '_database_PDB_matrix',
 '_refine',
 '_reflns_shell',
 '_pdbx_validate_torsion',
 '_struct_ref_seq_dif',
 '_pdbx_database_remark',
 '_exptl_crystal_grow',
 '_exptl_crystal',
 '_exptl',
 '_pdbx_struct_sheet_hbond',
 '_pdbx_database_related',
 '_cell',
 '_struct_sheet_order',
 '_pdbx_poly_seq_scheme',
 '_pdbx_struct_assembly_prop',
 '_struct_conf',
 '_pdbx_audit_revision_item',
 '_pdbx_nonpoly_scheme',
 '_audit_conform',
 '_diffrn_detector',
 '_struct_conn',
 '_refine_hist',
 '_struct_asym',
 '_entity_poly_seq',
 '_diffrn_source',
 '_struct_biol',
 '_struct_sheet_range',
 '_pdbx_entity_nonpoly',
 '_struct_site',
 '_atom_site',
 '_pdbx_database_status',
 '_diffrn_radiation',
 '_entity_src_gen',
 '_struct',
 '_pdbx_validate_rmsd_angle',
 '_pdbx_audit_revision_group',
 '_struct_conf_type',
 '_pdbx_SG_project',
 '_pdbx_refine_tls_group',
 '_struct_keywords',
 '_struct_ref',
 '_refine_ls_shell',
 '_database_2',
 '_struct_site_gen',
 '_pdbx_unobs_or_zero_occ_residues',
 '_struct_conn_type',
 '_diffrn_radiation_wavelength',
 '_refine_ls_restr',
 ]

import sys
with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)


def create_bioassembly_df(assembly_info,
                          pdb_id,
                          dfp,
                          cif_bioassembly_dir = qspaceDirs['bioassemblyStructuresDir'],
                          ignore_list = []):
    
    
    oligomeric_count = assembly_info['oligomeric_count']
    oligomeric_details = assembly_info['oligomeric_details']
    oligomeric_id = assembly_info['id']
    method_details = assembly_info['method_details']
    details = assembly_info['details']
    

    if type(oligomeric_id) != list:
        oligomeric_count = [oligomeric_count]
        oligomeric_details = [oligomeric_details]
        oligomeric_id = [oligomeric_id]
        method_details = [method_details]
        details = [details]
    
    
    info = []
    for i, assemblyId in enumerate(oligomeric_id):
#         print assemblyId, oligomeric_count[i],oligomeric_details[i],method_details[i],details[i]

        cif_bioassembly = op.join(cif_bioassembly_dir,'{}-assembly{}.cif'.format(pdb_id.lower(),assemblyId) ) 
        if not op.exists(cif_bioassembly):
            print ('missing : {}'.format(cif_bioassembly))
            continue
#         print ("reading bioassembly data")
        data_assembly = CifFileReader().read(cif_bioassembly , ignore=ignore_list)

        entity_dict = {}
#         print (i, assemblyId, data_assembly.keys())
        assembly_key = list(data_assembly.keys())[0]
        
        for k, entity_id in enumerate(data_assembly[assembly_key]['_entity']['id']):
            stoich =data_assembly[assembly_key]['_entity']['pdbx_number_of_molecules'][k]

            if entity_id not in dfp.entity_id.unique() and int(entity_id) not in dfp.entity_id.unique() :
#                 print pdb_id,entity_id, type(entity_id)
                continue

            df_entity = dfp[dfp.entity_id.isin([ str(entity_id),int(entity_id)])]
            entity_dict.update({"ent_{}".format(entity_id) : int(stoich)})

#             print entity_id , stoich

            
        info.append( [pdb_id,'{}-assembly{}.cif'.format(pdb_id.lower(),assemblyId), str(entity_dict),oligomeric_details[i] ,oligomeric_count[i],method_details[i],details[i]])
#         dfpbioassembply.loc[index,'pdb_entry'] = pdb_id
#         dfpbioassembply.loc[index,'assembly_id'] = '{}-assembly{}.cif'.format(pdb_id.lower(),assemblyId)
#         dfpbioassembply.loc[index,'entity_stoich'] = str(entity_dict)
#         dfpbioassembply.loc[index,'oligomeric_details'] = oligomeric_details[i]
#         dfpbioassembply.loc[index,'oligomeric_count'] = oligomeric_count[i]
#         dfpbioassembply.loc[index,'method_details'] = method_details[i]
#         dfpbioassembply.loc[index,'details'] = details[i]
#         index +=1

    return info


def get_BIO_gene_maps(pdbMaps, dfpdb_mapped, dfbioassembly, chain_source = 'auth_asym_ids'):

    outdict = {}
    for pdb_entry in tqdm(list(pdbMaps.keys())):
        dfp = dfpdb_mapped[dfpdb_mapped.pdb_entry == pdb_entry]
        dfb = dfbioassembly[dfbioassembly.pdb_entry == pdb_entry]

        for index, row in dfb.iterrows():
            entity_stoich = ast.literal_eval(row.entity_stoich)

            for ent, ent_stoich in entity_stoich.items():
                entity_id = ent.split('_')[-1]

                
                dfp_e = dfp[dfp.entity_id.isin([int(entity_id), entity_id])]
#                 print pdb_entry, index, ent, len(dfp_e)
                chains_in_entity = dfp_e.get(chain_source).tolist()[0]
                if type(chains_in_entity) == str:
                    chains_in_entity = ast.literal_eval(chains_in_entity)

                for i in range(int(ent_stoich / len(chains_in_entity))):

                    for chain in chains_in_entity:
                        chain_id = "{}-{}".format(chain, i)
                        for gene, entity_mmseq in pdbMaps[pdb_entry].items():
                            if entity_id in entity_mmseq:
                                score = entity_mmseq[entity_id]


                                if row.assembly_id not in outdict:
                                    outdict.update({row.assembly_id : {gene : {chain_id : [score >=0.5, score]}}})
                                elif gene not in outdict[row.assembly_id]:
                                    outdict[row.assembly_id].update({gene : {chain_id : [score >=0.5, score]}})
                                elif chain_id not in outdict[row.assembly_id][gene]:
                                    outdict[row.assembly_id][gene].update({chain_id : [score >=0.5, score]})
    return outdict
