# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def run_002B_downloadPDBs(list_of_pdbs,force_download = False,pdb_checked = [],outfolder = qspaceDirs['cifStructuresDir']):
#     print '-----------------------------------------'
    print ('Downloading .cif files for {} PDB entries...\n-------------------------------------'.format(len(list_of_pdbs)))
#     print '\t> {}'.format(outfolder)
    
    for pdb_id in tqdm(list_of_pdbs):
        if pdb_id in pdb_checked:
            continue


    #     if op.exists(op.join(dst_cif_folder_external, '{}.cif'.format(pdb_id.lower()))):
    #         pdb_checked  += [pdb_id]
    #         continue



        cif_file = download_structures.download_cif(pdb_id,
                                                     outfolder = outfolder,
                                                     force_download = force_download,
                                                     pdb_checked = pdb_checked)
        
        
        
        
def run_002B_downloadPDB_bioassemblies(dfpdb_mapped,
                                       force_download = False,
                                       outfolder = qspaceDirs['bioassemblyStructuresDir'],
                                       use_existing_data = True
                                      ):
    
    
    appender = []
    pdb_checked = []
    
            
    bioassembly_df_infile = op.join(qspaceDirs['DataOutput_dir'], '002B-ALL_mapped_BIOASSEMBLY.csv')
    if op.exists(bioassembly_df_infile):# and use_existing_data:
#         print '----------------------------------------------'
#         print 'Loading existing data...'
#         print '\t> {}'.format(bioassembly_df_infile)
        dfbioassembly = pd.read_csv(bioassembly_df_infile,index_col=0)
        pdb_checked = list(dfbioassembly.pdb_entry.unique())
        appender = list(dfbioassembly.to_records(index=False))


    pdbs_remaining = set(dfpdb_mapped.pdb_entry.unique().tolist()) - set(pdb_checked)
        
    print ('Downloading bioassemblies and determining oligomerization of entity/chains...\n------------------------------')
    run = 0 
    for pdb_id in tqdm(pdbs_remaining):
        cif_file = op.join(qspaceDirs['cifStructuresDir'],'{}.cif'.format(pdb_id.lower()))
        
        if not op.exists(cif_file): 
            log.warn("Error: missing file for {}".format(op.basename(cif_file)))
#             print "missing cif file : {}".format(cif_file)
            continue
            
        dfp = dfpdb_mapped[dfpdb_mapped.pdb_entry == pdb_id]
        if len(dfp) > 20:
            continue



#         print pdb_id ,

        bioassembly_list = download_structures.download_bioassemblies(pdb_id,
                                                                      outfolder = outfolder ,
                                                                      force_download = force_download,
                                                                     )

        data = CifFileReader().read(cif_file , ignore=handle_bioassemblies.ignore_list_bioassembly_mapping)
        try:
            assembly_info =  data[pdb_id]["_pdbx_struct_assembly"]
            info = handle_bioassemblies.create_bioassembly_df(assembly_info,
                                                          pdb_id=pdb_id,
                                                          dfp = dfp,
                                                          ignore_list = handle_bioassemblies.ignore_list_bioassembly_mapping,
                                                          cif_bioassembly_dir =outfolder
                                                         )


            appender +=info
        except (IndexError, KeyError) as e:
            log.warn("Error with determining oligomerization of bioassembly {}".format(op.basename(cif_file)))
#             print 'Error with bioassembly mapping ... {}'.format(cif_file)
            pass
        pdb_checked += [pdb_id]
        
        if run %200 == 0 :
            
            dfbioassembly = pd.DataFrame()
            dfcols = ['pdb_entry', 'assembly_id', 'entity_stoich', 'oligomeric_details', 'oligomeric_count',
                      'method_details', 'details']
            dfbioassembly = pd.DataFrame.from_records(appender, columns = dfcols)
            dfbioassembly.to_csv(bioassembly_df_infile)
        run +=1
    dfbioassembly = pd.DataFrame()
    dfcols = ['pdb_entry', 'assembly_id', 'entity_stoich', 'oligomeric_details', 'oligomeric_count',
              'method_details', 'details']
    dfbioassembly = pd.DataFrame.from_records(appender, columns = dfcols)
    dfbioassembly.to_csv(bioassembly_df_infile)
    log.info("Saving oligomerization state of bioassemblies in  dataframe...\n\t{}".format(bioassembly_df_infile))
#     print '----------------------------------------------'
#     print 'Saving bioassembly dataframe...'
#     print '\t> {}'.format(bioassembly_df_infile)

    return dfbioassembly