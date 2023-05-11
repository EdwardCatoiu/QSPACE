# -*- coding: utf-8 -*-

from .utils import *
log = logging.getLogger(__name__)

        
def run_001C(itasser_metadata,
             itasserRepoFolder, 
             blattnerUniprotMapping,
             needs_scoring = False,
             force_clean  = False,
             force_realign = False,
             force_remove_string = False,
            ):
    
    """
    ### Description
    The purpose of this notebook is to:
    - determine the TM-scores of itasser models (if not already found)
    - move the Itasser files from the ITASSER/REPOSITORY into our appropriate folder
    - clean/remove string-like regions of itasser files
    - determine the real chains in each ITASSER model
    - find missing alignment files (gene-chain-structure alignment)

    ### Input(s)
    1. ITASSER metadata
    2. ITASSER Model Repository

    ### Output(s)
    1. Find TM-scores for ITASSER models 
    2. Move ITASSER models into appropriate folder 
    3. Clean itasser models (ssbio)
    4. Remove string-like regions from itasser models
    5. Determine real protein chains in Itasser models
    6. Expected needle alignment files 

    ### Folder(s)
    ```
    -------------------------InPuts------------------------------------------------
    qspace
    └───inputs
    │   |   1. 001C-ITASSER_raw_metadata_hyperlink.csv
    |       ...
    //
    ItasserModelRepository <User Defined>
    │   2. {geneID}_ECOLI_model1.pdb
    |       ...
    -------------------------OutPuts------------------------------------------------
    qspace
    └───data
    │   |   1. 001C-ITASSER_raw_metadata_hyperlink_tmscores.csv
    │   |   4. 001C-itasser_stringlike_residues.csv
    │   |   5. 001C-ITASSER-model_chains_to_genes.json
    │   |   6. 001C-expected_alignment_files_ITASSER.json
    |       ...
    //
    externalHardDrive <User Defined>
    └───structures
    |   |
    |   └───all_itasser
    |   |   |   2. {geneID}_ECOLI_model1.pdb
    |   |       ...
    |   └───all_itasser_clean
    |   |   |   3. {geneID}_ECOLI_model1_clean.pdb
    |   |       ...
    |   └───all_itasser_clean_string_removed
    |   |   |   4. {geneID}_ECOLI_model1_clean_residues_removed.pdb
    ```
    """
    
    if needs_scoring:
#         print '------------------------------------------------'
        print ('Getting I-TASSER TM-score from... \t https://zhanggroup.org/Ecoli/ \n----------------------------------------------')        
        for index in tqdm(itasser_metadata.index.tolist()):
            row  = itasser_metadata.loc[index]
            if "TM-score" in row.dropna().index.tolist() and 'template' in row.dropna().index.tolist():
                continue


            url = row.hyperlink
        #     break

            page = requests.get(url)
            model1_tm_score_loc = page.text.find('TM-score')
            tmscore_txt = page.text[model1_tm_score_loc: model1_tm_score_loc + 30].split(' <br/')[0]
            tmscore_template = tmscore_txt.split(' to ')[-1]
            tmscore_score = float(tmscore_txt.split(' to ')[0].split('=')[-1])
            tmscore_template,tmscore_score
            itasser_metadata.loc[index,'TM-score'] = tmscore_score
            itasser_metadata.loc[index,'template'] = tmscore_template
        #     print index, tmscore_score, tmscore_template, url
            time.sleep(0.1)

        outfile = op.join(qspaceDirs['DataOutput_dir'], '001C-ITASSER_raw_metadata_hyperlink_tmscores.csv')
        log.info('Saving Quality of ITASSER models... \n\t{}'.format(outfile))

        itasser_metadata.to_csv(outfile)
    
    

    all_itasser_dir = qspaceDirs['itasserStructuresDir']
#     print '------------------------------------------------'
#     print 'Copying relevant models from the ITASSER repository...'
#     print '\t> {}'.format(itasserRepoFolder)
#     print '\t> {}'.format(all_itasser_dir)
    
    itasser_sfile = {}
    for index  in tqdm(itasser_metadata.index.tolist()):
        row = itasser_metadata.loc[index]
        newFile = op.join(all_itasser_dir,'{}_model1.pdb'.format(row.get('Entry Name')))
        itasser_metadata.loc[index,'sfile'] = newFile
    #     itasser_sfile.update({index : newFile})
        if not op.exists(newFile):    
            sfile = False
            for f in os.listdir(itasserRepoFolder):
                if '{}_model1.pdb'.format(row.get('Entry Name')) == f:
                    sfile = op.join(itasserRepoFolder, f)
                    break

            if sfile:
    #             print sfile
                itasser_metadata.loc[index,'sfile'] = newFile
                shutil.copy(src = sfile, dst = all_itasser_dir)
                continue
            else:
                pass
#                 print row.get("Entry Name")

                
   

    clean_itasser_dir = qspaceDirs['itasserCleanStructuresDir']
#     print '------------------------------------------------'
        
    itasser_outliers_outfile = op.join(qspaceDirs['DataOutput_dir'],'001C-ITASSER-stringOutliers.csv')
    
    
#     print '\t> {}'.format(clean_itasser_dir)
    if not op.exists(itasser_outliers_outfile) or force_clean:
        print ('Cleaning ITASSER MODELS...\n----------------------------------')

        for f in tqdm(os.listdir(all_itasser_dir)):
            if '.pdb' not in f:
                continue
            if f[0] =="E" or "ECOLI" in f:
                origFile = f
                s_id = f.split('_model')[0]
                origFilePath = op.join(all_itasser_dir, f)
                my_clean_pdb = cleanpdb.clean_pdb(pdb_file=origFilePath, outdir=clean_itasser_dir,force_rerun=force_clean)
    
    
    if not op.exists(itasser_outliers_outfile) or force_remove_string:

        print ('Removing string-like regions from ITASSER MODELS...\n---------------------------------')

        outlier_residues = {}
        folder = clean_itasser_dir
        for f in tqdm(os.listdir(folder)):
            filepath = op.join(folder, f)
            entry_name = f.split('_model')[0]
            if entry_name not in itasser_metadata.get("Entry Name").values:
                continue


            temp = {}
            df_atoms,inner_fence, view,total_residues_removed,total_residues  ,iteration_count= clean_itasser_string.visualize_outliers(filepath=filepath)

            temp.update({'residues_removed_IDs' : total_residues_removed})
            temp.update({'residues_removed' : float(len(total_residues_removed)) / float(total_residues)})
            temp.update({'model_similarity' : 1. - float(len(total_residues_removed)) / float(total_residues)})
            temp.update({'model_filepath' : filepath})
            temp.update({'iterations' :iteration_count })

            outlier_residues.update({f : temp})

        df_outliers = pd.DataFrame.from_dict(outlier_residues, orient = 'index')
        
    
        string_removed_folder = qspaceDirs['itasserCleanStringRemovedStructuresDir']
    #     print '\nWriting cleaned/string-removed ITASSER MODELS to...'
    #     print '\t> {}'.format(string_removed_folder)


        for index, row in tqdm(df_outliers.iterrows()):
            filepath = row.model_filepath
        #     if '/all_itasser_clean_string_removed/' in filepath:
        #         continue

            base_id = op.basename(filepath).rstrip('.pdb') + '_residues_removed.pdb'

            structure = StructureIO(filepath)
            outfile = op.join(string_removed_folder, base_id)
            if not op.exists or force_remove_string:
                with open (outfile, 'w' ) as f:
                    for atom in structure.structure[0].get_atoms():
                        if 'X_%s' %(atom.parent.id[1]) in row.residues_removed_IDs:
                            continue
                        line = clean_itasser_string.atom_line(atom)
                        f.write(line)
                    f.write('TER')
                    f.close()
            df_outliers.loc[index, 'new_model_filepath'] = outfile
            
        df_outliers.to_csv(itasser_outliers_outfile)
        log.info('Saving String-like regions of ITASSER models... \n\t{}'.format(itasser_outliers_outfile))


    
    df_outliers = pd.read_csv(itasser_outliers_outfile, index_col = 0 )
    
    print ("Getting total chains in ITASSER files....\n----------------------------------")
    total_chains = get_structure_chains_itasser(structureFolder=qspaceDirs['itasserCleanStringRemovedStructuresDir'] ,
                                                itasser_metadata=itasser_metadata,
                                                blattner_uniprot_mapping=blattnerUniprotMapping,)
    
#     print 'Saving all real protein chains in ITASSER MODELS...'
    outfile = op.join(qspaceDirs['DataOutput_dir'],'001C-ITASSER-model_chains_to_genes.json')
    with open(outfile, 'w') as f :
        json.dump(total_chains,f)
    log.info('Saving gene-chain maps for ITASSER models... \n\t{}'.format(outfile))
        

    outfilepath = op.join(qspaceDirs['DataOutput_dir'],'001C-expected_alignment_files_ITASSER.json' )
    find_missing_needle_files(total_chains_dict=total_chains, 
                              outfilepath = outfilepath,
                              blattnerUniprotMapping=blattnerUniprotMapping,
                              SequenceAlignmentDir= qspaceDirs['SequenceAlignmentDir'],
                              force_rerun= force_realign, 
                             )
#     print 'Saving...'
    

    return df_outliers