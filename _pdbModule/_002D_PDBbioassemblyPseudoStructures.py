# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

import Bio

def run_002D_pseudo_structures(bioMaps,dfrepseq, query_type = 'PDB', dfpdb_structure_info = pd.DataFrame()):
    structure_info = dfpdb_structure_info
    
    
    print ('Finding Pseudo-Structures for {} structures...\n----------------------------'.format(query_type))
    input_stoich_chains_filtered = pseudo.remake_chains_from_quality(bioMaps)
    pseudo_stoich_filtered,pseudo_stoich_chains_filtered  = pseudo.find_pseudo_matches(input_stoich_chains_filtered)

    print ('Getting gene-chain quality for Pseudo-Structures...\n----------------------------')
    pseudo_quality_filtered = pseudo.get_quals(pseudo_stoich_chains = pseudo_stoich_chains_filtered,
                                               input_quality = bioMaps)
    
    
    
    df,ssbio_errors = pseudo.make_dataframe_pseudo_structures(pseudo_stoich_filtered, pseudo_quality_filtered,dfrepseq)
    
    
    print ('Getting total structure quality for Pseudo-Structures...\n----------------------------')
    df3 = pseudo.assign_structure_quality(df,df_gene_length=dfrepseq)
    #normalized from 0-1 --> 0-100
    df3['total_quality'] = df3['total_quality'] * 100.
        
    print ('Ranking Pseudo-Structures by total quality...\n----------------------------')

    df4 = pseudo.assign_quality_ranks(df3)

    pdb_entry_dict = {}
    for index, row in tqdm(df4.iterrows()):
        pdb_entry = index.split('-assembly')[0]
        pdb_entry_dict.update({index : pdb_entry})

    df4['pdb_entry'] = pd.Series(pdb_entry_dict)
    
    #this line uses the structure's resolution to determine the 
    #best pseudo-structure from a list of ones that share the 
    #same normlized score from the PDB mmSeqs Seuqence Service
    print ('Determining best quality Pseudo-Structures...\n----------------------------')
    best_structures_list = pseudo.find_best_structures(df4,query_type=[query_type],structure_info = structure_info)

    dfbest =  df4[df4.index.isin(best_structures_list[query_type])]

    glist = []
    for gdi in dfbest.gene_stoichiometry:
        if type(gdi) == str:
            gdi = ast.literal_eval(gdi)
        glist += gdi.keys()
#     print "{} Genes accounted for in the {} structures".format(len(set(glist)),len(dfbest))

    outfile = op.join(qspaceDirs['DataOutput_dir'],'002D-all_structures_{}.csv'.format(query_type))
    df4.to_csv(outfile)
    log.info("Saving ALL {} pseudo-structures ...\n\t{}".format(query_type, outfile))
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'002D-best_structures_{}.csv'.format(query_type))
    dfbest.to_csv(outfile)
    log.info("Saving BEST {} pseudo-structures ...\n\t{}".format(query_type, outfile))
    return df4,dfbest

def run_002D_downloadPDBfiles(dfbest,
                              pdb_folder = qspaceDirs['pdbStructuresDir'],
                              force_rerun = False,
                              download_errors = []
                             ):
    print ('Downloading .pdb files for BEST quality Pseudo-Structures...\n----------------------------')
    for pdb_entry in tqdm(dfbest.pdb_entry.unique()):
        outfile = op.join(pdb_folder, "{}.pdb".format(pdb_entry.lower()))
        if op.exists(outfile) and not force_rerun:
            continue
        try:
            download_link = 'https://files.rcsb.org/download/{}.pdb'.format(pdb_entry)
            html = urllib.request.urlopen(download_link).read()
            with open(outfile, 'wb') as f:
                f.write(html)
                f.close()
        except (urllib.error.HTTPError) as e:
            download_errors += [pdb_entry]
            log.warn("Could not download .pdb for PDB entry :  {}".format(pdb_entry))

#             print '{} :  \t {}'.format(e, pdb_entry)   
    return download_errors

def run_002D_getAASeqInPDBFile(dfbest,
                               pdb_folder = qspaceDirs['pdbStructuresDir'],
                               cif_folder = qspaceDirs['cifStructuresDir'],
                               force_rerun = False,
                              ):
    errors_with_pdb = []
    checked_pdbs= []
    appender = []

    #use existing information for all the AA seqs of all chains in all PDBs
    outfile = op.join(qspaceDirs['DataOutput_dir'],"002D-chain_seqs_from_PDBs.csv")
    if op.exists(outfile) and not force_rerun:
        dfStructureSeqs = pd.read_csv(outfile, index_col = 0 )
        checked_pdbs = list(dfStructureSeqs.pdb_entry.unique())
        appender = list(dfStructureSeqs.to_records(index=False))
#         print "Using previous data..."
#         print "\t> {} ".format(outfile)
        
#     print '--------------------------------------------------'
    print ('Getting AA sequences in each chain of {} .pdb structure files......\n------------------'.format(len(dfbest.pdb_entry.unique())))
    for pdb_entry in tqdm(dfbest.pdb_entry.unique()):
        if pdb_entry in checked_pdbs or pdb_entry in errors_with_pdb:
            continue
        structureFile = op.join(pdb_folder, "{}.pdb".format(pdb_entry.lower()))
        print (structureFile)
        #check to see that we have a error-free structure file
        try:
            a_structure = StructProp(ident = pdb_entry,structure_path=structureFile,file_type ='pdb')
            parsed = a_structure.parse_structure()
        except (IOError,IndexError, Bio.PDB.PDBExceptions.PDBConstructionException, ValueError, KeyError) as e:
            errors_with_pdb +=[pdb_entry]
            log.warn("Error loading/parsing structure {}".format(op.basename(structureFile)))
#             print 'error : {}'.format(pdb_entry)
            continue

        #get the amino acid sequence for every chain in the PDB file
        chain_dict = {}
        for chain in parsed.first_model.get_chains():
            #     print chain.id
            chain_prop = a_structure.chains.get_by_id(chain.id)
            chain_prop_seq  = str(chain_prop.seq_record.seq)

            info = [pdb_entry, structureFile, chain.id, chain_prop_seq]
    #         chain_dict.update({chain.id: chain_prop_seq})
            appender.append(info)
        checked_pdbs += [pdb_entry]
        
    
    
        
    checked_cifs= []
    appender_cif = []
    errors_with_cif=[]
    
#     print '--------------------------------------------------'
    print ('\nThere were errors using the .pdb files for {} structures......\n-----------------------'.format(len(errors_with_pdb)))
    print ('Getting AA sequences in each chain of {} .cif structure files......\n-----------------------'.format(len(errors_with_pdb)))

    for pdb_entry in tqdm(errors_with_pdb):
        if pdb_entry in checked_cifs or pdb_entry in errors_with_cif:
            continue

        dfb = dfbest[dfbest.pdb_entry.isin([pdb_entry])]
        worth_the_computation = False
        for index, row in dfb.iterrows():
            pdb_quality = row.pdb_quality
            if type(pdb_quality) == str:
                pdb_quality = ast.literal_eval(pdb_quality)
            if len(pdb_quality.keys()) > 20:
#                 print len(pdb_quality.keys()), index 
                worth_the_computation = True

    #     if not worth_the_computation:
    #         continue

        structureFile = op.join(cif_folder, "{}.cif".format(pdb_entry.lower()))
        
#         print pdb_entry,
        #check to see that we have a error-free structure file
        try:
            a_structure = StructProp(ident = pdb_entry,structure_path=structureFile,file_type ='cif')
            parsed = a_structure.parse_structure()
        except (IOError,IndexError, Bio.PDB.PDBExceptions.PDBConstructionException,ValueError, KeyError) as e:
            errors_with_cif +=[pdb_entry]
            log.warn("Error loading/parsing structure {}".format(op.basename(structureFile)))
            continue

        #get the amino acid sequence for every chain in the PDB file
        chain_dict = {}
        for chain in parsed.first_model.get_chains():
            #     print chain.id
            chain_prop = a_structure.chains.get_by_id(chain.id)
            chain_prop_seq  = str(chain_prop.seq_record.seq)

            info = [pdb_entry, structureFile, chain.id, chain_prop_seq]
    #         chain_dict.update({chain.id: chain_prop_seq})
            appender_cif.append(info)
        checked_cifs += [pdb_entry]
    #     seq_dict.update({pdb_entry : chain_dict})
    
    
    
    dfStructureSeqs = pd.DataFrame()
    dfcols = ['pdb_entry', 'structure_file', 'auth_chain', 'chain_seq']
    dfStructureSeqs = pd.DataFrame.from_records(appender + appender_cif, columns = dfcols)
    dfStructureSeqs.to_csv(outfile )
    log.info("Saving Dataframe of all AA seqs for all chains in PDB/CIF files...\n\t{}".format(outfile))
#     print "\t> {}".format(outfile)
    return dfStructureSeqs


    
def run_002D_needle_alignment_quality(dfbest,
                                      dfStructureSeqs_all,
                                      dfrepseq,
                                      needle_dir = qspaceDirs['SequenceAlignmentDir'],
                                      checked = [],
                                      force_rerun = False
                                     ):

    
    
    pdb_quality_needle = {}
    
    print ("Using needle alignment to get sequence identify of all chains in {} structure files.....\n---------------------------------".format(len(dfbest)))
    for index in tqdm(dfbest.index.tolist()):
        if index in checked:
            continue
#         print index
        row = dfbest.loc[index]

    #     if row.pdb_entry !='7lo2':
    #         continue
        dfseqs_pdb = dfStructureSeqs_all[dfStructureSeqs_all.pdb_entry == row.pdb_entry]


        pdb_quality = row.pdb_quality
        if type(pdb_quality) == str:
            pdb_quality = ast.literal_eval(pdb_quality)

        for blattner, blattner_chains in (pdb_quality.items()):
            if 'UniProtId' in dfrepseq.loc[blattner].dropna().keys():
                seqId_gene = dfrepseq.loc[blattner,'UniProtId']
                seqAA_gene = dfrepseq.loc[blattner,'UniProtSeq']
            else:
                seqId_gene = "WT_consensus"
                seqAA_gene = dfrepseq.loc[blattner,'AlleleomeSeq']


            for chainId, chainInfo in blattner_chains.items():
                chainId_structure = chainId.split('-')[0]

                #force re-run essentially
                needle_file =  "{}_{}_{}-{}.needle".format(row.pdb_entry.lower(), blattner,seqId_gene,chainId_structure)
    #             if op.exists(op.join(needle_dir, needle_file)):
    #                 continue

                dfseqs_pdb_chain = dfseqs_pdb[dfseqs_pdb.auth_chain == chainId_structure]
                if len(dfseqs_pdb_chain) == 0:
                    log.warn("Error w/Chain-sequence for {}-{}".format(row.pdb_entry, chainId_structure))

#                     print 'Error w/ chain seqs : {}-{}'.format(row.pdb_entry, chainId_structure)
                    continue
                chain_seq = dfseqs_pdb_chain.chain_seq.values[0]




                needle_file = ssbioaln.run_needle_alignment(seq_a =chain_seq,
                                                         seq_b = seqAA_gene,
                                                         outfile = needle_file, 
                                                         outdir = needle_dir,
                                                         force_rerun=force_rerun)
#                 print needle_file
                try:
                    needle_statistics = ssbioaln.needle_statistics(needle_file)['asis_asis']
                except IOError:
                    log.warn("ERROR: could not fine needle file :{}".format(needle_file))
                    continue
                    
                try:
                    identity = needle_statistics['percent_identity']
                except KeyError:
                    identity = 0
                    log.warn("ERROR: needle alignment identity = 0% : {}".format(needle_file))

#                     print 'Statistics Error : {}'.format(needle_file)
               
                new_chain_info = [identity > 50, identity]
                
                #### fill out the new sequence identity quality of gene-chain map
                if index not in pdb_quality_needle:
                    pdb_quality_needle.update({index : {blattner : {chainId : new_chain_info}}})

                elif blattner not in pdb_quality_needle[index]:
                    pdb_quality_needle[index].update({blattner : {chainId : new_chain_info}})
                elif chainId not in pdb_quality_needle[index][blattner]:
                    pdb_quality_needle[index][blattner].update({chainId : new_chain_info})
        checked +=[index]
        
    dfbest['pdb_quality_needle'] = pd.Series(pdb_quality_needle)

    df3_b = pseudo.assign_structure_quality(dfbest,df_gene_length=dfrepseq,quality_key = 'pdb_quality_needle')
    df4_b = pseudo.assign_quality_ranks(df3_b)
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'002D-best_structures_PDB.csv')
    df4_b.to_csv(outfile)
    log.info("Saving BEST {} pseudo-structures w/needle alignment quality...\n\t{}".format("PDB", outfile))
#     print ">{} ".format(outfile)

    return df4_b
        

    