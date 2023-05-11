# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def run_002A_textService(blattner_uniprot_mapping,bestStructuresDir, force_ssbio_map = False,force_rerun_PDB_api= True):
    
    def getSsbioBestStructures(blattner_uniprot_mapping,bestStructuresDir,force_rerun=force_ssbio_map): 
        print ("Getting list of PDB structures through ssbio-UniProtId mapping...\n------------------------------------------------")
        c = 0 
        appender= []
        for gene, uniprot_id in tqdm(blattner_uniprot_mapping.items()): 
            outfile = op.join(bestStructuresDir, '{}_best_structures.json'.format(uniprot_id))
            if op.exists(outfile) and not force_rerun:
                with open(outfile, 'r') as f:
                    best_structures = json.load(f)
                    best_structures = best_structures[uniprot_id]
    #                 print gene

            else:    
                try:
                    uniprots_to_pdbs = bs_unip.mapping(fr='UniProtKB_AC-ID', to='PDB', query=[uniprot_id])
                    best_structures = ssbio_bs(uniprot_id,
                                           outname='{}_best_structures'.format(custom_slugify(uniprot_id)),
                                           outdir= bestStructuresDir,
                                           seq_ident_cutoff=0.0,
                                           force_rerun=force_rerun
                                          )

                except (TypeError) as e:
                    log.warn('Could not use SSBIO to find structures for {}'.format(uniprot_id))
                    continue
                    
            appender = PDB_search.get_info_from_best_structures(uniprot_id=uniprot_id, best_structures=best_structures, appender=appender)
            c +=1
        #     break

        df_uniprot_to_pdb = pd.DataFrame.from_records(appender, columns = ["uniprot_id", "pdb_id", "chain_id", "rank","coverage", "resolution", "start","end", "unp_start","unp_end", "experimental_method"])

        return df_uniprot_to_pdb

    
    
    
    df_uniprot_to_pdb = getSsbioBestStructures(blattner_uniprot_mapping,bestStructuresDir)
    outfile = op.join(qspaceDirs['DataOutput_dir'],'002A-uniprot_to_pdb.csv') 
    df_uniprot_to_pdb.to_csv(outfile)
    log.info("Saving ssbio functions uniprot-to-pdb dataframe...\n\t{}".format(outfile))
    
#     print "\t> {}".format(outfile)
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'002A-uniprot_to_pdb_full_info.csv') 
    print ("Getting structure-chain-uniprotIds for PDBs identfied via text service (uniprot ID)...\n------------------------------------------------")

    if op.exists(outfile) and not force_rerun_PDB_api:
        dfpdb_search_uniprotID = pd.read_csv(outfile,index_col = 0)
        
        
        log.info("Not using API. Using previous PDB-text service information...\n\t{}".format(outfile))
    else:
        list_of_PDB_entries = df_uniprot_to_pdb.pdb_id.unique()
        dfpdb_search_uniprotID = textService_getPdbInfo(list_of_PDB_entries)
    dfpdb_search_uniprotID.to_csv(outfile)
    log.info("Saving structure-chain-uniprotIds for PDBs identified via text service/ssbio...\n\t{}".format(outfile))
    
#     print "\t> {}".format(outfile)
    return df_uniprot_to_pdb,dfpdb_search_uniprotID


def run_002A_sequenceService(blattner_uniprot_mapping, dfseq, query_type = "UniProt",force_rerun_api = False, force_genes = set()):
    
    genome_results = {}
    error_genes = {}
    
    
    #only use sequences that are in the blattner/uniprot mapping
    dfseq = dfseq[dfseq.index.isin(blattner_uniprot_mapping.keys())]
    
    if query_type == "UniProt":
        dfseq = dfseq[dfseq.UniProtSeq.isna() == False]
        aa_sequences = dfseq.UniProtSeq.to_dict() 
    elif query_type =='Alleleome':
        dfseq = dfseq[dfseq.AlleleomeSeq.isna() == False]
        aa_sequences = dfseq.AlleleomeSeq.to_dict()
    else:
        raise (KeyError , 'query_type needs to be in ["UniProt","Alleleome"]')
    
    #use existing genome_results 
    json_file = op.join(qspaceDirs['DataOutput_dir'],'002A-BLAST_PDB_{}Seq.json'.format(query_type)) 
    if op.exists(json_file) and not force_rerun_api:
        with open(json_file, 'r') as f:
            genome_results = json.load(f)
        genome_results = ast.literal_eval(genome_results)
        log.info("Using Previous PDB search API sequence service results...\n\t{}".format(json_file))
#         print "\t> {}".format(json_file)

        
        
    #find remaining genes to check (for multiple runs)
    genes_to_check = set(aa_sequences.keys()) - set(error_genes.keys()) - set(genome_results.keys())
    genes_to_check = genes_to_check.union(force_genes)
#     print len(genes_to_check), len(genome_results), len(error_genes)
#     print '\n------------------------------------------------'
    print ("Running {} sequences using PDB sequenceService API for {} genes...\n------------------------------".format(query_type, len(genes_to_check)))

    for gene  in tqdm(genes_to_check):
        seq = aa_sequences[gene]

        url = blast_api.encode_seq_request_new(seq)
        req= requests.get(url)

        if req.status_code != 200:
            error_genes.update({gene : req.status_code})
            log.info("No PDBs available for {} sequence for gene {}".format(query_type,gene))
            continue
        response = req.text
        response = ast.literal_eval(response)
        results_dict = blast_api.get_results_dict_new(response)
        genome_results.update({gene :results_dict })
        #dump the genome_results with every iteration. yes its slower but its better to save inside the loop
        with open(json_file, 'w') as f:
            json.dump(json.dumps(genome_results), f)
#     print "\t> {}".format(json_file)

    return genome_results


def run_002A_structureInfo(PDB_files,force_rerun = False):
    
#     print '\n------------------------------------------------'
    print ("Getting experimental/resolution information for PDB structures...\n-----------------------------------------")
    appender = []
    pdb_checked_for_structure = []
    infile = op.join(qspaceDirs['DataOutput_dir'], '002A-PDBstructureInfo.csv')
    
    if op.exists(infile):# and not force_rerun:
        dfpdb_structure_info = pd.read_csv(infile,index_col = 0)
        pdb_checked_for_structure = list(dfpdb_structure_info.pdb_id.unique())
        appender = list(dfpdb_structure_info.to_records(index=False))
#         print "\t Using existing data..."
#         print "\t> {}".format(infile)

    for pdb_id in tqdm(PDB_files):
        if pdb_id in pdb_checked_for_structure:
            continue

        url_encoded = PDB_search.encode_query_structure_info(pdb_id)
        url = 'https://data.rcsb.org/graphql?query={}'.format(url_encoded)

        req = requests.get(url)
        if req.status_code != 200:
            log.info("No Structure information for PDB ID : {}".format(pdb_id))
            pass
        response = req.text
        response = response.replace('null','"null"')
        response_dict = ast.literal_eval(response)['data']['entries']

        resolution = []
        release_date = []
        structure_type = []
        for data in (response_dict):
            resolution += data['rcsb_entry_info']['resolution_combined']
            release_date += [data['rcsb_accession_info']['initial_release_date']]
            structure_type += [data['rcsb_entry_info']['structure_determination_methodology']]

    #     print resolution
        if len(resolution) == 1:
            resolution = resolution[0]
        if len(release_date) == 1:
            release_date = release_date[0]
        if len(structure_type) == 1:
            structure_type = structure_type[0]

        info = [pdb_id, str(release_date),str(resolution),str(structure_type)]
        appender.append(info)
        pdb_checked_for_structure+=[pdb_id]
    
    
    dfpdb_structure_info = pd.DataFrame.from_records(appender, columns=['pdb_id','release_date','resolution','structure_type'] )
    dfpdb_structure_info.to_csv(infile)
    log.info("Saving PDB experimental/resolution information...\n\t{}".format(infile))

#     print "Saving data..."
#     print "\t> {}".format(infile)
    return dfpdb_structure_info
    
    
    

    
