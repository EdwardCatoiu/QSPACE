# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)



def run_002C_mapAPI_DataToBioassemblies(sequenceService_api_result,
                                        dfpdb_mapped,
                                        dfbioassembly,
                                        query_type = 'Uniprot',
                                        chain_source = 'auth_asym_ids'):
    
#     print '\n---------------------------------------------'
    print ('Getting PDB-chain-gene mapping quality...\n---------------------------')
    pdbMaps = handle_pdbs.get_PDB_gene_maps(mmSeqs=sequenceService_api_result,dfpdb_mapped= dfpdb_mapped)

    outfile = op.join(qspaceDirs['DataOutput_dir'], '002C-{}-PDB_quality.json'.format(query_type))
    with open(outfile, 'w') as f:
        json.dump(pdbMaps, f)
    log.info("Saving PDB-chain-gene mapping...\n\t{}".format(outfile))
        
#     print '\n---------------------------------------------'
    print ('Transfering PDB-chain-gene mapping quality to PDB bioassembly files...\n---------------------------')
    
    bioMaps = handle_bioassemblies.get_BIO_gene_maps(pdbMaps=pdbMaps,
                                                     dfpdb_mapped= dfpdb_mapped, 
                                                     dfbioassembly = dfbioassembly,
                                                     chain_source =chain_source) 

    outfile = op.join(qspaceDirs['DataOutput_dir'], '002C-{}-BIO_quality.json'.format(query_type))
    with open(outfile,'w') as f:
        json.dump(bioMaps, f)
    log.info("Saving PDB Bioassembly-chain-gene mapping...\n\t{}".format(outfile))
   
        
        
    return pdbMaps,bioMaps