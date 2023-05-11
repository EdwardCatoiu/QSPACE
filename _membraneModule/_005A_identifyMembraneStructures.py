# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def run_005A_IdentifyMembrane(dfrepseq,
                              dfall,
                              dfiml1515= op.join(qspaceDirs['Input_dir'], "005A-iML1515_GP_subsystems_nbt.3956-S7.csv"),
                          dfecocyc_loc = op.join(qspaceDirs['Input_dir'], "005A-EcocycSmartTable-All-genes-of-E.-coli-K-12-substr.-MG1655.txt"),
                              uniprotSequenceDir = qspaceDirs['UniprotSeqsDir'],
                             ):
    
    print ('Finding all potential membrane structures...\n----------------------------------')
    
    if type(dfiml1515) == str:
        dfiml1515 = pd.read_csv(dfiml1515)
        
    if type(dfecocyc_loc) == str:
        dfecocyc_loc = pd.read_csv(dfecocyc_loc,sep = '\t')
        
    dfall = dfall[dfall.sfileChain_Residue.isna() == False]
    
    
    membraneGenesGO, membraneStructuresGO,GO_terms_data  = get_membrane_from_GO(dfrepseq = dfrepseq,
                                                                                      dfall=dfall ,
                                                                                      sequenceFolder=uniprotSequenceDir)
    
    
    
    membraneGenesIml, membraneStructuresIml = get_membrane_from_iML1515(dfiml1515=dfiml1515,
                                                                              dfall=dfall)
    
    membraneGenesUNI, membraneStructuresUNI = get_membrane_from_UniProt(dfrepseq=dfrepseq,
                                                                              dfall=dfall,
                                                                              sequenceFolder=uniprotSequenceDir)
    
    
    membraneGenesECO, membraneStructuresECO = get_membrane_from_Ecocyc(dfall=dfall,
                                                                             dfecocyc_loc = dfecocyc_loc)
    
    
    data = {"GO_genes" :list(membraneGenesGO) ,
            "GO_structures" : list(membraneStructuresGO),
            "IML_genes" : list(membraneGenesIml),
            "IML_structures" :list(membraneStructuresIml) ,
            "UNI_genes" :list(membraneGenesUNI) ,
            "UNI_structures" :list(membraneStructuresUNI) ,
            "ECO_genes" : list(membraneGenesECO),
            "ECO_structures" : list(membraneStructuresECO),
           }
    
    
    potential_membrane_structures = list(membraneStructuresGO.union(membraneStructuresIml).union(membraneStructuresUNI).union(membraneStructuresECO))
#     print '----------------------------------\nFound {} potential membrane structures...'.format(len(potential_membrane_structures))

    outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures.json")
    with open(outfile, 'w') as f:
        json.dump(potential_membrane_structures, f)
    log.info("Saving potential membrane structures......\n\t{}".format(outfile))

#     outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures_all_sources.json")
#     with open(outfile, 'w') as f:
#         json.dump(data, f)
#     log.info("Saving potential membrane structures All sources......\n\t{}".format(outfile))
    
#     outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures_ECOCYC.json")
#     with open(outfile, 'w') as f:
#         json.dump(membraneStructuresECO, f)
#     log.info("Saving potential membrane structures Ecocyc......\n\t{}".format(outfile))
    
#     outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures_GeneOntology.json")
#     with open(outfile, 'w') as f:
#         json.dump(membraneStructuresGO, f)
#     log.info("Saving potential membrane structures Gene Ontology......\n\t{}".format(outfile))
    
#     outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures_iML1515.json")
#     with open(outfile, 'w') as f:
#         json.dump(membraneStructuresIml, f)
#     log.info("Saving potential membrane structures iML1515......\n\t{}".format(outfile))
    
#     outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-potential_membrane_structures_UniProt.json")
#     with open(outfile, 'w') as f:
#         json.dump(membraneStructuresUNI, f)
#     log.info("Saving potential membrane structures UniProt......\n\t{}".format(outfile))
    
    return potential_membrane_structures, data                            
