# -*- coding: utf-8 -*-

from .utils import *
log = logging.getLogger(__name__)


def run_001E(input_quality,dfseq,query_type, dfswiss_model_metrics = pd.DataFrame(),alphafold_is_better = alphafold_is_better):

    
    if query_type in ['ITASSER','ALPHAFOLD']:
        structure_info = pd.DataFrame()
        
    elif query_type=='SWISS':
        structure_info = dfswiss_model_metrics
        if len(structure_info) == 0:
            raise IOError('Use the SWISS - Metrics')

        
    
    
#     print '\n----------------------------{}-------------------------------------'.format(query_type)
#     print 'Making Pseudo-Structures...'
    input_stoich_chains_filtered = pseudo.remake_chains_from_quality(input_quality)
    pseudo_stoich_filtered,pseudo_stoich_chains_filtered  = pseudo.find_pseudo_matches(input_stoich_chains_filtered)

#     print 'Getting Pseudo-Structure Quality...'
    pseudo_quality_filtered = pseudo.get_quals(pseudo_stoich_chains = pseudo_stoich_chains_filtered,
                                               input_quality = input_quality
                                              )

    df,ssbio_errors = pseudo.make_dataframe_pseudo_structures(pseudo_stoich_filtered, pseudo_quality_filtered,df_gene_length=dfseq)
    df3 = pseudo.assign_structure_quality(df,df_gene_length=dfseq)
    df4 = pseudo.assign_quality_ranks(df3)
    # dfpdb_structure_info = pd.read_csv('../data/raw/002-PDBstructureInfo.csv',index_col=0)

    best_structures_list = pseudo.find_best_structures(df4,query_type=[query_type],structure_info = structure_info)
    dfbest =  df4[df4.index.isin(best_structures_list[query_type])]

    ###### remove the monomers where alphafold is clearly better
    if query_type=='ITASSER':
        index_remove = []
        for index, row in dfbest.iterrows():
            if index.split('_&')[0] in alphafold_is_better.keys():
                index_remove +=[index]
#                 print index.split('_ECOLI')[0],
        dfbest=dfbest[dfbest.index.isin(index_remove) == False]

    glist = []
    for gdi in dfbest.gene_stoichiometry:
        if type(gdi) == str:
            gdi = ast.literal_eval(gdi)
        glist += gdi.keys()
#     print "{} Genes accounted for in the {} structures".format(len(set(glist)),len(dfbest))

    outfile = op.join(qspaceDirs['DataOutput_dir'],'001F-all_structures_{}.csv'.format(query_type))
    df4.to_csv(outfile)
#     print ">{} ".format(outfile)
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'001F-best_structures_{}.csv'.format(query_type))
    dfbest.to_csv(outfile)
#     print ">{} ".format(outfile)
    return df4,dfbest
