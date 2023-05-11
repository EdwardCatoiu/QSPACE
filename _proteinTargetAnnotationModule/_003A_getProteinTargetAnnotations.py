# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def filter_cplx_for_genes_not_in_query(proteinAnnotation, queryGeneList):
    remove_cplx = {}
    genes_found = []
    for index, row in proteinAnnotation.iterrows():
        gs = row.gene_stoich
        if type(gs) == str:
            gs = ast.literal_eval(gs)
        keep = True
        for gene in gs.keys():
            if gene not in queryGeneList:
                keep = False
                log.warn("{} outside of current query --> removed {}".format( gene,index))
                break
        if not keep:
            remove_cplx.update({index : gene})
            continue
        
        genes_found += list(gs.keys())
        
    print ('-----------------------------\nRemoved {} Protein complexes due insufficient query genes'.format(len(remove_cplx)))
    
    print ("-----------------------------\nFound   {} Genes".format(len(set(genes_found))))
    return  proteinAnnotation[proteinAnnotation.index.isin(remove_cplx) == False] , set(genes_found)


def run_003A_get_protein_annotations_ecocyc(df_ecocyc_recipe,queryGeneList):
    print ("\nEcoCyc annotations\n-------------------------------")
    
    gs = {}

    genes_covered = []
    ecocyc_recipe = {}
    for index, row in df_ecocyc_recipe.iterrows():
#         if 'cobraME-enzyme' in row.dropna().keys():
#             ecocyc_recipe.update({row.get('cobraME-enzyme') : ast.literal_eval(row.gene_stoich)})

#         else:
#             ecocyc_recipe.update({index : ast.literal_eval(row.gene_stoich)})
        gs_row  = row.gene_stoich
        if type(gs_row) == str:
            gs_row= ast.literal_eval(gs_row)
        gs.update({index :gs_row})
        genes_covered += list(gs_row.keys())
    genes_covered = set(genes_covered)
#     print "{} Genes in Ecocyc Enzymes".format(len(genes_covered))
#     print "{} Enzymes Used From EcoCyc".format(len(df_ecocyc_recipe))

    df_ecocyc_recipe['gene_stoich'] = pd.Series(gs)
    
    #filter out protein targets that contain genes outside the query list
    df_ecocyc_recipe,genes_covered = filter_cplx_for_genes_not_in_query(df_ecocyc_recipe, queryGeneList=queryGeneList)
    
    return df_ecocyc_recipe,genes_covered
    
    
def run_003A_get_protein_annotations_cobrame(proteinAnnotations_cobraME,queryGeneList):
    print ("\nCobraME iJL1678b-ME  annotations\n-------------------------------")

    dfcobrame = pd.DataFrame()

    genes_covered = []

    gs_series = {}
    for e, gs in proteinAnnotations_cobraME.items():

        if e in updated_me_recipe:
            gs = updated_me_recipe[e]

        gs_series.update({e : gs})


        for g in gs:
            genes_covered += [g]
    genes_covered = set(genes_covered)
#     print "\n{} Genes in ME Model".format(len(genes_covered))
    dfcobrame['gene_stoich'] = pd.Series(gs_series)
#     print "{} Enzymes Used From COBRAME".format(len(dfcobrame))
   
    #filter out protein targets that contain genes outside the query list
    dfcobrame,genes_covered = filter_cplx_for_genes_not_in_query(dfcobrame, queryGeneList=queryGeneList)

    return dfcobrame, genes_covered

def run_003A_get_protein_annotations_unmodeled(dfrepseq, genes_covered):
    print ("\nUnmodelled gene annotations\n-------------------------------")

    dfunmodelled = pd.DataFrame()
    unmodelled_recipe = {}
    gs = {}
    genes_unmodelled = []
    for g in dfrepseq.index.tolist():
        if g not in genes_covered: #these are the geenes in me model or in ecocyc 
            unmodelled_recipe.update({"Unmodelled_MONOMER_{}".format(g) : {g : 1}})
            genes_unmodelled += [g]
    dfunmodelled['gene_stoich'] = pd.Series(unmodelled_recipe)
#     print "\n{} Umodelled Genes".format(len(unmodelled_recipe))
    return dfunmodelled, genes_unmodelled




def run_003A_get_protein_annotations_final(proteinAnnotations_ecocyc,proteinAnnotations_cobrame,proteinAnnotations_unmodeled,dfrepseq):
    dfparent = pd.DataFrame(columns = ['gene_stoich', 'cobraME-enzyme', 'description', 'cobrame', 'unmodelled',
                                       'ecocyc', 'monomer_type', 'k_mer', 'gene_string', 'len_gene_string'])

#     print '----------------------------------------------'
    print ("Compiling Ecocyc Protein Complex Annotations\n------------------------------")
    data = {}
    for c in tqdm(proteinAnnotations_ecocyc.columns):
        data.update({c : proteinAnnotations_ecocyc[c].to_dict()})

#     print "Getting Unmodelled Data"
    print ("Compiling Unmodelled Protein Monomer Annotations\n------------------------------")

    gs = {}
    for index, row in tqdm(proteinAnnotations_unmodeled.iterrows()):
        gs.update({index: row.gene_stoich})
    #     dfparent.loc[index, 'gene_stoich'] = str(row.gene_stoich) 
    data['gene_stoich'].update(gs)

    print ("Compiling CobraME-iJL1678b Protein Complex Annotations\n------------------------------")
    gs = {}
    for index, row in tqdm(proteinAnnotations_cobrame.iterrows()):
    #     print index
        if row.gene_stoich in data['gene_stoich'].values():
            for cplx in proteinAnnotations_ecocyc[proteinAnnotations_ecocyc.gene_stoich ==row.gene_stoich].index:
                data['cobraME-enzyme'].update({cplx : index})
    #             dfparent.loc[cplx , 'cobraME-enzyme'] = index
    #             print cplx, index
            continue

    #         print index, row.gene_stoich
    #         continue
        gs.update({index: row.gene_stoich})
        data['cobraME-enzyme'].update({index : index})

    data['gene_stoich'].update(gs)

    dfparent['gene_stoich'] = pd.Series(data['gene_stoich'])
    for col , series in data.items():
        if col == 'gene_stoich':
            continue
        for index, value in series.items():
            dfparent.loc[index, col] = value

    print ("Determining source of protein complex annotations\n------------------------------")
    for index, row in tqdm(dfparent.iterrows()):
        if 'Unmodelled' in index:
            dfparent.loc[index,'unmodelled'] = True
        if 'cobraME-enzyme' in row.dropna().keys():
            dfparent.loc[index,'cobrame'] = True
        if index in proteinAnnotations_ecocyc.index:
            dfparent.loc[index,'ecocyc'] = True

    print ("Determining protein complex amino acid sequence\n------------------------------")
    genes_needed = []
    for index, row in tqdm(dfparent.iterrows()):
        genes = row.gene_stoich

        if len(genes) == 1:
            dfparent.loc[index,'monomer_type'] = 'mono'
        else:
            dfparent.loc[index,'monomer_type'] = 'hetero'

        dfparent.loc[index,'k_mer'] = np.sum(list(genes.values()))


        string = ''
        error = False
        for g , gs in genes.items():
            for i in range(int(gs)):


                try:
                    if str(dfrepseq.loc[g,'UniProtSeq']) == 'nan':
                        string += dfrepseq.loc[g,'AlleleomeSeq'] + ':' 

                    else:
                        string += dfrepseq.loc[g,'UniProtSeq'] + ':' 
                except KeyError:
                    error = True
                    print ('\t error: {} not in dfrepseq'.format(g))
                    genes_needed += [g]
                    break
    #         print g, gs, dfrepseq.loc[g,'rep_seq'],'\n'


        string = string.rstrip(':')
    #     print type(genes), genes
    #     dfparent.loc[cplx, 'gene_stoich'] = str(genes)

        if error:
            continue
        dfparent.loc[index, 'gene_string'] = string
        dfparent.loc[index, 'len_gene_string'] = len(string)
#     print "{} Enzymes Total".format(len(dfparent))
#     print "{} Enzymes w/error (could not determine aa string, genesnot in genelist)".format(len(dfparent[dfparent.gene_string.isna()]))
#     print "{} Enzymes that can be found from our genelist".format(len(dfparent[dfparent.gene_string.isna()==False]))
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'003A-enzyme_targets_prior_to_BFS.csv') 
    dfparent[dfparent.gene_string.isna()==False].to_csv(outfile)
    log.info("Saving Protein Target Annotations (prior to structure-guided re-annotation...\n\t{}".format(outfile))
#     print "Saving...\n\t> {}".format(outfile)

    return dfparent[dfparent.gene_string.isna()==False], dfparent[dfparent.gene_string.isna()], set(genes_needed)

