# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def run_003B_monomers(proteinTargets,df_structures, query = 'PDB'):

    quality_cutoffs = {'PDB':80,'SWISS': 70}
    quality_columns = {'PDB':'pdb_quality_needle','SWISS': 'pdb_quality'}
    
    cutoff_value = quality_cutoffs[query]
    col_check = quality_columns[query]
    
    
    
    monomers = proteinTargets[proteinTargets.k_mer.isin([1,'1'])]
    monomer_parent = monomers.gene_stoich.to_dict()
    for k, v in monomer_parent.items():
        if type(v) == str:
            v = ast.literal_eval(v)
        monomer_parent.update({k : v})
#     print (monomer_parent)

    print ('Found {} annotated or unmodelled monomers...'.format(len(monomer_parent)))

#     print "Mapping Monomers to {}...".format(query)
    results = prepBFS.get_data(df= df_structures, 
                                   parent_input = monomer_parent,
                                   stype = [query],
                                   quality_cutoff = cutoff_value,
                                   quality_source = col_check)
    [df, structure_stoich,structure_quality,output ,  to_bfs,  full, homo_oligo, hetero_oligo ] = results
#     print "-----{} Results-----".format(query)
#     print "Protein-Complexes Mapped to 1 structure : {}".format(len(full))
#     print "Protein-Complexes w/BFS subunit mapping available: {}".format(len(to_bfs))
#     print "Protein-Complexes w/structural evidence of oligomerization : {}".format(len(set(homo_oligo.keys() + hetero_oligo.keys())))
    
    missing = set(monomer_parent.keys()) - set(to_bfs) - set(full) - set(homo_oligo) - set(hetero_oligo) 
#     print "Protein-Complexes missing: {}".format(len(missing))
    
    return df, structure_stoich,structure_quality,output ,  to_bfs,  full, homo_oligo, hetero_oligo,monomer_parent 


def run_003B_HeteroOligos(structuralEvidence,monomer_parent,quality_pdb,blast_cutoff = 80):
#     print '---------------CASE V (all)---------------------\n1. Hetero - Oligomers from PBD .....'

    
    oa_dict_pdb = structuralEvidence['hetero_pdb']
    o_dict_pdb = structuralEvidence['homo_pdb']
    o_dict_swiss = structuralEvidence['homo_swiss']
    
    
    dfoa_hetero = oligomerization.get_dfhetero(parent = monomer_parent,
                                   real_structure_quality_pdb = quality_pdb ,
                                   oa_dict_pdb =oa_dict_pdb ,
                                   o_dict_pdb=o_dict_pdb,
                                   o_dict_swiss=o_dict_swiss, 
                                   blast_cutoff = blast_cutoff)

    
#     print len(dfoa_hetero), '--->',
    for index, row in dfoa_hetero.iterrows():
        dfoa_hetero.loc[index,'gene_stoichiometry'] = str({row.gene : 1})
    dfoa_hetero = dfoa_hetero[dfoa_hetero.Hetero_kmer_PDB.isna() == False]
#     print len(dfoa_hetero)
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003B-hetero_oligo_evidence_needs_curation.xlsx')
    dfoa_hetero.to_excel(outfile)
    
    log.info("Saving Hetero-oligomer evidence from PDB only for manual curation...\n\t{}".format(outfile))
    return dfoa_hetero
    
    
    
def run_003B_caseV_HeteroOligos(dfpdb_and_swiss_structures,
                                 new_recipes = {},    
                                 df_recipes_changed = pd.DataFrame(columns= ['finalGeneStoich', 'stoichChanged',
                                                                             'method', 'caseNum', 'evidenceSupported',
                                                                             'evidenceIgnored']),
                                 curation_confirmed = False,
                                 monomer_parent = {}
                                 
                                ):
#     
    print ("CASE V (passed)\n------------------------------------------")
    print ('Importing Manual Curation...')
    infile = op.join(qspaceDirs['Input_dir'], '003B-hetero_oligo_evidence_after_curation.xlsx')
    if not op.exists(infile) or not curation_confirmed:
        file_needed = op.join(qspaceDirs['DataOutput_dir'], '003B-hetero_oligo_evidence_needs_curation.xlsx')
        log.warn("Need to uploaded a manual curation file OR set curation_confirmed == False")
        raise  FileNotFoundError(
    errno.ENOENT, os.strerror(errno.ENOENT), infile)
        

    dfoa_hetero_manual = pd.read_excel(infile,index_col = 0 )
    
#     print 'Found manual curation information for {} protein complexes....'.format(len(dfoa_hetero_manual))
    #we need to remove the protein complexes if they are not in our monomer query....
    dfoa_for_query = dfoa_hetero_manual[dfoa_hetero_manual.index.isin(monomer_parent)]
#     print 'The monomer list contains {} of these protein complexes....'.format(len(dfoa_for_query))
    
    for cplx, gs in dfoa_for_query.manual_stoich.dropna().to_dict().items():
        
        
        if type(gs) != dict:
            gs = ast.literal_eval(gs)
#             print gs, cplx
        new_recipes.update({cplx : gs})
    

    #lots of effort for the supplement data...
    for index, row in  dfoa_for_query.iterrows():
        df_recipes_changed.loc[index, "finalGeneStoich"] = row.gene_stoichiometry
        gs = ast.literal_eval(row.gene_stoichiometry )
        df_recipes_changed.loc[index, "stoichChanged"] = np.sum(list(gs.values())) > 1
        df_recipes_changed.loc[index, "caseNum"] = 'V'
        df_recipes_changed.loc[index, "method"] = 'manual curation'

        if np.sum(list(gs.values())) > 1:
            structuralEvidence_dict = {}
            structuralEvidence = ast.literal_eval(row.Homo_kmer_PDB) +ast.literal_eval(row.Homo_kmer_SWISS)

            for structure in structuralEvidence:
                gene_stoich = str(dfpdb_and_swiss_structures.loc[structure,"gene_stoichiometry"])
                qual = dfpdb_and_swiss_structures.loc[structure,"total_quality"]
                if gene_stoich in structuralEvidence_dict:
                    structuralEvidence_dict[gene_stoich].update({structure : qual})
                else:
                    structuralEvidence_dict.update({gene_stoich : {structure : qual}})

            df_recipes_changed.loc[index, "evidenceSupported"] = str(structuralEvidence_dict)
            df_recipes_changed.loc[index, "evidenceQuality"] = str(structuralEvidence_dict)

            structuralEvidenceIgnored_dict = {}
            structuralEvidenceIgnored = list(ast.literal_eval(row.Hetero_kmer_PDB).keys())

            for structure in structuralEvidenceIgnored:
                try:
                    gene_stoich = str(dfpdb_and_swiss_structures.loc[structure,"gene_stoichiometry"])
                except KeyError:
                    log.warn("Unable to find Gene Stoich for {} --> structure evidence ignored".format(structure))
                    continue
                qual = dfpdb_and_swiss_structures.loc[structure,"total_quality"]
                if gene_stoich in structuralEvidenceIgnored_dict:
                    structuralEvidenceIgnored_dict[gene_stoich].update({structure : qual})
                else:
                    structuralEvidenceIgnored_dict.update({gene_stoich : {structure : qual}})

            df_recipes_changed.loc[index, "evidenceIgnored"] = str(structuralEvidenceIgnored_dict)

        else:
            structuralEvidenceIgnored_dict = {}
            structuralEvidenceIgnored = []
            try:

                for col in ['Hetero_kmer_PDB','Homo_kmer_PDB','Homo_kmer_SWISS']:
                    if col in row.dropna().keys() and col =='Hetero_kmer_PDB':
                        structuralEvidenceIgnored += ast.literal_eval(row.get(col)).keys()

                    elif col in row.dropna().keys() and col in ['Homo_kmer_PDB','Homo_kmer_SWISS']:
                        structuralEvidenceIgnored += ast.literal_eval(row.get(col))


                for structure in structuralEvidenceIgnored:
    #                 print structure
                    try:
                        gene_stoich = str(dfpdb_and_swiss_structures.loc[structure,"gene_stoichiometry"])
                        qual = dfpdb_and_swiss_structures.loc[structure,"total_quality"]

                    except KeyError:
                        if row.manual not in structuralEvidence_dict:
                            structuralEvidenceIgnored_dict.update({row.manual : [structure]})
                        else:
                            structuralEvidenceIgnored_dict[row.manual] += [structure]

                        continue
                    if gene_stoich in structuralEvidenceIgnored_dict:
                        structuralEvidenceIgnored_dict[gene_stoich].update({structure : qual})
                    else:
                        structuralEvidenceIgnored_dict.update({gene_stoich : {structure : qual}})

            except:
                structuralEvidenceIgnored_dict = row.manual


            df_recipes_changed.loc[index, "evidenceIgnored"] = str(structuralEvidenceIgnored_dict)   
    return new_recipes, df_recipes_changed


def run_003B_caseII_autoQCQA_SwissHomoOligos(structuralEvidence,
                                             proteinAnnotation,
                                      df_swiss,
                                      real_structure_stoich_swiss,
                                      dfpdb_and_swiss_structures , 
                                      new_recipes = {},    
                                      df_recipes_changed = pd.DataFrame(columns= ['finalGeneStoich', 'stoichChanged',
                                                                                  'method', 'caseNum', 'evidenceSupported',
                                                                                  'evidenceIgnored']),
                                      monomer_parent = {},
                                             
                                ):
    
    oa_dict_pdb = structuralEvidence['hetero_pdb']
    o_dict_pdb = structuralEvidence['homo_pdb']
    o_dict_swiss = structuralEvidence['homo_swiss']
    print ("CASE II (all)\n------------------------------------------")
    print('Auto QCQA for SWISS homo-oligomers...')
    
    df_swiss_results, swiss_only_dict,failed_swiss  = oligomerization.auto_qcqa_swiss(parent = monomer_parent,
                                            o_dict_swiss=o_dict_swiss,
                                            oa_dict_pdb=oa_dict_pdb,                   
                                            o_dict_pdb=o_dict_pdb,
                                            df_swiss=df_swiss,
                                            real_structure_stoich_swiss=real_structure_stoich_swiss,
                                           )
    log.info("Auto SWISS Curation : added new recipe for {} enzymes using only SWISS data".format(len(df_swiss_results.gene_stoich.to_dict())))
    new_recipes.update(df_swiss_results.gene_stoich.to_dict())

    for index, row in df_swiss_results.iterrows():
        df_recipes_changed.loc[index,'finalGeneStoich'] = str(row.gene_stoich)
        df_recipes_changed.loc[index,'stoichChanged'] =  True
        df_recipes_changed.loc[index,'method'] =  'automatic SWISS'
        df_recipes_changed.loc[index,'caseNum'] =  'II'
        gs = ast.literal_eval(str(row.gene_stoich))
        stoich = list(gs.values())[0]
        df_recipes_changed.loc[index,'evidenceSupported'] =  str( {str(gs) : {row.get(stoich) : row.get('total_quality')}})
        df_recipes_changed.loc[index,'evidenceQuality'] =  str( {str(gs) : {row.get(stoich) : row.get('total_quality')}})
#     for case in df_recipes_changed.caseNum.unique():
#         print case , '\t', len(df_recipes_changed[df_recipes_changed.caseNum == case])
#     len(df_recipes_changed)


    for index, info in failed_swiss.items():
        evidence_ignored_dict= {}
        for stoich, structurelist in info.items():

            for structure in structurelist:
                gs = str(dfpdb_and_swiss_structures.loc[structure,'gene_stoichiometry'])
                if gs not in evidence_ignored_dict:
                    evidence_ignored_dict.update({gs : {structure :dfpdb_and_swiss_structures.loc[structure,'total_quality'] }})
                elif structure not in evidence_ignored_dict[gs]:
                    evidence_ignored_dict[gs].update({structure :dfpdb_and_swiss_structures.loc[structure,'total_quality'] })

        df_recipes_changed.loc[index,'method'] =  'automatic SWISS'
        df_recipes_changed.loc[index,'caseNum'] =  'II'
        df_recipes_changed.loc[index,'evidenceIgnored'] =  str(evidence_ignored_dict)


        if 'evidenceSupported' in df_recipes_changed.loc[index].dropna().keys():
            continue
        else:
            df_recipes_changed.loc[index,'finalGeneStoich'] = str(proteinAnnotation.loc[index, 'gene_stoich'])
            df_recipes_changed.loc[index,'stoichChanged'] =  False


    return new_recipes,df_recipes_changed

def run_003B_caseIV_run_HomoOligos(structuralEvidence,
                            dfpdb_structures,
                            monomer_parent = {},
                            
                           ):
    print ("CASE IV (all)\n------------------------------------------")

    oa_dict_pdb = structuralEvidence['hetero_pdb']
    o_dict_pdb = structuralEvidence['homo_pdb']
    o_dict_swiss = structuralEvidence['homo_swiss']
    full_1_1_pdb = structuralEvidence['full_pdb']
    full_1_1_swiss= structuralEvidence['full_swiss']
    keys  = set(o_dict_pdb.keys()) - set(o_dict_swiss.keys())-set(oa_dict_pdb.keys())
    print ("Analyzing {} Structures with ONLY PDB Homo-k-mer Structures...".format(len(keys)))

    o_dict  = {}
    for cplx in keys:
        o_dict.update({cplx : o_dict_pdb[cplx]})

    structure_results = {}
    for k , model_list in o_dict.items():
        for s_id in model_list:
    # #         df_swiss.loc[swiss_id]

            kmer = dfpdb_structures.loc[s_id,'gene_stoichiometry']
            if type(kmer) == str:
                kmer = ast.literal_eval(kmer)
            kmer = np.sum(list(kmer.values()))

            if k not in structure_results:
                structure_results.update({k : {kmer : [s_id]}})
            elif kmer not in structure_results[k]:
                structure_results[k].update({kmer : [s_id]})
            elif s_id not in structure_results[k][kmer]:
                structure_results[k][kmer] += [s_id]

    total_quality = {}           
    df_structure_results = pd.DataFrame(    columns = ['gene','calculated','manual','total_quality'])
    for cplx, results in structure_results.items():
    #     if len(results) > 1 :
    #         print cplx
    #         print results
    #         break
    #     for kmer, swiss_id_list in results.items():
    #         if len(swiss_id_list) > 1:
    #             print cplx, kmer, swiss_id_list


        for kmer, swiss_id_list in results.items():
    #         swiss_id = swiss_id_list[0]
            if cplx not in total_quality:
                total_quality.update({cplx : dfpdb_structures[dfpdb_structures.index.isin(swiss_id_list)].total_quality.to_dict()})
            else:
                total_quality[cplx].update(dfpdb_structures[dfpdb_structures.index.isin(swiss_id_list)].total_quality.to_dict())

            df_structure_results.loc[cplx, kmer] = str(swiss_id_list)

        if cplx in full_1_1_pdb:
            if cplx in full_1_1_swiss:
                full11 = {1 : full_1_1_pdb[cplx] + full_1_1_swiss[cplx]}
            else:
                full11 = {1 : full_1_1_pdb[cplx]}
        elif cplx in full_1_1_swiss:
            full11 = {1 : full_1_1_swiss[cplx]}
        else:
            continue




        for kmer, swiss_id_list in full11.items():
    #         swiss_id = swiss_id_list[0]
            if cplx not in total_quality:
                total_quality.update({cplx : dfpdb_structures[dfpdb_structures.index.isin(swiss_id_list)].total_quality.to_dict()})
            else:
                total_quality[cplx].update(dfpdb_structures[dfpdb_structures.index.isin(swiss_id_list)].total_quality.to_dict())

            df_structure_results.loc[cplx, kmer] = str(swiss_id_list)#         is_good_model, metrics = manual_curation.check_swiss_model(df_swiss,swiss_id = swiss_id, total_quality_cutoff = 70)
    #         for metric, value in metrics.items():
    #             df_swiss_results.loc[cplx, metric] = value


    df_structure_results['total_quality'] = pd.Series(total_quality) 
    gs = {}
    parent = monomer_parent
    for index, row in df_structure_results.iterrows():
        row = row.drop('total_quality').dropna()
        if len(row) ==1:
            df_structure_results.loc[index, 'calculated'] = row.keys()[0]
            df_structure_results.loc[index, 'manual'] = row.keys()[0]
        df_structure_results.loc[index, 'gene'] = list(parent[index].keys())[0]

    col_stoich = []
    for col in df_structure_results.columns:
        try:
            int(col)
            col_stoich += [col]
        except ValueError:
            continue
    col_stoich.sort()
    df_structure_results = df_structure_results[["gene", 'calculated', 'manual'] + col_stoich + ["total_quality"]]

    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003B-df_homo_pdb_only_needs_curation.xlsx')
    df_structure_results.to_excel(outfile)
    log.info("Saving Homo-oligomer evidence from PDB only for manual curation...\n\t{}".format(outfile))
    
    return df_structure_results


def run_003B_caseIV_handle_HomoOligos(dfpdb_and_swiss_structures,
                                 new_recipes = {},    
                                 df_recipes_changed = pd.DataFrame(columns= ['finalGeneStoich', 'stoichChanged',
                                                                             'method', 'caseNum', 'evidenceSupported',
                                                                             'evidenceIgnored']),
                                 curation_confirmed = False,
                                 monomer_parent = {}
                                 
                                ):
    
    print ("CASE IV (passed)\n------------------------------------------")
#     print'Importing Manual Curation...'
    infile = op.join(qspaceDirs['Input_dir'], '003B-df_homo_pdb_only_after_curation.csv')
    if not op.exists(infile) or not curation_confirmed:
        file_needed = op.join(qspaceDirs['DataOutput_dir'], '003B-df_homo_pdb_only_needs_curation.xlsx')
        log.warn("Need to perform manual curation or set curation_confirmed == False")
        raise (FileNotFoundError, 'Please confirm manually the oligomerization states in...\n\t> {}\nAnd save them here\n\t> {}'.format(file_needed,infile))
#         return new_recipes,df_recipes_changed # this should be uncommented if not using manual curation, comment out the raise error statement in previous line

    
    df_structure_results_input = pd.read_csv(infile, index_col=0)

#     print 'Found manual curation information for {} protein complexes....'.format(len(df_structure_results_input))
    #we need to remove the protein complexes if they are not in our monomer query....
    df_structure_results_manual = df_structure_results_input[df_structure_results_input.index.isin(monomer_parent)]
#     print 'The monomer list contains {} of these protein complexes....'.format(len(df_structure_results_manual))
    
    #new recipes
    for index, row in df_structure_results_manual.iterrows():
        if int(row.manual) == 1:
            continue
        new_recipes.update({index : {row.gene : int(row.manual)}})
        
    #supplemental DF
    for index, row in df_structure_results_manual.iterrows():
        df_recipes_changed.loc[index,'finalGeneStoich'] =    str({row.gene : int(row.manual)})
        df_recipes_changed.loc[index,'stoichChanged'] =  int(row.manual) !=1
        df_recipes_changed.loc[index,'method'] =  'manual curation'
        df_recipes_changed.loc[index,'caseNum'] =  'IV'


        evidence_ignored_dict = {}
        for col, structureList in row.dropna().items():
            try:
                int(col)
            except ValueError:
                continue

            if str(col) == str(row.manual):
                continue    
    #         print col

            gs = str({row.gene : int(col)})
            qualDict = {}
            for structure in ast.literal_eval(structureList):
                try:
                    qualDict.update({structure: dfpdb_and_swiss_structures.loc[structure,'total_quality'] })
                except KeyError:
                    continue        
                evidence_ignored_dict.update({gs : qualDict})
        df_recipes_changed.loc[index,'evidenceIgnored'] =  str(evidence_ignored_dict)


        evidence_for = {}
        for col, structureList in row.dropna().items():
            try:
                int(col)
            except ValueError:
                continue

            if str(col) != str(row.manual):
                continue    
    #         print col, 'passed'
            gs = str({row.gene : int(col)})
            qualDict = {}
            for structure in ast.literal_eval(structureList):
                try:
                    qualDict.update({structure: dfpdb_and_swiss_structures.loc[structure,'total_quality'] })
                except KeyError:
                    continue
            evidence_for.update({gs : qualDict})
        df_recipes_changed.loc[index,'evidenceSupported'] =  str(evidence_for)
        df_recipes_changed.loc[index,'evidenceQuality'] =  str(evidence_for)

    return new_recipes,df_recipes_changed

def find_ignored_evidence(row, manual_key = 'manual_stoich', gene = False):

    evidence_ignored_dict = {}

    for key, value in row.dropna().items():
        try:
            int(key)
        except ValueError:
            continue

        if key in [str(row.get(manual_key)),int(row.get(manual_key))]:
            continue

        evidence_ignored = value
#         print key , evidence_ignored

        qualDict = {}
        
        if gene:
            gs = str({gene : int(key)})
        else:
            gs = str({list(ast.literal_eval(row.gene_stoich).keys())[0] : int(key)})
    
        geneStoichIgnored  = gs
        for structure in ast.literal_eval(evidence_ignored):
            row_qual = row.get('total_quality')
            if type(row_qual) == str:
                row_qual = ast.literal_eval(row_qual)
            
            qual = row_qual[structure]
            qualDict.update({structure :  qual})
    
        evidence_ignored_dict.update({str(geneStoichIgnored) : qualDict})
        
    return evidence_ignored_dict


def run_003B_caseI_homoOligos_PDB_and_SWISS(structuralEvidence,
                                      df_swiss,
                                      df_pdb,
                                      dfpdb_and_swiss_structures,
                                      real_structure_stoich_swiss ,
                                      oldProteinAnnotation ,
                                      monomer_parent= {},
                                      new_recipes = {},
                                      df_recipes_changed = pd.DataFrame(columns= ['finalGeneStoich','stoichChanged',
                                                                                  'method', 'caseNum','evidenceSupported',
                                                                                  'evidenceIgnored']),
                                     ):
    print ('Case I (passed) and Case III (all)\n-----------------------------')
    oa_dict_pdb = structuralEvidence['hetero_pdb']
    o_dict_pdb = structuralEvidence['homo_pdb']
    o_dict_swiss = structuralEvidence['homo_swiss']
    full_1_1_pdb = structuralEvidence['full_pdb']
    full_1_1_swiss= structuralEvidence['full_swiss']
    
    df_to_manual,df_calculated = oligomerization.auto_qcqa_swiss_and_pdb(parent = monomer_parent,                         
                                                                         o_dict_swiss=o_dict_swiss,
                                                                         oa_dict_pdb=oa_dict_pdb,                   
                                                                         o_dict_pdb=o_dict_pdb,
                                                                         df_swiss=df_swiss,
                                                                         df_pdb=df_pdb,
                                                                         dfbest_trim = dfpdb_and_swiss_structures,
                                                                         real_structure_stoich_swiss=real_structure_stoich_swiss,
                                                                        )
                                                                         
    ############df caclulcated ##############
    c = 0        
    for index, row in df_calculated.iterrows():
        try:
            old_gs  = oldProteinAnnotation.loc[index,'gene_stoich']
            if type(old_gs) == str:
                old_gs = ast.literal_eval(old_gs)
        except KeyError:
            continue

        k_mer = row.manual
        new_gs = {list(old_gs.keys())[0] : k_mer}
        new_recipes.update({index : new_gs})
        c+=1
        
        df_recipes_changed.loc[index, 'finalGeneStoich'] = str(new_gs)
        df_recipes_changed.loc[index, 'method'] = 'automatically PDB+SWISS'
        df_recipes_changed.loc[index, 'caseNum'] = 'I'

        evidence_ignored_dict = find_ignored_evidence(row, manual_key='manual', gene =list(old_gs.keys())[0] )
        df_recipes_changed.loc[index, 'evidenceIgnored'] = str(evidence_ignored_dict)
        df_recipes_changed.loc[index, 'stoichChanged'] = True


        if int(row.manual) == 1:
            df_recipes_changed.loc[index, 'stoichChanged'] = False
            print (index)
            continue
    #     df_recipes_changed.loc[index, 'evidenceSupported'] = row.get(row.manual)
        quals_for = {str(new_gs) : {}}
        for structure in ast.literal_eval(row.get(int(row.manual))):
            quals_for[str(new_gs)].update({structure : row.total_quality[structure]})


        df_recipes_changed.loc[index, 'evidenceQuality'] = str(quals_for)
        
    log.info("AUTO Curation added new recipe for {} enzymes using PDB & SWISS data".format(c))
    
    for index, row in df_to_manual.iterrows():
        df_to_manual.loc[index,'gene_stoich'] = str(oldProteinAnnotation.loc[index,"gene_stoich"])
        
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003B-df_homo_pdb_and_swiss_needs_curation.csv')
    df_to_manual.to_csv(outfile)
    log.info("Saving Homo-oligomer evidence from PDB and SWISS for manual curation...\n\t{}".format(outfile))
    
    return new_recipes,df_recipes_changed

def run_003B_caseIII_handle_HomoOligos_PDB_and_SWISS(new_recipes = {},
                                             df_recipes_changed = pd.DataFrame(columns= ['finalGeneStoich','stoichChanged',
                                                                                         'method', 'caseNum','evidenceSupported',
                                                                                         'evidenceIgnored']),
                                             curation_confirmed = False,
                                             monomer_parent = {}
                                             
                                            ):
    
    print ('Case III (passed)\n-----------------------------------')
#     print'Importing Manual Curation...'
    infile = op.join(qspaceDirs['Input_dir'], '003B-df_homo_pdb_and_swiss_after_curation.csv')
    if not op.exists(infile) or not curation_confirmed:
        file_needed = op.join(qspaceDirs['DataOutput_dir'], '003B-df_homo_pdb_and_swiss_manual_curation.csv')
        log.warn("Need to perform manual curation or set curation_confirmed == False")
        raise (FileNotFoundError, 'Please confirm manually the oligomerization states in...\n\t> {}\nAnd save them here\n\t> {}'.format(file_needed,infile))
#         return new_recipes,df_recipes_changed # this should be uncommented if not using manual curation, comment out the raise error statement in previous line
           
            
    df_manual_input = pd.read_csv(infile, index_col=0)
 
#     print 'Found manual curation information for {} protein complexes....'.format(len(df_manual_input))
    #we need to remove the protein complexes if they are not in our monomer query....
    df_to_manual = df_manual_input[df_manual_input.index.isin(monomer_parent)]
#     print 'The monomer list contains {} of these protein complexes....'.format(len(df_to_manual))
    
    for index, row in df_to_manual.iterrows():
    
        new_gs = str({list(ast.literal_eval(row.gene_stoich).keys())[0] : int(row.manual_stoich)})
        df_recipes_changed.loc[index, 'finalGeneStoich'] = new_gs
        df_recipes_changed.loc[index, 'method'] = 'manual curation'
        df_recipes_changed.loc[index, 'caseNum'] = 'III'

        evidence_ignored_dict = find_ignored_evidence(row)
        df_recipes_changed.loc[index, 'evidenceIgnored'] = str(evidence_ignored_dict)
        df_recipes_changed.loc[index, 'stoichChanged'] = True


        if int(row.manual_stoich) == 1:
            df_recipes_changed.loc[index, 'stoichChanged'] = False
            continue
        df_recipes_changed.loc[index, 'evidenceSupported'] = row.get(row.manual_stoich)


        quals_for = {new_gs : {}}
        for structure in ast.literal_eval(row.get(str(row.manual_stoich))):
            quals_for[new_gs].update({structure : ast.literal_eval(row.total_quality)[structure]})


        df_recipes_changed.loc[index, 'evidenceQuality'] = str(quals_for)


        new_recipes.update({index : {list(ast.literal_eval(row.gene_stoich).keys())[0] : int(row.manual_stoich)}})
    
    return new_recipes,df_recipes_changed

def run_003B_rename_enzymes(new_recipes,dfparent,dfrepseq):
    naming_prefix = {2 : 'DI',
                 3 : 'TRI',
                 4 : 'TETRA',
                 5 : 'PENTA',
                 6 : 'HEXA',
                 7 : 'HEPTA',
                 8 : 'OCTA',
                 9 : 'NONA',
                 10 : 'DECA',
                }
#     print '--------------------------------------\nRe-naming oligomer enzymes...'
    new_data = {}
    translation = {}
    count = 1
    for old_enzyme, new_gs in new_recipes.items():
        row_data = dfparent.loc[old_enzyme].to_dict()
        k_mer = np.sum(list(new_gs.values()))
    #     break
        row_data.update({'gene_stoich' :new_gs })
        row_data.update({'k_mer' :k_mer })

        string = ''
        error = False
        for g , gs in new_gs.items():
            for i in range(int(gs)):
                try:
                    string += dfrepseq.loc[g,'UniProtSeq'] + ':' 
                except KeyError:
                    error = True
                    print ('\t error: ', g)
        #         print g, gs, dfrepseq.loc[g,'rep_seq'],'\n'
        string = string.rstrip(':')
        row_data.update({'gene_string' : string })
        row_data.update({'len_gene_string' : len(string)  })

        if k_mer in naming_prefix:
            new_enzyme = old_enzyme.replace("MONO", naming_prefix[k_mer])
        else:
            new_enzyme = old_enzyme.replace("MONO", "{}-".format(k_mer))
    #     print "{}\t{}-->{}".format(count,old_enzyme,new_enzyme)
        new_data.update({new_enzyme : row_data})
        translation.update({new_enzyme : old_enzyme})
        count+=1
        
        
    print ("Original Enzymes : ", len(dfparent))
    dfparent2 = dfparent[dfparent.index.isin(new_recipes) == False]
    dfparent2_wrong = dfparent[dfparent.index.isin(new_recipes) == True]
    print ("Incorrect Monomers : ",len(dfparent2_wrong))
    print ("Identified New Oligomers : ",len(new_data))

    for index, data in tqdm(new_data.items()):
        dfparent2.loc[index] = data

#     print "Total Enzymes : ",len(dfparent2)
   
        
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003B-enzyme_targets_with_new_oligomers.csv')
    log.info("Saving NEW protein (target) complex annotations...\n\t{}".format(outfile))
    dfparent2.to_csv(outfile)
    return new_data,translation,dfparent2
    
    