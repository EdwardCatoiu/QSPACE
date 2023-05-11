# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def run_003C_oliogmers(proteinTargets,df_structures, query = 'PDB'):

    quality_cutoffs = {'PDB':80,'SWISS': 70}
    quality_columns = {'PDB':'pdb_quality_needle','SWISS': 'pdb_quality'}
    
    cutoff_value = quality_cutoffs[query]
    col_check = quality_columns[query]
    
    
    
    oligomer = proteinTargets[proteinTargets.k_mer.isin([1,'1']) == False]
    oligomer_parent = oligomer.gene_stoich.to_dict()

    for k, v in oligomer_parent.items():
        if type(v) == str:
            v = ast.literal_eval(v)
        oligomer_parent.update({k : v})

#     print 'Found {} annotated oligomers...'.format(len(oligomer_parent))

#     print "Mapping Oligomers to {}...".format(query)
    results = prepBFS.get_data(df= df_structures, 
                                   parent_input = oligomer_parent,
                                   stype = [query],
                                   quality_cutoff = cutoff_value,
                                   quality_source = col_check)
    [df, structure_stoich,structure_quality,output ,  to_bfs,  full, homo_oligo, hetero_oligo ] = results
#     print "-----{} Results-----".format(query)
#     print "Protein-Complexes Mapped to 1 structure : {}".format(len(full))
#     print "Protein-Complexes w/BFS subunit mapping available: {}".format(len(to_bfs))
#     print "Protein-Complexes w/structural evidence of oligomerization : {}".format(len(set(homo_oligo.keys() + hetero_oligo.keys())))
    
    missing = set(oligomer_parent.keys()) - set(to_bfs) - set(full) - set(homo_oligo) - set(hetero_oligo) 
#     print "Protein-Complexes missing: {}".format(len(missing))
    
    return df, structure_stoich,structure_quality,output ,  to_bfs,  full, homo_oligo, hetero_oligo,oligomer_parent 

def run_003C_find_missing_oligos(structuralEvidence,oligomer_parent,dfpdb_and_swiss_structures,df_swiss):
#     print'--------------------------------\nFinding oligomers with missing structural data.....'
    full_1_1_pdb = structuralEvidence['full_pdb']
    full_1_1_swiss = structuralEvidence['full_swiss']
#     print 'Removing SWISS-MODELS with low confidence of Oligomerization.....'
    dellist = []
    full_1_1_pdb = structuralEvidence['full_pdb']
    full_1_1_swiss = structuralEvidence['full_swiss']
#     print len(full_1_1_swiss)
    for cplx, swisslist in full_1_1_swiss.items():

        gs = dfpdb_and_swiss_structures.loc[swisslist[0],'gene_stoichiometry']
        if type(gs) == str:
            gs = ast.literal_eval(gs)

        is_monomer = np.sum(list(gs.values())) == 1

        is_good_model, metrics = oligomerization.check_swiss_model(df_swiss,swiss_id = swisslist[0], total_quality_cutoff = 70,is_monomer=is_monomer)

        if not is_good_model:
            dellist += [cplx]

    for k in dellist:
        del full_1_1_swiss[k]
    
#     print "Structures Found  for {} of {} oligomers".format(len(set(full_1_1_pdb.keys()+ full_1_1_swiss.keys())), len(oligomer_parent))
    df_oligo_missing = oligomer_parent[oligomer_parent.index.isin(list(full_1_1_pdb.keys())+ list(full_1_1_swiss.keys())) == False]
    df_oligo_missing= df_oligo_missing.sort_values(by = 'len_gene_string')
#     print "Structures Missing for {} oligomers ".format(len(df_oligo_missing))
    return df_oligo_missing


def run_003C_sendCplxSeqs_to_AFMultiFolder(df_oligo_missing, 
                                           outdir = qspaceDirs['AlphaMultiSeqDir'],
                                           alphafold_seq_len_cutoff = 2000,
                                           remove_old_seqs= True,
                                          ):
    
    to_alphafold = df_oligo_missing[df_oligo_missing.len_gene_string <= alphafold_seq_len_cutoff]

    print ("Sending {} of {} Protein complexes to Alphafold Multimer .....\n-----------------------------------------".format(len(to_alphafold), len(df_oligo_missing)))
#     print "\t> {}".format(outdir)
    if remove_old_seqs:
        for f in os.listdir(outdir):
            os.remove(op.join(outdir,f))

    # data= df_enzymes_to_alphafold_short#[df_enzymes_to_alphafold_short.source == 'ecocyc']
    for index, row in to_alphafold.iterrows():
    #     print index,
        cplx = SeqProp(id = index, seq = row.gene_string )
        outfile = op.join(outdir, cplx.id)
        cplx.write_fasta_file(outfile = outfile, force_rerun=True)
    #     print len(row.gene_string),'\t', index
    
    too_long_for_alphafold = df_oligo_missing[df_oligo_missing.len_gene_string > alphafold_seq_len_cutoff]
    print ("{} Protein complexes are large than the {} AA sequence cutoff".format(len(too_long_for_alphafold), alphafold_seq_len_cutoff))
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003C-too_long_for_alphafold.csv')
    log.info("Saving protein oligomers > 2000 AAs for manual curation and send to alphafold multimer...\n\t{}".format(outfile))
    too_long_for_alphafold.to_csv(outfile) 
    return to_alphafold
     
    
def run_003C_sendByPartsSeqs_to_AFMultiFolder(df_oligo_missing, 
                                              dfrepseq,
                                              outdir = qspaceDirs['AlphaMultiSeqDir'],
                                              alphafold_seq_len_cutoff = 2000,
                                              remove_old_seqs= False,

                                             ):
    
    infile  = op.join(qspaceDirs['Input_dir'], "003C-alphafoldMultimerByParts.csv")
    if op.exists(infile):
        alphafoldByPartsManual = pd.read_csv(infile, index_col=0) 
        log.info("Importing protein oligomers > 2000 AAs after manual curation ...\n\t{}".format(infile))

#         print'--------------------------------\nImporting partial-oligomers for {} protein-complex .....'.format(len(alphafoldByPartsManual))
#         print'\t> {}'.format(infile)
      
    else:
#         print 'No Partial Enzymes found in file {}'.format(infile)
        return pd.DataFrame()
        
    alphafoldByPartsManual_outsideOfQuery=alphafoldByPartsManual[alphafoldByPartsManual.index.isin(df_oligo_missing.index) == False]
#     print'Partial-oligomers for {} protein-complexes fall outside the current geneList query'.format(len(alphafoldByPartsManual_outsideOfQuery))


    alphafoldByPartsManual=alphafoldByPartsManual[alphafoldByPartsManual.index.isin(df_oligo_missing.index)]
    dfpartial = pd.DataFrame(columns = ['gene_stoich','sequence','seq_len'])
    
    for index, row in alphafoldByPartsManual.iterrows():
        for partial_s in ['s1','s2','s3','s4']:
            if partial_s in row.dropna().keys():


                partial_stoich = ast.literal_eval(row.get(partial_s).replace(' ',''))
                partial_ID = ""
                sequence = ""

                for g,s in partial_stoich.items():
                    partial_ID +='{}_{}__'.format(g,s)
                    seq = dfrepseq.loc[g,'UniProtSeq']
                    for i in range(s):
                        sequence += seq
                        sequence += ":"
    #             if 'b0617' in partial_stoich:
    #                 print g, s,sequence
                dfpartial.loc[partial_ID[0:-1], 'gene_stoich'] = str(partial_stoich)
                dfpartial.loc[partial_ID[0:-1], 'sequence'] = sequence.rstrip(':')
                dfpartial.loc[partial_ID[0:-1], 'seq_len'] = len(sequence.rstrip(':')            )
    dfpartial =dfpartial.sort_values(by = 'seq_len')
    if remove_old_seqs:
        for f in os.listdir(outdir):
            os.remove(op.join(outdir,f))
    
    c = 0 
    for index, row in dfpartial[dfpartial.seq_len <= alphafold_seq_len_cutoff].iterrows():
        cplx = SeqProp(id = index, seq = row.sequence)
        outfile = op.join(outdir, cplx.id)
        cplx.write_fasta_file(outfile = outfile, force_rerun=True)
        c+=1
#     print "Added sequences of {} partial protein-oligomers to {}".format(c, outdir)

    return dfpartial
