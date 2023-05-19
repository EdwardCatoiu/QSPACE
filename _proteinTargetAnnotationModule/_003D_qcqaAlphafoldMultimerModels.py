# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)
import ssbio.protein.sequence.utils.alignment as ssbioaln


def result_in_query(cplxfile, queryAF ):
        was_queried = []
        for query in queryAF:
            was_queried += [cplxfile.find(query)]
        
        if np.max(was_queried) == -1:
            return False
        
        else:
            return True
            
def run_003D_autoQCQA_AFMultimerModels(queryAF,
                                       dfseq,
                                       alphaFoldMultimer_dir = qspaceDirs['alphafoldMultimerDir'],
                                      ):
    
            
        
        
    print ("Running QCQA on Alphafold Multimers...\n---------------------------------")
    
    dfalpha2 = pd.DataFrame(columns = ['alphafold1_model_score', 'alphafold1_ipTM', 'alphafold1_pTM', 'plddt',
                                       'pdb_file', 'score_file', 'dfparent_index', 'sequence_length', 'QCQA_verdict',
                                       'PDB', 'SWISS', 'Comments'])
    c = 0
    for cplxfile in tqdm(os.listdir(alphaFoldMultimer_dir)):
        if '.result' not in cplxfile:
            continue
            
        # we should not be checking all of our AF results that fall outside of our query.
        if cplxfile.split('.result')[0].split('-MEM-')[0] not in queryAF:
            log.warn("Ignoring AF-Multimer model {}. Found results but not in current QSPACE query".format(cplxfile))
            continue
        cplx = cplxfile.split('.result')[0]
        cplxfolder = op.join(alphaFoldMultimer_dir, cplxfile)
    #     break

        found_scores = False
        found_pdb= False
        for f in os.listdir(cplxfolder):
            if 'scores' in f and 'rank_1' in f:
                found_scores = True
                score_file = op.join(cplxfolder,f)
            if 'scores' in f and 'rank_001' in f:
                found_scores = True
                score_file = op.join(cplxfolder,f)

            if '.pdb' in f and 'rank_1' in f:
                pdb_file = op.join(cplxfolder,f)
                found_pdb  = True
            if '.pdb' in f and 'rank_001' in f:
                pdb_file = op.join(cplxfolder,f)
                found_pdb  = True
            if found_pdb and found_scores:
                break

    #     print score_file
    #     print pdb_file,'\n'
        if not found_scores or not found_pdb:
            log.warn("Could not find AF-Multi .pdb or scores.json for {}".format(f))
            continue

        with open(score_file, 'rb') as f:
            scores= json.load(f)

        ipTM = scores['iptm']
        pTM = scores['ptm']
        model_score = 0.8*ipTM + 0.2*pTM
        plddt = np.mean(scores[u'plddt'])


        dfalpha2.loc[cplx, 'alphafold1_model_score'] = model_score
        dfalpha2.loc[cplx, 'alphafold1_ipTM'] = ipTM
        dfalpha2.loc[cplx, 'alphafold1_pTM'] = pTM
        dfalpha2.loc[cplx,'plddt'] =  plddt

        dfalpha2.loc[cplx, 'pdb_file'] = pdb_file
        dfalpha2.loc[cplx, 'score_file'] = score_file

        dfalpha2.loc[cplx,'dfparent_index'] =  cplx.split('-MEM')[0]
        c+=1
       
    
    #get the sequences from the AF Multimer Model files
    sequence_dict = {}
    chain_qual_dict = {}
    chain_stoich_dict = {}
    
    print ('Confirming gene-chain seq quality for {} AF multimer .pdb files\n---------------------------------'.format(len(dfalpha2)))
    
    indexlist = dfalpha2.index.tolist()
    for index in tqdm(indexlist):
        row = dfalpha2.loc[index]
        pdb_file = row.pdb_file

        s = StructProp(ident=index, structure_path= pdb_file, file_type = 'pdb')
        p= s.parse_structure()
        for chain in p.first_model.get_chains():
            chain_seq= ''
            for residue in chain.get_residues():
                chain_seq +=three_to_one(residue.resname)
            if index not in sequence_dict:
                sequence_dict.update({index : {chain.id : chain_seq}})
            elif chain.id not in sequence_dict[index]:
                sequence_dict[index].update({chain.id : chain_seq})
              
            
            ######### {AF-model : {gene : {chain : [plddt > 50, plddt]}}}
            for gene, row_gene in dfseq.iterrows():
                if chain_seq not in [row_gene.UniProtSeq, row_gene.AlleleomeSeq]:
                    continue
                    
                if index not in chain_qual_dict:
                    chain_qual_dict.update({index : {gene : {chain.id : [row.plddt > 50, row.plddt] }}})
                elif gene not in chain_qual_dict[index]:
                    chain_qual_dict[index].update({gene : {chain.id : [row.plddt > 50, row.plddt] }})
                elif chain.id not in chain_qual_dict[index]:
                    chain_qual_dict[index][gene].update( {chain.id : [row.plddt > 50, row.plddt]})
                
                if index not in chain_stoich_dict:
                    chain_stoich_dict.update({index : {gene : 1}})
                elif gene not in chain_stoich_dict[index]:
                    chain_stoich_dict[index].update({gene : 1})
                else:
                    chain_stoich_dict[index][gene] += 1
                    
    
    for index, row in tqdm(dfalpha2.iterrows()):
        total_length = 0
        for chain , seq in sequence_dict[index].items():
            total_length += len(seq)
        dfalpha2.loc[index,'sequence_length'] = total_length
    dfalpha2 = dfalpha2.sort_values(by = 'sequence_length')
    for index , row in dfalpha2[dfalpha2.alphafold1_model_score >= 0.8].iterrows():
        dfalpha2.loc[index, 'QCQA_verdict'] = 'Auto'

#     print "{} AF Multimer Models need QCQA".format(len(dfalpha2[dfalpha2.QCQA_verdict.isna()]))
#     print "{} of {} AF Multimer Models need manual QCQA".format(len(dfalpha2[dfalpha2.QCQA_verdict.isna()]),len(dfalpha2))
    outfile = op.join(qspaceDirs['DataOutput_dir'], '003D-alphafold_Multimer_Results_needs_manual_curation.csv')
    log.info('Saving AF multimer dataframe that needs manual curation of model w/scores <0.80 ....\n\t{}'.format(outfile))
#     print "Saving....\n\t> {}".format(outfile)
    dfalpha2.to_csv(outfile)
    return dfalpha2, chain_qual_dict, chain_stoich_dict

            
def run_003D_loadManualQCQA_AFMultimerModels(queryAF):
    
    infile_auto = op.join(qspaceDirs['DataOutput_dir'], '003D-alphafold_Multimer_Results_needs_manual_curation.csv')

    infile = op.join(qspaceDirs['Input_dir'], '003D-alphafold_Multimer_Results_after_manual_curation.csv')
    if not op.exists(infile):
        log.info("No manual curation QCQA dataframe found for AF multimers. Using only auto-qcqa results (AF-M score > 0.8)\n\t{}".format(infile_auto))

#         print 'No Manual QCQA file for AF multimers.\n\t> {}\nUsing Auto QCQA only...\n\t>{}'.format(infile,infile_auto)
        dfalpha = pd.read_csv(infile_auto,index_col = 0)
    else:
#         print 'Using Manual QCQA file for AF multimers.\n\t> {}\n'.format(infile)
        dfalpha = pd.read_csv(infile,index_col = 0)
        log.info("Using manual curation QCQA dataframe for AF multimers\n\t{}".format(infile))

        
    dfalpha2 = pd.DataFrame(columns = dfalpha.columns.tolist())
    
    for index, row in dfalpha.iterrows():
        if index.split('.result')[0].split('-MEM-')[0] not in queryAF:
            log.warn("Ignoring AF-Multimer model {}. Found manual QCQA results but not in current query".format(index))
            continue
        
        
#         if not result_in_query(cplxfile=index, queryAF= queryAF ):
#             continue
        dfalpha2.loc[index] = row
#     
#     print "{} AF Multimer Models fall outside of current query".format(len(dfalpha))
#     print "{} AF Multimer Models were found".format(len(dfalpha2))
    return dfalpha2

def run_003D_pseudoStructures_AF_best(dfalpha,
                                      geneStoich_dict,
                                      quality_dict):
    
    dfalpha_passed = dfalpha[dfalpha.QCQA_verdict.isin(['Good','Auto'])]

    dfalpa_best = pd.DataFrame()


    for index, row in dfalpha_passed.iterrows():
        key  = index + '_AlphaMulti'
        dfalpa_best.loc[key,'stype'] = 'ALPHAFOLD_MULTIMER'
        dfalpa_best.loc[key,'gene_stoichiometry'] = str(geneStoich_dict[index])
        dfalpa_best.loc[key,'pdb_quality'] = str(quality_dict[index])
        dfalpa_best.loc[key,'identical_structures'] = str([key])
        dfalpa_best.loc[key,'AA_length'] = row.sequence_length
        dfalpa_best.loc[key,'total_quality'] = row.plddt
        dfalpa_best.loc[key,'quality_rank'] = 1

        dfalpa_best.loc[key,'alphafold1_model_score'] = row.alphafold1_model_score
        dfalpa_best.loc[key,'alphafold1_ipTM'] = row.alphafold1_ipTM
        dfalpa_best.loc[key,'alphafold1_pTM'] = row.alphafold1_pTM
        dfalpa_best.loc[key,'plddt'] =  row.plddt

    outfile = op.join(qspaceDirs['DataOutput_dir'], '003D-best_structures_AlphafoldMultimer.csv')
    log.info('Saving AF multimer pseudo-structure information dataframe for {} models passing QCQA ....\n\t{}'.format(len(dfalpa_best),outfile))
    dfalpa_best.to_csv(outfile)
    
    return dfalpa_best


def run_003D_needle_alignment_AFMulti(AFchain_quality_dict,dfalphamulti,dfrepseq,force_rerun = True):
    
    structure_error = []
    error_needle = []
    
    expected_needle_files_dict ={}
    for structure, gene_chain_info in AFchain_quality_dict.items():
        expected_needle_files = []
        for gene , chaininfo in gene_chain_info.items():
            for chain in chaininfo.keys():
                expected_needle_files += ["{}_AlphaMulti_{}_{}-{}.needle".format(structure, gene, dfrepseq.loc[gene,'UniProtId'],chain)] 
        expected_needle_files_dict.update({structure : expected_needle_files})
       
    
    for structure, needle_file_list in expected_needle_files_dict.items():
        structure_id = "{}_AlphaMulti".format(structure)
        structure_file_path = dfalphamulti.loc[structure, 'pdb_file']
        a_structure = StructProp(ident = structure_id,structure_path=structure_file_path,file_type = 'pdb')
        try:
            parsed = a_structure.parse_structure()
        except (AttributeError, IndexError, PDBConstructionException):
            log.warn( "{} Structure Error".format(structure_id))
            structure_error +=[pdb_id]
            continue

        
        for needle_file in needle_file_list:
            if op.exists(op.join(qspaceDirs['SequenceAlignmentDir'], needle_file)) and not force_rerun:
                continue
                
            chain_id = needle_file.split('-')[-1].split('.needle')[0]
            gene = needle_file.split('AlphaMulti_')[-1].split('_')[0]
            seqId = needle_file.split('AlphaMulti_{}_'.format(gene))[-1].split('-')[0]
            
            if 'consensus' in needle_file:
                g_seq = dfrepseq.loc[gene,'AlleleomeSeq']
            else:
                g_seq = dfrepseq.loc[gene,'UniProtSeq']
            
            if g_seq is None or str(g_seq) == 'nan' :
#             print ( '\t Error : No Sequence >>>> {}' .format(needle_file))
                error_needle +=[needle_file]
                continue
        
        
            try:
                chain_prop = a_structure.chains.get_by_id(chain_id)
                chain_seq_record = chain_prop.seq_record
            except (KeyError,AttributeError):
    #             print ( '{}-{} Chain Does Not Exist'.format(pdb_id, chain_id))
                error_needle +=[needle_file]
                continue
            if not chain_seq_record:
    #             print  ( "Chain Sequence Not Parsed : {}-{}".format(pdb_id, chain_id))
                error_needle +=[needle_file]
                continue 
            if len(chain_seq_record) == 0:
    #             print ( "Not a real chain :  {}-{}".format(pdb_id, chain_id))
                error_needle +=[needle_file]
                continue       
            
            #SEQUENCES 
            struct_seq = str(chain_prop.seq_record.seq)
            if struct_seq.count('X') > 0:
                #eliminate the HOH X in the chain seq
                struct_seq = "".join([''.join(aa) for aa in str(chain_prop.seq_record.seq) if aa != 'X'])
                if len(struct_seq) == 0:
                    error_needle +=[needle_file]
                    continue
    #         print '\t Test >>>> {}'.format( needle_file)

            outfile_path = op.join(qspaceDirs['SequenceAlignmentDir'], needle_file)
            if not op.exists(outfile_path) or force_rerun:
                if type(g_seq) == str:
                    aln_file = ssbioaln.run_needle_alignment(seq_a = struct_seq,
                                                             seq_b = g_seq,
                                                             outfile = needle_file, 
                                                             outdir = qspaceDirs['SequenceAlignmentDir'],
                                                             force_rerun=force_rerun)
                    
                    log.info("Aligning...\t{}".format(needle_file))
    return  structure_error,error_needle