# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)



fixed_chains = {'1rc2-assembly3.cif_&_MATCH_1' : {'A-0' : 'A', 'A-1' : 'A-5', 'B-0' : 'A-6', 'B-1' : 'A-7' },
                'CHEA-CPLX_AlphaMulti' : {'B' : 'B', 'C' : 'A' },
#                 '3ulk-assembly1.cif_&_NOBFS_0' : {'A-0' : 'A','A-1' : 'A-2','B-0' : 'B-3','B-1' : 'B-4' },
                '6ag8-assembly2.cif_&_MATCH_0' : {'A-0' : 'A','C-0' : 'A-3'},
                '1n3s-assembly3.cif_&_NOBFS_0' : {'A-0' : 'A','B-0' : 'B','C-0' : 'C','D-0' : 'D','E-0' : 'E','F-0' : 'F-2','G-0' : 'G-2','H-0' : 'H-2','I-0' : 'I-2', 'J-0' : 'J-2', },
                '6erj-assembly2.cif_&_MATCH_2' : {'A-0' : 'A','B-0' : 'A'},
                '1xru-assembly2.cif_&_NOBFS_0' : {'A-0' : 'A','B-0' : 'B-2',},
                '1u60-assembly4.cif_&_MATCH_0' : {'A-0' : 'A','B-0' : 'B','C-0' : 'C-3','D-0' : 'D-3'},
                'CPLX-159_AlphaMulti' : {'B' : 'A','C' : 'B'},
                '3usy-assembly4.cif_&_NOBFS_0' : {'B-0' : 'B','A-0' : 'B-3'},
                '4yxb-assembly1.cif_&_MATCH_0' : {'A-0' : 'A','B-0' : 'B','D-0':'C'},
                '7wgu-assembly1.cif_&_NOBFS_0' : {'A-0' : 'A','B-0' : 'B-2'},                
                'CPLX0-841_AlphaMulti' : {'C' : 'A','B' : 'B'},
                'CPLX0-1669_AlphaMulti' : {'C' : 'A','B' : 'B'},
                'CPLX0-1668_AlphaMulti' : {'B' : 'A','C' : 'B',},
                '7k1r-assembly1.cif_&_NOBFS_0' : {'A-0' : 'A','B-0' : 'B-2','C-0' : 'C-3','D-0' : 'D',},
                'CHEZ-CPLX_AlphaMulti' : {'B' : 'B','C' : 'A',},
                '2vyc-assembly2.cif_&_MATCH_4' : {'A-0' : 'F','B-0' : 'G','C-0' : 'H','D-0' : 'I','E-0' : 'J','F-0' : 'F-3','G-0' : 'G-3','H-0' : 'H-3','I-0' : 'I-3', 'J-0' : 'J-3', },
                '2scu-assembly2.cif_&_MATCH_0' : {'B-0' : 'E-2','E-0' : 'E','A-0' : 'D-2','D-0' : 'D',},
                '6nhl-assembly1.cif_&_NOBFS_0' : {'A-0' : 'A','B-0' : 'B-2'},
                'CPLX-63_AlphaMulti' : {'C' : 'B','B' : 'A',},
                '5xn8-assembly2.cif_&_MATCH_2' : {'D-0': 'D', 'D-1': 'D-3', 'B-0': 'D-5', 'B-1': 'D-7', 'A-1': 'C-7', 'A-0': 'C-5', 'C-1': 'C-3', 'C-0': 'C'},
                '2qfi-assembly1.cif_&_MATCH_0' : {"A-0" : 'B-3'},
                '2qfi-assembly2.cif_&_MATCH_0' : {"A-0" : 'B-3'},
                "3ulk-assembly1.cif_&_NOBFS_0" : {'A-0':'A-3','A-1':'A-4'} ,
                "3ulk-assembly2.cif_&_NOBFS_0" : {'A-0':'A-3','A-1':'A-4'} ,
                "7lkl-assembly1.cif_&_NOBFS_0" : {'A-0':'A', "A-1" : 'A-2', 'B-0':'B', "B-1" : 'B-2', },
                "1qzt-assembly3.cif_&_MATCH_1" : {u'D-0':"D-2", u'A-0':"A", u'B-0':"B", u'C-0':"C-2"},
                "1jqj-assembly3.cif_&_NOBFS_0" : {u'A-0':"A",u'B-0':"B",u'C-0':"C",u'D-0':"D",u'E-0':"E",u'F-0':"A-2",u'G-0':"B-2",u'H-0':"C-2",u'I-0':"D-2",u'J-0':"E-2"},                

               }

def run_004B_QuaternaryProteomeSkeleton(dfqual,
                                              dfbest,
                                              dfrepseq,
                                             ):
    appenderDict = {}
    err_cplx=[]


    remaining_index = set(dfqual.index.tolist()) - set(appenderDict.keys())

    for cplx in tqdm(remaining_index ):
        if cplx in appenderDict:
            continue
    #     if cplx =='ribosome':
    #         continue
    # #     if cplx != 'CPLX0-7985':
    #         continue
#         print cplx
        appender = []

        row = dfqual.loc[cplx]
        match = row.final_match
        if type(match) == str:
            match = ast.literal_eval(match)
    #     dfb = dfbest.loc[structure]

        for pseudo_structure in (match.keys()):


            structure_stoich = match[pseudo_structure]

            try:
                structureType = dfbest.loc[pseudo_structure, "stype"]
            except KeyError:
                print (pseudo_structure)
                continue
            structureId = pseudo_structure.split('_&_')[0].split('.')[0]
    #         print structureId
            #get the structure file and structProp
            sfile = find_sfile(structureId=structureId, stype=structureType)
            if not sfile:
                print (cplx, pseudo_structure)
                raise IOError
            struct_prop = StructProp(ident=structureId, structure_path= sfile, file_type = op.basename(sfile).split('.')[-1])

            #get the base structure Id for the correct Needle File
            structureId_base = pseudo_structure.split('_&_')[0].split('-assembly')[0]
            structure_pdb_quality = dfbest.loc[pseudo_structure,'pdb_quality']
            if type(structure_pdb_quality) == str:
                structure_pdb_quality = ast.literal_eval(structure_pdb_quality)

            for gene in (structure_pdb_quality.keys()):
                chainInfo = structure_pdb_quality[gene]
#                 print (gene, chainInfo)
                geneSeqId = str(dfrepseq.loc[gene,'UniProtId'])
                geneSeq = dfrepseq.loc[gene,'UniProtSeq']
                fastaFile = op.join(qspaceDirs['UniprotSeqsDir'], '{}.fasta'.format(geneSeqId))

                if geneSeqId == 'nan':
                    geneSeqId = 'WT_consensus'
                    geneSeq = dfrepseq.loc[gene,'AlleleomeSeq']
                    fastaFile = op.join(qspaceDirs['AlleleomeSeqsDir'], '{}_{}.fasta'.format(gene,geneSeqId))

                #get the gene SeqProp  
                if not op.exists(fastaFile):
                    print ('no fasta file : {} {} {}'.format(gene, geneSeqId,fastaFile))
                    raise IOError
                    
                record = SeqIO.read(fastaFile, "fasta")
#                 print (geneSeqId)
                seqprop = SeqProp(id =geneSeqId,seq = str(record.seq), sequence_path =fastaFile)
                struct_prop.reference_seq = seqprop
                p = struct_prop.parse_structure()
#                 print (struct_prop.structure_path)
#                 print (fastaFile)

                for chainId in (chainInfo):
                    chainLabel = chainId.split('-')[0]
                    needle_file = op.join(qspaceDirs['SequenceAlignmentDir'], '{}_{}_{}-{}.needle'.format(structureId_base, gene, geneSeqId, chainLabel)
                                         )
#                     print (chainId)
#                     print (needle_file)
#                     needle_file = '../GEMPRO/needle_alignments/{}_{}_{}-{}.needle'.format(structureId_base,gene,geneSeqId,chainLabel)
                    if not op.exists(needle_file):
                        print ('needle file does not exists : {}' .format(needle_file))
                        raise IOError
    #                 print '\n',gene, structureId, chainLabel

                    dfaln = ssbioaln.get_alignment_df_from_file(needle_file)
                    dfaln_structure= dfaln[dfaln.id_a_pos.isna() == False]
    #                 dfaln_structure = dfaln_structure[dfaln_structure.id_a_aa !='X']

                    dfaln['structResNum'] = np.nan

                    try:
                        new_chainId = copy.deepcopy(chainId.split('-')[0])

                        if pseudo_structure in fixed_chains:
                            if chainId in fixed_chains[pseudo_structure]:
                                new_chainId = fixed_chains[pseudo_structure][chainId]
#                                 print 'Using Manual Chains for {} \t {}-->{}'.format(pseudo_structure, chainId, new_chainId)



                        structChain = struct_prop.chains.get_by_id(new_chainId) #use the A-0/B-0 ID provided
                        repchain_structure_mapping = structChain.seq_record.letter_annotations['structure_resnums']
#                         print (new_chainId,repchain_structure_mapping)
                        i = 0
                        for  index_dfaln in (dfaln_structure.index):
                            structaa1_dfaln =   dfaln.loc[index_dfaln,'id_a_aa']
                            if structaa1_dfaln == 'X':
                                i+=1
                                continue
                            try:
                                structResNum = str(repchain_structure_mapping[i])
                                while structResNum in ['inf','nan']:
                                    i +=1
                                    structResNum = str(repchain_structure_mapping[i])
    #                                 print index_dfaln, i

                            except IndexError:
                                print ('missing : {} {} {}'.format(pseudo_structure, new_chainId, gene))
                                continue


                            structaa1 = three_to_one(p.first_model[new_chainId][int(structResNum)].resname)
                            structaa1_dfaln =   dfaln.loc[index_dfaln,'id_a_aa']
                            if structaa1 != structaa1_dfaln:

                                print ('alignment issue: {} \t {}\t {}\t{} \t {}\t {}\t{} \t {}\t {}\t'.format(cplx,structureId, gene, chainId,new_chainId, structResNum,index_dfaln, structaa1, structaa1_dfaln ))
                                continue
    #                             raise AlignmentError
                            dfaln.loc[index_dfaln,'structResNum'] = int(structResNum)
                            i +=1

                    except KeyError:

                        print ('CHAIN DNE' , chainId, chainLabel, new_chainId, gene, cplx, '\t', pseudo_structure)

                        err_cplx +=[cplx]
#                         print "{}".format(sfile)
                        dfaln['structResNum'] = np.nan



                    for i , (index_dfaln, row_dfaln) in enumerate(dfaln.iterrows()):

                        match_type =  row_dfaln.get('type')

                        geneSeqNum = row_dfaln.id_b_pos
                        geneSeqAA = row_dfaln.id_b_aa


                        structChainAA = row_dfaln.id_a_aa
                        structResNum = row_dfaln.structResNum

                        genome_position =  "{}&{}&{}&{}".format(gene,geneSeqId,str(geneSeqNum),geneSeqAA)

                        
                        info = [cplx, gene, geneSeqId, geneSeqNum,geneSeqAA, genome_position,structureId, structure_stoich,match_type, structResNum, structChainAA, chainLabel, chainId, new_chainId]
                        appender.append(info)



        appenderDict.update({cplx : appender})

    global_appender = []
    for k , appender in appenderDict.items():
        global_appender += appender

    column_list=['cplx', 'gene', 'geneSeqId', 'seqNum', 'seqAA1', 'GenomeLocation', 'structureId', 'structureStoich', 'matchType', 
                 'structNum','structAA1', 'structChain','structChainId','structChainId_mod']
    df = pd.DataFrame.from_records(global_appender, columns=column_list )
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '004B-alldata_skeleton.csv') 
    df.to_csv(outfile)
#     print "Saving...\n\t> {}".format(outfile)
    
    
    resNumDict = df.structNum.to_dict()
    chainIdDict = df.structChainId_mod.to_dict()

    chain_residue = {}
    for index, resNum in resNumDict.items():
        if str(resNum) =='nan':
            continue
        chain_residue.update({index : "{}_{}".format(chainIdDict[index], int(resNum))})
    df['Chain_Residue'] = pd.Series(chain_residue)

    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '004B-alldata_skeleton.csv') 
    df.to_csv(outfile)
    log.info("Saving QSPACE proteome...\n\t{}".format(outfile))
    
    return df