# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)




def run_007C_mapLTEEmutants(dfalldata,
                           infile = op.join(qspaceDirs['Input_dir'], '007C-LTEE_mutations_input.csv'),
                           ):
    
    print ("Mapping LTEE mutations\n-----------------------------")

    dfmut = pd.read_csv(infile, index_col = 0 )

    global_mut_dict = {}
    checked = []

    for gene in tqdm(dfmut.get('Blattner Number').unique()):
        if gene in checked or gene not in dfalldata.gene.unique():
            continue    
        dfmut_gene = dfmut[dfmut.get('Blattner Number') == gene]
        dfg = dfalldata[dfalldata.gene== gene]
        if len(dfg) == 0:
            print ('No Gene : {}'.format( gene))

            continue

        for index, row in dfmut_gene.iterrows():
            aaRes = row.get(u'AA_residue')
            dfzindex = row.get('dfz_index')
    #         gene = row.get(u'Blattner Number')
    #         print aaRes, row.detail
            dfg_seq = dfg[dfg.seqNum.isin([aaRes,float(aaRes), int(aaRes), str(aaRes)])]

            if len(dfg_seq) == 0:
                log.info("Failed to map mutation | {} {} at dfz index {}".format(gene,row.detail, dfzindex ))

                continue
            found = True
            if  dfg_seq.seqAA1.values[0] != row.AA_orig:
                found = False
    #             dfg = dfalldata[dfalldata.gene== gene]
                dfg_seq = dfg[dfg.seqNum == dfzindex +1]
                if len(dfg_seq) == 0 :
                    print ('--26--No Residue : {} {} --{}'.format(aaRes, gene,dfzindex))

                    continue
                if dfg_seq.seqAA1.values[0] == row.AA_orig:
                    found = True
            if not found:
                log.info("Failed to map mutation | {} {} at index {} in input dataframe".format(gene,row.detail, index ))
                continue
        #     print index, '\t',row.detail,'\t', len(dfg),'\t', found



            for globIndex in dfg_seq.index.tolist():
                if globIndex not in global_mut_dict:
                    global_mut_dict.update({globIndex : [index]})
                else:
                    global_mut_dict[globIndex] += [index]

        checked += [gene]
    #     break

    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007C-StructuralProteome_LTEE_mutations.json') 
    with open(outfile,'w') as f:   
        json.dump(global_mut_dict, f)
    log.info("Saving QSPACE Proteome - LTEE mutations details...\n\t{}".format(outfile))
    
    
    
    glob_mut_gs = {}
    for index, mutIndexList in tqdm(global_mut_dict.items()):
        gsDict = {}
        for mutIndex in mutIndexList:
            aa1 = dfmut.loc[mutIndex,'AA_orig']
            aa2 = dfmut.loc[mutIndex,'AA_final']
            gs = grantham_score(aa1,aa2)[0]

            if "*" in [aa1, aa2]:
                gs = 'TERM'

                if "*" == aa1:
                    gs = 'ELONG'
    #             print index

            gsDict.update({"{}>{}".format(aa1,aa2) : gs})

        glob_mut_gs.update({index : gsDict})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007C-StructuralProteome_LTEE_Grantham.json') 
    with open(outfile,'w') as f:   
        json.dump(glob_mut_gs, f)
    log.info("Saving QSPACE Proteome - LTEE mutations Grantham Score...\n\t{}".format(outfile))
    
    return global_mut_dict,glob_mut_gs
        