
# -*- coding: utf-8 -*-
from .utils import *


def figS6_AB(data,fig = False, ax = False, save = False, query= 'genes'):
    """
    Data must be a list of length 4: [textService_UniToPDB, textService_PDBinfo, sequenceService_UniProt, sequenceService_Alleleome]
    """
    
    if not fig:
        fig, ax = plt.subplots()

    ##### Panels A-B ######
    # data for Uniprot text
    
    [textService_UniToPDB, textService_PDBinfo, sequenceService_UniProt, sequenceService_Alleleome] = data
    uniprotId_textService_genes_w_structures = set(textService_UniToPDB.uniprot_id.unique())
    uniprotId_textService_PDBs_mapped = set(textService_PDBinfo.pdb_entry.unique())

    #data for UniProt Seqs
    api_result = sequenceService_UniProt
    uniprot_sequenceService_genes_w_structures  = set(api_result.keys())
    uniprot_sequenceService_PDBs_mapped = []
    for gene, mapped_pdbs in api_result.items():
        uniprot_sequenceService_PDBs_mapped += mapped_pdbs.keys()
    uniprot_sequenceService_PDBs_mapped = set(uniprot_sequenceService_PDBs_mapped)

    #data for Alleleome Seqs
    api_result = sequenceService_Alleleome
    alleleome_sequenceService_genes_w_structures  = set(api_result.keys())
    alleleome_sequenceService_PDBs_mapped = []
    for gene, mapped_pdbs in api_result.items():
        alleleome_sequenceService_PDBs_mapped += mapped_pdbs.keys()
    alleleome_sequenceService_PDBs_mapped = set(alleleome_sequenceService_PDBs_mapped)


    if query =='genes':
        set1 = uniprotId_textService_genes_w_structures
        set2 = uniprot_sequenceService_genes_w_structures
        set3 = alleleome_sequenceService_genes_w_structures

        venn3_unweighted([set1, set2, set3], set_labels = ['Uniprot IDs','Uniprot Seqs','Alleleome Seqs'], ax = ax)
        ax.set_title('Genes w/PDB structures')
        outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS6A-002A-GenesWithStructuresFromAPIs.png')

    elif query =='structures':
        
        set1 = uniprotId_textService_PDBs_mapped
        set2 = uniprot_sequenceService_PDBs_mapped
        set3 = alleleome_sequenceService_PDBs_mapped

        venn3_unweighted([set1, set2, set3], set_labels = ['Uniprot IDs','Uniprot Seqs','Alleleome Seqs'], ax = ax)
        ax.set_title('PDB structures identified')
        outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS6B-002A-StructuresIdentifiedFromAPIs.png')

    else:
        raise (KeyError, 'query must be "genes" or "structures"')

    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax










def figS6_C(dfbioassembly,fig = False, ax = False, save = False):

    if not fig:
        fig, ax = plt.subplots()
    
    
    num_bios = {}
    for pdb_entry in tqdm(dfbioassembly.pdb_entry.unique()):
        dfb = dfbioassembly[dfbioassembly.pdb_entry == pdb_entry]
        num_bios.update({pdb_entry : len(dfb)})



    
    ax.hist(num_bios.values(), bins = np.linspace(1,15,15), color = 'grey', alpha = 0.5)
    ax.set_yscale('log')
    ax.set_ylim(0.5,100000)
    ax.set_xticks(np.array(range(16)) + 0.5)
    ax.set_xticklabels(range(16))

    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    ax.set_xlabel('Number of BioAssembly Files', size = 14)
    ax.set_ylabel('Number of PDB Entries', size = 14)
    
    if save:
        outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS6C-002B-NumberOfBioassembliesPerPDB.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax

def figS6_D(dfbioassembly,fig = False, ax = False, save = False):
    if not fig:
        fig, ax = plt.subplots()


    num_oligomer = {34 : 0}

    for i in range(1,500):
        if i <= 30:
        #     print i, '\t', len(dfbioassembly[dfbioassembly.oligomeric_count ==i])
            num_oligomer.update({i :len(dfbioassembly[dfbioassembly.oligomeric_count.isin([i, str(i)])]) })
        else:
            num_oligomer[34] += len(dfbioassembly[dfbioassembly.oligomeric_count.isin([i, str(i)])])

    # num_oligomer.update({34 :len(dfbioassembly[dfbioassembly.oligomeric_count >= 31]) })

    for i , num in num_oligomer.items():
        if i == 1:
            color = 'grey'
            alpha = 0.5
        else:
            color ='#243C68'
            alpha = .8
        ax.bar(i, num, width = 0.9, color = color, alpha = alpha)
    ax.set_yscale('log')

    xticks = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34]) + 0.5


    ax.set_xticks(xticks)

    xticklabels = ['', 1, '','', '', 5,'','','','', 10, '','','','', 15, '','','','', 20, '','','','', 25,'', '','','', 30, r'$\geq$31',]

    ax.set_xticklabels(xticklabels)
    # ax.set_xticklabels((range(35)))
    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    ax.set_xlabel('Oligomeric State', size = 14)
    ax.set_ylabel('BioAssembly Files', size = 14)
    ax.set_ylim(0.5,100000)
    ax.set_xlim(0,36)

    if save:
        outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS6D-002B-OligomericStateOfBioassemblies.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')

    return fig, ax
