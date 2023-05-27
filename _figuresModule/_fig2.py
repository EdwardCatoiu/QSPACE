# -*- coding: utf-8 -*-
from .utils import *


def fig2_B(dfqual, fig = False, ax = False, save = False,outfile = False):
    if not fig:
        fig,ax =plt.subplots()
        
    
    structure_list = []
    for index, row in dfqual.iterrows():
        structure_match = ast.literal_eval(row.final_match)
        for structure in structure_match:
            structure_list +=[structure.split('_&_')[0]]

    structure_dict = {}
    for structure in structure_list :
        if 'assembly' in structure:
            stype = 'PDB'
        elif 'AlphaMulti' in structure:
            stype = 'ALPHAFOLD_MULTIMER'
        elif 'AF-' in structure and 'model' in structure:
            stype= 'ALPHAFOLD'
        elif 'ECOLI' in structure and 'clean_residues' in structure:
            stype= 'ITASSER'
        else:
    #         print structure
            stype= 'SWISS'

        structure_dict.update({structure : stype})
    
    df = pd.DataFrame()
    df['stype'] = pd.Series(structure_dict)

    data_structure_source = {}
    for stype in ['PDB','ALPHAFOLD','ALPHAFOLD_MULTIMER','ITASSER','SWISS']:
        data_structure_source.update({stype : len(df[df.stype == stype])})
    
    ax.set_ylim(0,1.2*np.max(list(data_structure_source.values())))
    xticklabel = []
    for i , stype in enumerate(['PDB','SWISS','ITASSER', 'ALPHAFOLD', 'ALPHAFOLD_MULTIMER',  ]):
        ax.bar(i, data_structure_source[stype], color = 'k', alpha = 0.7, lw = 3)
        ax.annotate(data_structure_source[stype], xy = (i-0.2,data_structure_source[stype] *1.1), size = 14) 
        
        xticklabel += [stype]
    ax.set_xticks(np.array(list(range(5))))
    ax.set_xticklabels(xticklabel, rotation = 90)
    
#     ax.set_yticklabels('')

    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig2B-004A-SourcesOfStructuresUsed.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax



def fig2_C(dfqual, proteinTargetsOld, proteinTargetsNew, fig = False, ax = False, save = False,outfile = False):
   
    monomers = []
    oligomers = []
    subunits = []
    new_structure = []
    new_recipe = []
    old_recipe = []
    for index, row in tqdm(dfqual.iterrows()):
        gene_stoich = proteinTargetsNew.loc[index,'gene_stoich']
        if type(gene_stoich) == str:
            gene_stoich = ast.literal_eval(gene_stoich)
        if np.sum(list(gene_stoich.values())) == 1:
            monomers += [index]
            continue

        else: 
            oligomers += [index]

    for index, row in tqdm(dfqual.iterrows()):
        gene_stoich = proteinTargetsNew.loc[index,'gene_stoich']
        if type(gene_stoich) == str:
            gene_stoich = ast.literal_eval(gene_stoich)
           
        match_dict = row.final_match
        if type(match_dict) == str:
            match_dict = ast.literal_eval(match_dict)
            
        if np.sum(list(gene_stoich.values())) == 1:
            continue

        if np.sum(list(match_dict.values())) != 1:
            subunits += [index]
            continue 

        else:

            if 'AlphaMulti' in list(match_dict.keys())[0]:
                new_structure += [index]
                continue


            try:
#                 gene_stoich_old = ast.literal_eval(proteinTargetsOld.loc[index,'gene_stoich'])
                gene_stoich_old = proteinTargetsOld.loc[index,'gene_stoich']
                if type(gene_stoich_old) == str:
                    gene_stoich_old = ast.literal_eval(gene_stoich_old)
        #         if gene_stoich == gene_stoich_old:
                old_recipe +=[index]
#                 print index
    #         else:
    #             new_recipe += [index]
            except KeyError:
                new_recipe += [index]
    
    data_enzymetype = {"Monomer":len(monomers),
       "Oligomer":len(oligomers),
       "Old Recipe (2a.i-Complete)":len(old_recipe),
       'Old Recipe (2a.i-Subunit)':len(subunits),
       'New Recipe (2a.ii)':len(new_recipe),
       'New Structure (2a.iii)':len(new_structure),
       }
    
    #####begin figure
    if not fig:
        fig,ax =plt.subplots()
        
    xticklabels = []
    for i, label in enumerate(['Monomer', 'Oligomer', 'Old Recipe (2a.i-Complete)','Old Recipe (2a.i-Subunit)', 'New Recipe (2a.ii)','New Structure (2a.iii)' ]):
        value = data_enzymetype[label]
        ax.bar(i,value)
        ax.annotate(str(value), xy = (i-0.3, value * 1.1),size = 14)
        xticklabels += [label]
        
    ax.set_ylim(0,1.2*np.max(list(data_enzymetype.values())))

    ax.set_xticks(np.linspace(0, len(data_enzymetype)-1,len(data_enzymetype)))
    ax.set_xlim(-0.8,6)
    ax.set_xticklabels(xticklabels, rotation = -90, size = 12)#,horizontalalignment  = 'center')
    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    ax.set_ylabel('Proteins Found',size = 15)
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig2C-004A-TypesOfMatchesToStructures.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax 
    

    
    
def fig2_D(dfqual, proteinTargetsNew, fig = False, ax = False, save = False, outfile = False):
    #####begin figure
    if not fig:
        fig,ax =plt.subplots()

    homo_oligo = []
    hetero_oligo = []

    for index, row in tqdm(dfqual.iterrows()):
        gene_stoich = proteinTargetsNew.loc[index,'gene_stoich']
        if type(gene_stoich) == str:
            gene_stoich = ast.literal_eval(gene_stoich)
            
        if len(gene_stoich) == 1:
            homo_oligo += [np.min([np.sum(list(gene_stoich.values())), 33])]
        else:
            hetero_oligo += [np.min([np.sum(list(gene_stoich.values())), 33])]
    
    n1, bins1, p1 = ax.hist(homo_oligo, bins = np.linspace(2,35,34), alpha = 0.7, color = 'k', label = "Homo-oligomer")
    n2, bins2, p2 = ax.hist(hetero_oligo, bins = np.linspace(2,35,34), alpha = 0.7, color = 'orange', label = "Hetero-oligomer")
    ax.set_yscale('log')
    ax.legend(loc = 'upper right')
    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    ax.set_xlabel('Oligomer Size\n(number of genes)',size = 15)
    ax.set_ylabel('Protein Oligomers',size = 15)
    ax.set_xticks([0,2,5,10,15,20,25,30,33.6])
    ax.set_xticklabels([0,2,5,10,15,20,25,30,"  >33"])
    max_y = np.max([np.max(n1), np.max(n2)])
    ax.set_ylim(0.8,3*max_y)

    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig2D-004A-OligomerStateOfProteinComplexes.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')

    return fig, ax 