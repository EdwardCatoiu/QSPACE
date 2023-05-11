
# -*- coding: utf-8 -*-
from .utils import *

from matplotlib_venn import venn3,venn3_unweighted

def figS8_A(proteinTargets,fig = False, ax = False, save = False):
    if not fig:
        fig,ax =plt.subplots()
        
    set1 = set(proteinTargets[proteinTargets.cobrame == True].index.tolist())
    set2 = set(proteinTargets[proteinTargets.ecocyc == True].index.tolist())
    set3 = set(proteinTargets[proteinTargets.unmodelled == True].index.tolist())
    v = venn3_unweighted([set1, set2, set3], ('COBRAME','ECOCYC','UNMODELLED'), ax = ax)
    ax.set_title("Sources of Protein Complexes")

    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS8A-003A-SourcesOfProteinComplexes.png')
    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax


def figS8_B(proteinTargets,fig = False, ax = False, save = False):
    if not fig:
        fig,ax =plt.subplots()

    monomer = proteinTargets[proteinTargets.k_mer == 1].index.tolist()
    cobrame = proteinTargets[proteinTargets.get('cobrame') == True].index.tolist()
    ecocyc = proteinTargets[proteinTargets.get('ecocyc') == True].index.tolist()

    set1 = set(monomer)
    set2 = set(cobrame)
    set3 = set(ecocyc)

    v = venn3_unweighted([set1, set2, set3], ("Monomer",'COBRAME',  'ECOCYC'))
    ax.set_title("Annotated Oligomerization of Protein Complex ")
    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS8B-003A-AnnotatedOligomerizationfProteinComplexes.png')

    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    return fig,ax
    
    
def figS8_C(hetero_oligo_pdb, homo_oligo_pdb,homo_oligo_swiss,fig = False, ax = False, save = False):
    if not fig:
        fig,ax =plt.subplots()

    set1 = set(hetero_oligo_pdb.keys())
    set2 = set(homo_oligo_pdb.keys())
    set3 = set(homo_oligo_swiss.keys())

    venn3_unweighted([set1, set2, set3], ('Hetero-k-mer (PDB)','Homo-k-mer (PDB)',  'Homo-k-mer (SWISS)'))
    ax.set_title("Structural Evidence of Oligomerization\nfor {} 'MONOMER' Enzymes".format(len(set1.union(set2).union(set3))))
    # c=venn3_circles([set1, set2, set3], linestyle='-', linewidth=1, color="k")
    
    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS8C-003B-StructuralEvidenceOfOligomerization.png')
    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')

    return fig, ax


def figS8_D(df_recipes_changed,fig = False, ax = False, save = False):
    if not fig:
        fig,ax =plt.subplots()


    i = 0 
    xlabel = []
    for case in ['I','II','III','IV','V']:
        dfr = df_recipes_changed[df_recipes_changed.caseNum ==case]
        xlabel += [case]
        xlabel += [case]


        ax.bar(i, len(dfr[dfr.stoichChanged == False]),color ='grey')
        i+=1

        ax.bar(i, len(dfr[dfr.stoichChanged == True]), color = 'orange')
        xy = np.array((i, len(dfr[dfr.stoichChanged == True])))
        ax.annotate(len(dfr[dfr.stoichChanged == True]), xy = xy + np.array([0.1,10]))
        i+=1

    ax.set_xticklabels('')
    ax.set_xticks(np.array(range(10) )+ 0.4)
    ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
    ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    
    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS8D-003B-ResultsOfStoichReannotation.png')
    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')