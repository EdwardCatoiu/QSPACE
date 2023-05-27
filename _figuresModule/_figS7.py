
# -*- coding: utf-8 -*-
from .utils import *


def figS7_A(data,fig = False, ax = False, save = False, legend = False,label_text = True, outfile = False):
    """
    Data must be = [genes_swiss, genes_pdb,genes_itasser,genes_alpha]. Each genes_{source} must bet a set.
    """
    
    [genes_swiss, genes_pdb,genes_itasser,genes_alpha] = data
    if not fig:
        fig,ax =plt.subplots()


    labels = venn.get_labels((genes_swiss, genes_pdb,genes_itasser,genes_alpha))
    for k , v in labels.items():
        if v =='0':
            labels.update({k : ''})

    fig, ax = venn.venn4(labels, 
                          names = ["SWISS-MODEL", "PDB", "I-TASSER","AlphaFold"], 
                          ax = ax,
                          fig = fig, 
                          legend = legend,
                         label_text = label_text,)
    ax.set_title('Genes represented in protein structures')

    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7A-002D-GenesRepresentedInProteinStructures.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax

def figS7_BCDE(dfpseudo,fig = False, ax = False, save = False, legend = True, query = 'ITASSER',outfile = False):
    if not fig:
        fig,ax =plt.subplots()
    
    bins = np.linspace(0,100,51)

    quality = dfpseudo.total_quality.values
    
    if query =='ITASSER':
        color ='#FAF49B'
        label = 'I-TASSER ({})'.format(len(quality))
#         outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7B-001D-PseudoStructurQuality_ITASSER.png')

    elif query =='ALPHAFOLD':
        color ='#F9BEC0'
        label = 'AlphaFold ({})'.format(len(quality))
#         outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7C-001D-PseudoStructurQuality_ALPHAFOLD.png')
        
    elif query =='PDB':
        color ='#AECEEA'
        label = 'PDB ({})'.format(len(quality))
#         outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7E-002D-PseudoStructurQuality_PDB.png')

      
    elif query =='SWISS':
        color ='#AFE0B2'
        label = 'SWISS ({})'.format(len(quality))
#         outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7D-001D-PseudoStructurQuality_SWISS.png')

    ax.hist(quality,bins=bins,histtype='stepfilled',density=True,color=color, label = label )
    if legend:
        ax.legend(loc = 'upper left', fontsize = 11)
    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7-001D-PseudoStructurQuality_{}.png'.format(query))

        
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax

def figS7_F(data,fig = False, ax = False, save = False, legend = True, outfile = False):
    if not fig:
        fig,ax =plt.subplots()
    
    bins = np.linspace(0,100,51)
    a= 0.7

    for query in ['ITASSER','ALPHAFOLD','PDB','SWISS']:
        data_to_plot  = data[query].total_quality.tolist()
        if query =='ITASSER':
            color ='#FAF49B'
            label = 'I-TASSER ({})'.format(len(data_to_plot))

        elif query =='ALPHAFOLD':
            color ='#F9BEC0'
            label = 'AlphaFold ({})'.format(len(data_to_plot))

        elif query =='PDB':
            color ='#AECEEA'
            label = 'PDB ({})'.format(len(data_to_plot))

        elif query =='SWISS':
            color ='#AFE0B2'
            label = 'SWISS ({})'.format(len(data_to_plot))
            
        ax.hist(data_to_plot,bins=bins,histtype='stepfilled',alpha = a, density=True,color=color, label = label)
        
    if legend:
        ax.legend(loc = "upper left")
    
   
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS7D-002D-PseudoStructurQuality_All_best.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax