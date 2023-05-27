
# -*- coding: utf-8 -*-
from .utils import *



def figS2(dfqual,fig = False, ax = False, save = False, outfile = False):
    if not fig:
        fig,ax =plt.subplots()


    cplx_stype = {}
    for stype in ['ALPHAFOLD','ALPHAFOLD_MULTIMER','PDB','SWISS','ITASSER']:
        dfqual_s  = dfqual[dfqual.get('{}'.format(stype)).isna() == False]
        cplx_stype.update({stype : set(dfqual_s.index.tolist())})
    
    labels = venn.get_labels(list(cplx_stype.values()))
    for k ,v in labels.items():
        if v =='0':
            labels.update({k : ''})

    fig, ax = venn.venn5(labels  , names=list(cplx_stype.keys()), fig = fig, ax = ax)
    ax.set_title('Proteins are mapped to structures from various databases')

    fig.set_figheight(6)
    fig.set_figwidth(6)

    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS2-004A-ProteinComplexesMappedToStructures.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax