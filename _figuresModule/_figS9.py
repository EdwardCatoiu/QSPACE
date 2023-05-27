
# -*- coding: utf-8 -*-
from .utils import *


def figS9_A(dfalpha_manual,fig = False, ax = False, save = False,histtype= 'stepfilled',    alpha = 0.5, outfile = False):
    if not fig:
        fig,ax =plt.subplots()
        
    bins = np.linspace(0,2000,21)
    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good','Auto']) == False]
    ax.hist(data.sequence_length.values, 
               bins = bins,
               color = 'red',
               histtype= histtype,
               alpha = alpha, 
               label ='Fail : {}'.format(len(data))
              )   
    
    
    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Auto'])]
    ax.hist(data.sequence_length.values,
               bins = bins,
               color = u'#274f37',
               histtype= histtype,
               alpha = alpha, 
               label ='AutoPass : {}'.format(len(data))
              )

    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good'])]
    ax.hist(data.sequence_length.values, 
               bins = bins,
               edgecolor =u'#274f37',
               hatch = '//',
               color = 'white',
               histtype= histtype,
               alpha = alpha, 
               label ='ManualPass : {}'.format(len(data))
              )

    ax.set_xlabel('Sequence Length (AAs)')

    ###########
    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS9A-003D-AFMultiModel_vs_seqLength.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax


def figS9_B(dfalpha_manual,fig = False, ax = False, save = False,histtype= 'stepfilled',alpha = 0.5, outfile = False):
    if not fig:
        fig,ax =plt.subplots()
        
    bins = np.linspace(0,100,51)
    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good','Auto']) == False]
    ax.hist(data.plddt.values, 
               bins = bins,
               color = 'red',
               histtype= histtype,
               alpha = alpha, 
               label ='Fail : {}'.format(len(data))
              )

    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Auto'])]
    ax.hist(data.plddt.values,
               bins = bins,
               color = u'#274f37',
               histtype= histtype,
               alpha = alpha, 
               label ='AutoPass : {}'.format(len(data))
              )

    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good'])]
    ax.hist(data.plddt.values, 
               bins = bins,
               edgecolor =u'#274f37',
               hatch = '//',
               color = 'white',
               histtype= histtype,
               alpha = alpha, 
               label ='ManualPass : {}'.format(len(data))
              )

    ax.set_xlabel('pLLDT')

    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS9B-003D-AFMultiModel_vs_ModelScore.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax


def figS9_C(dfalpha_manual,fig = False, ax = False, save = False,histtype= 'stepfilled',alpha = 0.5, outfile = False):
    if not fig:
        fig,ax =plt.subplots()
 
    bins = np.linspace(0,1,51)
    ######
    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good','Auto']) == False]
    ax.hist(data.alphafold1_model_score.values, 
               bins = bins,
               color = 'red',
               histtype= histtype,
               alpha = alpha,  
               label ='Fail : {}'.format(len(data))
              )

    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Auto'])]
    ax.hist(data.alphafold1_model_score.values, 
               bins = bins,
               color = u'#274f37',
               histtype= histtype,
               alpha = alpha, 
               label ='AutoPass : {}'.format(len(data))
              )


    data = dfalpha_manual[dfalpha_manual.QCQA_verdict.isin(['Good'])]
    ax.hist(data.alphafold1_model_score.values,
               bins = bins,
               edgecolor = u'#274f37',
               hatch = '//',
               color = 'white',
               histtype= histtype,
               alpha = alpha,  
               label ='ManualPass : {}'.format(len(data))
              )

    ax.set_xlabel('Alphafold Multimer Score')
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS9C-003D-AFMultiModel_vs_pLLDT.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax
    


