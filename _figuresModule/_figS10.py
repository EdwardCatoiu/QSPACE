
# -*- coding: utf-8 -*-
from .utils import *


def figS10_A(data,fig = False, ax = False, save = False):
    if not fig:
        fig,ax =plt.subplots()
    
    set1 = set(data['ProteinTarget'].keys())-set(data['FullMatch'].keys())
    set2 = set(data['FullMatch'].keys())
    set3 = set(data['BFSAvailable'].keys())

    v = venn3_unweighted([set1, set2, set3], ('Single Stucture (Low Quality)','Single Stucture (High Quality)','Multi-Structure'), ax = ax)

    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS10A-003D-AFMultiModel_vs_seqLength.png')
    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax


def figS10_BE(dfmatch_for_supp,dfqual_for_supp, save = False):
    color_palette = sns.color_palette('RdYlGn', 101
                                 ).as_hex()
    sns.palplot(color_palette)
    
    
    final_numlist = [[],[],[],[]]
    final_qualist = [[],[],[],[]]
    final_cplxlist = [[],[],[],[]]

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig,ax =plt.subplots()
    fig.set_figheight(5) 
    fig.set_figwidth(15)
    gs1 = gridspec.GridSpec(1,3 ,height_ratios= [1], width_ratios=[1,1,1])
    gs1.update( wspace=0.35, hspace = 0.2)
    ax0 = plt.subplot(gs1[0])
    ax1 = plt.subplot(gs1[1])
    ax2 = plt.subplot(gs1[2])


    for i, (cplx, row_final_match) in enumerate(dfqual_for_supp.iterrows()):
        final_qual = row_final_match.final_qual                                                      
        final_num_structures = row_final_match.num_pdbs

        dfm = dfmatch_for_supp[dfmatch_for_supp.cplx == cplx]
        dfm  = dfm.sort_values(by = ['num_pdbs','match_quality'], ascending=[True,True])
        index_keep = []
        min_quality = 0
        for num_pdbs in dfm.num_pdbs.unique():
            dfm_n = dfm[dfm.num_pdbs == num_pdbs]

            best_at_level = dfm_n.last_valid_index()
            if dfm_n.loc[best_at_level, 'match_quality'] >min_quality:
                index_keep +=[best_at_level]
                min_quality = dfm_n.loc[best_at_level, 'match_quality']


        dfm = dfm[dfm.index.isin(index_keep)]
        c = color_palette[int(final_qual)]

        greater = dfm[dfm.match_quality > final_qual]
        less = dfm[dfm.match_quality < final_qual]

        if len(less) > 0:
            ax1.plot( less.match_quality.values.tolist()+[final_qual],less.num_pdbs.values.tolist()+ [final_num_structures], color = c, alpha = 0.7)#, marker ='.', markeredgecolor=c, markerfacecolor ='white')

            final_qualist[1] +=[final_qual]
            final_numlist[1] +=[final_num_structures]
            final_cplxlist[1] +=[cplx]

        if len(greater) > 0:
            ax2.plot([final_qual] +greater.match_quality.values.tolist(),[final_num_structures]+greater.num_pdbs.values.tolist(), color = c, alpha = 0.7)#, marker ='.', markeredgecolor=c, markerfacecolor ='white')
    #         if final_qual < 40:
    #             ax2.plot([final_qual] +greater.match_quality.values.tolist(),[final_num_structures]+greater.num_pdbs.values.tolist(), color = 'k',lw=3, alpha = 1)

            final_qualist[2] +=[final_qual]
            final_numlist[2] +=[final_num_structures]
            final_cplxlist[2] +=[cplx]

        if len(less) == 0 and len(greater) == 0:
            dfm = dfmatch_for_supp[dfmatch_for_supp.cplx == cplx]
            dfm  = dfm.sort_values(by = ['num_pdbs','match_quality'], ascending=[True,True])
            ax0.plot(dfm.match_quality.values.tolist(),dfm.num_pdbs.values.tolist(), color = c, alpha = 0.7)#, marker ='.', markeredgecolor=c, markerfacecolor ='white')

            if len(dfm) == 1:
                final_qualist[0] +=[final_qual]
                final_numlist[0] +=[final_num_structures]
                final_cplxlist[0] +=[cplx]

            else:
                final_qualist[3] +=[final_qual]
                final_numlist[3] +=[final_num_structures]
                final_cplxlist[3] +=[cplx]

    #         print cplx

    k = 0
    for i , q in enumerate(final_qualist[k]):
        n = final_numlist[k][i]
        c = color_palette[int(q)]
        cplx = final_cplxlist[k][i]
        if i ==0:
            label = 'No Choice'
        else:
            label = ''
        ax0.scatter(q, n, color = 'white', edgecolor = c,marker = 'o', s= 150, label =label)

    k = 3
    for i , q in enumerate(final_qualist[k]):
        n = final_numlist[k][i]
        c = color_palette[int(q)]
        cplx = final_cplxlist[k][i]
        if i ==0:
            label = 'No Tradeoff'
        else:
            label = ''
        ax0.scatter(q, n, color = c,marker = 'o', s= 150, label  =label)

    k = 1
    l1 = False
    l2 = False
    for i , q in enumerate(final_qualist[k]):
        n = final_numlist[k][i]
        c = color_palette[int(q)]
        cplx = final_cplxlist[k][i]

        if cplx not in final_cplxlist[2]:
            if not l1:
                ax1.scatter(q, n, color = c,marker = 'o', s= 150, label ='Choose Quality')
                l1 = True
            else:
                ax1.scatter(q, n, color = c,marker = 'o', s= 150, label ='')
        else:
            if not l2:
                ax1.scatter(q, n, color = 'white' , edgecolor=c,marker = 's', s= 150, alpha = 0.5, label ='Mixed Tradeoff')
                l2 = True
            else:
                ax1.scatter(q, n, color = 'white' , edgecolor=c,marker = 's', s= 150, alpha = 0.5, label ='')

    s = 150
    color = 'white'

    k = 2
    l1 = False
    l2= False
    for i , q in enumerate(final_qualist[k]):
        n = final_numlist[k][i]
        c = color_palette[int(q)]
        cplx = final_cplxlist[k][i]

        if cplx not in final_cplxlist[1]:
            if not l1:
                ax2.scatter(q, n, color = c,marker = 'o', s= 150, label ='Fewest Structures')
                l1 = True
            else:
                ax2.scatter(q, n, color = c,marker = 'o', s= 150, label ='')
        else:
            if not l2:
                ax2.scatter(q, n, color = 'white' , edgecolor=c,marker = 's', s= 150, alpha = 0.5, label ='Mixed Tradeoff')
                l2 = True
            else:
                ax2.scatter(q, n, color = 'white' , edgecolor=c,marker = 's', s= 150, alpha = 0.5, label ='')
    #     if i >40:break
    #     print cplx, 
    
    for ax in [ax1,ax2,ax0]:
    #     ax.set_yticks(np.linspace(0,22,23)+0.4)
        ax.set_ylim(0,30)
        ax.set_xlim(0,100)
        ax.legend(loc = 'upper left')
        ax.set_xlabel('Sequence Identity')
        ax.set_ylabel('# of Structures used in protein-complex representation')
    #     ax.set_yticklabels('')
    #     ax.set_xticklabels('')
        ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
        ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    
    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS10CDE-004A-tradeoffs_in_matches.png')
    if save:
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
        
    #########################second figure (panel B)
    
    set0_only_one_match = set(final_cplxlist[0])
    set1_increase_in_quality = set(final_cplxlist[1])
    set2_decrease_in_quality = set(final_cplxlist[2])
    set3_no_valid_tradeoff = set(final_cplxlist[3])
    
    fig2, ax2 = plt.subplots()
    fig2.set_figheight(6)
    labels = venn.get_labels(data = [set0_only_one_match,set3_no_valid_tradeoff,set1_increase_in_quality,set2_decrease_in_quality],
                            )
    for k , v in labels.items():
        if v =='0':
            labels.update({k : ''})

    fig2, ax2 = venn.venn4(labels, names = ["No Choice (Single Choice)",'No Tradeoff (Easy Choice)','Choose Quality',"Choose Quaternary (Fewer) Structures"], fig = fig2,ax = ax2, legend=False)
    
    outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'FigS10B-004A-TradeOffVenn.png')
    if save:
        fig2.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')


    
    # fig.savefig('../figures/raw/figs10-004A-tradeoffs_in_matches.png',dpi = 900,transparent = True, bbox_inches = 'tight')
