
# -*- coding: utf-8 -*-

from .utils import *


def fig1_E(dfskeleton):
    cplxList = dfskeleton.cplx.unique()
    genesList = dfskeleton.gene.unique()
    structuresList = dfskeleton.structureId.unique()
    num_AAs = len(dfskeleton)
    unique_genome_positions =  dfskeleton.GenomeLocation.unique()

    print ("Number of Genes       : {}".format(len(genesList)))
    print ("Number of Enzymes     : {}".format(len(cplxList)))

    print ("Number of Structures  : {}".format(len(structuresList)))
    print ("Number of Amino Acids : {}".format(num_AAs))
    print ("Number of Genome Positions : {}".format(len(unique_genome_positions)))
    
    
def fig1_F(dfskeleton,fig = False, ax = False, save = False,save_outfile = False,figureOutfile =False):
    if not fig:
        fig,ax =plt.subplots()
    
    def getStructuralCoverage(dfskeleton):
        dfstructure_coverage = pd.DataFrame()
        for cplx in tqdm(dfskeleton.cplx.unique()):
        #     if cplx != 'ABC-34-CPLX':
        #         continue

            dfc = dfskeleton[dfskeleton.cplx == cplx]
            total_aas = 0 
            total_aas_covered = 0

            for structureId in dfc.structureId.unique():
                dfs = dfc[dfc.structureId == structureId]
                dfs = dfs[dfs.geneSeqId.isna() == False]
                structureStoich = dfs.structureStoich.values[0]

                total_aas += structureStoich*len(dfs)
                total_aas_covered += structureStoich*len(dfs[dfs.structNum.isna() == False])

        #     print total_aas, total_aas_covered
            dfstructure_coverage.loc[cplx, 'total_aas'] = total_aas
            dfstructure_coverage.loc[cplx, 'total_aas_covered'] = total_aas_covered

        dfstructure_coverage['coverage'] = dfstructure_coverage.total_aas_covered / dfstructure_coverage.total_aas
        dfstructure_coverage = dfstructure_coverage.sort_values(by = ['total_aas','coverage'], ascending= [ False, False])
        dfstructure_coverage = dfstructure_coverage.sort_values(by = 'coverage')

        if save_outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'DataS2-structural_coverage_of_enzymes.csv')
#         print 'Saving....\n\t> {}'.format(outfile)
            dfstructure_coverage.to_csv(outfile)
        return dfstructure_coverage
    
    
    
    ###get data
    dfstructure_coverage = getStructuralCoverage(dfskeleton)
    global_average = float(len(dfskeleton.structNum.dropna())) / float(len(dfskeleton))
    ### begin figure
    
    ax.hist(dfstructure_coverage.coverage, bins = np.linspace(0,1,51), color = u"grey", alpha = 1, orientation = 'horizontal' )
    ax.set_xscale('log')
    ax.set_xlabel('Enzyme Count',size = 14)
    ax.set_ylim(0., 1.02)
    ax.set_xlim(0.5, 10000)
    ax.set_ylabel('Structural Coverage',size = 14)


    for ax in [ax]:
        ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
        ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    #     ax.set_xticklabels('')
    #     ax.set_yticklabels('')
    ax.plot((0,100000), (global_average,global_average), '--k')
    ax.annotate('Global Avg. = {:.2f}'.format(global_average), xy = (120,0.88),)
    

    if save:
        if not figureOutfile:
            figureOutfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig1F-004C-StructuralCoverageOfProteome.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax
    
    