
# -*- coding: utf-8 -*-
from .utils import *



def fig4_A(potentialMembraneData,fig = False, ax = False, save = False,legend = False, outfile = False):
    
    if not fig:
        fig,ax =plt.subplots()
        fig.set_figheight(6)
        fig.set_figwidth(6)

    uni = potentialMembraneData['UNI_structures']
    iml = potentialMembraneData['IML_structures']
    go  = potentialMembraneData['GO_structures']
    eco = potentialMembraneData['ECO_structures']
    
    labels = venn.get_labels([uni,iml,go,eco] )
    for k ,v in labels.items():
        if v =='0':
            labels.update({k : ''})
    
    fig, ax = venn.venn4(labels  , names=['UNIPROT','iML1515','Gene Ontology (GO)', "EcoCyc"], fig = fig, ax = ax, legend = legend,)
    ax.set_title('Identify membrane structures')

   

    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig4A-005A-PotentialMembraneStructures.png')

        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')
    
    return fig, ax


def fig4E_aa(dfQuatProteome,fig = False, ax = False, save = False, outfile = False):
    if not fig:
        fig,ax =plt.subplots()
        fig.set_figheight(6)
        fig.set_figwidth(6)
    

    compartments_AA = {18:['Extracellular (secreted)'],
                       17:['OM-Extracellular-Bulb'],
                       16:['OM-Extracellular-Embedded','OM-Periplasmic_E-Embedded'],
                       15:['OM-Periplasmic_E-Bulb'],
                       14:['Periplasm (OM-Associated)'],
                       13:['Periplasm'],
                       12:['Periplasm (IM-Associated)'],
                       11:['IM-Periplasmic_C-Bulb'],
                       10:['IM-Periplasmic_C-Embedded','IM-Cytoplasmic-Embedded'],
                       9:['IM-Cytoplasmic-Bulb'],
                       8:['Cytosol (IM-Associated)'],
                       7:['Cytosol'],
                       6:['OM-OM_1-Embedded', 'OM-OM_2-Embedded',],
                       5:['OM-OM_1-Bulb', 'OM-OM_2-Bulb',],
                       4:['Periplasm (MEM-Associated)'],
                       3:[ 'IM-IM_1-Embedded', 'IM-IM_2-Embedded',],
                       2:['IM-IM_1-Bulb', 'IM-IM_2-Bulb',],
                       1:['Inner Membrane (No Area)','Membrane Unknown (No Area)', 'Outer Membrane (No Area)'],
                       0:['Unknown'],
                      }
    yticklabel = []
    for index, aaComps in compartments_AA.items():
        ax.barh(y= index, width = len(dfQuatProteome[dfQuatProteome.get('005F_AA_Compartment').isin(aaComps)]) )
        yticklabel += [aaComps[0]]
    yticklabel.reverse()
    ax.set_yticks(range(0,19))
    ax.set_yticklabels(yticklabel)
    ax.set_xscale('log')
    for ax in [ax]:
        ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
        ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    ax.set_xlim(10,10000000)
    ax.set_xlabel('Amino Acids', size = 15)
    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig4E_aa-005A-aa_level_compartments.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')

    return fig, ax

def fig4E_prot(dfQuatProteome,fig = False, ax = False, save = False,legend = True, outfile = False):
    compartments_AA = {18:['Extracellular (secreted)'],
                       17:['OM-Extracellular-Bulb'],
                       16:['OM-Extracellular-Embedded','OM-Periplasmic_E-Embedded'],
                       15:['OM-Periplasmic_E-Bulb'],
                       14:['Periplasm (OM-Associated)'],
                       13:['Periplasm'],
                       12:['Periplasm (IM-Associated)'],
                       11:['IM-Periplasmic_C-Bulb'],
                       10:['IM-Periplasmic_C-Embedded','IM-Cytoplasmic-Embedded'],
                       9:['IM-Cytoplasmic-Bulb'],
                       8:['Cytosol (IM-Associated)'],
                       7:['Cytosol'],
                       6:['OM-OM_1-Embedded', 'OM-OM_2-Embedded',],
                       5:['OM-OM_1-Bulb', 'OM-OM_2-Bulb',],
                       4:['Periplasm (MEM-Associated)'],
                       3:[ 'IM-IM_1-Embedded', 'IM-IM_2-Embedded',],
                       2:['IM-IM_1-Bulb', 'IM-IM_2-Bulb',],
                       1:['Inner Membrane (No Area)','Membrane Unknown (No Area)', 'Outer Membrane (No Area)'],
                       0:['Unknown'],
                      }
    
    yticklabel = []
    for index, aaComps in compartments_AA.items():
        yticklabel += [aaComps[0]]
    yticklabel.reverse()

    protein_locations = {}

    for s in tqdm(dfQuatProteome.structureId.unique()):
        if s in protein_locations:
            continue
        dfs = dfQuatProteome[dfQuatProteome.structureId == s]
        protein_compartments = dfs.get('005F_Protein_Compartment').dropna().unique()

        if len(protein_compartments) > 1:
            print (s)
            break
        if len(protein_compartments) == 1:
            protein_locations.update({s : protein_compartments[0]})
        else:
            protein_locations.update({s : 'Unknown'})

    
    if not fig:
        fig,ax =plt.subplots()
        fig.set_figheight(6)
        fig.set_figwidth(6)
        
        
        
    
    dfproteinLocations = pd.DataFrame()
    dfproteinLocations['Compartment'] = pd.Series(protein_locations)
    for compartment in dfproteinLocations.Compartment.unique():
        dfc = dfproteinLocations[dfproteinLocations.Compartment == compartment]
        found_yval = False
        for y_val, aaComps in compartments_AA.items():
            if compartment in aaComps:
                found_yval = True
                break
        if found_yval:
            ax.scatter(len(dfc), y_val, marker = 's',s=200, label = compartment)
        else:
            if compartment == 'TransMembrane':
                ax.scatter(len(dfc), 5, marker = 's',s=200, label = compartment)
            elif compartment == 'Inner Membrane':
                ax.scatter(len(dfc), 10, marker = 's',s=200, label = compartment)
            elif compartment == 'Outer Membrane':
                ax.scatter(len(dfc), 16, marker = 's',s=200, label = compartment)
            else:
                print (compartment)
    ax.set_xlim(0.5,10000)
    ax.set_yticks(range(0,19))
    
    ax.set_yticklabels(yticklabel) #from prev fig
    ax.set_xscale('log')
    if legend:
        ax.legend(loc=(1,0))
    ax.set_xlabel('Protein Structures',size = 15)
    ax.set_ylabel('Cellular Compartment',size = 15)
    for ax in [ax]:
        ax.tick_params(axis='both', which='major', labelsize=14,length = 8, width = 1.5,)
        ax.tick_params(axis='both', which='minor', labelsize=14,length = 4, width = 1.5)
    
    if save:
        if not outfile:
            outfile = op.join(qspaceDirs['FiguresOutput_dir'], 'Fig4E_prot-005A-prot_level_compartments.png')
        fig.savefig(outfile,dpi = 900,transparent = True, bbox_inches = 'tight')

    return fig, ax