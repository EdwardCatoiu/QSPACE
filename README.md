# QSPACE

# Description:
This package generates the <b>Q</b>uaternary <b>S</b>tructural <b>P</b>roteome <b>A</b>tlas of a <b>Ce</b>ll (QSPACE) as described in ___link___.  

For any list of genes and protein-complexes with annotated gene-stoichiometries, QPSACE identifies the protein structures that best represent the oligomeric state of the protein complexes that can be re-created from the gene list input. For genes outside the annotated oligomeric protein complexes, QSPACE will use structure-guided analysis to semi-automatically confirm oligomerization of monomer/non-annotated genes in the gene list. When no structures from existing databses can re-create the oligomeric proteins, QSPACE will genereate multi-gene fasta files that can be modelled with ColabFold/Alphafold Multimer. Once the oligomeric structures are determined, protein structural properties, functionally important enzyme domains, wild-type sequence variation, and laboratory mutants are mapped in 3D. QSPACE structures are also oriented across the cell membrane using data from UniProt/OPM/deepTMHMM. 

This repository contains data generated for the QSPACE of <i>E. coli</i> K-12 MG1655 that can be used to re-create the figures in the main text of our manuscript.  As an example, QSPACE is set up to run for a set of 24 genes (QSPACE-24) that showcase the various modules. Data is also provided for the set of genes in <i>E. coli</i> metabolic model iML1515 (QSPACE-1515) and for <i>E. coli</i> K-12 MG1655 at genome-scale (QSPACE-GS). 

# Package Dependencies:
QSPACE was successfully tested using:
- python v3.7.9
- biopython v1.81 (https://github.com/biopython/biopython)
- nglview v2.7.7 (https://github.com/nglviewer/nglview)
- numpy v1.19.5 (https://github.com/numpy/numpy)
- matplotlib v3.3.3 (https://github.com/matplotlib/matplotlib)
- matplotlib-venn v0.11.9 (https://github.com/konstantint/matplotlib-venn)
- pandas v1.1.5 (https://github.com/pandas-dev/pandas)
- SCRATCH-1D v1.1 (https://scratch.proteomics.ics.uci.edu/)
- scipy v1.5.4 (https://github.com/scipy/scipy)
- seaborn v0.10.1 (https://github.com/mwaskom/seaborn)
- selenium v3.14.0 (https://github.com/SeleniumHQ/selenium)
- ssbio (https://github.com/SBRG/ssbio)

# External Servers / Software:
- ColabFold w/AFmultimer (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=ADDuaolKmjGW)
- DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)
- Orientations of Proteins in Membranes (OPM) (https://opm.phar.umich.edu/ppm_server2_cgopm)
- ScanNet (https://github.com/jertubiana/ScanNet)

# Install Instructions

