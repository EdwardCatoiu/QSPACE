# QSPACE
This package generates the <b>Q</b>uaternary <b>S</b>tructural <b>P</b>roteome <b>A</b>tlas of a <b>Ce</b>ll (QSPACE) as described in ___link___. 

Direct your questions to ecatoiu@ucsd.edu!

# Contents
- [Overview](#overview)
- [Package Dependencies](#package-dependencies)
- [External Servers / Software](#external-servers--software)
- [Install Instructions](#install-instructions)
- [Run times](#a-note-on-run-times)
- [Graphical Abstract](#graphical-abstract)
- [Description of modules](#description-of-select-modules)
- [License](https://github.com/EdwardCatoiu/QSPACE/blob/main/LICENSE)

# Overview
For any list of genes and protein-complexes with annotated gene-stoichiometries, QPSACE identifies the protein structures that best represent the oligomeric state of the protein complexes that can be re-created from the gene list input. For genes outside the annotated oligomeric protein complexes, QSPACE will use structure-guided analysis to semi-automatically confirm oligomerization of monomer/non-annotated genes in the gene list. When no structures from existing databses can re-create the oligomeric proteins, QSPACE will genereate multi-gene fasta files that can be modelled with ColabFold/Alphafold Multimer. Once the oligomeric structures are determined, protein structural properties, functionally important enzyme domains, wild-type sequence variation, and laboratory mutants are mapped in 3D. QSPACE structures are also oriented across the cell membrane using data from UniProt / OPM / DeepTMHMM. 

This repository contains data generated for the QSPACE of <i>E. coli</i> K-12 MG1655 that can be used to re-create the figures in the main text of our manuscript.  As an example, QSPACE is set up to run for a set of 24 genes (QSPACE-24) that showcase the various modules. Data is also provided for the set of genes in <i>E. coli</i> metabolic model iML1515 (QSPACE-1515) and for <i>E. coli</i> K-12 MG1655 at genome-scale (QSPACE-GS). 

# Package Dependencies
QSPACE was successfully tested using:
- python v3.7.9
- biopython v1.81 (https://github.com/biopython/biopython)
- DisEMBL v1.4 (http://dis.embl.de/)
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

# External Servers / Software
- ColabFold w/AFmultimer (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=ADDuaolKmjGW)
- DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)
- Orientations of Proteins in Membranes (OPM) (https://opm.phar.umich.edu/ppm_server2_cgopm)
- ScanNet (https://github.com/jertubiana/ScanNet)

# Install Instructions
- install package dependecies 
- git clone https://github.com/EdwardCatoiu/QSPACE.git
- follow the tutorial in [demo_QSPACE.ipynb](https://github.com/EdwardCatoiu/QSPACE/blob/main/demo_QSPACE.ipynb) 
- ... or generate the figures in the [manuscript](https://github.com/EdwardCatoiu/QSPACE/blob/main/Generate_Manuscript_Files.ipynb)

# A note on run-times

We provide all QSPACE-generated results information in the Demo_Results/data folder for QSPACE-24_example. We provided QSPACE-generated results for Module 002 for QSPACE-1515 and QSPACE-GS. These are the PDB/API results so that users can avoid sending identical queries to PDB for E. coli proteins while testing this repository.      

The QSPACE-generated results for E. coli QSPACE-GS are 8 GB large.  The data downloaded (structures from all sources, sequences, membrane calculations, opm structures, structural properties, etc.) to return the QSPACE-GS results dataframe is 260 GB large. A link will be added in this repository once this dataset and the results datasets are hosted elsewhere. Using these datasets will avoid additional runtimes associated with generating the QSPACEs.

The following modules require additional run-times
- Module 000B-C (UniProt download sequences / Alleleome download sequences)
- Module 001A-C (Download Alphafold structures / Swiss-Prot structures / I-tasser structures)
- Module 002A-B (PDB API/download structures)
- Module 003B/D (Manual curation of oligomeric structures, calculate AlphaFold Multimer structures)
- Module 005C (OPM server calculations)
- Module 006B/D (SCRATCH, ScanNet)

<p align="left">
  <img width="400" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Expected_run_times.png">
</p>

# Graphical Abstract
This package generates the <b>Q</b>uaternary <b>S</b>tructural <b>P</b>roteome <b>A</b>tlas of a <b>Ce</b>ll (QSPACE). This oligomeric structural representation of protein complexes includes the calculation of protein structural properties, residue-level membrane integration and subcellular compartmentalization, and mapping of functionally important regions and mutational databases.

<p align="left">
  <img width="500" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/From_Manuscript/fig1_qspace_low_qual.jpg">
</p>

# Description of Select Modules

#### Module 3-4
<i>Module 3A</i> determines the protein structures that can be mapped to the protein complexes (Panel 2a.i).
<i>Module 3B</i> performs a structure-guided re-annotation of falsely annotated protein monomers (Panel 2a.ii).
<i>Module 3C/D</i> generates protein complex sequences and for integration with alphafold multimer/colabfold and QCQA of multimer predictions (Panel 2a.iii).
<i>Module 4A</i> determines the final protein represenation of all protein complexes (Panel 2b-d).
<p align="left">
  <img width="500" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/From_Manuscript/fig2_qspace_low_qual.jpg">
</p>

#### Module 5
<i>Module 5A</i> finds all potential membrane structures (Panel 4a).
<i>Module 5B-D</i> calculates membrane embeddededness of proteins (Panel 4b).
<i>Module 5E</i> performs QCQA analysis of the membrane calculations and identifies viable membrane planes (Panel 4c).
<i>Module 5F</i> orients each protein across the membrane with angstrom-level precision and maps the residue-level information back to the QSPACE (Panel 4d-e).
<p align="left">
  <img width="500" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/From_Manuscript/fig4_qspace_low_qual.jpg">
</p>

### Module 7
<i>Module 7A</i> maps all enzymatically important regions to the QSPACE (Panel 3a, left).
<i>Module 7B/C</i> maps laboratory-acquired mutations from [ALE](https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/Inputs/007B-ALE_mutations_input.csv) and [LTEE](https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/Inputs/007C-LTEE_mutations_input.csv) datasets to the QSPACE (Panel 3a, center).
<i>Module 7D</i> maps [the alleleome](https://www.pnas.org/doi/10.1073/pnas.2218835120) wild-type sequence variation [data](https://github.com/EdwardCatoiu/Alleleome) to the QSPACE (Panel 3a,right).

<p align="left">
  <img width="500" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/From_Manuscript/fig3_module7A-D_low_qual.jpg">
</p>

