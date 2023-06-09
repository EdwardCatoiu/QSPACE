# QSPACE
This package generates the <b>Q</b>uaternary <b>S</b>tructural <b>P</b>roteome <b>A</b>tlas of a <b>Ce</b>ll (QSPACE). 

Direct your questions to ecatoiu@ucsd.edu!

# Contents
- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Run times](#a-note-on-run-times)
- [Graphical Abstract](#graphical-abstract)
- [Description of modules](#description-of-select-modules)
- [License](https://github.com/EdwardCatoiu/QSPACE/blob/main/LICENSE)

# Overview
For any list of genes and protein-complexes with annotated gene-stoichiometries, QPSACE identifies the protein structures that best represent the oligomeric state of the protein complexes that can be re-created from the gene list input. For genes outside the annotated oligomeric protein complexes, QSPACE will use structure-guided analysis to semi-automatically confirm oligomerization of monomer/non-annotated genes in the gene list. When no structures from existing databses can re-create the oligomeric proteins, QSPACE will genereate multi-gene fasta files that can be modelled with ColabFold/Alphafold Multimer. Once the oligomeric structures are determined, protein structural properties, functionally important enzyme domains, wild-type sequence variation, and laboratory mutants are mapped in 3D. QSPACE structures are also oriented across the cell membrane using data from UniProt / OPM / DeepTMHMM. 

This repository contains data generated for the QSPACE of <i>E. coli</i> K-12 MG1655 that can be used to re-create the figures in the main text of our manuscript.  As an example, QSPACE is set up to run for a set of 24 genes (QSPACE-24) that showcase the various modules. Data is also provided for the set of genes in <i>E. coli</i> metabolic model iML1515 (QSPACE-1515) and for <i>E. coli</i> K-12 MG1655 at genome-scale (QSPACE-GS). 

# System Requirements
## Hardware Requirements
The QSPACE package requireds a standard computer with ~8GB RAM to support compiling the final QSPACE dataframe and 300GB (for a query of 4000+ genes) of free disk space. For this reason, we suggest using an external hard drive to save the QSPACE results. 

The runtimes provided below are generated using a computer with the following specs:
- 32 GB RAM
- 4 cores @ 3.60GHz
- 750Mbps internet speed
- Seagate External Hard Drive (STEA4000400), we reccomend a SSD with faster read/write speeds than this one

# Software Requirements
## OS Requirements

The development version of this package has been tested on the following systems:

- Linux: Ubuntu 16.04.3 LTS
- Mac OSX:
- Windows:

# Package Dependencies
QSPACE requires the following packages:
- python v3.7.9
- biopython v1.81 (https://github.com/biopython/biopython)
- nglview v2.7.7 (https://github.com/nglviewer/nglview)
- numpy v1.19.5 (https://github.com/numpy/numpy)
- matplotlib v3.3.3 (https://github.com/matplotlib/matplotlib)
- matplotlib-venn v0.11.9 (https://github.com/konstantint/matplotlib-venn)
- pandas v1.1.5 (https://github.com/pandas-dev/pandas)
- scipy v1.5.4 (https://github.com/scipy/scipy)
- seaborn v0.10.1 (https://github.com/mwaskom/seaborn)
- selenium v3.14.0 (https://github.com/SeleniumHQ/selenium)
- ssbio (https://github.com/SBRG/ssbio)
- scikit-learn v.24.2 (https://scikit-learn.org/stable/install.html)
- pdbecif v1.5 (https://pypi.org/project/PDBeCif/)
- cloudpickle v1.6.0 (https://github.com/cloudpipe/cloudpickle)
- openpyxl v3.0.10 (https://pypi.org/project/openpyxl/)

## Local Software
QSPACE interfaces with following software packages that must be downloaded (details in tutorial):
- ScanNet (https://github.com/jertubiana/ScanNet)
- DisEMBL v1.4 (http://dis.embl.de/)
- SCRATCH-1D v1.1 (https://scratch.proteomics.ics.uci.edu/)
- EMBOSS (https://emboss.sourceforge.net/download/#Stable)

## External Software/Servers
QSPACE interfaces with the following external software/servers (details in tutorial):
- ColabFold w/AFmultimer (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=ADDuaolKmjGW)
- DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)
- Orientations of Proteins in Membranes (OPM) (https://opm.phar.umich.edu/ppm_server2_cgopm)

# Installation Guide
Installation should take < 5 minutes after installing dependencies
1. install package dependecies and local software 
2. git clone https://github.com/EdwardCatoiu/QSPACE.git
3. copy/replace the four files in your ssbio package with their corresponding files in QSPACE/ssbio_files/
4. run the tutorial/demo

# Demo
- Option 1: [follow the tutorial](https://github.com/EdwardCatoiu/QSPACE/blob/main/demo_QSPACE.ipynb) for QSPACE-24 (24 genes that showcase the different modules in this package)
- Option 2: [re-create the figures](https://github.com/EdwardCatoiu/QSPACE/blob/main/Generate_Manuscript_Files.ipynb) in the manuscript from data provided in the supplementary material and/or from data generated by QSPACE-GS

# A note on run-times

We provide <b>all</b> QSPACE-generated results information in the [Demo_Results/data folder](https://github.com/EdwardCatoiu/QSPACE/tree/main/Demo_Results/data) for QSPACE-24_example and only Module 002 results for QSPACE-1515 and QSPACE-GS. These are the PDB/API results so that users can avoid sending identical queries to PDB for <i>E. coli</i> proteins while using the tutorial.      

The QSPACE-generated results for <i>E. coli</i> QSPACE-GS are 8 GB large.  The data downloaded, calculated, and saved (e.g structures from all sources, sequences, membrane calculations, opm structures, structural properties, etc.) when running QSPACE-GS (4000+ genes of <i>E. coli</i>) requires 260 GB of disk space. A link will be added in this repository once this dataset and the results datasets are hosted elsewhere. Using these datasets will avoid additional runtimes associated with generating genome-scale QSPACEs for <i>E. coli</i>.

The following modules require additional run-times
- Module 000B-C (Download gene sequences from UniProt and the wild-type <i>E. coli</i> Alleleome)
- Module 001A-C (Download homology models from Alphafold DB / Swiss-Prot / I-TASSER)
- Module 002A-B (Download structures from PDB using APIs)
- Module 003B/D (Manual curation of oligomeric structures, generate AlphaFold Multimer structures)
- Module 005C (OPM server calculations)
- Module 006B/D (calculations with local software SCRATCH and ScanNet)

The following files are [hosted](https://drive.google.com/drive/folders/1OkXnPK2YP3WAk62Mmu1p00z2dKiS1HQN?usp=share_link) on Google Drive that can be used to reduce run-times for <i>E. coli</i> projects
- All Alphafold Multimer models calculated for <i>E. coli</i> QSPACE-GS (Module 003C/D)
- All OPM server calculations for <i>E. coli</i> QSPACE-GS (Module 005C)
- Folders of images used to manually oriented membrane proteins from Dataset S6-7 (Module 5E)
- [Project-specific results](https://drive.google.com/drive/folders/1_NiNF8I07x632h2vMkk8MO5jyJQDf4b5?usp=share_link) generated in QSPACE-1515 and -GS

<p align="left">
  <img width="400" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Expected_run_times.png">
</p>

# Graphical Abstract
This package generates the <b>Q</b>uaternary <b>S</b>tructural <b>P</b>roteome <b>A</b>tlas of a <b>Ce</b>ll (QSPACE). This oligomeric structural representation of protein complexes includes the calculation of protein structural properties, residue-level membrane integration and subcellular compartmentalization, and mapping of functionally important regions and mutational databases.

<p align="left">
  <img width="500" src="https://github.com/EdwardCatoiu/QSPACE/blob/main/Manuscript/From_Manuscript/fig1_qspace_low_qual.jpg">
</p>

# Description of Select Modules
Additional details are provided in Table S1 of the Supplentary Material of our manuscript.
#### Module 3-4
<i>Module 3A</i> determines the protein structures that can be mapped to the protein complexes (Panel 2a.i).
<i>Module 3B</i> performs a structure-guided re-annotation of falsely annotated protein monomers (Panel 2a.ii).
<i>Module 3C/D</i> generates protein complex sequences for integration with alphafold multimer/colabfold and provides QCQA of multimer predictions (Panel 2a.iii).
<i>Module 4A</i> determines the final structural represenation of all protein complexes (Panel 2b-d).
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

