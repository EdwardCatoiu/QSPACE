import logging
import os
import os.path as op
import urllib3
import sys
import json
import copy
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logging.getLogger("requests").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.ERROR)
log = logging.getLogger(__name__)



def is_ipynb():
    """Return True if the module is running in IPython kernel,
    False if in IPython shell or other Python shell.

    Copied from: http://stackoverflow.com/a/37661854/1592810
    There are other methods there too

    >>> is_ipynb()
    False

    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole?
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal running IPython?
            return False
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter
    
if is_ipynb():
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm    
    
class QSPACE():
    """Generic class containing all directories necessary 
    """

    def __init__(self,
                 ProjectName,
                 create_dirs=True, 
                 root_dir=None,
                 external_hardDrive_dir = None,
                 results_dir = None,
                 input_dir = None,
                demo = False):
        """Initialize the QSPACE project.

        Specify the name of your project, 
        the root directory, 
        the external hardDrive directory where you wish to save all structure-related files, 
        and the path to your ssbio directory.
        
        ssbio can be downloaded at https://github.com/SBRG/ssbio

        """
        self.id = ProjectName


        # Create directories
        self.root_dir = None
        self.externalHardDrive_dir = external_hardDrive_dir
        
        # Project Related directories
        self.Input_dir = None
        self.DataOutput_dir = None
        self.Results_dir = None
#         self.FinalDataOutput_dir = None
        self.FiguresOutput_dir = None
        
        #Global Directories on Local Disk
#         self.ssbioDir = ssbio_dir
        self.FunctionsDir = None
#         self.NotebookDir = None
        self.UniprotSeqsDir = None
        self.AlleleomeSeqsDir = None
        self.SequenceAlignmentDir = None
        self.AlphaMultiSeqDir = None
        self.featuresUniprotDir = None
        self.tmhmmResultsDir = None

        
        #Global Directories on External HardDrive
        # Protein Structures
        self.alphafoldStructuresDir = None
        self.alphafoldMetricsDir = None
        self.alphafoldMultimerDir = None
        self.swissStructuresDir = None
        self.itasserStructuresDir = None
        self.itasserCleanStructuresDir = None
        self.itasserCleanStringRemovedStructuresDir = None
        self.cifStructuresDir = None
        self.pdbStructuresDir = None
        self.bestStructuresDir = None
        self.bioassemblyStructuresDir = None
        self.cifToPdbFolder = None
        self.scannetStructuresDir = None
        #Protein Properties
        self.scratchResultsDir = None
        self.msmsResultsDir = None
        self.disulfideDir = None
        self.proteinInterfaceDir = None
        self.scannetResultsDir = None
        self.dsspResultsDir = None
        self.AlleleomeWTcsvDir = None
        #Membrane Module
        self.opmStructuresToSendDir = None
        self.opmOutputStructuresDir = None
        self.opmOutputDataDir = None
        self.opmCsvDir = None
        self.opmManualOrientDir = None
        self.UniprotCsvDir = None
        self.TMHMMcsvDir = None
#         self.____dir = None
#         self.____dir = None
#         self.____dir = None
#         self.____dir = None
        
        
        
        if create_dirs:
            self._create_dirs(root_dir, input_dir, results_dir, self.id, demo)
        
        directories = copy.deepcopy(self.__dict__)
        for key in ['id']:
            del directories[key]

        with open('qspace_directories.json','w') as f:
            json.dump(directories,f)
#         print (directories)
        
    def _create_dirs(self, root_dir, input_dir, results_dir, ProjectName,demo):
        """Internal method to create QSPACE directories.

        data_dir is created at data/ProjectName
        input_dir is determined by user
        figures_dir is created at figures/ProjectName
        
        all other directories are shared (i.e. you can re-use PDBs you have downloaded in the shared folder for a different project
        (just make sure that you re-align the gene-sequences to the PDB structures in Module #002D).
        
        """
        list_of_dirs = []

        if not root_dir:
            root_dir = os.getcwd()

        if not op.exists(root_dir):
            raise ValueError('{}: folder does not exist'.format(root_dir))

        self.root_dir = root_dir       
        ############################## Project-related on local disk ###############################################
        # Input_dir - directory where all relevant inputs are stored
        self.Input_dir = input_dir
        list_of_dirs.append(self.Input_dir)
        
        self.Results_dir = results_dir
        list_of_dirs.append(self.Results_dir)

        # DataOutput_dir - directory where all data (dataframes and more) will be stored
        self.DataOutput_dir = op.join(self.Results_dir, 'data', ProjectName)
        list_of_dirs.append(self.DataOutput_dir)

        # FinalDataOutput_dir - directory where all final data (dataframes and more) will be stored
#         self.FinalDataOutput_dir = op.join(self.root_dir, 'finalData', ProjectName)
#         list_of_dirs.append(self.FinalDataOutput_dir)
        
        
        # FiguresOutput_dir - directory where all figures will be stored
        self.FiguresOutput_dir = op.join(self.Results_dir, 'figures', ProjectName)
        list_of_dirs.append(self.FiguresOutput_dir)

        ############################### Global Folders on local disk ##############################################
        
        # FunctionsDir - directory where all global functions can be found
        self.FunctionsDir = op.join(self.root_dir, 'functions')
        list_of_dirs.append(self.FunctionsDir)
        
        # tmhmmResultsDir - directory where all membrane-residues from TMHMM will be stored
        self.tmhmmResultsDir = op.join(self.root_dir, '005D3-TMHMM_results')
        list_of_dirs.append(self.tmhmmResultsDir)
       

        ############################### Global Folders on hard drive ##############################################
        
        ### Inputs
        
        # AlleleomeWTcsvDir - directory where all AA-level WT sequence variation .csv files will be stored
        self.AlleleomeWTcsvDir = op.join(self.externalHardDrive_dir,'inputs', 'WTAlleleomeCsv')
        list_of_dirs.append(self.AlleleomeWTcsvDir)
        
        
        
        ### Sequences ##########################################
        # UniprotSeqsDir - directory where all uniprot sequence files will be stored
        self.UniprotSeqsDir = op.join(self.externalHardDrive_dir, 'sequences','uniprot')
        list_of_dirs.append(self.UniprotSeqsDir)
        # AlleleomeSeqsDir - directory where all alleleome sequence files will be stored
        self.AlleleomeSeqsDir = op.join(self.externalHardDrive_dir, 'sequences','alleleome')
        list_of_dirs.append(self.AlleleomeSeqsDir)        
        # AlphaMultiSeqDir - directory where all sequences that need to be sent to Alphafold Multimer  will be stored
        self.AlphaMultiSeqDir = op.join(self.externalHardDrive_dir, 'sequences','to_alphafold_multimer')
        list_of_dirs.append(self.AlphaMultiSeqDir)
        

        ### Needle Alignments ##################################
        # SequenceAlignmentDir - directory where all needle_alignments will be stored
        self.SequenceAlignmentDir = op.join(self.externalHardDrive_dir,'needle_alignments')
        list_of_dirs.append(self.SequenceAlignmentDir)
        
        
        ### Structures ##########################################
        # alphafoldStructuresDir - directory where all alphafold database models will be stored
        self.alphafoldStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_alphafold')
        list_of_dirs.append(self.alphafoldStructuresDir)
        # alphafoldMetricsDir - directory where all alphafold database model metrics will be stored
        self.alphafoldMetricsDir = op.join(self.externalHardDrive_dir, 'structures','all_alphafold-metrics')
        list_of_dirs.append(self.alphafoldMetricsDir)
        # alphafoldMultimerDir - directory where all alphafold multimer models will be stored
        if not demo:
            self.alphafoldMultimerDir = op.join(self.externalHardDrive_dir, 'structures','all_alphafoldMultimer')
        else:
            self.alphafoldMultimerDir = op.join(self.root_dir, "Demo_Results" , 'AFMstructures')
        list_of_dirs.append(self.alphafoldMultimerDir)
        # swissStructuresDir - directory where all SWISS MODELS will be stored
        self.swissStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_swiss')
        list_of_dirs.append(self.swissStructuresDir)
        # itasserStructuresDir - directory where all ITASSER models will be stored
        self.itasserStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_itasser')
        list_of_dirs.append(self.itasserStructuresDir)
        # itasserCleanStructuresDir - directory where all ITASSER cleaned models will be stored
        self.itasserCleanStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_itasser_clean')
        list_of_dirs.append(self.itasserCleanStructuresDir)
        # itasserCleanStringRemovedStructuresDir - directory where all ITASSER cleaned and string-like regions removed will be stored
        self.itasserCleanStringRemovedStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_itasser_clean_string_removed')
        list_of_dirs.append(self.itasserCleanStringRemovedStructuresDir)
        # cifStructuresDir - directory where all .cif files from PDB will be stored
        self.cifStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_cifs')
        list_of_dirs.append(self.cifStructuresDir)
        # pdbStructuresDir - directory where all .pdb files from PDB will be stored
        self.pdbStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_pdbs')
        list_of_dirs.append(self.pdbStructuresDir)
        # bioassemblyStructuresDir - directory where all .cif bioassembly files from PDB will be stored
        self.bioassemblyStructuresDir = op.join(self.externalHardDrive_dir, 'structures','all_bioassemblies')
        list_of_dirs.append(self.bioassemblyStructuresDir)
        # bestStructuresDir - directory where all UniProt-PDB best structures using SSbio/bioservices text API will be stored
        self.bestStructuresDir = op.join(self.externalHardDrive_dir, 'structures','ssbio_best_structures')
        list_of_dirs.append(self.bestStructuresDir)
        
        ##### Modified structures for various property calculations
        # cifToPdbFolder - directory where .cif files are converted to .pdb
        self.cifToPdbFolder = op.join(self.externalHardDrive_dir, 'structures','cif_to_pdb')
        list_of_dirs.append(self.cifToPdbFolder)
        # scannetStructuresDir - directory where all .cif files will be cleaned to send to Scannet
        self.scannetStructuresDir = op.join(self.externalHardDrive_dir, 'structures','Scannet_chain_removed')
        list_of_dirs.append(self.scannetStructuresDir)
        
        
        ### Results   ##########################################
        # featuresUniprotDir - directory where all uniprot feature dataframes will be stored
        self.featuresUniprotDir = op.join(self.externalHardDrive_dir, 'results','results_UniProtFeatures')
        list_of_dirs.append(self.featuresUniprotDir)
        
        # scratchResultsDir - directory where all SCRATCH RESULTS will be stored
        self.scratchResultsDir = op.join(self.externalHardDrive_dir,'results', 'results_SCRATCH')
        list_of_dirs.append(self.scratchResultsDir)
        
        # msmsResultsDir - directory where all MSMS results will be stored
        self.msmsResultsDir = op.join(self.externalHardDrive_dir,'results', 'results_MSMS')
        list_of_dirs.append(self.msmsResultsDir)
        
        # disulfideDir - directory where all Disfulfide-bridge information will be stored
        self.disulfideDir = op.join(self.externalHardDrive_dir,'results','results_DISULFIDE_BRIDGES')
        list_of_dirs.append(self.disulfideDir)
        
        # proteinInterfaceDir - directory where all geometrically calculated protein protein interfaces will be stored
        self.proteinInterfaceDir = op.join(self.externalHardDrive_dir,'results','results_Interface_Regions')
        list_of_dirs.append(self.proteinInterfaceDir)
        
        # scannetResultsDir - directory where all SCANNET protein-protein interfaces will be stored
        self.scannetResultsDir = op.join(self.externalHardDrive_dir,'results','results_SCANNET')
        list_of_dirs.append(self.scannetResultsDir)
        
        # dsspResultsDir - directory where all DSSP results will be stored
        self.dsspResultsDir = op.join(self.externalHardDrive_dir,'results','results_DSSP')
        list_of_dirs.append(self.dsspResultsDir)
        
       
        
        #### Membrane Module
                
        # opmStructuresToSendDir - directory where all membrane structures are converted to PDB and uploaded to OPM will be stored
        self.opmStructuresToSendDir = op.join(self.externalHardDrive_dir, 'results','results_OPM', 'to_opm_pdbs')
        list_of_dirs.append(self.opmStructuresToSendDir)
        
        # opmOutputStructuresDir - directory where all OPM-generated .pdb structures will be stored
        if not demo:
            self.opmOutputStructuresDir = op.join(self.externalHardDrive_dir,'results','results_OPM', 'OPMstructures_new')
        else: #byass the opm server for demo.
            self.opmOutputStructuresDir = op.join(self.root_dir,'Demo_Results','fromOPMstructures')
        list_of_dirs.append(self.opmOutputStructuresDir)
        
        # opmOutputDataDir - directory where all OPM-generated data will be stored
        self.opmOutputDataDir = op.join(self.externalHardDrive_dir,'results','results_OPM', 'OPMdata_new')
        list_of_dirs.append(self.opmOutputDataDir)
        
        # opmCsvDir - directory where all .csv files for OPM-generated .pdb structures will be stored (w/membrane annotation)
        self.opmCsvDir = op.join(self.externalHardDrive_dir,'results','results_OPM', 'OPMcsv_new')
        list_of_dirs.append(self.opmCsvDir)
        
        # opmManualOrientDir - directory where OPM-generated structures (.png) will be stored (these are manually oriented inputs)
        self.opmManualOrientDir = op.join(self.externalHardDrive_dir,'results','results_OPM', 'OPMmanualPNGs')
        list_of_dirs.append(self.opmManualOrientDir)
        
        # UniprotCsvDir - directory where all .csv files for Uniprot-generated membrane calculations will be stored
        self.UniprotCsvDir = op.join(self.externalHardDrive_dir,'results','results_UniProtMEM', 'UniProtCsv')
        list_of_dirs.append(self.UniprotCsvDir)
        
        # TMHMMcsvDir - directory where all .csv files for DeepTMHMM-generated membrane calculations will be stored
        self.TMHMMcsvDir = op.join(self.externalHardDrive_dir,'results','results_tmhmmMEM')
        list_of_dirs.append(self.TMHMMcsvDir)
        
        
        for directory in list_of_dirs:
            if not op.exists(directory):
                os.makedirs(directory)
                log.info('{}: created directory'.format(directory))
            else:
                log.debug('{}: directory already exists'.format(directory))
#         list_of_dirs += [self.ssbioDir]
        
    
    
            
            
if __name__ == '__main__':
    pass