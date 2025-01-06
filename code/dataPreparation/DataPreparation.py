

# Handles the loading and preparation of genome data (for training purposes and for solving classification for queries)
# Responsible for reading input genome sequences and organizing them into a usable format for feature extraction


import os


class DataPreparation:
    
    def __init__(self):
        #Initializes this class
        pass
  
    def loadTrainingGenomeSequences(self, filepath: str) -> list[tuple[str, str, str]]:
        """
        Description: Loads annotated genome sequences for training throught executing a given script
        Args: filepath (str): Path to the script that loads the annotated genome sequence files
        Returns: List/Dictonary of genome sequences with their name (taxonomic ID) and corresponding host name(s)
        """
        # Not too sure how to do: needs to run script "virusHostQuery.py" no args


        # Following section has not been tested, just a skeleton idea

        dataPath = "../../viral_genomes/"

        genome_data = []
    
        # Traverse the host folders
        for host_folder in os.listdir(dataPath):
            host_path = os.path.join(dataPath, host_folder)
            
            # Ensure it's a folder
            if not os.path.isdir(host_path):
                continue
            
            # Process each file in the folder
            for virus_file in os.listdir(host_path):
                virus_path = os.path.join(host_path, virus_file)
                
                # Ensure the file is a fasta file
                if not (virus_file.endswith(".fasta") or virus_file.endswith(".fna")):
                    continue
                
                # Extract host taxid (folder name) and virus taxid (file name without extension)
                host_taxid = host_folder
                virus_taxid = os.path.splitext(virus_file)[0]
                
                # Read the file and extract the sequence
                try:
                    with open(virus_path, "r") as file:
                        lines = file.readlines()
                        sequence = "".join(lines[1:]).replace("\n", "")  # Skip header and join the rest
                        
                    # Create a tuple and add to the list
                    genome_data.append((host_taxid, virus_taxid, sequence))
                
                except Exception as e:
                    print(f"Error reading file {virus_path}: {e}")
        
            return genome_data
        
    def loadQueryGenomeSequence(self, filepath: str) -> list[tuple[str, str]]:
        """
        Description: Loads genome sequence from a file (in FASTA format)
        Args: filepath (str): Path to the FASTA file of the query
        Returns: Tupel of genome sequence of query (and the name (taxonomic ID) of the queried virus)
        """
        pass

    def organizeTrainingDataByHost(self, genomeData: list[tuple[str, str, str]]) -> list[tuple[str, str, str]]:
        """
        Description: Organizes the training sequences for the genome data by the  host type
        Args: genomeData (list[tuple[str, str, str]]): List of genome sequences (as FASTA files annotated by the virusname and the corresponding host name)
        Returns: Mapping of host types to genome sequences for preprocessing of the training set (individual training of all PU classifier)
        """
        pass
