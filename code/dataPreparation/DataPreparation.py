

# Handles the loading and preparation of genome data (for training purposes and for solving classification for queries)
# Responsible for reading input genome sequences and organizing them into a usable format for feature extraction


class DataPreparation:
    
    def __init__(self):
        #Initializes this class
        pass
  
    def loadTrainingGenomeSequences(self, filepath: str):
        """
        Description: Loads annotated genome sequences for training throught executing a given script
        Args: filepath (str): Path to the script that loads the annotated genome sequence files
        Returns: List/Dictonary of genome sequences with their name, (taxonomic ID) and corresponding host name(s)
        """
        pass
    
    def loadQueryGenomeSequence(self, filepath: str):
        """
        Description: Loads genome sequence from a file (in FASTA format)
        Args: filepath (str): Path to the FASTA file of the query
        Returns: Tupel of genome sequence of query and the name of the queried virus
        """
        pass

    def organizeTrainingDataByHost(self, genomeData: list):
        """
        Description: Organizes the training sequences for the genome data by the  host type
        Args: genomeData (list): List of genome sequences (as FASTA files annotated by the virusname and the corresponding host name)
        Returns: Mapping of host types to genome sequences for preprocessing of the training set (individual training of all PU classifier)
        """
        pass
