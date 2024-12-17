
#  Preprocessing of dataobjects (extraction of valueable features out of given genome sequence)

class DataPreprocessing:

    def __init__(self):
        #Initializes this class
        pass

    def extractKmerFrequencies(self, genomeSequence: str, k: int):
        """
        Description: Extracts k-mer frequencies from a genome sequence for a given k
        Args: genomeSequence (str): DNA sequence and k (int): Size of the k-mers (hyperparameter -> should be individually optimized)
        Returns: K-mer frequencies as a (indexed) list/dictonary
        """
        pass

    def calculateCodonUsage(self, genomeSequence: str):
        """
        Description: Calculates codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Codon usage  usage frequencies as a (indexed) list/dictonary
        """
        pass
    

    def calculateDicodonUsage(self, genomeSequence: str):
        """
        Description: Calculates di-codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Di-Codon usage  usage frequencies as a (indexed) list/dictonary
        """
        pass

    def computeDinucleotideFrequencies(self, genomeSequence: str):
        """
        Description: Computes dinucleotide frequencies from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Dinucleotide frequencies as a (indexed) list/dictonary
        """
        pass
    

    def bindFeaturesToMatrix(self, chosenFeatures: list):
        """
        Description: Binds all chosen features (list of lists) to a feature matrix that will be the input of the classifiers
        Args: chosenFeatures (list): List of Feature Lists for all chosen features (or names of chosen features -> depends on actual implementation)
        Returns: Matrix with datapoints (genomes) as rows and the relevant features as columns
        """
        pass
