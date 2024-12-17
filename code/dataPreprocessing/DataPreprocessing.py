
#  Preprocessing of dataobjects (extraction of valueable features out of given genome sequence)

class DataPreprocessing:

    def __init__(self):
        #Initializes this class
        pass

    def extractKmerFrequencies(self, genomeSequence: str, k: int) -> dict[str, int]:
        """
        Description: Extracts k-mer frequencies from a genome sequence for a given k
        Args: genomeSequence (str): DNA sequence and k (int): Size of the k-mers (hyperparameter -> should be individually optimized)
        Returns: K-mer frequencies as a (indexed) list/dictonary
        """
        pass

    def calculateCodonUsage(self, genomeSequence: str) -> dict[str, int]:
        """
        Description: Calculates codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Codon usage  usage frequencies as a (indexed) list/dictonary
        """
        pass
    

    def calculateDicodonUsage(self, genomeSequence: str) -> dict[str, int]:
        """
        Description: Calculates di-codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Di-Codon usage  usage frequencies as a (indexed) list/dictonary
        """
        pass

    def computeDinucleotideFrequencies(self, genomeSequence: str) -> dict[str, int]:
        """
        Description: Computes dinucleotide frequencies from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Dinucleotide frequencies as a (indexed) list/dictonary
        """
        pass
    

    def bindFeaturesToMatrix(self, chosenFeatures: dict[list]) -> list[list]:
        """
        Description: Binds all chosen features of all training genomes (dict of lists) to a feature matrix that will be the input of the classifiers
        Args: chosenFeatures (dict[list]): Dictonary of Feature Lists for all chosen features (or names of chosen features -> depends on actual implementation) for each genome
        Returns: Matrix with datapoints (genomes) as rows and the relevant features as columns (-> Structure: (virusnames, input features, hostnames))
        """
        pass
