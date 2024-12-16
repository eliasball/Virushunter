
#  Preprocessing of dataobjects (extraction of valueable features out of given genome sequence)

class DataPreprocessing:

    def __init__(self):
        #Initializes this class
        pass

    def extractKmerFrequencies(self, genomeSequence: str, k: int):
        """
        Description: Extracts k-mer frequencies from a genome sequence.
        Args: genomeSequence (str): DNA sequence and k (int): Size of the k-mers (hyperparameter -> should be individually optimized)
        Returns: K-mer frequencies as a (indexed) list 
        """
        pass

    def calculateCodonUsage(self, genomeSequence: str):
        """
        Description: Calculates codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Codon usage  usage frequencies as a (indexed) list
        """
        pass
    

    def calculateDicodonUsage(self, genomeSequence: str):
        """
        Description: Calculates di-codon usage from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Di-Codon usage  usage frequencies as a (indexed) list
        """
        pass

    def computeDinucleotideFrequencies(self, genomeSequence: str):
        """
        Description: Computes dinucleotide frequencies from the genome sequence.
        Args: genomeSequence (str): DNA sequence
        Returns: Dinucleotide frequencies as a (indexed) list
        """
        pass
