from classification import PUClassifier

#  Performs the application of the PU binary classifiers and composes then to one model

class HostClassification:

    def __init__(self):
        #Initializes this class
        pass


    def classify(self, genomeSequence: str, classifier: PUClassifier) -> float:
        """
        Description: Classifies given query genome sequence based on the given classifier
        Args: genomeSequence (str): Genome sequence to classify and classifier (PUClassifier): The model that is used for classification
        Returns: Confidence score for the host type that corresponds to the given classifier
        """
        pass
    

    def classifyWithMultipleModels(self, genomeSequence: str, classifiers: dict[str, PUClassifier]) -> tuple[str, float]:
        """
        Description: Composes a set of PU classifiers to predict the host with the highest confidence score
        Args: genomeSequence (str): The genome sequence to classify and classifiers (dict[str, PUClassifier]): A dictionary (keys = host names and values = trained PUClassifier instances)
        Returns: A tupel of the host with the highest confidence score and the score itself (e.g. (hostName, confidenceScore))
        """
    

