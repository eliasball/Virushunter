from classification import PUClassifier

#  Performs the training and application of PU binary classifiers and composes then to one model

class Classification:

    def __init__(self):
        #Initializes this class
        pass

    def trainPUClassifier(self, trainingDatapoints: list):
        """
        Description: Trains a positive and unlabeled learning classifier (with bootstrapping and confidence scoring)
        Args: trainingDatapoints (list): List of tuples (input features and binary labels (1 for positive, 0 for unlabeled) for data points)
        Returns: Trained PU classifier (maybe just weights or the whole trained clasifier if we use library functions and classes)
        """
        pass

    def classify(self, genomeSequence: str, classifier: PUClassifier):
        """
        Description: Classifies given query genome sequence based on the given classifier
        Args: genomeSequence (str): Genome sequence to classify and classifier and (PUClassifier): THe model that is used for classification
        Returns: Confidence score for the host type that corresponds to the given classifier
        """
        pass
    

    def classifyWithMultipleModels(self, genomeSequence: str, classifiers: dict):
        """
        Description: Composes a set of PU classifiers to predict the host with the highest confidence score
        Args: genomeSequence (str): The genome sequence to classify and classifiers (dict): A dictionary (keys = host names and values = trained PUClassifier instances)
        Returns: A tupel of the host with the highest confidence score and the score itself (e.g. (hostName, confidenceScore))
        """
    

