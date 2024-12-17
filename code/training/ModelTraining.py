from classification import PUClassifier

#  Performs the training of the individual PU binary classifiers

class ModelTraining:

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