from classification import PUClassifier

#  Performs the training of the individual PU binary classifiers

class ModelTraining:

    def __init__(self):
        #Initializes this class
        pass

    def trainPUClassifier(self, trainingDatapoints: list[tuple[list[float], int]]) -> PUClassifier:
        """
        Description: Trains a positive and unlabeled learning classifier (with bootstrapping and confidence scoring)
        Args: trainingDatapoints (list[tuple[list[float], int]]): List of tuples (input features and binary labels (1 for positive, 0 for unlabeled) for data points)
        Returns: Trained PU classifier (maybe just weights or the whole trained clasifier if we use library functions and classes)
        """
        pass
    

    def extractRelevantTrainingData(self, dataMatrix: list[list], hostID: str) -> list[tuple[list[float], int]]:
        """
        Description: Extracts for a specific classifier relevant datapoints out of all training data and give them labels (positive, if hostID matches or unlabeled, esle)
        Args: dataMatrix (list[list]): Matrix of all datapoints (each row in style: virusname, input features, hostname(s)) and hostID (str): ID of the host for which the data is preperated
        Returns: List of tuples (x,y) where x represent the input features and y the label (1 = positive, 0 = unlabeled)
        """
        pass