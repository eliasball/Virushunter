# A helper class for representing (or encapsulating) Positive-Unlabeled learning classifiers

class PUClassifier:

    self.host = None

    def __init__(self):
        #Initializes this class
        pass

    def initializeModel(self, params: list) -> None:
        """
        Description: Initializes the PU learning model with the hard parameters
        Args: params (list): List of hard parameters that are necessary to set up before the training process starts (e.g. learning rate)
        """
        pass

    def trainModel(self, trainingDatapoints: list[tuple[list[float], int]]) -> None:
        """
        Description: Trains the PU learning model with the given training data
        Args: trainingDatapoints (list): List of tuples (input features and binary labels (1 for positive, 0 for unlabeled) for data points)
        """
        pass

    def classify(self, features: list[float]) -> int:
        """
        Description: Predicts the class of the given query data point
        Args: queryDatapoint (list): List of input features for the query data point
        Returns: The predicted class of the query data point (1 for positive, 0 for negative)
        """
        pass