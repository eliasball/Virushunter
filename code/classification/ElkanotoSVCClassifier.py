from pulearn import ElkanotoPuClassifier
from sklearn.svm import SVC

class ElkanotoSVCClassifier(PUClassifier):
    """
    Based on the unweighted Elkan & Noto method (https://pulearn.github.io/pulearn/doc/pulearn/elkanoto.html)
    with Support Vector Classification (https://scikit-learn.org/1.5/modules/generated/sklearn.svm.SVC.html)
    as a predictor .
    """

    def __init__(self):
        self.model = None
        self.internalModel = None

    def initializeModel(self, params: list) -> None:
        # use radial basis function and enable estimation of probabilities when classifying
        self.internalModel = SVC(kernel='rbf', probability=True)

        # use basic Elkanoto PU classifier with 0.1 hold-out
        
        # the postive and unlabeled part basically works like this (according to the source code):
        # A hold-out is taken from the positive training datapoints and the model is trained on the rest as positive
        # and all unlabeled datapoints as negatives. The trained model then predicts the probabilites for the hold-out
        # and saves the mean probability.
        # Whenever predicting a new datapoint, the positive probability is predicted using the underlying predictor,
        # and is then divided by the mean probability of the hold-out.
        # The prediction is labeled positive if the ratio is bigger than the given threshold.
        self.model = ElkanotoPuClassifier(self.internalModel, hold_out_ratio=0.1)

    def trainModel(self, trainingDatapoints: list[tuple[list[float], int]]) -> None:
        # get feature vectors and labels from training data
        X = [datapoint[0] for datapoint in trainingDatapoints]
        y = [datapoint[1] for datapoint in trainingDatapoints]
        # make sure to use -1 for negative labels
        y = [1 if label == 1 else -1 for label in y]
        self.model.fit(X, y)

    def classify(self, features: list[float]) -> int:
        # predict single sample and return its classification
        prediction = self.model.predict([features], threshold=0.5)[0]
        # make sure to return 0 for negative labels
        return prediction if prediction == 1 else 0