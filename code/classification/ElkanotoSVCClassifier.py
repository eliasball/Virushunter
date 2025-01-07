from pulearn import ElkanotoPuClassifier
from sklearn.svm import SVC
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.multiclass import unique_labels

class ElkanotoSVCClassifier(ClassifierMixin, BaseEstimator):
    """
    Based on the unweighted Elkan & Noto method (https://pulearn.github.io/pulearn/doc/pulearn/elkanoto.html)
    with Support Vector Classification (https://scikit-learn.org/1.5/modules/generated/sklearn.svm.SVC.html)
    as a predictor.
    """

    def __init__(self, threshold=0.5, hold_out_ratio=0.1):
        self.threshold = threshold
        self.hold_out_ratio = hold_out_ratio

    def fit(self, X, y):
        ## init models
        # use radial basis function and enable estimation of probabilities when classifying
        self.internalModel = SVC(kernel='rbf', probability=True)

        # use basic Elkanoto PU classifier with hold-out
        
        # the postive and unlabeled part basically works like this (according to the source code):
        # A hold-out is taken from the training datapoints and the model is trained on the rest as positive
        # and all unlabeled datapoints as negatives. The trained model then predicts the probabilites for the
        # positives in the hold-out and saves the mean probability.
        # Whenever predicting a new datapoint, the positive probability is predicted using the underlying predictor,
        # and is then divided by the mean probability of the hold-out.
        # The prediction is labeled positive if the ratio is bigger than the given threshold.
        self.model = ElkanotoPuClassifier(self.internalModel, hold_out_ratio=self.hold_out_ratio)

        ## fit
        self.classes_ = unique_labels(y)
        self.model.fit(X, y)

        return self

    def predict(self, X):
        return self.model.predict(X, threshold=self.threshold)