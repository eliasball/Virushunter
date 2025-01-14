import numpy as np
import pandas as pd
import ast
import csv
import json
import os
from sklearn.model_selection import KFold, cross_validate
from sklearn.metrics import recall_score, make_scorer

from classification.ElkanotoSVCClassifier import ElkanotoSVCClassifier
from dataPreparation import DataPreparation
from dataPreprocessing import DataPreprocessing


featureFilePath = "genomesWithFeatures.csv"

# Training of the PU-Classifier

# Define hosts used to read files and hosts to build a classifier for
hosts = [10090, 9534]
hostsToClassify = [10090]

classifiers = {}
for host in hostsToClassify:
    classifiers[host] = ElkanotoSVCClassifier()

# Create dictionary of HostTaxIDs (as key) and empty lists (of feature lists for each datapoint) for each host (as value)
featuresPerHost = {}
for host in hosts:
    featuresPerHost[host] = []

# Extract features of virus sequences for each host in seperate list and convert them to a np array
print("Reading data")
with open(featureFilePath) as file:
    for line in file:
        values = line.rstrip().split(",")
        host = int(values[0])
        if host in hosts:
            featuresPerHost[host].append(list(map(float, values[2:])))

for host in hosts:
    featuresPerHost[host] = np.array(featuresPerHost[host])

# Transform features (X) and host classes (y) into expected format of PU Classifier and run k-fold cross validation for each host
X = np.concatenate([featuresPerHost[host] for host in hosts])

def getYForHost(h):
    y = []
    for host in hosts:
        y += [1 if host == h else 0] * len(featuresPerHost[host])
    return np.array(y)

scores = {}
specificity = make_scorer(recall_score, pos_label=0)
for host in hostsToClassify:
    print("Running cross validation for host", host)
    # load the training data
    X = X
    y = getYForHost(host)

    # run k-fold cross validation
    cv = KFold(n_splits=4, shuffle=True)
    #cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=42)
    scores[host] = cross_validate(classifiers[host], X, y, cv=cv, scoring={'balanced_accuracy': 'balanced_accuracy', 'recall': 'recall', 'specificity': specificity, 'precision': 'precision', 'f1': 'f1'})

print(scores)