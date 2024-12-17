from dataPreparation import DataPreparation
from dataPreprocessing import DataPreprocessing
from classification import HostClassification, PUClassifier
from training import ModelTraining

# Main method to control the flow of the whole pipeline 
# TODO: We should outsource the training process into a extra main when we are leaving the sandbox level

def main():

    # Data Preparation of TRaining Data
    dataPreperator = DataPreparation.DataPreparation()
    trainingGenomeData = list()
    organizedData = dataPreperator.organizeTrainingDataByHost(trainingGenomeData)

    # Preprocessing of training data
    dataPreprocessor = DataPreprocessing.DataPreprocessing()
    features = []
    for host, sequences in organizedData.items():
        for seq in sequences:
            kmerFrequencies = dataPreprocessor.extractKmerFrequencies(seq, k=4)
            codonUsage = dataPreprocessor.calculateCodonUsage(seq)
            dicodonUsage = dataPreprocessor.calculateDicodonUsage(seq)
            dinucleotideFrequencies = dataPreprocessor.computeDinucleotideFrequencies(seq)
            featureVector = {**kmerFrequencies, **codonUsage, **dicodonUsage, **dinucleotideFrequencies}
            features.append((featureVector, host))
            # TODO: make use of dataPreprocessor.bindFeaturesToMatrix(featureVector)

    # Classification (Training and Testing)
    trainer = ModelTraining.ModelTraining()
    classifier = HostClassification.HostClassification()
    puModels = {}

    # Train all models
    for host in organizedData.keys():
        # TODO: Make use of trainer.extractRelevantTrainingData(features, host) for easier handling
        positiveSamples = list()  # TODO: Collect positive samples
        unlabeledSamples = list() # TODO: Collect unlabeled samples
        labels = [1] * len(positiveSamples) + [0] * len(unlabeledSamples)
        allFeatures = positiveSamples + unlabeledSamples
        puModel = trainer.trainPUClassifier(allFeatures, labels)
        puModels[host] = puModel

    # Execute the classifiers for all query sequences
    querySequences = []  # TODO: Add first query genome sequences for little testing here

    for query in querySequences:
        # Call the method to classify the query sequence with the trained classifiers
        hostPrediction, confidenceScore = classifier.classifyWithMultipleModels(query, puModels)
        
        print(f"Query Sequence: {query}")
        print(f"Predicted Host: {hostPrediction}, Confidence Score: {confidenceScore}")

if __name__ == "__main__":
    main()