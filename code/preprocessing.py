import os
from dataPreprocessing import DataPreprocessing

featureFilePath = "genomesWithFeatures.csv"
chosenFeatures = { "kmer": [2, 3, 6, 8]}

hosts = [f.name for f in os.scandir("../viral_genomes") if f.is_dir()]
filesForHosts = {host: [f.name for f in os.scandir("../viral_genomes/" + host) if f.is_file()] for host in hosts}

with open(featureFilePath, "a") as featureFile:
    dataPreprocessor = DataPreprocessing.DataPreprocessing()
    for host in hosts:
        print ("Processing host: " + host)
        for filename in filesForHosts[host]:
            with open("../viral_genomes/" + host + "/" + filename, "r") as file:
                lines = file.readlines()
                sequence = "".join(lines[1:]).replace("\n", "")
                # flatten feature dict to list
                features = [x for v in dataPreprocessor.extractFeaturesFromGenome(sequence, chosenFeatures).values() for x in v]
                # write as csv
                featureFile.write(host + "," + filename.split(".")[0] + "," + ",".join(map(str, features)) + "\n")