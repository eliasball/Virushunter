# Virushunter
## Main goal
Given a nucleotide sequence, determine the probable hosts.

The implemented pipeline can be found in `code/trainingProcess.ipynb`. It includes the documentation of the exectuted steps. The remaining code for data download, preprocessing, and training can be found in the relevant folders in `code`.

The requirements for running the code are Python version `3.12` and the packages specified in `requirements.txt`.

The used genomic data in `viral_genomes` was downloaded from the NCBI Genome Dataset. The generated feature files used for training will appear in `features` when running the pipeline.