from itertools import product
import pandas as pd

#  Preprocessing of dataobjects (extraction of valueable features out of given genome sequence)


class DataPreprocessing:

    def __init__(self):
        #Initializes this class
        pass

    def extractKmerFrequencies(self, genomeSequence: str, k: int) -> tuple[list[float], float]:
        """
        Description: Extracts relative k-mer frequencies from a genome sequence for a given k. Handles invalid k-mers separately in counter.
        Args: DNA sequence and size of the k-mers
        Returns: Array of relative K-mer frequencies (indexed) and the relative frequency of invalid k-mers.
        """

        # Change RNA to DNA if it is necessary
        if 'U' in genomeSequence:
            genomeSequence = genomeSequence.replace('U', 'T')  
        
        # Calculate the reverse complement of the given sequence
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_genome_sequence = ''
        for base in reversed(genomeSequence):
            if base in complement:
                reverse_genome_sequence += complement[base]
            else:
                reverse_genome_sequence += base   # N is complemetary to N 
        
        # Instanciate Dict to store the kmer Frequencies and initialize counter for invalid (within N) and total kmers
        bases = ['A', 'G', 'C', 'T']    
        kmer_combinations = [''.join(kmer) for kmer in product(bases, repeat=k)]    #list of all possible combinations of bases
        kmer_frequencies = [0] * len(kmer_combinations)     #list of 0 to initiate ocurrence of each bases combination
        kmer_index = {kmer: idx for idx, kmer in enumerate(kmer_combinations)}   #dict with kmer as key and index as value
        total_kmers_count = 0       
        invalid_kmer_count = 0

        # count k-mers of the genome string in 6-frame read style 
        for i in range(len(genomeSequence) - k + 1):
            
            # Count Forward String
            kmer = genomeSequence[i:i + k]
            if kmer in kmer_index:  
                kmer_frequencies[kmer_index[kmer]] += 1
                total_kmers_count += 1
            else:  
                invalid_kmer_count += 1
                total_kmers_count += 1

            # Count Backward String
            kmerReverse = reverse_genome_sequence[i:i + k]
            if kmerReverse in kmer_index:  
                kmer_frequencies[kmer_index[kmerReverse]] += 1
                total_kmers_count += 1
            else:  
                invalid_kmer_count += 1
                total_kmers_count += 1

        # Convert absolute to relative frequencies
        if total_kmers_count > 0:
            kmer_frequencies = [freq / total_kmers_count for freq in kmer_frequencies]    #calculate ratio for each entry in kmer_frequencies
            invalid_kmer = invalid_kmer_count / total_kmers_count
        else:
            kmer_frequencies = [0.0] * len(kmer_frequencies)
            invalid_kmer = 0.0

        # Return array of relative K-mer frequencies (lexiographically ordered) and the relative frequency of invalid k-mers
        return kmer_frequencies, invalid_kmer


    def calculateCodonTranslation(self, genomeSequence: str) -> tuple[list[int], int]:
        """
        Description: Calculates translation of codon frequencies.
        Args: DNA or RNA sequence
        Returns: Array of size 22 where each index corresponds to a specific translation of a codon (aminoacid or start/stop)
                 in alphabetical order (but start (M) and stop (Z) in the end) and the relative frequency of invalid codons.
        """
        
        # Calculate codon combinations
        bases = ['A', 'G', 'C', 'T']
        codon_combinations = [''.join(kmer) for kmer in product(bases, repeat=3)]

        # Calculate codon frequencies with k=3 
        codon_frequencies, invalid_codon = self.extractKmerFrequencies(genomeSequence, 3)


        # Translate codon frquencies into their translated frequencies (Define Stop Codon as Z and rest Codons as usual)
        codons_translated = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                                'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                                'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                                'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                                'TAA': 'Z', 'TAC': 'Y', 'TAG': 'Z', 'TAT': 'Y',
                                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                                'TGA': 'Z', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                                'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
        
        # Translate codons into their translation in ribosomes
        translated_frequencies = {codon: 0 for codon in 'ACDEFGHIKLNPQRSTVWXYMZ'}
        for codon, frequency in zip(codon_combinations, codon_frequencies):
            codon_translation = codons_translated.get(codon, None)
            if codon_translation:
                translated_frequencies[codon_translation] += frequency
        
        # Normalize again (translation messed normalization of kmer function)
        total_frequencies = sum(translated_frequencies.values()) + invalid_codon
        if total_frequencies > 0:
            normalized_frequencies = [freq / total_frequencies for freq in translated_frequencies.values()]
        else:
            normalized_frequencies = [0.0] * len(translated_frequencies)
        invalid_codon = invalid_codon / total_frequencies

        # Return array where each index corresponds to a specific translation of a codon and the relative frequency of invalid codons
        return normalized_frequencies, invalid_codon
      

    def calculateDicodonTranslation(self, genomeSequence: str) -> dict[str, int]:
        """
        Description: Calculates di-codons from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Returns: Array of size 2^22 where each index corresponds to a specific translation of a dicodon (aminoacid or start/stop)
                 in alphabetical order (but start (M) and stop (Z) in the end) and the relative frequency of invalid dicodons.
        """

        # Calculate codon combinations
        bases = ['A', 'G', 'C', 'T']
        codon_combinations = [''.join(kmer) for kmer in product(bases, repeat=6)]

        # Calculate dicodon frequencies with k=6
        dicodon_frequencies, invalid_dicodon = self.extractKmerFrequencies(genomeSequence, 6)


        # Translate codon frquencies into their translated frequencies (Define Stop Codon as Z and rest Codons as usual)
        codons_translated = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                                'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                                'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                                'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                                'TAA': 'Z', 'TAC': 'Y', 'TAG': 'Z', 'TAT': 'Y',
                                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                                'TGA': 'Z', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                                'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
        
        # Initialize dicodon frequency dictionary (combinations of two consecutive codons) and translate dicodons
        translated_dicodon_frequencies = {f"{codon1}{codon2}": 0 for codon1 in 'ACDEFGHIKLNPQRSTVWXYMZ' for codon2 in 'ACDEFGHIKLNPQRSTVWXYMZ'}
        for dicodon, frequency in zip(codon_combinations, dicodon_frequencies):
            mid = len(dicodon) // 2
            codon1_translation = codons_translated.get(dicodon[:mid], None)
            codon2_translation = codons_translated.get(dicodon[mid:], None)
            dicodon_translation = codon1_translation + codon2_translation
            if codon1_translation and codon2_translation:
                translated_dicodon_frequencies[dicodon_translation] += frequency
        
        # Normalize again (translation messed normalization of kmer function)
        total_frequencies = sum(translated_dicodon_frequencies.values()) + invalid_dicodon
        if total_frequencies > 0:
            normalized_frequencies = [freq / total_frequencies for freq in translated_dicodon_frequencies.values()]
        else:
            normalized_frequencies = [0.0] * len(translated_dicodon_frequencies)
        invalid_dicodon = invalid_dicodon / total_frequencies

        # Return array where each index corresponds to a specific translation of a dicodon and the relative frequency of invalid dicodons
        return normalized_frequencies, invalid_dicodon


  
    def extractFeaturesFromGenome(self, genomeSequence: str, chosenFeatures: list[str]) -> dict[str, list[float]]:
        """
        Description: Extracts selected features for one genome sequence based on chosenFeatures.
        Args: genomeSequence (str): A DNA sequence. chosenFeatures (List[str]): Dict of feature types to extract, e.g., {"kmer": [2,3,5], "codon": [], "dicodon": []}.
        Returns: Dict[str, List[float]]: A dictionary where each key is a feature type (e.g., "Kmer5", "Codon") and the value is a list of corresponding feature values.
        """
        feature_dict = {}

        # Check and extract k-mers for all specified k values
        if "kmer" in chosenFeatures:
            for k in chosenFeatures["kmer"]:
                kmer_label = f"Kmer{k}"
                kmer_frequencies, invalid_kmer_ratio = self.extractKmerFrequencies(genomeSequence, k)
                feature_dict[kmer_label] = kmer_frequencies + [invalid_kmer_ratio]

        # Extract codon frequencies
        if "codon_translation" in chosenFeatures:
            codon_frequencies, invalid_codon_ratio = self.calculateCodonTranslation(genomeSequence)
            feature_dict["Codon"] = codon_frequencies + [invalid_codon_ratio]

        # Extract dicodon frequencies
        if "dicodon_translation" in chosenFeatures:
            dicodon_frequencies, invalid_dicodon_ratio = self.calculateDicodonTranslation(genomeSequence)
            feature_dict["Dicodon"] = dicodon_frequencies + [invalid_dicodon_ratio]

        return feature_dict

    
    
    def calculateFeatureMatrix(self, genomes: list[tuple[str, str, str]], chosenFeatures: dict[str, list]) -> pd.DataFrame:
        """
        Description: Calculates features of multiple genomes and structures them in a feature matrix based on chosenFeatures that will be the input of the classifiers.
        Args: genomes (List[tuple[str, str, str]]): List of tuples where each tuple contains (hostname, virusname, sequence).
            chosenFeatures (Dict[str, list]): Dictionary specifying which features to extract and their parameters.
        Returns: pd.DataFrame: Matrix with datapoints (genomes) as rows and columns 'virusname', 'feature_vector', 'feature_vector_flat' and 'hostname'.
        """
        feature_matrix = []

        for genome_tuple in genomes:
            hostname, virusname, genome_sequence = genome_tuple

            # Calculation of the features
            features_dict = self.extractFeaturesFromGenome(genome_sequence, chosenFeatures)
            features_flat = []
            for key, value in features_dict.items():
              if isinstance(value, list):  
                  if all(isinstance(v, float) for v in value):  # Make sure List contains just floats
                      features_flat.extend(value)
                  else:
                      raise ValueError(f"Unexpected type in list for key '{key}': {value}")
              elif isinstance(value, float):
                  features_flat.append(value)
              else:
                  raise ValueError(f"Unexpected value type for key '{key}': {value}")
            
            # Add Data to DataFrame
            feature_matrix.append({
                'HostTaxID': hostname,
                'VirusTaxID': virusname,
                'FeaturesDictionary': features_dict,
                'FeaturesFlat': features_flat
            })

        # Erstelle das DataFrame
        feature_matrix_df = pd.DataFrame(feature_matrix)
        
        return feature_matrix_df


    def saveFeatureMatrixAsCSV(self, feature_matrix_df: pd.DataFrame, file_path: str):
        """
        Description: Saves a feature matrix DataFrame to a CSV file
        Args: feature_matrix_df (pd.DataFrame): The feature matrix to be saved.
              file_path (str): The path to the CSV file where the data will be saved.
        Returns: None
        """
        feature_matrix_df.to_csv(file_path, index=False)
        print(f"Feature matrix successfully saved to {file_path}")

    
if __name__ == "__main__":

    # Generate Mock data
    chosenFeatures = { "kmer": [2, 3, 5], "codon_translation": [], "dicodon_translation" : []}

    genomes = [("hosttaxID1", "virustaxID1", "AAAAAAAAAAAAAA"), ("hosttaxID1", "virustaxID2", "TGCATGCATGCAT"),
               ("hosttaxID1", "virustaxID3", "AAAAAAAAATTTTTTTTTTTTGGGGGGGGGGGGGCCCCCCCCCCCATATATATATATATATA"),
               ("hosttaxID2", "virustaxID4", "TGCATGCATGCATTGTGACGACTGAACGTACGTGCACGATGCAATGTTGCCGAATTGGCAATGC"),
               ("hosttaxID2", "virustaxID5", "AUGAUGAUGAUGAUGAUGAUGAUG")]

    preprocessor = DataPreprocessing()
    feature_dataframe = preprocessor.calculateFeatureMatrix(genomes,chosenFeatures)

    #Print the results for the mock data
    for index, row in feature_dataframe.iterrows():

      genome = genomes[index]
      features = row['FeaturesDictionary']

      print(f"Host: {row['HostTaxID']}, Virus: {row['VirusTaxID']}, Sequence: {genome[2]}")
      for feature_name, feature_values in features.items():
        print(f"{feature_name}: {feature_values}")
      print(f"Features (Flatten): {row['FeaturesFlat']}")
      print(len(row['FeaturesFlat'])) # Make sure that every flat list has the same length
      print("-" * 100) #Split the different genomes visually

    #Save features in a csv
    preprocessor.saveFeatureMatrixAsCSV(feature_dataframe, "code/mock_features_preprocessing.csv")

