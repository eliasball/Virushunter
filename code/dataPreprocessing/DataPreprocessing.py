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
        from itertools import product

        bases = ['A', 'G', 'C', 'T']    
        kmer_combinations = [''.join(kmer) for kmer in product(bases, repeat=k)]    #list of all possible combinations of bases

        
        kmer_frequencies = [0] * len(kmer_combinations)     #list of 0 to initiate ocurrence of each bases combination

        kmer_index = {kmer: idx for idx, kmer in enumerate(kmer_combinations)}   #dict with kmer as key and index as value

        total_kmers_count = 0         #initialise counter to track kmers to calculate ratio 
        invalid_kmer_count = 0

        # count k-mers of the genome string in 3-frame read style 
        for i in range(len(genomeSequence) - k + 1):
            kmer = genomeSequence[i:i + k]  
            if kmer in kmer_index:  
                kmer_frequencies[kmer_index[kmer]] += 1
                total_kmers_count += 1
            else:  
                invalid_kmer_count += 1
                total_kmers_count += 1

        if total_kmers_count > 0:
            kmer_frequencies = [freq / total_kmers_count for freq in kmer_frequencies]    #calculate ratio for each entry in kmer_frequencies
            invalid_kmer = invalid_kmer_count / total_kmers_count
        else:
            kmer_frequencies = [0.0] * len(kmer_frequencies)
            invalid_kmer = 0.0

        return kmer_frequencies, invalid_kmer


    def calculateCodon(self, genomeSequence: str) -> tuple[list[int], int]:
        """
        Description: Calculates codon frequencies for all possible combinations of 3 bases (A, G, C, T).
                 If the input is RNA (contains 'U'), it is converted to DNA.
                 If the sequence length is not divisible by 3, it adds extra bases ('N') to make it valid.
        Args: DNA or RNA sequence
        Returns: Array of size 64 where each index corresponds to a specific codon in alphabetical order and a count of invalid codons.
        """
        
        from itertools import product

    
        bases = ['A', 'G', 'C', 'T']
        codon_combinations = [''.join(codon) for codon in product(bases, repeat=3)]

        # Initialise with 0
        codon_frequencies = [0] * len(codon_combinations)

        # enumerates creates index for each entry in codon_combinations and with comprehension and codon as key creates dictonary
        codon_to_index = {codon: idx for idx, codon in enumerate(codon_combinations)}

        total_codons_count = 0
        invalid_codon_count = 0

        # change RNA to DNA 
        if 'U' in genomeSequence:
            genomeSequence = genomeSequence.replace('U', 'T')  


        # count codons in genome string 
        for i in range(len(genomeSequence) - 3 + 1):
            codon = genomeSequence[i:i + 3]  
            if codon in codon_to_index:
                codon_frequencies[codon_to_index[codon]] += 1
                total_codons_count += 1
            else:
                invalid_codon_count += 1
                total_codons_count += 1

        if total_codons_count > 0:
            codon_frequencies = [freq / total_codons_count for freq in codon_frequencies]
            invalid_codon = invalid_codon_count / total_codons_count
        else:
            codon_frequencies = [0.0] * len(codon_frequencies)
            invalid_codon = 0.0

        return codon_frequencies, invalid_codon

    def specificCodonFrequency(self, genomeSequence: str, targetCodon: str) -> tuple[str,float]:
        """
        Description: Extracts the frequency of a specific codon from the genome sequence.
        Args: genomeSequence (str): DNA sequence
        targetCodon (str): The codon whose frequency is to be calculated.
        Returns: int: Frequency of the specified codon.
        """
        #get all frequencies of all codons
        codon_frequencies = self.calculateCodon(genomeSequence)
        
        #extract the total frequency of target codon
        target_frequency = codon_frequencies.get(targetCodon, 0)
    
        # calculate realtive frequency of target codon
        total_codons_count = sum(codon_frequencies.values())

        if total_codons_count > 0:
            relative_ratio = target_frequency / total_codons_count
        else:
            relative_ratio = 0.0  

        #return targetCodon, relative_ratio
        pass
        
    

    def calculateDicodon(self, genomeSequence: str) -> dict[str, int]:
        """
        Description: Calculates di-codons from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Di-Codon usage  usage frequencies as a (indexed) list/dictonary
        """

        stable_k = 6

        #if len(genomeSequence) % stable_k != 0:     #check that sequence can be divided by the codon
            #print ("Some genome Sequence at the end of the string can not me matched due to length of Dicodon and genome string .")
            #rest = len(genomeSequence) % stable_k
            #print("The last", rest, "bases will not the checked for all possible codons in 3-frame read manner.")

        return self.extractKmerFrequencies(genomeSequence, stable_k)

    def computeDinucleotideFrequencies(self, genomeSequence: str) -> dict[str, int]:
        """ 
        Description: Computes dinucleotide frequencies from the genome sequence
        Args: genomeSequence (str): DNA sequence
        Returns: Dinucleotide frequencies as a (indexed) list/dictonary
        """
        stable_k = 2

        #if len(genomeSequence) % stable_k != 0:     #check that sequence can be divided by the kmer
            #print ("Some genome Sequence at the end of the string can not me matched due to length of 2-mer and genome string.")
            #rest = len(genomeSequence) % stable_k
            #print("The last", rest, "bases will not the checked for all possible 2-mers in 3-frame read manner.")

        return self.extractKmerFrequencies(genomeSequence, stable_k)
    


    def choosenFeaturesoneGenome(self, genomeSequence: str, chosenFeatures: dict[str, list]) -> list[float]:
        """
        Description: Extracts selected features for one genome sequence based on chosenFeatures.
        Args: A DNA sequence, chosenFeatures specifying which features to extract and their parameters.
            E.g.: {"kmer": [2], "codon": [], "dicodon": [], "dinucleotide": []}
        Returns: A list of numerical features extracted from the genome, including invalid k-mer and codon frequencies.
        """
        feature_vector = []

        if "kmer" in chosenFeatures:
            k = chosenFeatures["kmer"][0] if chosenFeatures["kmer"] else 2  
            kmer_frequencies, invalid_kmer_ratio = self.extractKmerFrequencies(genomeSequence, k)
            feature_vector += kmer_frequencies  
            feature_vector.append(invalid_kmer_ratio)  


        if "codon" in chosenFeatures:
            codon_frequencies, invalid_codon_ratio = self.calculateCodon(genomeSequence)
            feature_vector += codon_frequencies  
            feature_vector.append(invalid_codon_ratio)  

        # pass
        if "specific_codon" in chosenFeatures:
            target_codon = chosenFeatures["specific_codon"][0]  
            _, specific_codon_ratio = self.specificCodonFrequency(genomeSequence, target_codon)
            feature_vector.append(specific_codon_ratio)

      
        if "dicodon" in chosenFeatures:
            dicodon_frequencies, invalid_dicodon_ratio = self.calculateDicodon(genomeSequence)  
            feature_vector += dicodon_frequencies  
            feature_vector.append(invalid_dicodon_ratio)  

    
        if "dinucleotide" in chosenFeatures:
            dinucleotide_frequencies, invalid_dinucleotide_ratio = self.computeDinucleotideFrequencies(genomeSequence) 
            feature_vector += dinucleotide_frequencies  
            feature_vector.append(invalid_dinucleotide_ratio)  

        return feature_vector

    
    
    def bindFeaturesToMatrix(self, genomes: dict[str, str], chosenFeatures: dict[str, list]) -> list[list[float]]:
        """
        Description: Binds features of multiple genomes into a feature matrix based on chosenFeatures that will be the input of the classifiers
        Args: Dictionary of genome names (keys) and their sequences (values).
              Dictionary specifying which features to extract and their parameters.
        Returns: Matrix with datapoints (genomes) as rows and the relevant features as columns (-> Structure: (virusnames, input features, hostnames))
        """
        feature_matrix = []

        for genome_name, genome_sequence in genomes.items():
            # for one genome
            feature_vector = self.choosenFeaturesoneGenome(genome_sequence, chosenFeatures)
            # add feature_vector to matrix
            feature_matrix.append(feature_vector)

        return feature_matrix



if __name__ == "__main__":

    genomeSequence = "agctaaggcct"
    k = 2
    chosenFeatures = { "kmer": [3], "codon": [], "dinucleotide": []}

    genomes = genomes = {
    "genome1": "agctaaggcct",
    "genome2": "ttgcaaggtcc",
    "genome3": "ccgtaaggcat"
    }



    preprocessor = DataPreprocessing()
    kmer_frequencies = preprocessor.extractKmerFrequencies(genomeSequence, k)
    codons = preprocessor.calculateCodon(genomeSequence)
    #specific_codon = preprocessor.specificCodonFrequency(genomeSequence, targetCodon="agc")
    dicodons = preprocessor.calculateDicodon(genomeSequence)
    dinucleotide_frequencies = preprocessor.computeDinucleotideFrequencies(genomeSequence)
    choosen_feature = preprocessor.choosenFeaturesoneGenome(genomeSequence,chosenFeatures)
    feature_matrix = preprocessor.bindFeaturesToMatrix(genomes,chosenFeatures)

    #Print returns
    print("K-mer frequencies:", kmer_frequencies)
    print("Codon:", codons)
    #print("Specific codon:",specific_codon)
    print("Dicodon:",dicodons)
    print("DinucleotideFrequencies :",dinucleotide_frequencies)
    print("choosen_feature_for_genome:",choosen_feature)
    print("matrix_for_choosen_feature_for_genomes:",feature_matrix)


