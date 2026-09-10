"""Greedy motif finding algorithm.

Reads DNA sequences and motif length from an input file and finds the
most likely set of motifs using Greedy Motif Search.

Author: Jimmy Capecci
"""

import argparse
import math


def profile(motifs: list, k: int) -> dict[str, dict[int, float]]:
    """Builds a k-mer profile in a nested dictionary.
    {base: {k-mer position: probability}}

    Example:
    motifs = ['GGC', 'AGG']

    profile = {
        'A': {0: 0.5, 1: 0.0, 2: 0.0},
        'T': {0: 0.0, 1: 0.0, 2: 0.0},
        'C': {0: 0.0, 1: 0.0, 2: 0.5},
        'G': {0: 0.5, 1: 1.0, 2: 0.5}
    }

    :param motifs: list containing the motifs
    :param k: length of each motif
    :return: profile containing the probability of each base at each k-mer position
    """

    # initalize empty profile
    prof = {
        'A': {i: 1 for i in range(k)},
        'T': {i: 1 for i in range(k)},
        'C': {i: 1 for i in range(k)},
        'G': {i: 1 for i in range(k)}
    }

    # for each motif in the list of motifs
    # add count to profile
    for motif in motifs:
        for pos, base in enumerate(motif):
            prof[base][pos] += 1

    # find length of profile
    total = len(motifs) + 4

    # for each base in profile find percentage of occurence
    for base in prof:
        for pos in prof[base]:
            prof[base][pos] /= total

    return prof


def profile_most_probable(prof: dict[str, dict[int, float]], seq: str, k:int) -> str:
    """Finds the motif with the highest probability according to a profile.

    :param prof: nested dictionary containing the motif profile
    :param seq: sequence to find the best motif
    :param k: size of motifs
    :return: motif with the highest probability
    """

    # intialize best percentage at -1
    largest_percentage = -1
    best_motif = ''

    # form each motif in sequence find the best motif
    for i in range(len(seq) - k + 1):

        # find motif
        motif = seq[i:i+k]
        
        # intialize percentage = to 1
        percentage =  1

        # find the likelhood of current motif
        # that it matches the profile
        for pos, base in enumerate(motif):
            percentage *= prof[base][pos]

        # if percentage is greater than the current largest 
        # it is the best motif
        if percentage > largest_percentage:
            largest_percentage = percentage 
            best_motif = motif
    return best_motif

    
def consensus(prof: dict[str, dict[int, float]], k: int) -> str:
    """Builds the consensus motif from a profile

    :param prof: motif profile
    :param k: size of motifs
    :return: consensus motif
    """

    # initialize consensus motif
    consensus_string = ''

    # for each position in k
    for pos in range(k):

        # initalize best base
        best_probability = 0

        # for each base in profile
        for base in prof:

            # find the base with the best probability
            if prof[base][pos] > best_probability:
                best_probability = prof[base][pos]
                best_base = base

        # build consesnsus motif
        consensus_string += best_base

    return consensus_string


def hamming(motif1: str, motif2: str) -> int:
    """Find hamming distance, the number of times two sequences differ

    :param motif1: first sequence
    :param motif2: second sequence
    :return: hamming distance
    """

    # find hamming distance
    # the number of times two sequences differ from one another
    h_distance = 0
    for base1, base2 in zip(motif1, motif2):
        if base1 != base2:
            h_distance += 1
    return h_distance


def entropy(profile: dict[str, dict[int, float]], k) -> int:
    """Entropy scoring function that just uses profile to determine 
    how uncertain 

    :param profile: profile of motifs
    :return: Entropy score
    """
    e_total = 0
    for pos in range(k):
        for base in ('A', 'T', 'C', 'G'):
            prob = profile[base][pos]
            e_total += -prob * math.log2(prob)


    return e_total
    


def score(prof: dict[str, dict[int, float]], k: int, motifs: list, scoring_function: str) -> int:
    """find score of set of motifs using the consensus motif 
    and hamming distance

    :param prof: profile of motifs
    :param k: size of motifs
    :param motifs: set of motifs
    :param scoring_function: determines either consensus or entropy scoring function
    :return: score
    """

    if scoring_function == 'Hamming':
        # find the consensus motif
        consensus_string = consensus(prof, k)
        final_score = 0

        #  for each motif in set of motifs find the final score
        for motif in motifs:
            final_score += hamming(consensus_string, motif)
        return final_score

    else:
        return entropy(prof, k)
            

def main():

    # set up parser and add input and output arguements
    parser = argparse.ArgumentParser(description='Greedy Motif Search')
    parser.add_argument('-i', '--input', required=True, help='input file')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-s', '--scoring_function', required=False, default='Entropy', 
                        help='Defines scoring function used')

    # parse arguements
    args = parser.parse_args()

    # open input file in read mode and find 
    # k and the sequences
    with open(args.input, 'r') as f:
        k = int(f.readline().strip().split(' ')[0])
        dna_seq = []
        for line in f:
            dna_seq.append(line.strip())

    # initalize list of best motifs and best score
    best_motifs = []
    best_score = float('inf')

    # find first sequence
    dna_1 = dna_seq[0]

    # for each kmer in the first sequence
    for index in range(len(dna_1) - k + 1):

        # build inital kmer and add it to current motifs
        initial_kmer = dna_1[index:index+k]
        current_motifs = [initial_kmer]

        # for each sequence except the first sequence
        for seq in dna_seq[1:]:

            # find the current profile
            current_profile = profile(current_motifs, k)

            # find the best motif according to profile and add it to current motifs
            current_motifs.append(profile_most_probable(current_profile, seq, k))

        # find final profile and score 
        final_profile = profile(current_motifs, k)
        current_score = score(final_profile, k, current_motifs, args.scoring_function)

        # if score is better, its the best kmer profile
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs
        
    # write best motifs to ouput file
    with open(args.output, 'w') as o:
        for motif in best_motifs:
            o.write(f'{motif}\n')
    

if __name__ == '__main__':
    main() 
