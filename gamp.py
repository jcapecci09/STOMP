"""Gibbs-based Algorithm for Motif Prediction (GAMP)

Reads DNA sequences and motif length from an input file and finds the
most likely set of motifs using a gibbs sampler, motif profiles, 

Author: Jimmy Capecci
"""

import argparse
import math
from random import randint, random
import logomaker
import matplotlib.pyplot as plt


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
    count_prof = {
        'A': {i: 1 for i in range(k)},
        'T': {i: 1 for i in range(k)},
        'C': {i: 1 for i in range(k)},
        'G': {i: 1 for i in range(k)}
    }

    # for each motif in the list of motifs
    # add count to profile
    for motif in motifs:
        for pos, base in enumerate(motif):
            count_prof[base][pos] += 1

    # find length of profile
    total = len(motifs) + 4
    perc_prof = {
        base: count_prof[base].copy()
        for base in count_prof
    }

    # for each base in profile find percentage of occurence
    for base in perc_prof:
        for pos in perc_prof[base]:
            perc_prof[base][pos] /= total

    return perc_prof, count_prof


def update_profile(count_profile, motif, num_of_motifs, action):
    if action == "remove":
        change = -1
        num_of_motifs -= 1

    elif action == "add":
        change = 1
        num_of_motifs += 1

    else:
        raise ValueError("action must be 'add' or 'remove'")

    for position, base in enumerate(motif):
        count_profile[base][position] += change

    perc_profile = {}

    for base in count_profile:
        perc_profile[base] = {}

        for position in count_profile[base]:
            perc_profile[base][position] = (
                count_profile[base][position] /
                (num_of_motifs + 4)
            )

    return perc_profile, count_profile


def gibbs_sampler(prof: dict[str, dict[int, float]], seq: str, k:int) -> str:
    """Finds the motif with the highest probability according to a profile.

    :param prof: nested dictionary containing the motif profile
    :param seq: sequence to find the best motif
    :param k: size of motifs
    :return: motif with the highest probability
    """

    # 
    kmer_prob_dict = {}
    prob_sum = 0

    # form each motif in sequence find the best motif
    for i in range(len(seq) - k + 1):

        # grab kmer and add to list
        kmer = seq[i:i+k]

        # intialize prob = to 1
        prob =  1

        # find the likelhood of current kmer
        # that it matches the profile
        for pos, base in enumerate(kmer):
            prob *= prof[base][pos]
        kmer_prob_dict[kmer] = prob
        prob_sum += prob

    # normalize the probabilities so they are between [0, 1]
    kmer_prob_normalized = {key: value / prob_sum for key, value in kmer_prob_dict.items()}

    # sort probabilities in ascending order
    sorted_kmer_prob = dict(sorted(kmer_prob_normalized.items(), key=lambda x: x[1]))

    # draw random number between [0, 1]
    draw = random()
    cumulative = 0.0

    # probabilistically pick a kmer
    for kmer, prob in sorted_kmer_prob.items():
        cumulative += prob
        if draw <= cumulative:
            return kmer

    # fallback if rounding makes the draw fall just above the end
    return max(sorted_kmer_prob, key=sorted_kmer_prob.get)


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
    """find score of set of motifs

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
            

def parse_k(value):
    if "-" in value:
        start, end = map(int, value.split("-"))
        return range(start, end + 1)
    return [int(value)]


def make_logo(best_motifs, output):
    information = logomaker.alignment_to_matrix(
        best_motifs,
        to_type="information"
    )

    logo = logomaker.Logo(
        information,
        color_scheme="classic"
    )

    logo.ax.set_ylabel("Bits")
    logo.ax.set_ylim(0, 2)

    plt.savefig(
        output,
        format="svg",
        bbox_inches="tight"
    )

    plt.close()


def main():

    # set up parser and add input and output arguements
    parser = argparse.ArgumentParser(description='Greedy Motif Search')
    parser.add_argument('-i', '--input', required=True, help='input file')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-s', '--scoring_function', required=False, default='Entropy', 
                        help='Defines scoring function used')
    parser.add_argument('-k', '--kmer_size', required=False, default='8',
                        help='Size of motifs')
    parser.add_argument('--iterations', required=False, help='Number of iterations loop makes')

    # parse arguements
    args = parser.parse_args()

    dna_seq = []
    current_seq = []

    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if current_seq:
                    dna_seq.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            dna_seq.append("".join(current_seq))

    if args.iterations is None:
        iterations = len(dna_seq) * 10
    else:
        iterations = int(args.iterations)

    # initalize list of best motifs and best score
  
    k_list = parse_k(args.kmer_size)

    for k in k_list:
        best_motifs = []
        best_score = float('inf')

        # initalize starting motifs
        for seq in dna_seq:
            random_start = randint(0, len(seq) - k)
            best_motifs.append(seq[random_start:random_start+k])

        current_motifs = best_motifs.copy()
        current_profile, c_profile = profile(current_motifs, k)

        for i in range(iterations):
            random_index = randint(0, len(current_motifs) - 1)

            progress = (i + 1) / iterations
            bar_width = 40
            filled = int(progress * bar_width)
            bar = '#' * filled + '-' * (bar_width - filled)
            print(f'\r[k={k}] [{bar}] {i + 1}/{iterations} {progress * 100:.1f}%', end='', flush=True)

            # Remove old motif
            remove_motif = current_motifs.pop(random_index)

            current_profile, c_profile = update_profile(
                c_profile,
                remove_motif,
                len(current_motifs) + 1,
                'remove'
            )

            # Sequence corresponding to removed motif
            seq = dna_seq[random_index]

            # Choose new motif using profile without old motif
            new_motif = gibbs_sampler(current_profile, seq, k)

            # Add new motif
            current_motifs.insert(random_index, new_motif)

            current_profile, c_profile = update_profile(
                c_profile,
                new_motif,
                len(current_motifs) - 1,
                'add'
            )

            # Calculate score
            new_score = score(
                current_profile,
                k,
                current_motifs,
                args.scoring_function
            )

            if best_score > new_score:
                best_motifs = current_motifs.copy()
                best_score = new_score

        print('\r' + ' ' * 80 + '\r', end='')

        output = args.output.split('.')
        output.insert(1, f'_k{str(k)}.')
        output_string = ''.join(output)
        # write best motifs to ouput file
        with open(output_string, 'w') as o:
            for motif in best_motifs:
                o.write(f'{motif}\n')
        make_logo(best_motifs=best_motifs, output=output_string.replace('.txt', '.svg'))        

if __name__ == '__main__':
    main() 
