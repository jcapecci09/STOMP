"""Gibbs-based Algorithm for Motif Prediction (GAMP)

Reads DNA sequences and motif length from an input file and finds the
most likely set of motifs using a gibbs sampler, motif profiles, 

Author: Jimmy Capecci
"""

import argparse
import os
import math
from random import randint, random
import logomaker
import matplotlib.pyplot as plt
from typing import Literal
import math
from copy import deepcopy
import time


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


def update_profile(
    count_profile: dict[str, dict[int, int]],
    motif: str,
    num_of_motifs: int,
    action: Literal["add", "remove"],
) -> tuple[dict[str, dict[int, float]], dict[str, dict[int, int]]]:
    """Update a motif count profile and rebuild the probability profile.

    :param count_profile: Count table keyed by base then position.
    :param motif: The motif being added or removed from the profile.
    :param num_of_motifs: Number of motifs currently in the set.
    :param action: Either "add" or "remove" to indicate the update direction.
    :raises ValueError: If action is not "add" or "remove".
    :return: A tuple containing the updated probability profile and updated count profile.
    """

    # if action is remove set variables to remove from profile
    if action == "remove":
        change = -1
        num_of_motifs -= 1

    # if action is add set variables to add to profile
    elif action == "add":
        change = 1
        num_of_motifs += 1

    # raise error if action isnt add or remove
    else:
        raise ValueError("action must be 'add' or 'remove'")

    # perform change to count profile
    for position, base in enumerate(motif):
        count_profile[base][position] += change

    # intialize empty perc_profile
    perc_profile = {}

    # caluclate perc profile for each base
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

    # intialize empty probability dict and sum
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

    # normalize the probabilities so they are floats between [0, 1]
    kmer_prob_normalized = {key: value / prob_sum for key, value in kmer_prob_dict.items()}

    # sort probabilities in ascending order
    sorted_kmer_prob = dict(sorted(kmer_prob_normalized.items(), key=lambda x: x[1]))

    # draw random float between [0, 1]
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

    # for each position calculate shannon entropy for each base
    # add to total
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
            

def parse_k(value: str) -> list[int]:
    """Parse user input K to return a list of K to explore

    :param value: range of k in string
    :return: list of integers k's
    """
    if "-" in value:
        start, end = map(int, value.split("-"))
        return range(start, end + 1)
    return [int(value)]


def make_logo(best_motifs: list[str], output: str):
    """makes motif logo and saves it to output

    :param best_motifs: motifs to make logo
    :param output: output of graph
    """

    # Convert the motif alignment to an information-content matrix in bits
    information = logomaker.alignment_to_matrix(
        best_motifs,
        to_type="information"
    )

    # Create a sequence logo from the information-content matrix
    logo = logomaker.Logo(
        information,
        color_scheme="classic"
    )

    # set y label to bits and y range 0 to 2
    logo.ax.set_ylabel("Bits")
    logo.ax.set_ylim(0, 2)

    # save figure to output
    plt.savefig(
        output,
        format="svg",
        bbox_inches="tight"
    )

    plt.close()


def main():
    print('Finding motifs')
    

    # set up parser and add input and output arguements
    parser = argparse.ArgumentParser(description='Gibbs-based Algorithm for Motif Prediction')
    parser.add_argument('-i', '--input', required=True, help='input file')
    parser.add_argument('-o', '--output', required=False, help='output file')
    parser.add_argument('-s', '--scoring_function', required=False, default='Entropy', 
                        help='Defines scoring function used')
    parser.add_argument('-k', '--kmer_size', required=False, default='8',
                        help='Size of motifs')
    parser.add_argument('--iterations', required=False, help='Number of iterations loop makes')
    parser.add_argument('-t', '--temp', required=False, type=int, default=10, help='Temperature of simulated annealing algorithm')

    # parse arguements
    args = parser.parse_args()

    # intialize empty list to hold sequences
    dna_seq = []
    current_seq = []
    headers = []

    # with input file do the following
    with open(args.input, "r") as f:

        # for each line strip new lines
        for line in f:
            line = line.strip()

            # if line is a header grab
            if line.startswith(">"):
                headers.append(line) #  grab the header
                if current_seq:
                    dna_seq.append("".join(current_seq)) # join line and add to dna_seq
                    current_seq = [] # reintialize empty current_seq

            # if line isnt a header add to current seq list
            else:
                current_seq.append(line)

        # account for edge case, the last sequence
        if current_seq:
            dna_seq.append("".join(current_seq))

    # set default iterartions if not given
    if args.iterations is None:
        iterations = len(dna_seq) * 10 # ten times the size of totla sequences
    else:
        iterations = int(args.iterations) # turn iterartions into a integer

    # set default output if not given
    if args.output is None:
        os.makedirs('Results', exist_ok=True) 
        output = 'Results/motifs'
    else:
        output = args.output

    if args.temp <= 0:
        raise ValueError("--temp must be greater than 0")

    # standardize cooling rate based off of temp and iterations
    # ensures final temp is 0.1
    cool = (0.1 / args.temp) ** (1 / iterations)

    # parse input k to  make it a list
    k_list = parse_k(args.kmer_size)

    # use every value of k and do the following
    for k in k_list:
        start = time.perf_counter()

        # intialize empty list of best_motifs
        best_motifs = []

        # initalize starting motifs
        for seq in dna_seq:
            random_start = randint(0, len(seq) - k)
            best_motifs.append(seq[random_start:random_start+k])

        # set up starting profile and starting motifs
        current_motifs = best_motifs.copy()
        current_profile, c_profile = profile(current_motifs, k)

        # initalize current and best scores
        current_score = score(current_profile, k, current_motifs, args.scoring_function)
        best_score = current_score
        best_motifs = current_motifs.copy()

        # for each iteration do the following
        for i in range(iterations):

            # keep a snapshot of the current state for rollback on rejected moves
            prev_motifs = current_motifs.copy()
            prev_profile = deepcopy(current_profile)
            prev_count_profile = deepcopy(c_profile)

            # set random index
            random_index = randint(0, len(current_motifs) - 1)

            # create progress bar
            progress = (i + 1) / iterations
            bar_width = 40
            filled = int(progress * bar_width)
            bar = '#' * filled + '-' * (bar_width - filled)
            print(f'\r[k={k}] [{bar}] {i + 1}/{iterations} {progress * 100:.1f}%', end='', flush=True)

            # Remove old motif
            remove_motif = current_motifs.pop(random_index)

            # set new profiles without motif
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

            # update profile with new motif
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

            # Simulated annealing algorithm:
            # always accept better scores
            if new_score < current_score:
                current_score = new_score
                if new_score < best_score:
                    best_motifs = current_motifs.copy()
                    best_score = new_score

            # sometimes accept worse scores
            else:
                delta = new_score - current_score # take delta score
                prob = math.exp(-delta / args.temp) # find acceptance probability
                if random() < prob: # draw random floats and sometiems accept worse scores
                    current_score = new_score
                    if new_score < best_score:
                        best_motifs = current_motifs.copy()
                        best_score = new_score

                # else keep previous motifs
                else:
                    current_motifs = prev_motifs
                    current_profile = prev_profile
                    c_profile = prev_count_profile

            # lower temperature to be less exploratory 
            args.temp *= cool

        print()

        # create output string wiuth k value attached to end
        output_string = output + f'_k{str(k)}'

        # write best motifs to ouput file
        with open(output_string + '.txt', 'w') as o:
            for motif in best_motifs:
                o.write(f'{motif}\n')
        make_logo(best_motifs=best_motifs, output=output_string + '.svg')   

        end = time.perf_counter()
        print(f'Time to run: {end-start:.2f} seconds')     

if __name__ == '__main__':
    main() 
