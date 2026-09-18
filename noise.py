"""Create noisy datasets filled with random sequences to test
STOMP

Author: Jimmy Capecci
"""


from multifasta_collector import mfc
from random import shuffle, randint, choice
import os


def main():
    # List of input motif site files to process. Each file contains a set of
    # sequences for a specific transcription factor binding site family.
    infiles = [
        'Data/MA0001.1.sites',
        'Data/MA0007.2.sites',
        'Data/MA0073.1.sites',
        'Data/MA0106.1.sites',
        'Data/MA0106.2.sites',
    ]

    # Make sure the output folder exists before writing FASTA files.
    os.makedirs('Noise', exist_ok=True)

    # Output names match the motif IDs used in the filenames.
    outputs = ['MA0001.1', 'MA0007.2', 'MA0073.1', 'MA0106.1', 'MA0106.2']

    # Different levels of noisy sequence addition to generate.
    # 10%, 25%, 50%, 75%, and 95% noise fractions.
    noise_percentage = [0.1, 0.25, 0.5, 0.75, 0.95]

    # Store the sequence lists for each input dataset and the min/max lengths.
    list_of_seq = []
    mins = []
    maxs = []

    # Loop over each motif dataset, read its sequences, and compute the range of
    # sequence lengths found in that dataset.
    for infile in infiles:
        # mfc reads the file and returns a dictionary-like collection of IDs to
        # sequences. We convert the values to a list for easier manipulation.
        d = mfc(infile)

        min_length = float('inf')
        max_length = 0
        sequence_list = list(d.values())
        list_of_seq.append(sequence_list)

        # Determine the shortest and longest sequence in this dataset.
        for value in sequence_list:
            value_length = len(value)
            if value_length < min_length:
                min_length = value_length
            if value_length > max_length:
                max_length = value_length

        mins.append(min_length)
        maxs.append(max_length)

    # The synthetic noise sequences are built from a 4-letter alphabet: A, C, G, T.
    aa_list = ['A', 'C', 'G', 'T']

    # For each motif family, generate noisy sequence sets at each noise level.
    for mini, maxi, factor, output in zip(mins, maxs, list_of_seq, outputs):
        for noise in noise_percentage:
            # Start from a copy of the real sequences so we can append synthetic
            # noise without altering the original dataset.
            sequence_list = factor.copy()

            # Number of real sequences currently in the dataset.
            num_of_seqs = len(sequence_list)

            # Estimate how many random noisy sequences to add so that the final set
            # ends up with approximately the requested noise fraction.
            # Formula: new_noise = noise * N / (1 - noise), where N is the number
            # of original sequences.
            num_noise = round((noise * num_of_seqs) / (1 - noise))

            # Generate random sequence lengths within the observed range for this
            # motif family, so the noise sequences have similar size distribution.
            range_list = list(range(mini, maxi + 1))

            noise_collector = []

            # Create the synthetic noisy sequences.
            for _ in range(num_noise):
                random_length = choice(range_list)

                random_sequence = ''.join(
                    choice(aa_list)
                    for _ in range(random_length)
                )

                noise_collector.append(random_sequence)

            # Add all generated synthetic sequences to the real dataset and then
            # shuffle them so the noisy sequences are not clustered at the end.
            sequence_list.extend(noise_collector)
            shuffle(sequence_list)

            # Convert noise percentage like 0.10 into a filename-safe integer such
            # as "10" to produce files like MA0001.1_10.fasta.
            string_noise = str(round(noise * 100))

            # Write the final FASTA file for this motif and noise level.
            with open(f'Noise/{output}_{string_noise}.fasta', 'w') as o:
                for indx, seq in enumerate(sequence_list):
                    o.write(f'>{indx}\n')
                    o.write(f'{seq}\n')


if __name__ == '__main__':
    main()
