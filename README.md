# 🥾 STOMP 🧬
> **Stochastic Technique for Optimization of Motifs (STOMP)** — A command-line motif discovery tool that uses Gibbs sampling and simulated annealing to identify recurring patterns in DNA sequences.

---

![](assets/09182-ezgif.com-video-to-gif-converter.gif)

---
## 🧪 Methodology 

This tool identifies conserved sequence motifs within unaligned DNA sequences by combining **Gibbs Sampling** with **Simulated Annealing** to avoid local optima. 

Here is a step-by-step breakdown of how the algorithm operates:

*   **Initialization:** The tool reads the input sequences and randomly selects a starting $k$-mer (of user-defined size $k$) from each sequence to build an initial probability profile (Position Weight Matrix) using pseudocounts.
*   **Gibbs Sampling Iteration:** During each iteration, one motif is randomly removed from the active set. The profile is recalculated without it, and a new replacement motif is probabilistically drawn from that sequence based on the updated profile.
*   **Simulated Annealing Optimization:** The tool evaluates the newly selected motif set using either a **Shannon Entropy** or **Hamming Distance** scoring function. 
    *   Changes that *improve* (lower) the score are always accepted.
    *   Changes that *worsen* the score may still be accepted based on a probability equation ($e^{-\Delta / T}$) and a cooling "temperature" schedule. This allows the algorithm to escape local minimums early in the run.
*   **Multiple Restarts:** To ensure robust results, the algorithm performs full restarts from scratch (default is 3), independently searching different random starting states and saving the best overall motif alignment.
---

## ⚙️ Usage Instructions

1. Install the tool locally:
```bash
pip install git+https://github.com/jcapecci09/STOMP.git
```
2. Run STOMP on a FASTA file of your choice. You can use an example FASTA file in the `Data` folder:
```bash
stomp -i Data/MA0001.1.sites
```

*Note: There are additional flags you can use. They are listed below.*

## 🛠️ Command-Line Arguments

| Flag | Description | Default |
|---|---|---|
| `-i`, `--input` | Input FASTA file | Required |
| `-o`, `--output_directory` | Output directory | `Results` |
| `-s`, `--scoring_function` | Scoring method (`entropy` or `hamming`) | `entropy` |
| `-k`, `--kmer_size` | Size of motifs | `8` |
| `--iterations` | Number of iterations | `50 × number of sequences` |
| `--restart` | Number of full algorithm restarts | `3` |
| `-t`, `--temp` | Simulated annealing temperature | `10` |

---

## 📁 Output Files

Each run creates a timestamped results folder containing:

| File | Description |
|---|---|
| `motifs_k#.txt` | Best motif found for each input sequence |
| `motifs_k#.svg` | Sequence logo showing the discovered motif |
| `motifs_k#_pwm.csv` | Position weight matrix (PWM) for the discovered motif |
| `motifs_k#_scores.txt` | Run statistics including consensus, scores, iterations, restarts, and runtime |

---

**[📊 Benchmarks](benchmarks.md)**

---
## 👤 Author

I'm Jimmy Capecci, a Bioinformatics master's student at Loyola University Chicago. I was introduced to algorithms used in motif discovery in my advanced bioinformatics course, and I wanted to make a reusable tool to really understand the algorithms. STOMP was created as a way to apply these concepts in practice while building a tool that can identify motifs in DNA sequences.

---

## 🙏 Acknowledgements

This project is largely based on the concepts discussed in Chapter 3 of *Bioinformatics Algorithms: An Active-Learning Approach* by Phillip Compeau and Pavel Pevzner.
