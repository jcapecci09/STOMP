# 🥾 STOMP 🧬
> Stochastic Technique for Optimization of Motifs (STOMP). A motif discovery tool that uses Gibbs sampling and simulated annealing to identify recurring patterns in DNA sequence.

---
![](assets/09182-ezgif.com-video-to-gif-converter.gif)
---
## ⚙️ Usage Instructions


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

## 👤 Author
I'm Jimmy Capecci, a Bioinformatics master's student at Loyola University Chicago. I was introduced to algorithms used in motif discovery in my advanced bioinformatics course, and I wanted to make a reusable tool to really understand the algorithms. STOMP was created as a way to apply these concepts in practice while building a tool that can identify motifs in DNA sequences.
---
## 🙏 Acknowledgements
This project is largely based on the concepts discussed in Chapter 3 of *Bioinformatics Algorithms: An Active-Learning Approach* by Phillip Compeau and Pavel Pevzner.
