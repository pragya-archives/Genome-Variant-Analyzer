# Genome Variant Analyzer 🧬

A small Python project I made to compare two DNA sequences and find the differences between them.

The program takes a reference FASTA file and a mutant FASTA file, aligns the sequences, detects variants and checks whether the changes affect the resulting protein.

## What it can do

* Read DNA sequences from FASTA files
* Check whether the sequences contain valid DNA bases
* Align the reference and mutant sequences
* Find SNPs, insertions and deletions
* Translate DNA sequences into proteins
* Check the effect of a mutation on the protein
* Save the detected variants in a CSV file

## Tools used

* Python
* Biopython
* VS Code

## Example

For the example sequences included in this repository:

**Reference:** `ATGAAACCCGGGTAA`
**Mutant:** `ATGAAATCCGGGTAA`

The program detects:

`Position 7: C → T (SNP)`

The protein changes from:

`MKPG* → MKSG*`

So the mutation is classified as **Missense**.

## How to run

Make sure Biopython is installed, then run:

```bash
python main.py
```

When prompted, enter:

```text
reference.fasta
mutant.fasta
```

The detected variants are also saved in `variant_report.csv`.
