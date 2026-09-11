from Bio import Align, SeqIO
from Bio.Seq import Seq
import csv


def read_fasta(filename):
    try:
        return str(SeqIO.read(filename, "fasta").seq)
    except (FileNotFoundError, ValueError):
        print(f"Error reading {filename}.")
        exit()


def validate_sequence(sequence):
    return all(base in "ATGC" for base in sequence)


ref_file = input("Enter reference FASTA filename: ")
mut_file = input("Enter mutant FASTA filename: ")

reference = read_fasta(ref_file)
mutant = read_fasta(mut_file)

if not validate_sequence(reference) or not validate_sequence(mutant):
    print("Error: FASTA contains invalid DNA bases.")
    exit()


aligner = Align.PairwiseAligner()
alignment = aligner.align(reference, mutant)[0]

aligned_ref = alignment[0]
aligned_mut = alignment[1]

variants = []

for i, (ref, mut) in enumerate(zip(aligned_ref, aligned_mut), 1):
    if ref != mut:
        if ref == "-":
            variant_type = "Insertion"
        elif mut == "-":
            variant_type = "Deletion"
        else:
            variant_type = "SNP"

        variants.append((i, ref, mut, variant_type))


ref_protein = Seq(reference).translate()
mut_protein = Seq(mutant).translate()

if "*" in mut_protein and "*" not in ref_protein:
    effect = "Nonsense"
elif ref_protein != mut_protein:
    effect = "Missense"
else:
    effect = "Synonymous"


mutation_count = len(variants)
mutation_rate = (mutation_count / len(reference)) * 100


with open("variant_report.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Position", "Reference", "Mutant", "Type"])

    writer.writerows(variants)


print("\n" + "=" * 40)
print("       GENOME VARIANT ANALYZER")
print("=" * 40)

print(f"\nReference : {ref_file}")
print(f"Mutant    : {mut_file}")

print(f"\nSequence length : {len(reference)} bp")
print(f"Variants found  : {mutation_count}")
print(f"Mutation rate   : {mutation_rate:.2f}%")

print("\n" + "-" * 40)
print("VARIANTS")
print("-" * 40)

for position, ref, mut, variant_type in variants:
    print(f"Position {position}: {ref} → {mut} ({variant_type})")

print("\n" + "-" * 40)
print("PROTEIN")
print("-" * 40)

print(f"Reference : {ref_protein}")
print(f"Mutant    : {mut_protein}")
print(f"Effect    : {effect}")

print("\n" + "=" * 40)
print("Analysis complete!")
print("CSV report saved: variant_report.csv")
print("=" * 40)