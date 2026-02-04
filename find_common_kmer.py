import logging
from pathlib import Path
from collections import Counter

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")

K = 11
INPUT_DIR = Path("/Users/gopal/Desktop/Research/pathogens/EnteroBac/GENOMES")
OUTPUT_DIR = Path("/Users/gopal/Desktop/Research/pathogens/EnteroBac/OUTPUT")
OUTPUT_FILE = OUTPUT_DIR / "common_kmers_enterobac.csv"

def read_fasta_sequence(fasta_file):
    """Read FASTA and return concatenated sequence."""
    seq = []
    with open(fasta_file) as f:
        for line in f:
            seq.append(line.strip().upper())
    return "".join(seq)


def extract_kmers(sequence, k):
    """Return Counter of kmers (excluding N)."""
    return Counter(
        sequence[i:i+k]
        for i in range(len(sequence) - k + 1)
        if "N" not in sequence[i:i+k]
    )


def main():
    fasta_files = sorted(INPUT_DIR.glob("*.fna"))

    if not fasta_files:
        logging.error("No FASTA files found")
        return

    common_kmers = None
    kmer_frequencies = {}

    for idx, fasta in enumerate(fasta_files, start=1):
        genome_name = fasta.stem
        logging.info(f"Processing {genome_name}")

        sequence = read_fasta_sequence(fasta)
        kmer_counts = extract_kmers(sequence, K)

        if idx == 1:
            # Initialize with first genome
            common_kmers = set(kmer_counts.keys())
            kmer_frequencies = dict(kmer_counts)
        else:
            genome_kmers = set(kmer_counts.keys())

            # Intersect
            common_kmers &= genome_kmers

            # Update frequencies ONLY for surviving kmers
            kmer_frequencies = {
                k: kmer_frequencies[k] + kmer_counts[k]
                for k in common_kmers
            }

        logging.info(
            f"{genome_name}: common kmers now {len(common_kmers)}"
        )

        if not common_kmers:
            logging.warning("No common kmers remain — stopping early")
            break

    # Write output
    with open(OUTPUT_FILE, "w") as out:
        out.write("kmer\tfrequency\n")
        for kmer, freq in sorted(kmer_frequencies.items()):
            out.write(f"{kmer}\t{freq}\n")

    logging.info(f"Saved common k=11 kmers to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()