# Count the kmers of a given genome and populate the frequency
import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

wd = '/Users/gopal/Desktop/Research/KlebPneum'
os.chdir(wd)

def convert_fasta_to_single_line(input_file, output_file):
    """Convert a multi-line FASTA file to a single-line format."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header, sequence = None, []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    outfile.write(header + '\n')
                    outfile.write(''.join(sequence) + '\n')
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            outfile.write(header + '\n')
            outfile.write(''.join(sequence) + '\n')

def read_genome_length(fasta_file):
    """Read the genome sequence and return its length."""
    length = 0
    try:
        with open(fasta_file, 'r') as infile:
            for line in infile:
                if not line.startswith(">"):  # Skip header lines
                    length += len(line.strip())
        print(f"Genome length: {length}")
    except Exception as e:
        print(f"Error reading genome length: {e}")
    return length

def count_kmers(genome_length, kmer_size):
    """Count the number of k-mers possible in the genome."""
    if genome_length < kmer_size:
        return 0
    no_of_kmers = (genome_length - kmer_size) + 1
    print(f"Number of kmers for size {kmer_size}: {no_of_kmers}")
    return no_of_kmers

def read_kmer_counts(file_path):
    """Read k-mer counts from a Jellyfish dump file."""
    kmers, counts = [], []
    if not os.path.isfile(file_path):
        print(f"File not found: {file_path}")
        return kmers, counts
    with open(file_path, 'r') as file:
        current_count = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                try:
                    current_count = int(line[1:])  # Extract count after '>'
                except ValueError:
                    print(f"Skipping malformed header line: {line}")
                    current_count = None
            elif line:
                if current_count is not None:
                    kmers.append(line)
                    counts.append(current_count)
                else:
                    print(f"Skipping k-mer line '{line}' without a preceding count.")
    return kmers, counts

def create_dataframe(kmers, counts):
    """Create a DataFrame from k-mer lists."""
    return pd.DataFrame({'kmer': kmers, 'count': counts})

def save_to_csv(df, filename):
    """Save DataFrame to CSV file."""
    df.to_csv(filename, index=False)
    print(f"DataFrame saved to {filename}")

def save_to_txt(df, filename):
    """Save k-mer frequencies to a text file."""
    with open(filename, 'w') as file:
        for _, row in df.iterrows():
            file.write(f"{row['kmer']}\t{row['count']}\n")
    print(f"K-mer frequencies saved to {filename}")

def plot_top_kmers(df, n=100):
    """Plot the top n k-mers by count."""
    if df.empty:
        print("DataFrame is empty. No data to plot.")
        return
    df_sorted = df.sort_values(by='count', ascending=False)
    top_kmers = df_sorted.head(n)
    plt.figure(figsize=(10, 6))
    plt.bar(top_kmers['kmer'], top_kmers['count'], color='skyblue')
    plt.xlabel('K-mer')
    plt.ylabel('Count')
    plt.title(f'Top {n} K-mers by Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def main():
    input_file = 'kpkmergenome.fasta'
    output_fasta_file = 'output.fasta'
    convert_fasta_to_single_line(input_file, output_fasta_file)
    # Read the length of the genome
    genome_length = read_genome_length(output_fasta_file)
    # Specify kmer lengths
    kmer_lengths = [10, 11, 12, 13, 14]
    # To store k-mer counts for each size
    kmer_counts_data = []  
    for kmer_length in kmer_lengths:
        no_of_kmers = count_kmers(genome_length, kmer_length)
        kmer_counts_data.append({'kmer_size': kmer_length, 'number_of_kmers': no_of_kmers})
        mer_counts_file = f'mer_counts_k{kmer_length}.jf'
        mer_dump_file = f'mer_counts_dumps_k{kmer_length}.fa'
        frequencies_file_csv = f'kmer_frequencies_k{kmer_length}.csv'
        frequencies_file_txt = f'kmer_frequencies_k{kmer_length}.txt'
        # Counting k-mers with Jellyfish using subprocess
        subprocess.run(f'jellyfish count -m {kmer_length} -s 100M -t 10 -C {output_fasta_file} -o {mer_counts_file}', shell=True)
        subprocess.run(f'jellyfish dump {mer_counts_file} > {mer_dump_file}', shell=True)
        kmers, counts = read_kmer_counts(mer_dump_file)
        df = create_dataframe(kmers, counts)
        print(f"DataFrame head for k={kmer_length}:")
        print(df.head())
        print(df.groupby("count").size())
        if not df.empty:
            save_to_csv(df, frequencies_file_csv)
            save_to_txt(df, frequencies_file_txt)
            print("Basic statistics:")
            print(df.describe())
            plot_top_kmers(df)
    # Creates a dataframefor kmer_counts_data
    kmer_counts_df = pd.DataFrame(kmer_counts_data)
    print("K-mer counts based on genome length:")
    print("Hello")
    print(kmer_counts_df)

if __name__ == "__main__":
    main()
