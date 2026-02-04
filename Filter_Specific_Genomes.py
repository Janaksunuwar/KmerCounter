# Import information from assembly_summary_genbank.txt and download specific genome FASTA files
import os
import subprocess
import logging
import pandas as pd
import shutil
from pathlib import Path
import gzip

# --------------------------- SETUP LOGGING ---------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# --------------------------- LOAD ASSEMBLY SUMMARY ---------------------------
def load_assembly_summary(file_path):
    """
    Load assembly_summary_genbank.txt into a DataFrame.
    Handles lines starting with '#'.
    """
    # Read file manually to handle the header line
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#assembly_accession"):
                header_line = line.lstrip("#").strip()
                break
        else:
            raise ValueError("No header line found starting with '#assembly_accession'")

    # Read the data using the correct header line
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        names=header_line.split("\t"),
        dtype=str,
        low_memory=False
    )

    # Strip any spaces in column names
    df.columns = df.columns.str.strip()

    # Verify required columns exist
    if not {"assembly_accession", "ftp_path"}.issubset(df.columns):
        raise KeyError(f"Expected columns not found. Found: {list(df.columns)}")

    df = df[["assembly_accession", "ftp_path"]].dropna()
    logging.info(f"Loaded {len(df)} assemblies from {file_path}")
    return df

# DOWNLOAD FUNCTION
def download_with_wget(url, out_path):
    """Download a file from URL using wget."""
    try:
        subprocess.run(["wget", "-q", "-O", str(out_path), url], check=True)
        logging.info(f"Downloaded {url}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"wget failed for {url}: {e}")
        return False

# MAIN DOWNLOAD HANDLER
def download_genomes(df, download_dir):
    """
    Download genome FASTA files from the FTP paths provided in the dataframe.
    Adds status information and saves a summary file.
    """
    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)

    records = []

    for _, row in df.iterrows():
        acc = row["assembly_accession"]
        ftp_path = row["ftp_path"]
        asm_name = ftp_path.split("/")[-1]

        file_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        local_fname = f"{acc}_{asm_name}_genomic.fna.gz"
        local_path = download_dir / local_fname

        # Skip existing files
        if local_path.exists() and local_path.stat().st_size > 0:
            logging.info(f"⏭️ Skipping existing file: {local_fname}")
            status = "Already exists"
        else:
            success = download_with_wget(file_url, local_path)
            status = "Downloaded" if success else "Failed"

        records.append({
            "assembly_accession": acc,
            "asm_name": asm_name,
            "ftp_url": file_url,
            "local_file": str(local_path),
            "status": status
        })

    summary_df = pd.DataFrame(records)
    summary_path = download_dir / "download_summary_enter.xlsx"
    summary_df.to_excel(summary_path, index=False)
    logging.info(f"📄 Download summary saved to {summary_path}")
    return summary_df

# UNZIP FUNCTION
def unzip_gz_files(download_dir):
    """Unzip all .gz files in the directory."""
    for gz_file in Path(download_dir).glob("*.gz"):
        unzipped = gz_file.with_suffix("")  # Removes .gz
        with gzip.open(gz_file, "rb") as f_in:
            with open(unzipped, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        logging.info(f"🗜️ Unzipped {gz_file.name} -> {unzipped.name}")

# MAIN SCRIPT
def main():
    assembly_summary_path = Path("INSERT PATH/assembly_summary_genbank.txt")
    download_dir = Path("INSERT PATH")

    # Load assembly summary file
    assembly_df = load_assembly_summary(assembly_summary_path)

    # FILTER SPECIFIC ACCESSIONS
    accessions = [
    "INSERT LIST,"
    ]

    assembly_df = assembly_df[assembly_df["assembly_accession"].isin(accessions)]

    logging.info(f"Filtered to {len(assembly_df)} selected assemblies")

    # ------------------ DOWNLOAD & UNZIP ------------------
    logging.info("🚀 Starting genome download...")
    download_genomes(assembly_df, download_dir)
    unzip_gz_files(download_dir)
    logging.info("All genomes downloaded and unzipped successfully!")

if __name__ == "__main__":
    main()
