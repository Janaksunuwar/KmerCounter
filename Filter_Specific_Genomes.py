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
                header_line = line.lstrip("#").strip()  # remove leading '#' and trailing '\n'
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

# --------------------------- DOWNLOAD FUNCTION ---------------------------
def download_with_wget(url, out_path):
    """Download a file from URL using wget."""
    try:
        subprocess.run(["wget", "-q", "-O", str(out_path), url], check=True)
        logging.info(f"Downloaded {url}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"wget failed for {url}: {e}")
        return False

# --------------------------- MAIN DOWNLOAD HANDLER ---------------------------
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

# --------------------------- UNZIP FUNCTION ---------------------------
def unzip_gz_files(download_dir):
    """Unzip all .gz files in the directory."""
    for gz_file in Path(download_dir).glob("*.gz"):
        unzipped = gz_file.with_suffix("")  # Removes .gz
        with gzip.open(gz_file, "rb") as f_in:
            with open(unzipped, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        logging.info(f"🗜️ Unzipped {gz_file.name} -> {unzipped.name}")

# --------------------------- MAIN SCRIPT ---------------------------
def main():
    assembly_summary_path = Path("/Users/gopal/Desktop/Research/assembly_summary_genbank.txt")
    download_dir = Path("/Users/gopal/Desktop/Research/pathogens/EnteroBac/GENOMES")

    # Load assembly summary file
    assembly_df = load_assembly_summary(assembly_summary_path)

    # ------------------ FILTER SPECIFIC ACCESSIONS ------------------
    accessions = [
    "GCA_002968455.1", "GCA_003073995.1", "GCA_003186415.1", "GCA_003288475.1",
    "GCA_020540925.1", "GCA_024599655.1", "GCA_032298955.1", "GCA_049600895.1",
    "GCA_049600905.1", "GCA_049605525.1", "GCA_051549645.1", "GCA_051549745.1",
    "GCA_052975065.1", "GCA_031379055.2", "GCA_043905735.1", "GCA_050957065.1",
    "GCA_050957085.1", "GCA_050957105.1", "GCA_050957685.1", "GCA_050957345.1",
    "GCA_003324425.1", "GCA_024706625.1", "GCA_024706725.1", "GCA_017328485.1",
    "GCA_004343625.1", "GCA_017328285.1", "GCA_038080455.1", "GCA_031378785.2",
    "GCA_038080135.1", "GCA_032336445.2", "GCA_031377915.2", "GCA_031377625.2",
    "GCA_031378175.2", "GCA_031378355.2", "GCA_030868435.2", "GCA_033200655.2",
    "GCA_031378305.2", "GCA_031378745.1", "GCA_031377665.2", "GCA_031378025.2",
    "GCA_031377845.2", "GCA_031378285.2", "GCA_031378085.2", "GCA_015685515.2",
    "GCA_017328425.1", "GCA_017328505.1", "GCA_017328205.1", "GCA_017328105.1",
    "GCA_017345915.1", "GCA_032182145.2", "GCA_033726955.1", "GCA_031378545.1",
    "GCA_031378465.2", "GCA_033726975.1", "GCA_031378915.2", "GCA_031377495.2",
    "GCA_018019835.1", "GCA_017328585.1", "GCA_017328385.1", "GCA_038080795.1",
    "GCA_036669255.1", "GCA_032331115.2", "GCA_036669195.1", "GCA_032331815.1",
    "GCA_030868415.2", "GCA_031835335.2", "GCA_017328265.1", "GCA_038080175.1",
    "GCA_031377565.2", "GCA_033726515.1", "GCA_031377725.2", "GCA_031377585.2",
    "GCA_017345855.1", "GCA_032337615.2", "GCA_032330635.2", "GCA_031378525.2",
    "GCA_017328125.1", "GCA_030868395.2", "GCA_031378985.2", "GCA_017328365.1",
    "GCA_031376755.2", "GCA_038510485.1", "GCA_033726575.1", "GCA_032329195.2",
    "GCA_032329995.2", "GCA_031776635.2", "GCA_032332155.2", "GCA_031378105.2",
    "GCA_032329475.2", "GCA_032335295.2", "GCA_032330795.2", "GCA_017328165.1",
    "GCA_031378935.2", "GCA_031378505.2", "GCA_031378235.1", "GCA_017345835.1",
    "GCA_031378695.2", "GCA_031378455.2", "GCA_033726715.1", "GCA_031378725.2",
    "GCA_031378615.2", "GCA_033726775.1", "GCA_038510615.1", "GCA_031377705.2",
    "GCA_031379095.2", "GCA_017328225.1", "GCA_031378855.2", "GCA_031377765.2",
    "GCA_017345875.1", "GCA_032330595.2", "GCA_031377865.2", "GCA_031378065.2",
    "GCA_030868455.2", "GCA_031378145.2", "GCA_031378045.2", "GCA_029686665.1",
    "GCA_031378665.2", "GCA_033200455.2", "GCA_031378825.2", "GCA_031377955.2",
    "GCA_003260505.1", "GCA_014840875.1", "GCA_014842795.1", "GCA_003000605.1",
    "GCA_029686345.1", "GCA_014842835.1", "GCA_036691285.1", "GCA_045490375.1",
    "GCA_044164185.1", "GCA_044164045.1", "GCA_044164205.1", "GCA_044164265.1",
    "GCA_044164085.1", "GCA_044164095.1", "GCA_045490365.1", "GCA_045490295.1",
    "GCA_044164135.1", "GCA_044164245.1", "GCA_044164035.1", "GCA_044164025.1",
    "GCA_015673375.1", "GCA_015673435.1", "GCA_015673495.1", "GCA_049102015.1",
    "GCA_049102445.1", "GCA_049383405.1", "GCA_049102585.1", "GCA_049102565.1",
    "GCA_049102545.1", "GCA_049383605.1", "GCA_017328605.1", "GCA_002246615.1",
    "GCA_020552745.2", "GCA_020552775.1", "GCA_029675465.1"
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