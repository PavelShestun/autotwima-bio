"""Download reference genome (phiX174) for benchmarking."""

import urllib.request
from pathlib import Path

NCBI_PHIX174_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
OUTPUT_PATH = Path(__file__).resolve().parent / "phix174.fasta"


def download_benchmark_data() -> Path:
    if OUTPUT_PATH.exists():
        print(f"[OK] File already exists: {OUTPUT_PATH}")
        return OUTPUT_PATH

    print(f"Downloading phiX174 reference genome from NCBI to {OUTPUT_PATH}...")
    req = urllib.request.Request(
        NCBI_PHIX174_URL,
        headers={"User-Agent": "Autotwima-Bio/1.0"},
    )
    with (
        urllib.request.urlopen(req) as response,
        open(OUTPUT_PATH, "wb") as out_file,
    ):
        out_file.write(response.read())

    print(f"[DONE] Saved {OUTPUT_PATH.stat().st_size} bytes.")
    return OUTPUT_PATH


if __name__ == "__main__":
    download_benchmark_data()
