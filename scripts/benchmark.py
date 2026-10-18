"""Benchmark: Phi-group folding & quandle hashing on real DNA data."""

import sys
import time
from pathlib import Path

from autotwima_bio.stranded_dna import dna_quandle_hash, fold_dna_strands

DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "phix174.fasta"


def load_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main() -> None:
    if not DATA_FILE.exists():
        print("Data file not found. Running downloader first...")
        from data.download_data import download_benchmark_data

        download_benchmark_data()

    seq = load_fasta(DATA_FILE)
    length = len(seq)
    print(f"=== Autotwima-Bio Benchmark on phiX174 ({length} bp) ===")

    # 1. Quandle Hash Benchmark
    iterations = 200
    t0 = time.perf_counter()
    for _ in range(iterations):
        _ = dna_quandle_hash(seq)
    t1 = time.perf_counter()

    elapsed = t1 - t0
    m_bp_per_sec = (length * iterations) / (elapsed * 1e6)
    print(f"DNA Quandle Hashing: {iterations} iterations in {elapsed:.3f} s")
    print(f"Throughput: {m_bp_per_sec:.2f} Mbp/s")

    # 2. Stranded Coordinate Folding Benchmark
    intervals = [(i % 150, i % 3 == 0) for i in range(10_000)]
    t0 = time.perf_counter()
    for _ in range(200):
        _ = fold_dna_strands(intervals)
    t1 = time.perf_counter()

    fold_throughput = (10_000 * 200) / (t1 - t0)
    print(f"Stranded Coordinate Folding: {fold_throughput:,.0f} intervals/s")
    print("=== Benchmark Finished Successfully ===")


if __name__ == "__main__":
    sys.exit(main())
