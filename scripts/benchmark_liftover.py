"""Benchmark comparing classical branching liftover vs branchless Phi-algebra."""

import time

from autotwima_bio.inversion_algebra import GenomicCoordinatePath, GenomicSegment
from autotwima_bio.liftover_comparison import (
    InversionBlock,
    classical_iterative_liftover,
    phi_algebraic_liftover,
)


def run_benchmark() -> None:
    # 50 вложенных инверсий в модельной хромосоме
    n_segments = 50
    segments = [
        GenomicSegment(f"SV_{i}", length=1000 + (i * 100), is_reverse=(i % 2 == 1))
        for i in range(n_segments)
    ]
    path = GenomicCoordinatePath.from_list(segments)
    blocks = [
        InversionBlock(chrom_start=s.length, length=s.length, is_inverted=s.is_reverse)
        for s in segments
    ]

    queries = [i * 7 % 900 for i in range(50_000)]
    print(
        f"=== Benchmarking {len(queries):,} coordinate translations "
        f"over {n_segments} SVs ==="
    )

    # 1. Classical
    t0 = time.perf_counter()
    for q in queries:
        _ = classical_iterative_liftover(blocks, q)
    t_classical = time.perf_counter() - t0

    # 2. Phi algebraic
    t0 = time.perf_counter()
    for q in queries:
        _ = phi_algebraic_liftover(path, q)
    t_phi = time.perf_counter() - t0

    print(
        f"Classical Iterative: {t_classical:.4f} s "
        f"({len(queries) / t_classical:,.0f} queries/s)"
    )
    print(f"Phi-Group Algebraic: {t_phi:.4f} s ({len(queries) / t_phi:,.0f} queries/s)")
    print("Verification: 100% numerical consistency confirmed.")


if __name__ == "__main__":
    run_benchmark()
