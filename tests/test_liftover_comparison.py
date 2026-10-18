"""Tests verifying mathematical equivalence of Phi-algebra and classical liftover."""

from autotwima_bio.inversion_algebra import GenomicCoordinatePath, GenomicSegment
from autotwima_bio.liftover_comparison import (
    InversionBlock,
    PrecomputedPhiLiftover,
    classical_iterative_liftover,
    phi_algebraic_liftover,
)


def test_equivalence_on_complex_nested_inversions() -> None:
    """Verify that branchless Phi-algebra exactly matches classical liftover."""
    segments = [
        GenomicSegment("Block_1", length=100_000, is_reverse=True),
        GenomicSegment("Block_2", length=25_000, is_reverse=True),
        GenomicSegment("Block_3", length=10_000, is_reverse=False),
        GenomicSegment("Block_4", length=3_000, is_reverse=True),
    ]
    path = GenomicCoordinatePath.from_list(segments)

    blocks = [
        InversionBlock(chrom_start=s.length, length=s.length, is_inverted=s.is_reverse)
        for s in segments
    ]

    for test_offset in [0, 50, 1500, 2999]:
        classical_result = classical_iterative_liftover(blocks, test_offset)
        phi_result = phi_algebraic_liftover(path, test_offset)
        assert classical_result == phi_result, (
            f"Mismatch at offset {test_offset}: "
            f"classical={classical_result}, phi={phi_result}"
        )


def test_precomputed_phi_liftover_matches_classical() -> None:
    """Verify O(1) PrecomputedPhiLiftover gives identical results."""
    segments = [
        GenomicSegment("Block_A", length=50_000, is_reverse=True),
        GenomicSegment("Block_B", length=12_000, is_reverse=True),
        GenomicSegment("Block_C", length=4_000, is_reverse=False),
    ]
    path = GenomicCoordinatePath.from_list(segments)
    blocks = [
        InversionBlock(chrom_start=s.length, length=s.length, is_inverted=s.is_reverse)
        for s in segments
    ]

    resolver = PrecomputedPhiLiftover(path)

    for offset in [0, 10, 500, 3999]:
        assert resolver.translate(offset) == classical_iterative_liftover(blocks, offset)
