"""Comparison between classical iterative liftover and Autotwima Phi-group algebra."""

from dataclasses import dataclass

from autotwima_bio.inversion_algebra import GenomicCoordinatePath, GenomicSegment


@dataclass(frozen=True, slots=True)
class InversionBlock:
    """Classical genomic block for interval lookup (e.g., chain/BED format)."""

    chrom_start: int
    length: int
    is_inverted: bool


def classical_iterative_liftover(blocks: list[InversionBlock], local_offset: int) -> int:
    """Classical genomics algorithm: iterative boundary translation.

    Uses conditional branching on strand orientation.
    """
    current_pos = local_offset
    for block in blocks:
        if block.is_inverted:
            current_pos = block.chrom_start - current_pos
        else:
            current_pos = block.chrom_start + current_pos
    return current_pos


def phi_algebraic_liftover(path: GenomicCoordinatePath, local_offset: int) -> int:
    """Autotwima algebraic liftover: pure Phi-multiplication without branching."""
    query = GenomicSegment(name="Query_Variant", length=local_offset, is_reverse=False)
    full_path = GenomicCoordinatePath.from_list([*path.segments, query])
    return full_path.evaluate_phi().coord
