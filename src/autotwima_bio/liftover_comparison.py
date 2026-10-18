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

    Iterates through all blocks in O(N), updating cumulative coordinate
    and strand orientation.
    """
    if not blocks:
        return local_offset

    current_coord = blocks[0].chrom_start
    current_inverted = blocks[0].is_inverted

    for block in blocks[1:]:
        sgn = -1 if current_inverted else 1
        current_coord += sgn * block.chrom_start
        current_inverted = block.is_inverted

    final_sgn = -1 if current_inverted else 1
    return current_coord + final_sgn * local_offset


def phi_algebraic_liftover(path: GenomicCoordinatePath, local_offset: int) -> int:
    """Autotwima algebraic liftover: pure Phi-multiplication without branching."""
    query = GenomicSegment(name="Query_Variant", length=local_offset, is_reverse=False)
    full_path = GenomicCoordinatePath.from_list([*path.segments, query])
    return full_path.evaluate_phi().coord


class PrecomputedPhiLiftover:
    """O(1) coordinate resolver via Phi-group state collapse.

    Precomputes cumulative path element P = A_1 * A_2 * ... * A_N once in O(N).
    Translates any query coordinate in strict O(1) arithmetic without branching.
    """

    __slots__ = ("_cached_state",)

    def __init__(self, path: GenomicCoordinatePath) -> None:
        self._cached_state = path.evaluate_phi()

    def translate(self, local_offset: int) -> int:
        """O(1) branchless coordinate translation: P * (local_offset, False)."""
        return self._cached_state.coord + self._cached_state.sgn * local_offset
