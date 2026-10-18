"""Bioinformatics application: Stranded DNA coordinate folding & Knot hashing."""

from autotwima_bio.phi_group import PhiInfElement


def fold_dna_strands(intervals: list[tuple[int, bool]]) -> PhiInfElement:
    """Folds a sequence of stranded DNA reads/intervals into M_inf.

    Non-associative folding models directional coordinate drift across inversions.
    Each element: (length, is_reverse_strand).
    """
    if not intervals:
        return PhiInfElement.identity()

    acc = PhiInfElement(coord=intervals[0][0], is_reverse=intervals[0][1])
    for length, is_rev in intervals[1:]:
        elem = PhiInfElement(coord=length, is_reverse=is_rev)
        acc = acc.mul(elem)
    return acc


def dna_quandle_hash(sequence: str) -> int:
    """Generates a quandle-based topological knot hash from a DNA string.

    Captures hairpins, palindromes, and strand-symmetry via non-associative phase shifts.
    """
    valid_bases = {"A", "C", "G", "T", "N"}
    clean_seq = [b for b in sequence.upper() if b in valid_bases]
    if not clean_seq:
        return 0

    # Encoding: purines (A, G) = Forward, pyrimidines (C, T) = Reverse
    state = PhiInfElement.identity()
    for base in clean_seq:
        val = ord(base)
        is_pyr = base in ("C", "T")
        elem = PhiInfElement(coord=val, is_reverse=is_pyr)
        state = state.mul(elem)

    return state.coord
