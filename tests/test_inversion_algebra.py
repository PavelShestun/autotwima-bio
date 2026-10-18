"""Tests for genomic segments and their Phi-group representation."""

import pytest

from autotwima_bio.inversion_algebra import GenomicCoordinatePath, GenomicSegment


def test_genomic_segment_forward_strand_mapping() -> None:
    """Verify forward strand segment maps to sgn = +1 in M_inf."""
    seg = GenomicSegment(name="Exon_1", length=1200, is_reverse=False)
    phi_elem = seg.to_phi()

    assert phi_elem.coord == 1200
    assert not phi_elem.is_reverse
    assert phi_elem.sgn == 1


def test_genomic_segment_reverse_strand_mapping() -> None:
    """Verify inverted/reverse segment maps to sgn = -1 in M_inf."""
    seg = GenomicSegment(name="Inversion_Breakpoint_A", length=4500, is_reverse=True)
    phi_elem = seg.to_phi()

    assert phi_elem.coord == 4500
    assert phi_elem.is_reverse
    assert phi_elem.sgn == -1


def test_genomic_segment_rejects_negative_length() -> None:
    """Verify domain constraint: segment length must be non-negative."""
    with pytest.raises(ValueError, match="must be non-negative"):
        GenomicSegment(name="Invalid_Region", length=-10, is_reverse=False)


def test_associativity_defect_forward_strand_is_zero() -> None:
    """Verify Theorem 5.4: when segment A is on forward strand (eps1 = 1), defect is 0."""
    path = GenomicCoordinatePath.from_list(
        [
            GenomicSegment(name="Region_A", length=1000, is_reverse=False),
            GenomicSegment(name="Region_B", length=500, is_reverse=True),
            GenomicSegment(name="Region_C", length=250, is_reverse=False),
        ]
    )

    defect = path.associativity_defect(0, 1, 2)
    assert defect == 0


def test_associativity_defect_inverted_strand_matches_lean_theorem_5_4() -> None:
    """Verify Theorem 5.4: eps2 * (1 - eps1) * x3 under inversion."""
    # eps1 = -1 (reverse), eps2 = -1 (reverse), x3 = 300
    # Expected defect: (-1) * (1 - (-1)) * 300 = -1 * 2 * 300 = -600
    path = GenomicCoordinatePath.from_list(
        [
            GenomicSegment(name="Inv_1", length=5000, is_reverse=True),
            GenomicSegment(name="Inv_2", length=2000, is_reverse=True),
            GenomicSegment(name="Marker_3", length=300, is_reverse=False),
        ]
    )

    defect = path.associativity_defect(0, 1, 2)
    assert defect == -600


def test_associativity_defect_mixed_strand_orientation() -> None:
    """Verify Theorem 5.4 with eps1 = -1 and eps2 = +1."""
    # eps1 = -1 (reverse), eps2 = +1 (forward), x3 = 450
    # Expected defect: (+1) * (1 - (-1)) * 450 = 2 * 450 = 900
    path = GenomicCoordinatePath.from_list(
        [
            GenomicSegment(name="Inv_1", length=1500, is_reverse=True),
            GenomicSegment(name="Norm_2", length=800, is_reverse=False),
            GenomicSegment(name="Marker_3", length=450, is_reverse=False),
        ]
    )

    defect = path.associativity_defect(0, 1, 2)
    assert defect == 900
