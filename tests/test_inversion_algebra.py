"""Tests for genomic segments and their Phi-group representation."""

import pytest

from autotwima_bio.inversion_algebra import GenomicSegment


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
