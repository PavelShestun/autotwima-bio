"""Tests for bioinformatics application."""

from autotwima_bio.stranded_dna import dna_quandle_hash, fold_dna_strands


def test_strand_fold_empty() -> None:
    result = fold_dna_strands([])
    assert result.coord == 0
    assert not result.is_reverse


def test_strand_fold_inversion_preserves_recovery() -> None:
    intervals = [(100, False), (50, True), (30, False)]
    result = fold_dna_strands(intervals)
    assert result.coord == (100 + 50) + (-30)


def test_quandle_hash_deterministic() -> None:
    seq = "ATGCGATCGATCGATCGATC"
    h1 = dna_quandle_hash(seq)
    h2 = dna_quandle_hash(seq)
    assert h1 == h2
