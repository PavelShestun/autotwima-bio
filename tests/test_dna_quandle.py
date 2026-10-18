"""Tests for DNA quandle topological invariants and Reidemeister moves."""

from autotwima_bio.dna_quandle import CircularDNAKnot, DNAStrandPhase


def test_reidemeister_move_i_idempotence() -> None:
    """Verify Reidemeister I: strand twist untwisting corresponds to a |> a = a."""
    strand = DNAStrandPhase(phase_id=1, sign=1)
    assert strand.right_action(strand).phase_id == strand.phase_id


def test_reidemeister_move_ii_invertibility() -> None:
    """Verify Reidemeister II: passing one strand over another is reversible."""
    a = DNAStrandPhase(phase_id=0, sign=1)
    b = DNAStrandPhase(phase_id=2, sign=-1)

    # (a |> b) <| a == b
    crossed = a.right_action(b)
    restored = a.left_action(crossed)
    assert restored.phase_id == b.phase_id


def test_reidemeister_move_iii_self_distributivity() -> None:
    """Verify Reidemeister III: Theorem 7.2 self-distributivity.

    a |> (b |> c) = (a |> b) |> (a |> c).
    """
    a = DNAStrandPhase(phase_id=1)
    b = DNAStrandPhase(phase_id=0)
    c = DNAStrandPhase(phase_id=2)

    lhs = a.right_action(b.right_action(c))
    rhs = a.right_action(b).right_action(a.right_action(c))

    assert lhs.phase_id == rhs.phase_id


def test_topoisomerase_knot_distinction() -> None:
    """Verify quandle coloring distinguishes relaxed circular DNA from Trefoil knot."""
    # Relaxed circle (Unknot 0_1): 1 strand, no crossings
    unknot = CircularDNAKnot(name="Relaxed_Plasmid", crossings=())
    unknot_colorings = unknot.count_quandle_colorings(domain_size=3)

    # Trefoil knot (3_1): standard 3-crossing alternating supercoil
    # Crossings: (0 over 1->2), (1 over 2->0), (2 over 0->1)
    trefoil = CircularDNAKnot(
        name="Supercoiled_Trefoil_Plasmid",
        crossings=((0, 1, 2), (1, 2, 0), (2, 0, 1)),
    )
    trefoil_colorings = trefoil.count_quandle_colorings(domain_size=3)

    # Unknot has 3 trivial colorings, Trefoil has 9 valid colorings
    assert unknot_colorings == 3
    assert trefoil_colorings == 9
    assert trefoil_colorings != unknot_colorings
