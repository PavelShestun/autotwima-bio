"""Tests for Lean-verified Autotwima Phi-group theorems."""

from autotwima_bio.phi_group import (
    PhiInfElement,
    composite_invert,
)


def test_strict_nonassociativity() -> None:
    """Verify non-associativity in M_inf (Section 5.2)."""
    a = PhiInfElement(coord=0, is_reverse=True)
    b = PhiInfElement(coord=0, is_reverse=False)
    c = PhiInfElement(coord=1, is_reverse=False)

    left_assoc = a.mul(b).mul(c)
    right_assoc = a.mul(b.mul(c))

    assert left_assoc.coord == 1
    assert right_assoc.coord == -1
    assert left_assoc != right_assoc


def test_four_term_coherence() -> None:
    """Theorem 2.2: ((a * b) * c) * d = a * (b * (Phi_a(c) * Phi_b(Phi_a(d))))."""
    a = PhiInfElement(coord=5, is_reverse=True)
    b = PhiInfElement(coord=3, is_reverse=False)
    c = PhiInfElement(coord=7, is_reverse=True)
    d = PhiInfElement(coord=2, is_reverse=False)

    lhs = a.mul(b).mul(c).mul(d)

    # Canonical RHS reduction
    phi_a_c = a.phi(c)
    phi_b_phi_a_d = b.phi(a.phi(d))
    inner = phi_a_c.mul(phi_b_phi_a_d)
    rhs = a.mul(b.mul(inner))

    assert lhs == rhs


def test_left_cancellation_and_action() -> None:
    """Theorem 4.1 & Equation 10: a^(-1) * (a * x) = x."""
    a = PhiInfElement(coord=12, is_reverse=True)
    x = PhiInfElement(coord=42, is_reverse=False)

    assert a.inv().mul(a.mul(x)) == x
    assert a.mul(a.inv().mul(x)) == x


def test_composite_inversion_theorem() -> None:
    """Theorem 4.3: Exact recovery of x from (a * b) * x = y."""
    a = PhiInfElement(coord=10, is_reverse=True)
    b = PhiInfElement(coord=-4, is_reverse=True)
    x = PhiInfElement(coord=100, is_reverse=False)

    y = a.mul(b).mul(x)
    recovered_x = composite_invert(a, b, y)

    assert recovered_x == x


def test_right_cancellation_fails() -> None:
    """Theorem 4.4: Right cancellation a * c = b * c does NOT imply a = b."""
    # Counterexample: a and b have different coords, but sgn cancels out
    a = PhiInfElement(coord=10, is_reverse=True)
    b = PhiInfElement(coord=10, is_reverse=False)
    c = PhiInfElement(coord=0, is_reverse=False)

    assert a.mul(c) == b.mul(c)
    assert a != b
