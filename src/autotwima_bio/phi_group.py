"""Autotwima Phi-Magma and Phi-Group theoretical foundation.

Formalized in Lean 4 (Paper sections 2, 3, 4, 5, 7).
"""

from dataclasses import dataclass
from typing import Self


@dataclass(frozen=True, slots=True)
class PhiInfElement:
    """Discrete-infinite Phi-Group model M_inf = Z x Bool (Section 5.2).

    coord: genomic coordinate / translation offset in Z
    is_reverse: strand orientation (False = Forward / +1, True = Reverse / -1)
    """

    coord: int
    is_reverse: bool = False

    @property
    def sgn(self) -> int:
        """Orientation factor sgn(F) = 1, sgn(T) = -1."""
        return -1 if self.is_reverse else 1

    def mul(self, other: Self) -> "PhiInfElement":
        """Phi-multiplication (Equation 17).

        (x1, u1) * (x2, u2) = (x1 + sgn(u1)*x2, u2)
        """
        return PhiInfElement(
            coord=self.coord + self.sgn * other.coord,
            is_reverse=other.is_reverse,
        )

    def phi(self, target: Self) -> "PhiInfElement":
        """Phase twist automorphism Phi_a(c) (Equation 18).

        Phi_(x,u)(y, v) = (sgn(u)*y, v)
        """
        return PhiInfElement(
            coord=self.sgn * target.coord,
            is_reverse=target.is_reverse,
        )

    def inv(self) -> "PhiInfElement":
        """Inverse under quasigroup action (Equation 19).

        (x, u)^(-1) = (-sgn(u)*x, u)
        """
        return PhiInfElement(
            coord=-self.sgn * self.coord,
            is_reverse=self.is_reverse,
        )

    @classmethod
    def identity(cls) -> "PhiInfElement":
        """Neutral element e = (0, False) (Equation 20)."""
        return cls(coord=0, is_reverse=False)


def composite_invert(
    a: PhiInfElement, b: PhiInfElement, y: PhiInfElement
) -> PhiInfElement:
    """Theorem 4.3 (Equation 16): Solving (a * b) * x = y for x.

    x = Phi_a^(-1) ( b^(-1) * ( a^(-1) * y ) )
    Note: For M_inf, Phi_a^(-1) == Phi_a since sgn(u) in {+1, -1}.
    """
    step1 = a.inv().mul(y)
    step2 = b.inv().mul(step1)
    # Apply inverse twist Phi_a^(-1)
    return a.phi(step2)


def quandle_conjugate(a: PhiInfElement, b: PhiInfElement) -> PhiInfElement:
    """Theorem 7.3: Quandle of conjugation f |> g = f^(-1) o g o f.

    In Aut(M), mapped through Phi: Phi(a * b) = Phi(a) |> Phi(b).
    """
    return PhiInfElement(
        coord=a.phi(b).coord,
        is_reverse=b.is_reverse,
    )
