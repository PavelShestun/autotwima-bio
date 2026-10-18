"""Topological Quandle invariants for circular DNA and topoisomerase action.

Based on Section 7 of the Autotwima paper (Theorems 7.2 and 7.3)
and Joyce Knot Quandle theory for supercoiled DNA.
"""

from dataclasses import dataclass
from typing import Self


@dataclass(frozen=True, slots=True)
class DNAStrandPhase:
    """Represents a topological strand state in the automorphism quandle Aut(M).

    Attributes:
        phase_id: Discrete identifier of the strand or supercoil domain.
        sign: Crossing orientation (+1 for right-handed, -1 for left-handed).

    """

    phase_id: int
    sign: int = 1

    def right_action(self, other: Self) -> "DNAStrandPhase":
        """Quandle binary operation: self |> other = f^(-1) o g o f (Eq. 24).

        In modular strand crossing algebra: (2 * self - other) mod N.
        """
        return DNAStrandPhase(
            phase_id=(2 * self.phase_id - other.phase_id) % 3,
            sign=self.sign * other.sign,
        )

    def left_action(self, other: Self) -> "DNAStrandPhase":
        """Inverse quandle operation: self <| other = f o g o f^(-1).

        Guarantees invertibility of crossings (Reidemeister II).
        """
        return self.right_action(other)


@dataclass(frozen=True, slots=True)
class CircularDNAKnot:
    """Model of knotted circular plasmid or viral DNA (e.g. phiX174).

    Represented as a sequence of crossings under topoisomerase activity.
    """

    name: str
    crossings: tuple[tuple[int, int, int], ...]  # (over_strand, under_in, under_out)

    def count_quandle_colorings(self, domain_size: int = 3) -> int:
        """Computes the topological Joyce quandle invariant.

        Counts the number of distinct valid quandle colorings of the knot.
        - Unknot (relaxed circle): returns exactly domain_size colorings (trivial).
        - Trefoil knot (3_1 supercoil): returns > domain_size colorings.
        """
        valid_count = 0
        num_strands = len(self.crossings)
        if num_strands == 0:
            return domain_size

        # Exhaustive search over all strand colorings in Z_domain
        import itertools

        for coloring in itertools.product(range(domain_size), repeat=num_strands):
            is_valid = True
            for over, under_in, under_out in self.crossings:
                c_over = DNAStrandPhase(phase_id=coloring[over])
                c_in = DNAStrandPhase(phase_id=coloring[under_in])
                expected_out = c_over.right_action(c_in).phase_id
                if coloring[under_out] != expected_out:
                    is_valid = False
                    break
            if is_valid:
                valid_count += 1

        return valid_count
