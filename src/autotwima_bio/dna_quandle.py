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
        if not self.crossings:
            return domain_size

        # Dynamically map only active participating strand IDs to dense indices
        active_strands = sorted({idx for c in self.crossings for idx in c})
        strand_to_idx = {strand: i for i, strand in enumerate(active_strands)}
        num_strands = len(active_strands)

        valid_count = 0
        import itertools

        for coloring in itertools.product(range(domain_size), repeat=num_strands):
            is_valid = True
            for over, under_in, under_out in self.crossings:
                c_over = DNAStrandPhase(phase_id=coloring[strand_to_idx[over]])
                c_in = DNAStrandPhase(phase_id=coloring[strand_to_idx[under_in]])
                expected_out = c_over.right_action(c_in).phase_id
                if coloring[strand_to_idx[under_out]] != expected_out:
                    is_valid = False
                    break
            if is_valid:
                valid_count += 1

        return valid_count


@dataclass(frozen=True, slots=True)
class TopoisomerState:
    """State of a circular plasmid topoisomer during enzymatic relaxation.

    Attributes:
        knot: Topological knot structure.
        quandle_invariant: Count of valid quandle colorings.
        relative_mobility: Predicted electrophoretic mobility in agarose gel (0.0 to 1.0).

    """

    knot: CircularDNAKnot
    quandle_invariant: int
    relative_mobility: float


class TopoisomeraseII:
    """Enzymatic model of Type II Topoisomerase (DNA Gyrase / Topo II).

    Relaxes supercoiled circular DNA via ATP-dependent double-strand passage.
    """

    @staticmethod
    def resolve_crossing(knot: CircularDNAKnot, crossing_idx: int) -> CircularDNAKnot:
        """Simulates strand passage through a double-strand break.

        Resolves or removes the specified crossing from the circular DNA.
        """
        if not knot.crossings or crossing_idx >= len(knot.crossings):
            return knot

        new_crossings = tuple(
            c for i, c in enumerate(knot.crossings) if i != crossing_idx
        )
        # Any knot diagram with < 3 crossings collapses to Unknot by Reidemeister I & II
        if len(new_crossings) < 3:
            new_crossings = ()

        return CircularDNAKnot(
            name=f"{knot.name}_relaxed_step",
            crossings=new_crossings,
        )

    @classmethod
    def simulate_relaxation_pathway(
        cls, initial_knot: CircularDNAKnot
    ) -> list[TopoisomerState]:
        """Simulates full enzymatic relaxation from knotted DNA down to unknot.

        Tracks discrete quandle invariants and electrophoretic mobility bands.
        """
        trajectory: list[TopoisomerState] = []
        current = initial_knot

        while True:
            q_inv = current.count_quandle_colorings(domain_size=3)
            # Mobility model: compact knots (high Q) migrate faster
            mobility = 0.5 + 0.5 * min(1.0, (q_inv - 3) / 6.0)
            trajectory.append(
                TopoisomerState(
                    knot=current,
                    quandle_invariant=q_inv,
                    relative_mobility=round(mobility, 3),
                )
            )
            if not current.crossings or q_inv == 3:
                break
            # Topo II unties one crossing per catalytic cycle
            current = cls.resolve_crossing(current, crossing_idx=0)

        return trajectory
