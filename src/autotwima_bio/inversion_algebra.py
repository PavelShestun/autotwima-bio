"""Inversion-aware genomic coordinate components based on Autotwima Phi-groups."""

from dataclasses import dataclass

from autotwima_bio.phi_group import PhiInfElement, composite_invert


@dataclass(frozen=True, slots=True)
class GenomicSegment:
    """A genomic interval with length and strand orientation.

    Attributes:
        name: Chromosome region, gene or structural variant identifier.
        length: Non-negative segment length in base pairs.
        is_reverse: False for forward strand (+1), True for reverse strand (-1).

    """

    name: str
    length: int
    is_reverse: bool = False

    def __post_init__(self) -> None:
        if self.length < 0:
            msg = (
                f"Segment length must be non-negative, got {self.length} "
                f"for '{self.name}'"
            )
            raise ValueError(msg)

    def to_phi(self) -> PhiInfElement:
        """Converts genomic segment to M_inf element (length in Z, strand in Bool)."""
        return PhiInfElement(coord=self.length, is_reverse=self.is_reverse)


@dataclass(frozen=True, slots=True)
class GenomicCoordinatePath:
    """An ordered path of genomic segments representing chromosomal structure.

    Evaluates coordinate transitions through non-associative Phi-algebra.
    """

    segments: tuple[GenomicSegment, ...]

    @classmethod
    def from_list(cls, segments: list[GenomicSegment]) -> "GenomicCoordinatePath":
        """Constructs a path from an ordered list of genomic segments."""
        return cls(tuple(segments))

    def evaluate_phi(self) -> PhiInfElement:
        """Evaluates cumulative coordinates using left-to-right Phi-multiplication."""
        if not self.segments:
            return PhiInfElement.identity()
        acc = self.segments[0].to_phi()
        for seg in self.segments[1:]:
            acc = acc.mul(seg.to_phi())
        return acc

    def solve_relative_coordinate(
        self, seg_a: GenomicSegment, seg_b: GenomicSegment, target_coord: int
    ) -> int:
        """Finds the local coordinate x inside nested segments using Theorem 4.3.

        Given (A * B) * X = Y in global coordinates, recovers local feature coordinate X:
            X = Phi_a^(-1) ( b^(-1) * ( a^(-1) * Y ) )

        Args:
            seg_a: First structural segment (outer environment).
            seg_b: Second structural segment (inner nested environment).
            target_coord: Observed global coordinate Y.

        Returns:
            Recovered exact local coordinate x.

        """
        a = seg_a.to_phi()
        b = seg_b.to_phi()
        y = PhiInfElement(coord=target_coord, is_reverse=False)

        x_elem = composite_invert(a, b, y)
        return x_elem.coord

    def associativity_defect(self, a_idx: int, b_idx: int, c_idx: int) -> int:
        """Computes the associativity defect between three consecutive segments.

        Verifies Theorem 5.4 from the paper:
            Defect(a, b, c) = ((a * b) * c)_x - (a * (b * c))_x = eps2 * (1 - eps1) * x3.
        """
        a = self.segments[a_idx].to_phi()
        b = self.segments[b_idx].to_phi()
        c = self.segments[c_idx].to_phi()

        left_bracket = a.mul(b).mul(c).coord
        right_bracket = a.mul(b.mul(c)).coord
        defect_measured = left_bracket - right_bracket

        eps1 = a.sgn
        eps2 = b.sgn
        x3 = c.coord
        defect_theoretical = eps2 * (1 - eps1) * x3

        if defect_measured != defect_theoretical:
            msg = (
                f"Associativity defect mismatch! Measured {defect_measured}, "
                f"expected theoretical {defect_theoretical}"
            )
            raise AssertionError(msg)

        return defect_measured
