"""Inversion-aware genomic coordinate components based on Autotwima Phi-groups."""

from dataclasses import dataclass

from autotwima_bio.phi_group import PhiInfElement


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
