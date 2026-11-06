"""Locate exon-exon junctions by aligning an mRNA to its genomic sequence.

PerlPrimer used NCBI's Spidey for this; Spidey is no longer distributed.
This module performs a simple exact-match alignment: the mRNA is split into
the longest segments found, in order, in the genomic sequence, and each
junction is shifted (within the ambiguity of repeated bases) to the
canonical GT...AG intron boundary where possible.
"""

from __future__ import annotations

from dataclasses import dataclass

from .sequence import clean_sequence

MIN_EXON = 12


class AlignmentError(ValueError):
    """Raised when the mRNA cannot be aligned to the genomic sequence."""


@dataclass(frozen=True)
class Exon:
    """An exon in mRNA and genomic coordinates (0-based, half-open)."""

    mrna_start: int
    mrna_end: int
    genomic_start: int
    genomic_end: int

    @property
    def length(self) -> int:
        return self.mrna_end - self.mrna_start


@dataclass(frozen=True)
class ExonMap:
    exons: list[Exon]
    unaligned_tail: int  # bases at the 3' end of the mRNA not in the genome (e.g. poly-A)

    @property
    def boundaries(self) -> list[int]:
        """mRNA positions at which each exon after the first begins."""
        return [exon.mrna_start for exon in self.exons[1:]]


def map_exons(mrna: str, genomic: str, *, max_tail: int = 50) -> ExonMap:
    """Align an mRNA to genomic DNA and return its exon structure.

    Args:
        mrna: Spliced transcript sequence (mRNA or cDNA).
        genomic: Genomic sequence containing the whole gene.
        max_tail: Unaligned 3' bases tolerated (e.g. a poly-A tail).

    Raises:
        AlignmentError: if the mRNA does not align as exact exon segments.
    """
    mrna = clean_sequence(mrna).upper()
    genomic = clean_sequence(genomic).upper()
    if not mrna or not genomic:
        raise AlignmentError("Both an mRNA and a genomic sequence are required")

    exons: list[Exon] = []
    m_pos = g_pos = 0
    while m_pos < len(mrna):
        length, g_start = _longest_prefix_match(mrna, m_pos, genomic, g_pos)
        if length < MIN_EXON:
            remaining = len(mrna) - m_pos
            if exons and remaining <= max_tail:
                return ExonMap(_refine_junctions(exons, genomic), remaining)
            raise AlignmentError(
                f"Could not align the mRNA from position {m_pos + 1}: the sequences must "
                "match exactly (check that the genomic sequence covers the whole gene)"
            )
        exons.append(Exon(m_pos, m_pos + length, g_start, g_start + length))
        m_pos += length
        g_pos = g_start + length
    return ExonMap(_refine_junctions(exons, genomic), 0)


def _longest_prefix_match(mrna: str, m_pos: int, genomic: str, g_pos: int) -> tuple[int, int]:
    """Longest prefix of ``mrna[m_pos:]`` occurring in ``genomic[g_pos:]``."""
    lo, hi = 0, len(mrna) - m_pos
    best_start = -1
    while lo < hi:  # binary search: if a prefix matches, so do all shorter ones
        mid = (lo + hi + 1) // 2
        found = genomic.find(mrna[m_pos : m_pos + mid], g_pos)
        if found >= 0:
            lo, best_start = mid, found
        else:
            hi = mid - 1
    return lo, best_start


def _refine_junctions(exons: list[Exon], genomic: str) -> list[Exon]:
    """Shift each junction left, where the sequence allows, to a GT...AG intron."""
    refined = list(exons)
    for i in range(len(refined) - 1):
        left, right = refined[i], refined[i + 1]
        if left.genomic_end == right.genomic_start:
            continue  # no intron between these segments
        shift = 0
        while (
            shift < left.length - 1
            and genomic[left.genomic_end - shift - 1] == genomic[right.genomic_start - shift - 1]
        ):
            shift += 1
        for s in range(shift + 1):
            intron = genomic[left.genomic_end - s : right.genomic_start - s]
            if intron.startswith("GT") and intron.endswith("AG"):
                refined[i] = Exon(
                    left.mrna_start, left.mrna_end - s, left.genomic_start, left.genomic_end - s
                )
                refined[i + 1] = Exon(
                    right.mrna_start - s, right.mrna_end, right.genomic_start - s, right.genomic_end
                )
                break
    return refined
