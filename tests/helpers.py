"""Deterministic test sequence generators."""

from __future__ import annotations

import random

from pyprimer.sequence import GENETIC_CODE

SENSE_CODONS = sorted(codon for codon, aa in GENETIC_CODE.items() if aa != "*")


def random_dna(rng: random.Random, length: int, gc: float = 0.5) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def orf_sequence(rng: random.Random, upstream: int, codons: int, downstream: int) -> str:
    """Random sequence with an ATG...stop ORF of ``codons`` sense codons."""
    orf = "ATG" + "".join(rng.choices(SENSE_CODONS, k=codons)) + rng.choice(["TAA", "TAG", "TGA"])
    return random_dna(rng, upstream) + orf + random_dna(rng, downstream)


def cpg_rich(rng: random.Random, length: int) -> str:
    """GC-rich sequence with frequent CpG dinucleotides."""
    out = []
    while sum(map(len, out)) < length:
        out.append(rng.choice(["CG", "CG", "GC", "C", "G", "A", "T", "CCG", "GCG"]))
    return "".join(out)[:length]
