"""Basic DNA sequence utilities: parsing, complementation and base composition."""

from __future__ import annotations

import re

_COMPLEMENT = str.maketrans("AGCTagct", "TCGAtcga")
_NON_DNA = re.compile(r"[^ATGC]", re.IGNORECASE)
_STRIP = re.compile(r"[\s\d-]")

#: Standard genetic code (stop codons are ``*``).
GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}  # fmt: skip


def parse_sequence(text: str) -> tuple[str | None, str]:
    """Parse raw or FASTA-formatted input into ``(name, sequence)``.

    Only the first FASTA record is used. Whitespace, digits and gap
    characters (``-``) are removed; letter case is preserved because
    PerlPrimer uses case to mark features (ORFs, CpG islands).
    """
    name: str | None = None
    lines = text.lstrip().splitlines()
    if lines and lines[0].startswith(">"):
        name = lines[0][1:].strip() or None
        body: list[str] = []
        for line in lines[1:]:
            if line.startswith(">"):
                break
            body.append(line)
        text = "\n".join(body)
    return name, _STRIP.sub("", text)


def clean_sequence(text: str) -> str:
    """Return the bare sequence from raw or FASTA-formatted input."""
    return parse_sequence(text)[1]


def complement(seq: str) -> str:
    """Complement a DNA sequence (case preserved, other characters untouched)."""
    return seq.translate(_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA sequence."""
    return complement(seq)[::-1]


def bisulphite_convert(seq: str) -> str:
    """Convert every (upper-case) C to T, as for fully unmethylated DNA."""
    return seq.replace("C", "T")


def is_degenerate(seq: str) -> bool:
    """True if the sequence contains anything other than A, C, G or T."""
    return bool(_NON_DNA.search(seq))


def gc_percent(seq: str) -> float:
    """Percentage of G and C bases (case-insensitive)."""
    upper = seq.upper()
    return (upper.count("G") + upper.count("C")) / len(seq) * 100


def c_percent(seq: str) -> float:
    """Percentage of upper-case C bases (as used for bisulphite primers)."""
    return seq.count("C") / len(seq) * 100


def translate(seq: str, start: int, end: int) -> str:
    """Translate codons from ``start`` up to ``end`` (unknown codons give X)."""
    return "".join(GENETIC_CODE.get(seq[i : i + 3].upper(), "X") for i in range(start, end, 3))
