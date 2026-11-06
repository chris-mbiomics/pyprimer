"""Restriction enzyme sites (REBASE) and cloning tails for primers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cache
from importlib import resources

_GCG_LINE = re.compile(r"^;*([A-Za-z0-9;]+)\s+\d+\s+([A-z'_]+)\s+")
_AMBIGUOUS = re.compile(r"[rymkswbdhvn]", re.IGNORECASE)

DEFAULT_ANCHOR = "GCGCGC"


@dataclass(frozen=True)
class Enzyme:
    """A restriction enzyme.

    Attributes:
        name: Enzyme name, e.g. ``EcoRI``.
        pattern: REBASE notation (``'`` top-strand cut, ``_`` bottom-strand cut).
        site: Recognition sequence without cut marks.
    """

    name: str
    pattern: str
    site: str

    @property
    def label(self) -> str:
        return f"{self.name} - {self.pattern}"


def parse_rebase(text: str) -> list[Enzyme]:
    """Parse a REBASE file in GCG format."""
    enzymes = []
    for line in text.splitlines():
        match = _GCG_LINE.match(line)
        if not match:
            continue
        name, pattern = match.groups()
        enzymes.append(Enzyme(name, pattern, pattern.replace("'", "").replace("_", "").upper()))
    return enzymes


@cache
def rebase_enzymes() -> tuple[Enzyme, ...]:
    """Enzymes from the bundled REBASE database (version 603)."""
    data = resources.files("pyprimer").joinpath("data/gcg.603").read_text(encoding="latin-1")
    return tuple(parse_rebase(data))


def cloning_enzymes(region: str = "", exclude_found: bool = True) -> list[Enzyme]:
    """Non-degenerate six-cutters suitable for cloning, sorted by name.

    Args:
        region: The amplified sequence.
        exclude_found: Omit enzymes that cut within ``region``.
    """
    region = region.upper()
    chosen = {}
    for enzyme in rebase_enzymes():
        if len(enzyme.site) != 6 or _AMBIGUOUS.search(enzyme.site):
            continue
        if exclude_found and enzyme.site in region:
            continue
        chosen[enzyme.name] = enzyme
    return sorted(chosen.values(), key=lambda e: e.name.lower())


def find_enzyme(name: str) -> Enzyme | None:
    for enzyme in rebase_enzymes():
        if enzyme.name.lower() == name.lower():
            return enzyme
    return None


def cloning_tail(enzyme: Enzyme, anchor: str = DEFAULT_ANCHOR) -> str:
    """Return a primer tail ``anchor_site`` (the underscore marks the site)."""
    return f"{anchor.upper()}_{enzyme.site}"


@dataclass(frozen=True)
class ClonedPrimers:
    forward: str | None
    reverse: str | None


def add_cloning_tails(
    forward: str | None,
    forward_pos: int,
    reverse: str | None,
    reverse_end: int,
    *,
    gene_frame: int = 0,
    forward_tail: str | None = None,
    reverse_tail: str | None = None,
    forward_frame: int | None = None,
    reverse_frame: int | None = None,
) -> ClonedPrimers:
    """Add 5' cloning tails to a primer pair, optionally keeping the ORF in frame.

    Tails are written as ``anchor_site`` (e.g. ``GCGCGC_GAATTC``). When a
    frame is given (the frame of the site in the destination vector: 0, 1 or
    2), A spacers are inserted so that the fusion stays in frame. The result
    joins tail and primer with ``_`` to make the primer boundary visible.

    Args:
        forward / reverse: Primer sequences (5'->3'); None to skip.
        forward_pos: Top-strand position of the forward primer.
        reverse_end: Top-strand position of the reverse primer's 5' end.
        gene_frame: Frame of the ORF (its start position modulo 3).
    """
    cloned_f = cloned_r = None
    if forward and forward_tail:
        site = forward_tail.split("_")[1] if "_" in forward_tail else ""
        site_offset = 3 - len(site)
        tail = _strip_marks(forward_tail)
        spacer = ""
        if forward_frame is not None:
            join_frame = 3 - forward_frame
            spacer = "A" * (((forward_pos - gene_frame) % 3 + join_frame + site_offset) % 3)
        cloned_f = f"{tail}{spacer}_{forward}".upper()
    if reverse and reverse_tail:
        tail = _strip_marks(reverse_tail)
        spacer = ""
        if reverse_frame is not None:
            spacer = "A" * ((2 - (reverse_end - gene_frame) % 3 + reverse_frame) % 3)
        cloned_r = f"{tail}{spacer}_{reverse}".upper()
    return ClonedPrimers(cloned_f, cloned_r)


def _strip_marks(tail: str) -> str:
    return tail.replace("|", "").replace("_", "")
