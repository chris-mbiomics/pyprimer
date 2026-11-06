"""Primer-dimer prediction.

Each primer pair is aligned at every offset using a binding matrix. Single
matched bases flanked by mismatches are discarded, and the free energy of
the remaining duplex is summed with nearest-neighbour (including mismatch)
parameters. The most stable (most negative dG) alignment is reported.

By default only *extensible* 3' dimers are considered, i.e. alignments in
which both 3' terminal bases are paired; ``full=True`` scores every offset.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

from .sequence import is_degenerate, reverse_complement
from .thermo import DEFAULT_CONDITIONS, Conditions, free_energy_table

_UNSTABLE_SINGLE_MATCH = re.compile(r"01(?=[^1])")
# XOR of two aligned byte strings is 0 where the bases match.
_MATCH_TABLE = bytes(ord("1") if i == 0 else ord("0") for i in range(256))


@dataclass(frozen=True)
class Dimer:
    """The most stable dimer found between two primers.

    Attributes:
        primer_f: First primer (5'->3').
        primer_r: Second primer (5'->3').
        dg: Free energy of the most stable dimer (kcal/mol); 0 if none.
        offset: Alignment offset of that dimer, or None if there is none.
        binding: Binding string at that offset ('1' paired, '0' unpaired,
            '2' overhang), or None.
    """

    primer_f: str
    primer_r: str
    dg: float
    offset: int | None
    binding: str | None

    @property
    def found(self) -> bool:
        return self.offset is not None and self.dg != 0

    def diagram(self) -> str:
        """Render the dimer alignment as three lines of text."""
        if self.offset is None or self.binding is None:
            return ""
        pos = self.offset
        pfl, prl = len(self.primer_f), len(self.primer_r)
        pl = max(pfl, prl)

        fspace = " " * (max(pl - pfl, 0) + (pos - pl + 1 if pos > pl - 1 else 0))
        rspace = " " * (pl - pos - 1) if pos + 1 >= pfl else fspace
        pr_space = " " * (pfl - pos - 1)

        marks = []
        for j in range(min(pos + 1, prl)):
            base = self.binding[j] if j < len(self.binding) else "0"
            marks.append({"1": "|", "0": "."}.get(base, " "))

        return (
            f"{fspace}5' {self.primer_f} 3'\n"
            f"{rspace}   {pr_space}{''.join(marks)}\n"
            f"{rspace}{pr_space}3' {self.primer_r[::-1]} 5'\n"
        )


def primer_dimer(
    primer_f: str,
    primer_r: str,
    conditions: Conditions = DEFAULT_CONDITIONS,
    *,
    full: bool = False,
    extensible_only: bool = True,
) -> Dimer:
    """Find the most stable dimer formed between two primers.

    Args:
        primer_f: First primer, 5'->3'.
        primer_r: Second primer, 5'->3' (may equal ``primer_f`` for self-dimers).
        conditions: Reaction conditions used for the dG parameters.
        full: Score every alignment, not just 3' overlaps.
        extensible_only: When not ``full``, only score alignments whose 3'
            ends are both paired (i.e. that DNA polymerase can extend).

    Raises:
        ValueError: if either primer is empty or degenerate.
    """
    primer_f, primer_r = primer_f.upper(), primer_r.upper()
    for primer in (primer_f, primer_r):
        if not primer or is_degenerate(primer):
            raise ValueError(f"Invalid primer sequence for dimer calculation: {primer!r}")
    return _primer_dimer(primer_f, primer_r, conditions, full, extensible_only)


@lru_cache(maxsize=200_000)
def _primer_dimer(
    primer_f: str, primer_r: str, conditions: Conditions, full: bool, extensible_only: bool
) -> Dimer:
    dg_table = free_energy_table(conditions)
    pfl, prl = len(primer_f), len(primer_r)
    pl = max(pfl, prl)

    rprimer_r = primer_r[::-1]
    # At offset k, base j of primer_r (read 3'->5') faces base pfl-1-k+j of
    # primer_f; they pair exactly when that base equals rcomp_r[j].
    rcomp_r = reverse_complement(primer_r)
    f_bytes, r_bytes = primer_f.encode(), rcomp_r.encode()

    def pairs(j: int, i: int) -> bool:
        """Binding matrix entry: base j of rcomp_r vs base i of reversed primer_f."""
        return primer_f[pfl - 1 - i] == rcomp_r[j]

    scores: list[float] = []
    best: dict[float, tuple[int, str]] = {}
    last_offset = pfl + prl - 1 if full else pl - 2
    for k in range(last_offset + 1):
        scores.append(0.0)

        if extensible_only and not full:
            start = min(k, pfl - 1)
            end = min(k, prl - 1)
            if not pairs(0, start):
                continue
            col = start - k  # <= 0: PerlPrimer indexes from the end of the row
            if -col > pfl or not pairs(end, col % pfl):
                continue

        # Leading overhang ("2"), then match/mismatch bits for the aligned bases.
        overhang = max(0, min(prl, k - pfl + 1))
        last = min(k, prl - 1)
        bits = ""
        if last >= overhang:
            a = f_bytes[pfl - 1 - k + overhang : pfl - k + last]
            b = r_bytes[overhang : last + 1]
            diff = int.from_bytes(a, "big") ^ int.from_bytes(b, "big")
            bits = diff.to_bytes(len(a), "big").translate(_MATCH_TABLE).decode("ascii")
        binding = _UNSTABLE_SINGLE_MATCH.sub("00", "2" * overhang + bits)

        if "1" not in binding:
            continue
        first, last = binding.index("1"), binding.rindex("1")

        dg = 0.0
        for j in range(first, last):
            if binding[j] == "2" or binding[j : j + 2] == "00":
                continue
            f_index = pfl - k + j - 1
            dg += dg_table.get(primer_f[f_index : f_index + 2] + rprimer_r[j : j + 2], 0.0)
        dg += dg_table["init" + rprimer_r[first]]
        dg += dg_table["init" + rprimer_r[last]]

        score = round(dg, 2)
        scores[k] = score
        best[score] = (k, binding)  # later offsets win ties, as in PerlPrimer

    min_score = min(scores)
    if min_score in best:
        offset, binding = best[min_score]
        return Dimer(primer_f, primer_r, min_score, offset, binding)
    return Dimer(primer_f, primer_r, min_score, None, None)
