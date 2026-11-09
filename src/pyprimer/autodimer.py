"""AutoDimer: screen a set of oligos for primer-dimers and hairpins.

A port of NIST AutoDimer (Vallone & Butler, BioTechniques 37:226-231, 2004),
as served by the NIST DNA Analysis web tools. The original Visual Basic
program is only distributed compiled, so this module was reconstructed from
the paper and calibrated against the NIST web service; the recorded
reference results live in ``tests/data/autodimer_nist.json``.

Every pair of oligos, including each oligo with itself, is slid against each
other, and every hairpin fold (4 or 5 base loops) is tried. Each alignment
gets a score of Watson-Crick matches minus mismatches; alignments scoring at
least the threshold (one for dimers, one for hairpins) are reported with a
rough Tm and dG from SantaLucia (1998) nearest-neighbour parameters.

Quirks of the NIST implementation are kept deliberately (see README.md,
"AutoDimer details"): when the oligos differ in length,
only alignments anchored at either end of the shorter oligo are tried, plus a
few truncated extra shifts; unknown bases count as mismatches; and the
mismatch penalty is not applied to the last position of an alignment.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal

from .sequence import reverse_complement

#: SantaLucia (1998) unified nearest-neighbour parameters: dH (kcal/mol), dS (eu).
NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}  # fmt: skip

#: Initiation with a terminal G-C or A-T pair: dH (kcal/mol), dS (eu).
INIT: dict[str, tuple[float, float]] = {
    "G": (0.1, -2.8), "C": (0.1, -2.8), "A": (2.3, 4.1), "T": (2.3, 4.1),
}  # fmt: skip

#: Symmetry correction for self-complementary duplexes (eu).
SYMMETRY_DS = -1.4
#: Average penalty for a mismatched base pair, added to dH (kcal/mol).
MISMATCH_DH = 0.438
#: Destabilisation of a 4 or 5 base hairpin loop at 37 C (kcal/mol).
LOOP_DG = 5.6
#: Hairpin loop sizes screened.
LOOPS = (4, 5)
#: Upper limit on the number of oligos in one screen (as on the NIST web site).
MAX_OLIGOS = 100

_R = 1.987  # cal/(K mol)
_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
_UNKNOWN = re.compile(r"[^ACGT]")


@dataclass(frozen=True)
class AutoDimerSettings:
    """AutoDimer screening parameters.

    Attributes:
        sodium_m: Monovalent salt concentration (M).
        strand_um: Total strand concentration (uM); affects Tm only.
        dg_temp_c: Temperature at which dG is evaluated (C).
        min_score: Report primer-dimers scoring at least this (matches - mismatches).
        hairpin_min_score: Report hairpins scoring at least this.
    """

    sodium_m: float = 0.085
    strand_um: float = 1.0
    dg_temp_c: float = 37.0
    min_score: int = 7
    hairpin_min_score: int = 4

    def __post_init__(self) -> None:
        if self.sodium_m <= 0:
            raise ValueError("Salt concentration must be greater than zero")
        if self.strand_um <= 0:
            raise ValueError("Total strand concentration must be greater than zero")
        if self.dg_temp_c <= -273.15:
            raise ValueError("dG temperature must be above absolute zero")


DEFAULT_SETTINGS = AutoDimerSettings()


@dataclass(frozen=True)
class Oligo:
    label: str
    sequence: str


@dataclass(frozen=True)
class Dimer:
    """One primer-primer alignment (``self_dimer``: an oligo with itself).

    ``lineup`` draws the longer oligo on top. ``tm`` is ``None`` below 0 C
    and ``dg`` is ``None`` when not negative, which AutoDimer reports as
    "< 0" and "> 0".
    """

    label_1: str
    label_2: str
    sequence_1: str
    sequence_2: str
    matches: int
    score: int
    match_sequence: str
    tm: float | None
    dg: float | None
    lineup: str
    self_dimer: bool = False


@dataclass(frozen=True)
class Hairpin:
    """One hairpin fold (see :class:`Dimer` for ``tm`` and ``dg``)."""

    label: str
    sequence: str
    loop: int
    stem: int
    matches: int
    score: int
    match_sequence: str
    tm: float | None
    dg: float | None
    lineup: str


@dataclass(frozen=True)
class DimerScreen:
    hits: list[Dimer]
    comparisons: int
    settings: AutoDimerSettings


@dataclass(frozen=True)
class HairpinScreen:
    hits: list[Hairpin]
    comparisons: int
    settings: AutoDimerSettings


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------


def normalise(seq: str) -> str:
    """Upper-case, drop whitespace, read U as T and anything else as N."""
    seq = "".join(seq.split()).upper().replace("U", "T")
    return _UNKNOWN.sub("N", seq)


def parse_oligos(text: str) -> list[Oligo]:
    """Parse AutoDimer input: ``>label`` lines, each followed by a sequence.

    Sequence lines are joined until the next label; blank lines are ignored.
    Unlabelled sequences (one per line) are named ``Seq 1``, ``Seq 2``, ...
    """
    oligos: list[Oligo] = []
    label: str | None = None
    body: list[str] = []

    def flush() -> None:
        seq = normalise("".join(body))
        if seq:
            oligos.append(Oligo(label or f"Seq {len(oligos) + 1}", seq))

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            flush()
            label, body = line[1:].strip() or None, []
        elif line:
            body.append(line)
            if label is None:  # unlabelled input: one sequence per line
                flush()
                body = []
    flush()
    if len(oligos) > MAX_OLIGOS:
        raise ValueError(f"Too many sequences ({len(oligos)}; the limit is {MAX_OLIGOS})")
    return oligos


# --------------------------------------------------------------------------
# Thermodynamics
# --------------------------------------------------------------------------


def _duplex(strand: str, paired: list[bool], sodium_m: float) -> tuple[float, float]:
    """dH (kcal/mol) and salt-corrected dS (eu) of an alignment.

    ``strand`` is read 5'->3' and ``paired`` flags its Watson-Crick pairs.
    Stacks are counted only between two adjacent pairs, initiation only for
    paired ends, and the mismatch penalty for every mismatch but the last
    position.
    """
    dh = ds = 0.0
    for i in range(len(strand) - 1):
        if paired[i] and paired[i + 1]:
            h, s = NN[strand[i : i + 2]]
            dh += h
            ds += s
    for end in (0, len(strand) - 1):
        if paired[end]:
            h, s = INIT[strand[end]]
            dh += h
            ds += s
    dh += MISMATCH_DH * paired[:-1].count(False)
    ds += 0.368 * len(strand) * math.log(sodium_m)
    return dh, ds


def _report(dg: float, tm: float) -> tuple[float | None, float | None]:
    return (dg if dg < 0 else None), (tm if tm >= 0 else None)


def dimer_thermo(
    strand: str, paired: list[bool], settings: AutoDimerSettings = DEFAULT_SETTINGS
) -> tuple[float | None, float | None]:
    """``(tm, dg)`` of a primer-dimer alignment (``None``: below 0 C / not negative)."""
    dh, ds = _duplex(strand, paired, settings.sodium_m)
    ct = settings.strand_um * 1e-6 / 4
    if all(paired) and strand == reverse_complement(strand):  # self-complementary duplex
        ds += SYMMETRY_DS
        ct *= 4
    dg = dh - (settings.dg_temp_c + 273.15) * ds / 1000
    denominator = ds + _R * math.log(ct)
    tm = dh * 1000 / denominator - 273.15 if denominator else -273.15
    dg_out, tm_out = _report(dg, tm)
    return tm_out, dg_out


def hairpin_thermo(
    strand: str, paired: list[bool], settings: AutoDimerSettings = DEFAULT_SETTINGS
) -> tuple[float | None, float | None]:
    """``(tm, dg)`` of a hairpin stem plus a 5.6 kcal/mol (entropic) loop penalty.

    Follows the AutoDimer paper; the NIST web service reports different
    hairpin values that could not be reproduced (see README.md).
    """
    dh, ds = _duplex(strand, paired, settings.sodium_m)
    ds -= LOOP_DG * 1000 / 310.15
    dg = dh - (settings.dg_temp_c + 273.15) * ds / 1000
    tm = dh * 1000 / ds - 273.15 if ds else -273.15
    dg_out, tm_out = _report(dg, tm)
    return tm_out, dg_out


# --------------------------------------------------------------------------
# Screening
# --------------------------------------------------------------------------


def _compare(top: str, bottom: str) -> list[bool]:
    return [(a, b) in _PAIRS for a, b in zip(top, bottom, strict=True)]


def _score(paired: list[bool]) -> int:
    matches = sum(paired)
    return matches - (len(paired) - matches)


def _match_sequence(strand: str, paired: list[bool]) -> str:
    return "".join(b if p else "N" for b, p in zip(strand, paired, strict=True))


def _marks(paired: list[bool], mismatch: str) -> str:
    return "".join("|" if p else mismatch for p in paired)


def _dimer_alignments(top: str, bottom: str) -> list[tuple[str, str, str, list[bool]]]:
    """Every alignment tried for two oligos (``top`` at least as long).

    Returns ``(lineup, strand, partner, paired)``, where ``strand`` is the
    5'->3' oligo segment used for thermodynamics and ``partner`` the bases
    it pairs with.
    """
    n, m = len(top), len(bottom)
    rev_top, rev_bottom = top[::-1], bottom[::-1]
    out = []
    for k in range(1, m + 1):
        # 5' ends overlapping by k: top 5'->3' over the bottom reversed.
        pad = m - k
        seg, partner = top[:k], rev_bottom[m - k :]
        paired = _compare(seg, partner)
        lineup = "\n".join(
            [
                " " * (1 + pad) + f"5'-{top}-3'",
                " " * (4 + pad) + _marks(paired, "x"),
                f" 3'-{rev_bottom}-5'",
            ]
        )
        out.append((lineup, seg, partner, paired))
        # 3' ends overlapping by k: top reversed over the bottom 5'->3'.
        seg, partner = bottom[m - k :], rev_top[:k]
        paired = _compare(partner, seg)
        lineup = "\n".join(
            [
                " " * (1 + pad) + f"3'-{rev_top}-5'",
                " " * (4 + pad) + _marks(paired, "x"),
                f" 5'-{bottom}-3'",
            ]
        )
        out.append((lineup, seg, partner, paired))
    # NIST quirk: a few further shifts of the 5'-anchored frame. They are
    # drawn shifted by d but compare the top from base d - 1, and only within
    # its first m - 1 bases. (NIST fails once these would run out of bases;
    # they are capped here instead.)
    for d in range(1, min(n - m - 3, m - 1) + 1):
        seg, partner = top[d - 1 : m - 1], rev_bottom[: m - d]
        paired = _compare(seg, partner)
        lineup = "\n".join(
            [
                f" 5'-{top}-3'",
                " " * (4 + d) + _marks(paired, "x"),
                " " * (1 + d) + f"3'-{rev_bottom}-5'",
            ]
        )
        out.append((lineup, seg, partner, paired))
    return out


def screen_dimers(
    oligos: list[Oligo],
    settings: AutoDimerSettings = DEFAULT_SETTINGS,
    *,
    self_dimers: bool = True,
) -> DimerScreen:
    """Screen every pair of oligos for primer-dimers.

    Pairs are visited as NIST does (last oligo first, each against all
    earlier ones), each oligo followed by its self-dimer; hits keep that
    order. NIST does not screen self-dimers (``self_dimers=False``).
    """
    hits: list[Dimer] = []
    for j in range(len(oligos) - 1, -1, -1):
        one = oligos[j]
        for i, two in enumerate(oligos[: j + 1]):
            is_self = i == j
            if is_self and not self_dimers:
                continue
            top, bottom = (one, two) if len(one.sequence) > len(two.sequence) else (two, one)
            alignments = _dimer_alignments(top.sequence, bottom.sequence)
            if is_self:  # the full overlap is drawn from both ends; keep one
                del alignments[2 * len(one.sequence) - 1]
            for lineup, strand, _, paired in alignments:
                score = _score(paired)
                if score < settings.min_score:
                    continue
                tm, dg = dimer_thermo(strand, paired, settings)
                hits.append(
                    Dimer(
                        label_1=one.label,
                        label_2=two.label,
                        sequence_1=one.sequence,
                        sequence_2=two.sequence,
                        matches=sum(paired),
                        score=score,
                        match_sequence=_match_sequence(strand, paired),
                        tm=tm,
                        dg=dg,
                        lineup=lineup,
                        self_dimer=is_self,
                    )
                )
    n = len(oligos)
    pairs = n * (n + 1) // 2 if self_dimers else n * (n - 1) // 2
    return DimerScreen(hits, pairs, settings)


def _hairpin_folds(seq: str) -> list[tuple[int, int, str, str, list[bool], str]]:
    """Every hairpin fold tried: ``(loop, stem, arm5, arm3, paired, lineup)``.

    The 5' part (``a`` bases, the last two in the loop) is drawn reversed
    above the 3' part; ``arm3`` is the 3' stem strand read 5'->3'.
    """
    n = len(seq)
    out = []
    for a in range(4, n - 5):
        top = seq[:a][::-1]
        for loop in LOOPS:
            bottom = seq[a + loop - 4 :]
            stem = min(a, len(bottom)) - 2
            arm5, arm3 = top[2 : 2 + stem], bottom[2 : 2 + stem]
            paired = _compare(arm5, arm3)
            marker = "{" if loop == 4 else seq[a]
            lineup = "\n".join([f" {top}-5'", f"{marker}  {_marks(paired, 'X')}", f" {bottom}-3'"])
            out.append((loop, stem, arm5, arm3, paired, lineup))
    return out


def screen_hairpins(
    oligos: list[Oligo], settings: AutoDimerSettings = DEFAULT_SETTINGS
) -> HairpinScreen:
    """Screen each oligo for hairpins with 4 or 5 base loops."""
    hits: list[Hairpin] = []
    for oligo in oligos:
        for loop, stem, _, arm3, paired, lineup in _hairpin_folds(oligo.sequence):
            score = _score(paired)
            if score < settings.hairpin_min_score:
                continue
            tm, dg = hairpin_thermo(arm3, paired, settings)
            hits.append(
                Hairpin(
                    label=oligo.label,
                    sequence=oligo.sequence,
                    loop=loop,
                    stem=stem,
                    matches=sum(paired),
                    score=score,
                    match_sequence=_match_sequence(arm3, paired),
                    tm=tm,
                    dg=dg,
                    lineup=lineup,
                )
            )
    return HairpinScreen(hits, len(oligos), settings)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------


def format_value(value: float | None, missing: str) -> str:
    """Format to one decimal place, rounding halves up as NIST does."""
    if value is None:
        return missing
    return str(Decimal(value).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


def report_text(screen: DimerScreen | HairpinScreen) -> str:
    """Plain-text results in the style of the AutoDimer output file."""
    s = screen.settings
    hairpins = isinstance(screen, HairpinScreen)
    kind = "Hairpin" if hairpins else "Primer-Primer"
    threshold = s.hairpin_min_score if hairpins else s.min_score
    temp = f"{s.dg_temp_c:g}"
    lines = [
        f"Na+ = {s.sodium_m:g} M   Ct = {s.strand_um:g} uM   Threshold Score = {threshold}",
        f"Total Number of {kind} Comparisons = {screen.comparisons}",
        f"Number of Hits = {len(screen.hits)}",
        "-" * 60,
    ]
    for hit in screen.hits:
        if isinstance(hit, Hairpin):
            names = hit.label
        elif hit.self_dimer:
            names = f"{hit.label_1} (self-dimer)"
        else:
            names = f"{hit.label_1} versus {hit.label_2}"
        lines += [
            "",
            names,
            hit.lineup,
            f"Matches = {hit.matches}   Score = {hit.score}",
            f"est. Tm = {format_value(hit.tm, 'less than 0')} C"
            f"   dG @ {temp} C = {format_value(hit.dg, '> 0')} kcal/mol",
        ]
    if not screen.hits:
        what = "HAIRPIN" if hairpins else "PRIMER-DIMER"
        lines.append(f"\nNO {what} HITS FOR THIS SCORE THRESHOLD")
    return "\n".join(lines) + "\n"


def full_report(dimers: DimerScreen, hairpins: HairpinScreen) -> str:
    """Primer-dimer and hairpin results in one text report."""
    return "PRIMER-DIMERS\n\n" + report_text(dimers) + "\n\nHAIRPINS\n\n" + report_text(hairpins)
