"""Primer analysis and plain-text reports."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime

from .design import RepeatFilter
from .dimer import Dimer, primer_dimer
from .features import Region
from .sequence import GENETIC_CODE, gc_percent, is_degenerate
from .thermo import DEFAULT_CONDITIONS, Conditions, thermodynamics

MAP_WIDTH = 80


@dataclass(frozen=True)
class PrimerAnalysis:
    seq: str
    length: int
    tm: float
    gc: float
    dh: float
    ds: float
    dg: float
    run: bool
    repeat: str | None


@dataclass(frozen=True)
class LabelledDimer:
    label: str
    dimer: Dimer


@dataclass
class PairAnalysis:
    forward: PrimerAnalysis | None
    reverse: PrimerAnalysis | None
    extensible: list[LabelledDimer] = field(default_factory=list)
    non_extensible: list[LabelledDimer] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    dimer_temp_c: float = 37.0

    def dimer_text(self) -> str:
        """Primer warnings and dimer diagrams, formatted as in PerlPrimer."""
        lines = []
        for name, primer in (("forward", self.forward), ("reverse", self.reverse)):
            if primer and primer.run:
                lines.append(f"Warning: {name} primer run found")
            if primer and primer.repeat:
                lines.append(f"Warning: {name} primer repeat found : {primer.repeat}")
        temp = f"{self.dimer_temp_c:g}"
        lines.append(f"Most stable 3' extensible primer-dimers (at {temp}°C), if any\n")
        lines.extend(_format_dimer(d) for d in self.extensible)
        lines.append(f"\nMore stable non-extensible primer-dimers (at {temp}°C), if any\n")
        lines.extend(_format_dimer(d) for d in self.non_extensible)
        return "\n".join(lines) + "\n"


def _format_dimer(d: LabelledDimer) -> str:
    return f"{d.label}: {d.dimer.dg:.2f} kcal/mol\n\n{d.dimer.diagram()}"


def analyse_primer(
    seq: str,
    conditions: Conditions = DEFAULT_CONDITIONS,
    repeats: RepeatFilter | None = None,
) -> PrimerAnalysis:
    """Calculate Tm, composition and thermodynamics of a single primer."""
    repeats = repeats or RepeatFilter()
    seq = seq.upper().strip()
    t = thermodynamics(seq, conditions)
    return PrimerAnalysis(
        seq=seq,
        length=len(seq),
        tm=round(t.tm, 2),
        gc=round(gc_percent(seq), 1),
        dh=round(t.dh, 2),
        ds=round(t.ds, 2),
        dg=round(t.dg, 2),
        run=repeats.run_found(seq),
        repeat=repeats.repeat_found(seq),
    )


def analyse_pair(
    forward: str | None,
    reverse: str | None,
    conditions: Conditions = DEFAULT_CONDITIONS,
    repeats: RepeatFilter | None = None,
) -> PairAnalysis:
    """Analyse one or two primers and the dimers they can form.

    Lists extensible 3' dimers with dG < 0, then any non-extensible dimers
    that are more stable than the corresponding extensible dimer.
    """
    result = PairAnalysis(None, None, dimer_temp_c=conditions.dimer_temp_c)
    primers: dict[str, str] = {}
    for name, seq in (("Forward", forward), ("Reverse", reverse)):
        seq = (seq or "").upper().strip()
        if not seq:
            continue
        if is_degenerate(seq):
            result.errors.append(
                f"The {name.lower()} primer contains a degenerate or non-DNA character; "
                "the Tm of degenerate sequences cannot be calculated"
            )
            continue
        primers[name] = seq
        analysis = analyse_primer(seq, conditions, repeats)
        if name == "Forward":
            result.forward = analysis
        else:
            result.reverse = analysis

    combos = [
        (f"{a} vs. {b}", primers[a], primers[b])
        for a, b in (("Forward", "Forward"), ("Forward", "Reverse"), ("Reverse", "Reverse"))
        if a in primers and b in primers
    ]
    for label, a, b in combos:
        extensible = primer_dimer(a, b, conditions)
        if extensible.dg != 0:
            result.extensible.append(LabelledDimer(label, extensible))
        full = primer_dimer(a, b, conditions, full=True)
        if full.dg < extensible.dg:
            result.non_extensible.append(LabelledDimer(label, full))
    return result


@dataclass(frozen=True)
class MapPrimers:
    """Primer placement for a sequence map (top-strand coordinates)."""

    forward: str | None = None
    forward_pos: int = 0
    reverse: str | None = None
    reverse_end: int = 0
    reverse_len: int = 0
    forward_len: int = 0


@dataclass(frozen=True)
class MapLines:
    """Unwrapped lines of a sequence map, aligned column by column."""

    ruler: str
    primers: str
    sequence: str
    features: str


def map_lines(
    seq: str,
    primers: MapPrimers | None = None,
    features: list[Region] | None = None,
    *,
    translate: bool = True,
) -> MapLines:
    """Build the ruler, primer, sequence and translation lines of a map.

    With ``translate=False`` features are drawn as dashes (e.g. CpG islands).
    """
    return MapLines(
        ruler=_ruler(len(seq)),
        primers=_primer_line(primers),
        sequence=seq,
        features=_feature_line(seq, primers, features or [], translate),
    )


def sequence_map(
    seq: str,
    primers: MapPrimers | None = None,
    features: list[Region] | None = None,
    *,
    translate: bool = True,
    width: int = MAP_WIDTH,
) -> str:
    """Render a text map wrapped to ``width`` columns."""
    m = map_lines(seq, primers, features, translate=translate)
    lines = [m.ruler, m.primers, m.sequence, m.features]
    blocks = [
        "\n".join(line[start : start + width].rstrip() for line in lines)
        for start in range(0, len(seq), width)
    ]
    return "\n\n".join(blocks) + "\n"


def _ruler(length: int) -> str:
    ruler = "0........."
    for i in range(1, length // 10 + 1):
        ruler += str(i * 10) + "." * (9 - len(str(i)))
    return ruler


def _primer_line(p: MapPrimers | None) -> str:
    if not p or not p.forward:
        return ""
    line = " " * (p.forward_pos - 3) + f"5' {p.forward} 3'"
    if p.reverse:
        gap = p.reverse_end - p.reverse_len - p.forward_pos - p.forward_len - 5
        line += " " * gap + f"3' {p.reverse[::-1]} 5'"
    return line


def _feature_line(
    seq: str, primers: MapPrimers | None, features: list[Region], translate: bool
) -> str:
    out = ""
    prev_end = 0
    for region in features:
        start, end = region.start, region.end
        if not translate:
            out += " " * (start - prev_end) + "-" * (end - start)
            prev_end = end
            continue
        # Translate from the forward primer if it lies upstream of the ORF.
        fpos = primers.forward_pos if primers and primers.forward else 0
        if fpos and fpos < start:
            start = fpos + (region.start % 3 - fpos % 3)
        peptide = "".join(
            GENETIC_CODE.get(seq[i : i + 3].upper(), "X") + "  " for i in range(start, end, 3)
        )
        out += " " * (start - prev_end) + peptide
        prev_end = end
    return out


def text_report(
    title: str,
    analysis: PairAnalysis,
    seq_map: str,
    *,
    forward_label: str | None = None,
    reverse_label: str | None = None,
    amplicon: str = "",
) -> str:
    """Assemble a full plain-text report for a primer (pair)."""
    stamp = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    parts = [f"[{title} - Report generated at {stamp}]", "", "PRIMERS", "-------", ""]
    if analysis.forward:
        parts.append(f"Forward: 5' {forward_label or analysis.forward.seq} 3'")
    if analysis.reverse:
        parts.append(f"Reverse: 5' {reverse_label or analysis.reverse.seq} 3'")
    parts += ["", amplicon, "", "PRIMER DETAILS", "--------------", ""]
    temp = f"{analysis.dimer_temp_c:g}"
    for name, p in (("Forward", analysis.forward), ("Reverse", analysis.reverse)):
        if p:
            parts += [
                f"{name} primer: {p.seq}",
                "",
                f"Tm: {p.tm:.2f}°C\t\tLength: {p.length} bases",
                f"dS°: {p.ds:.2f} eu\t\tdH°: {p.dh:.2f} kcal/mol",
                f"dG°{temp}: {p.dg:.2f} kcal/mol",
                "",
                "",
            ]
    parts += [analysis.dimer_text(), "", "SEQUENCE MAP", "------------", "", seq_map]
    return "\n".join(parts)
