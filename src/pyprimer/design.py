"""Primer design for standard PCR, bisulphite PCR, real-time PCR and sequencing.

Positions are 0-based. Forward primers are described by their 5' position on
the given (top) strand. Reverse primers are described both by their position
on the reverse-complement strand (``reverse_pos``) and by the top-strand
coordinate of their 5' end (``reverse_end``), so that a reverse primer covers
``[reverse_end - reverse_len + 1, reverse_end]`` of the top strand.

Target ranges follow PerlPrimer: the forward primer must start between
``outer_start`` and ``inner_start``, and the reverse primer must end between
``inner_end`` and ``outer_end``.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field, replace
from enum import Enum

from .dimer import primer_dimer
from .features import Region, find_cpg_islands, find_orf
from .sequence import (
    bisulphite_convert,
    c_percent,
    clean_sequence,
    complement,
    gc_percent,
    is_degenerate,
    reverse_complement,
)
from .thermo import DEFAULT_CONDITIONS, Conditions, melting_temperature

#: How far each automatic search step moves the range (bp).
SEARCH_STEP = 10


class DesignError(ValueError):
    """Raised when primer design cannot run with the given input."""


class Search(str, Enum):
    """How to adjust the target range when no primer pairs are found."""

    FIXED = "fixed"  # use the ranges as given
    OUTWARDS = "outwards"  # widen the outer range 10 bp at a time
    INWARDS = "inwards"  # narrow the inner range 10 bp at a time


# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RepeatFilter:
    """Exclusion of mononucleotide runs and dinucleotide (or longer) repeats.

    Attributes:
        enabled: Apply the filter.
        run: Exclude runs of this many identical bases or more.
        repeat: Exclude motifs of two or more bases repeated this many times.
    """

    enabled: bool = True
    run: int = 4
    repeat: int = 3

    def __post_init__(self) -> None:
        if self.run < 2 or self.repeat < 2:
            raise ValueError("Run and repeat lengths must be at least 2")

    def _patterns(self) -> tuple[re.Pattern[str], re.Pattern[str]]:
        n = self.run
        return (
            re.compile(f"A{{{n},}}|C{{{n},}}|G{{{n},}}|T{{{n},}}"),
            re.compile(f"(.{{2,}})\\1{{{self.repeat - 1},}}"),
        )

    def run_found(self, seq: str) -> bool:
        return bool(self._patterns()[0].search(seq))

    def repeat_found(self, seq: str) -> str | None:
        """Return the repeated motif, if any."""
        match = self._patterns()[1].search(seq)
        return match.group(1) if match else None

    def rejects(self, seq: str) -> bool:
        return self.enabled and (self.run_found(seq) or self.repeat_found(seq) is not None)


@dataclass(frozen=True)
class GCFilter:
    """%GC limits applied when GC exclusion is enabled."""

    min_gc: float = 40.0
    max_gc: float = 60.0


@dataclass(frozen=True)
class PCRSettings:
    """Settings for standard PCR (and the base for real-time PCR)."""

    min_tm: float = 57.0
    max_tm: float = 63.0
    max_tm_diff: float = 3.0
    min_length: int = 20
    max_length: int = 24
    exclude_gc: bool = True
    gc_clamp: bool = True
    gc: GCFilter = field(default_factory=GCFilter)
    repeats: RepeatFilter = field(default_factory=RepeatFilter)

    def __post_init__(self) -> None:
        _check_lengths(self.min_length, self.max_length)
        _check_tm(self.min_tm, self.max_tm)


@dataclass(frozen=True)
class QPCRSettings(PCRSettings):
    """Settings for real-time PCR primers on a spliced (mRNA/cDNA) sequence.

    Attributes:
        min_amplicon / max_amplicon: Amplicon size limits (bp).
        span_boundary: Amplicon must contain an exon-exon junction.
        primer_on_boundary: At least one primer must straddle a junction ...
        boundary_overlap: ... with at least this many bases on either side.
        first_exon / last_exon: Optionally restrict primers to these exons
            (1-based, inclusive).
    """

    min_tm: float = 58.0
    max_tm: float = 62.0
    max_tm_diff: float = 2.0
    min_amplicon: int = 100
    max_amplicon: int = 300
    span_boundary: bool = True
    primer_on_boundary: bool = True
    boundary_overlap: int = 7
    first_exon: int | None = None
    last_exon: int | None = None


@dataclass(frozen=True)
class BisulphiteSettings:
    """Settings for bisulphite PCR primers.

    Attributes:
        min_c_percent: Minimum %C of the unconverted primer (gives
            specificity for the converted template).
        exclude_cpg: Reject primers containing any CpG (otherwise up to two
            are allowed and made degenerate).
        three_prime_c: Require a 3' C, or two Cs in the last three bases.
        repeats_after_conversion: Check runs/repeats on the converted
            sequence rather than the original.
    """

    min_tm: float = 55.0
    max_tm: float = 65.0
    max_tm_diff: float = 5.0
    min_length: int = 25
    max_length: int = 30
    min_c_percent: float = 30.0
    exclude_cpg: bool = False
    three_prime_c: bool = True
    repeats_after_conversion: bool = True
    repeats: RepeatFilter = field(default_factory=RepeatFilter)

    def __post_init__(self) -> None:
        _check_lengths(self.min_length, self.max_length)
        _check_tm(self.min_tm, self.max_tm)


@dataclass(frozen=True)
class SequencingSettings:
    """Settings for a tiled set of sequencing primers.

    Attributes:
        min_spacing / max_spacing: Distance between consecutive primers (bp).
        exclude_dimers: Reject primers whose self-dimer dG is below
            ``-max_dimer_dg`` kcal/mol.
    """

    min_tm: float = 56.0
    max_tm: float = 64.0
    min_length: int = 20
    max_length: int = 24
    min_spacing: int = 500
    max_spacing: int = 700
    exclude_gc: bool = True
    gc_clamp: bool = True
    exclude_dimers: bool = True
    max_dimer_dg: float = 5.0
    gc: GCFilter = field(default_factory=GCFilter)
    repeats: RepeatFilter = field(default_factory=RepeatFilter)

    def __post_init__(self) -> None:
        _check_lengths(self.min_length, self.max_length)
        _check_tm(self.min_tm, self.max_tm)
        if not 0 < self.min_spacing <= self.max_spacing:
            raise ValueError("Primer spacing must satisfy 0 < min <= max")


def _check_lengths(min_length: int, max_length: int) -> None:
    if not 1 < min_length <= max_length:
        raise ValueError("Primer length must satisfy 1 < min <= max")


def _check_tm(min_tm: float, max_tm: float) -> None:
    if min_tm >= max_tm:
        raise ValueError("Minimum Tm must be below maximum Tm")


# --------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class TargetRange:
    """Target region (see module docstring). Any field may be None (unset)."""

    outer_start: int | None = None
    inner_start: int | None = None
    inner_end: int | None = None
    outer_end: int | None = None

    @property
    def is_set(self) -> bool:
        return self.inner_start is not None and self.inner_end is not None


@dataclass(frozen=True)
class Candidate:
    """A single primer that passed all filters."""

    pos: int
    length: int
    seq: str
    tm: float
    original: str | None = None  # bisulphite: unconverted template window


@dataclass(frozen=True)
class PrimerPair:
    """A forward/reverse primer pair."""

    forward: str
    forward_pos: int
    forward_len: int
    forward_tm: float
    reverse: str
    reverse_pos: int
    reverse_len: int
    reverse_tm: float
    reverse_end: int
    amplicon: int
    dimer_dg: float  # most stable extensible 3' dimer
    full_dimer_dg: float  # most stable dimer of any kind
    forward_degenerate: str | None = None  # bisulphite: CpG-degenerate forms
    reverse_degenerate: str | None = None


@dataclass(frozen=True)
class SequencingPrimer:
    seq: str
    pos: int
    length: int
    tm: float
    dimer_dg: float


@dataclass
class DesignResult:
    pairs: list[PrimerPair]
    target: TargetRange
    sequence: str
    messages: list[str] = field(default_factory=list)


@dataclass
class SequencingResult:
    primers: list[SequencingPrimer]
    target: TargetRange
    sequence: str
    messages: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------
# Window scanning
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class _Window:
    """Everything that determines where primers are searched for."""

    min_len: int
    max_len: int
    target: TargetRange = TargetRange()
    min_amplicon: int | None = None
    max_amplicon: int | None = None
    exon_bounds: tuple[int, ...] = ()
    first_exon: int | None = None
    last_exon: int | None = None


def _scan_bounds(length: int, reverse: bool, w: _Window) -> tuple[int, int] | None:
    """Return the inclusive range of start positions to scan, or None.

    Mirrors PerlPrimer, which treats unset *and zero* values as "not set" here.
    """
    t = w.target
    inner_start, inner_end = t.inner_start, t.inner_end
    max_amp, min_amp = w.max_amplicon or 0, w.min_amplicon or 0
    lo = inner_start or 0
    hi = inner_end

    if not reverse and inner_end and max_amp and inner_end - max_amp > 0:
        lo = inner_end - max_amp
    if reverse and inner_start and length - inner_start - max_amp > 0:
        lo = length - inner_start - max_amp

    if max_amp and min_amp:
        if not reverse:
            hi = min(inner_start, length - min_amp) if inner_start else length - min_amp
        else:
            hi = length - inner_end + w.max_len if inner_end else length

    # Real-time PCR: restrict to a range of exons.
    bounds = w.exon_bounds
    first = (w.first_exon or 0) - 2  # boundary preceding the first exon
    last = max((w.last_exon or 0) - 1, 0)  # boundary following the last exon
    if not reverse:
        if w.first_exon and 0 <= first < len(bounds):
            lo = bounds[first] - w.max_len
        if w.last_exon and last < len(bounds):
            hi = bounds[last] + w.max_len
    else:
        if w.last_exon and last < len(bounds):
            lo = length - (bounds[last] + w.max_len)
        if w.first_exon and 0 <= first < len(bounds):
            hi = length - (bounds[first] - w.max_len)

    # Explicit outer ranges take precedence.
    if t.outer_start is not None and t.outer_end is not None:
        inner_start, inner_end = t.inner_start or 0, t.inner_end or 0
        if not (t.outer_start == inner_start and t.outer_end == inner_end):
            if not reverse:
                lo, hi = t.outer_start, inner_start + w.max_len
            else:
                lo, hi = length - t.outer_end - 1, length - inner_end + w.max_len

    if hi is None:
        return None
    return lo, hi


def _windows(seq: str, reverse: bool, w: _Window) -> Iterator[tuple[int, str]]:
    """Yield ``(position, window)`` for every candidate primer location."""
    length = len(seq)
    bounds = _scan_bounds(length, reverse, w)
    if bounds is None:
        return
    lo, hi = bounds
    for i in range(max(lo, 0), hi + 1):
        for size in range(w.min_len, w.max_len + 1):
            if i + size > length:
                continue
            window = seq[i : i + size]
            if not is_degenerate(window):
                yield i, window


def _pcr_candidates(
    seq: str,
    reverse: bool,
    w: _Window,
    settings: PCRSettings | SequencingSettings,
    conditions: Conditions,
) -> list[Candidate]:
    found = []
    for pos, window in _windows(seq, reverse, w):
        if settings.exclude_gc:
            gc = gc_percent(window)
            if gc < settings.gc.min_gc or gc > settings.gc.max_gc:
                continue
        # 3' GC clamp: two of the last three bases G or C.
        if settings.gc_clamp and not gc_percent(window[-3:]) > 50:
            continue
        if settings.repeats.rejects(window):
            continue
        tm = melting_temperature(window, conditions)
        if settings.min_tm < tm < settings.max_tm:
            found.append(Candidate(pos, len(window), window, tm))
    return found


def _count_cpg(window: str, reverse: bool) -> int:
    # Reverse windows are read from the reversed top strand, where a CpG
    # appears as "GC". (PerlPrimer searched for "CG" here, i.e. GpC sites.)
    return window.count("GC" if reverse else "CG")


def _bisulphite_candidates(
    seq: str, reverse: bool, w: _Window, settings: BisulphiteSettings, conditions: Conditions
) -> list[Candidate]:
    found = []
    for pos, window in _windows(seq, reverse, w):
        if c_percent(window) < settings.min_c_percent:
            continue
        cpgs = _count_cpg(window, reverse)
        if cpgs > (0 if settings.exclude_cpg else 2):
            continue
        if settings.three_prime_c and not (window[-1] == "C" or c_percent(window[-3:]) > 50):
            continue
        converted = bisulphite_convert(window)
        checked = converted if settings.repeats_after_conversion else window
        if settings.repeats.rejects(checked):
            continue
        # The converted strands are no longer complementary, so reverse
        # primers are designed against the (reversed) converted top strand.
        primer = complement(converted) if reverse else converted
        tm = melting_temperature(primer, conditions)
        if settings.min_tm < tm < settings.max_tm:
            found.append(Candidate(pos, len(window), primer, tm, original=window))
    return found


def _degenerate_bisulphite(primer: str, original: str, reverse: bool) -> str:
    """Mark CpG cytosines as Y (forward) or R (reverse) to tolerate methylation."""
    out = list(primer)
    for i, base in enumerate(original):
        if reverse:
            if base == "C" and i > 0 and original[i - 1] == "G":
                out[i] = "R"
        elif base == "C" and original[i + 1 : i + 2] == "G":
            out[i] = "Y"
    return "".join(out)


# --------------------------------------------------------------------------
# Pairing
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class _PairRules:
    max_tm_diff: float
    target: TargetRange
    qpcr: QPCRSettings | None = None
    exon_bounds: tuple[int, ...] = ()


def _pair_primers(
    seq_len: int,
    forwards: list[Candidate],
    reverses: list[Candidate],
    rules: _PairRules,
    conditions: Conditions,
    bisulphite: bool = False,
) -> list[PrimerPair]:
    t = rules.target
    pairs = []
    for f in forwards:
        for r in reverses:
            reverse_end = seq_len - r.pos - 1
            amplicon = reverse_end - f.pos

            if not reverse_end - r.length > f.pos + f.length:
                continue
            if (
                t.inner_start
                and t.inner_end
                and not (f.pos < t.inner_start and t.inner_end < reverse_end)
            ):
                continue
            if not (r.tm > f.tm - rules.max_tm_diff and f.tm > r.tm - rules.max_tm_diff):
                continue
            if rules.qpcr and not _qpcr_ok(f, r, reverse_end, amplicon, rules):
                continue

            dimer = min(
                primer_dimer(f.seq, f.seq, conditions).dg,
                primer_dimer(f.seq, r.seq, conditions).dg,
                primer_dimer(r.seq, r.seq, conditions).dg,
            )
            full_dimer = min(
                primer_dimer(f.seq, f.seq, conditions, full=True).dg,
                primer_dimer(f.seq, r.seq, conditions, full=True).dg,
                primer_dimer(r.seq, r.seq, conditions, full=True).dg,
            )
            pairs.append(
                PrimerPair(
                    forward=f.seq,
                    forward_pos=f.pos,
                    forward_len=f.length,
                    forward_tm=round(f.tm, 2),
                    reverse=r.seq,
                    reverse_pos=r.pos,
                    reverse_len=r.length,
                    reverse_tm=round(r.tm, 2),
                    reverse_end=reverse_end,
                    amplicon=amplicon,
                    dimer_dg=dimer,
                    full_dimer_dg=full_dimer,
                    forward_degenerate=(
                        _degenerate_bisulphite(f.seq, f.original, False)
                        if bisulphite and f.original
                        else None
                    ),
                    reverse_degenerate=(
                        _degenerate_bisulphite(r.seq, r.original, True)
                        if bisulphite and r.original
                        else None
                    ),
                )
            )
    # Least stable dimers first: by full dimer dG, then (stably) extensible dG.
    pairs.sort(key=lambda p: p.full_dimer_dg, reverse=True)
    pairs.sort(key=lambda p: p.dimer_dg, reverse=True)
    return pairs


def _qpcr_ok(
    f: Candidate, r: Candidate, reverse_end: int, amplicon: int, rules: _PairRules
) -> bool:
    q = rules.qpcr
    assert q is not None
    if amplicon < q.min_amplicon or amplicon > q.max_amplicon:
        return False
    bounds = rules.exon_bounds
    if q.span_boundary and not any(f.pos < b < reverse_end for b in bounds):
        return False
    if q.primer_on_boundary and q.boundary_overlap:
        ov = q.boundary_overlap
        if not any(
            (f.pos < b - ov and b + ov < f.pos + f.length)
            or (reverse_end - r.length < b - ov and b + ov < reverse_end)
            for b in bounds
        ):
            return False
    return True


# --------------------------------------------------------------------------
# Target ranges
# --------------------------------------------------------------------------


def _resolve_target(
    seq: str, target: TargetRange | None, find_region: Callable[[str], list[Region]]
) -> TargetRange:
    """Fill in a target range, defaulting to a detected feature (ORF/island)."""
    if target is None or not target.is_set:
        regions = find_region(seq)
        if not regions:
            raise DesignError("No target range given and no suitable feature was found")
        start, end = regions[0].start, regions[-1].end - 1
        outer_start = target.outer_start if target and target.outer_start else start
        outer_end = target.outer_end if target and target.outer_end else end
        target = TargetRange(min(outer_start, start), start, end, max(outer_end, end))
    inner_start, inner_end = target.inner_start, target.inner_end
    if inner_start is None or inner_end is None:
        raise DesignError("Both ends of the inner target range are required")
    outer_start = inner_start if target.outer_start is None else target.outer_start
    outer_end = inner_end if target.outer_end is None else target.outer_end

    length = len(seq)
    if not 0 <= outer_start <= inner_start < inner_end <= outer_end <= length:
        raise DesignError(
            "Target range must satisfy 0 <= outer start <= inner start < inner end "
            f"<= outer end <= sequence length ({length})"
        )
    return TargetRange(outer_start, inner_start, inner_end, outer_end)


def _amplicon_window(target: TargetRange, min_len: int, max_len: int) -> _Window:
    assert target.outer_start is not None and target.outer_end is not None
    assert target.inner_start is not None and target.inner_end is not None
    return _Window(
        min_len=min_len,
        max_len=max_len,
        target=target,
        min_amplicon=target.inner_end - target.inner_start,
        max_amplicon=max(
            target.outer_end - target.outer_start, target.inner_end - target.inner_start
        ),
    )


def _step(target: TargetRange, search: Search, length: int) -> TargetRange | None:
    """Move the target one step; None when it cannot move any further."""
    assert target.outer_start is not None and target.outer_end is not None
    assert target.inner_start is not None and target.inner_end is not None
    if search is Search.INWARDS:
        start, end = target.inner_start + SEARCH_STEP, target.inner_end - SEARCH_STEP
        if start >= end:
            return None
        return replace(target, inner_start=start, inner_end=end)
    if target.outer_start == 0 and target.outer_end == length:
        return None
    return replace(
        target,
        outer_start=max(target.outer_start - SEARCH_STEP, 0),
        outer_end=min(target.outer_end + SEARCH_STEP, length),
    )


def _prepare(sequence: str) -> str:
    seq = clean_sequence(sequence).upper()
    if not seq:
        raise DesignError("No DNA sequence")
    return seq


# --------------------------------------------------------------------------
# Public design functions
# --------------------------------------------------------------------------


def design_pcr(
    sequence: str,
    target: TargetRange | None = None,
    settings: PCRSettings | None = None,
    conditions: Conditions = DEFAULT_CONDITIONS,
    search: Search = Search.FIXED,
) -> DesignResult:
    """Design primer pairs for standard PCR.

    If no inner range is given, the longest ORF is used as the target.
    """
    settings = settings or PCRSettings()
    seq = _prepare(sequence)
    target = _resolve_target(seq, target, _orf_regions)

    def run(t: TargetRange) -> list[PrimerPair]:
        w = _amplicon_window(t, settings.min_length, settings.max_length)
        forwards = _pcr_candidates(seq, False, w, settings, conditions)
        reverses = _pcr_candidates(reverse_complement(seq), True, w, settings, conditions)
        rules = _PairRules(settings.max_tm_diff, t)
        return _pair_primers(len(seq), forwards, reverses, rules, conditions)

    return _search(seq, target, search, run)


def design_bisulphite(
    sequence: str,
    target: TargetRange | None = None,
    settings: BisulphiteSettings | None = None,
    conditions: Conditions = DEFAULT_CONDITIONS,
    search: Search = Search.FIXED,
    cpg: Callable[[str], list[Region]] = find_cpg_islands,
) -> DesignResult:
    """Design primer pairs for bisulphite-converted DNA.

    Primers are designed against the converted top strand. If no inner range
    is given, the predicted CpG island(s) are used as the target.
    """
    settings = settings or BisulphiteSettings()
    seq = _prepare(sequence)
    target = _resolve_target(seq, target, cpg)

    def run(t: TargetRange) -> list[PrimerPair]:
        w = _amplicon_window(t, settings.min_length, settings.max_length)
        forwards = _bisulphite_candidates(seq, False, w, settings, conditions)
        reverses = _bisulphite_candidates(seq[::-1], True, w, settings, conditions)
        rules = _PairRules(settings.max_tm_diff, t)
        return _pair_primers(len(seq), forwards, reverses, rules, conditions, bisulphite=True)

    return _search(seq, target, search, run)


def design_qpcr(
    mrna: str,
    exon_boundaries: list[int],
    settings: QPCRSettings | None = None,
    conditions: Conditions = DEFAULT_CONDITIONS,
) -> DesignResult:
    """Design real-time PCR primer pairs on an mRNA/cDNA sequence.

    Args:
        mrna: The spliced transcript sequence.
        exon_boundaries: Positions in ``mrna`` where each exon after the
            first begins (see :mod:`pyprimer.exons`).
    """
    settings = settings or QPCRSettings()
    seq = _prepare(mrna)
    bounds = tuple(sorted(exon_boundaries))
    if any(not 0 < b < len(seq) for b in bounds):
        raise DesignError("Exon boundaries must lie within the mRNA sequence")
    if not bounds and (settings.span_boundary or settings.primer_on_boundary):
        raise DesignError("Exon boundaries are required to design real-time PCR primers")

    w = _Window(
        min_len=settings.min_length,
        max_len=settings.max_length,
        min_amplicon=settings.min_amplicon,
        max_amplicon=settings.max_amplicon,
        exon_bounds=bounds,
        first_exon=settings.first_exon,
        last_exon=settings.last_exon,
    )
    forwards = _pcr_candidates(seq, False, w, settings, conditions)
    reverses = _pcr_candidates(reverse_complement(seq), True, w, settings, conditions)
    rules = _PairRules(settings.max_tm_diff, TargetRange(), settings, bounds)
    pairs = _pair_primers(len(seq), forwards, reverses, rules, conditions)
    return DesignResult(pairs, TargetRange(), seq, _found_message(pairs))


def design_sequencing(
    sequence: str,
    start: int | None = None,
    end: int | None = None,
    settings: SequencingSettings | None = None,
    conditions: Conditions = DEFAULT_CONDITIONS,
) -> SequencingResult:
    """Choose a set of evenly spaced forward sequencing primers.

    Starting from ``start``, each primer is the first acceptable candidate
    between ``min_spacing`` and ``max_spacing`` bases after the previous one.
    """
    settings = settings or SequencingSettings()
    seq = _prepare(sequence)
    if not end:
        start, end = 0, len(seq)
    start = start or 0
    end = min(end, len(seq))
    if start >= end:
        raise DesignError("Sequencing range start must be before its end")
    target = TargetRange(inner_start=start, inner_end=end)

    w = _Window(settings.min_length, settings.max_length, target)
    candidates = _pcr_candidates(seq, False, w, settings, conditions)

    messages = []
    primers: list[SequencingPrimer] = []
    last: int | None = None
    for c in candidates:
        if last is None:
            if c.pos < start:
                continue
            if c.pos > start + settings.max_spacing:
                messages.append("Contiguous primers could not be found within the set parameters")
                break
        else:
            if c.pos > end or c.pos < last + settings.min_spacing:
                continue
            if c.pos > last + settings.max_spacing:
                messages.append("Contiguous primers could not be found within the set parameters")
                break
        dg = primer_dimer(c.seq, c.seq, conditions, full=True).dg
        if settings.exclude_dimers and dg < -settings.max_dimer_dg:
            continue
        primers.append(SequencingPrimer(c.seq, c.pos, c.length, round(c.tm, 2), dg))
        last = c.pos

    messages.append(f"Found {len(primers)} primers")
    return SequencingResult(primers, target, seq, messages)


def _orf_regions(seq: str) -> list[Region]:
    orf = find_orf(seq)
    return [orf] if orf else []


def _found_message(pairs: list[PrimerPair]) -> list[str]:
    return [f"Found {len(pairs)} primer pairs"]


def _search(
    seq: str,
    target: TargetRange,
    search: Search,
    run: Callable[[TargetRange], list[PrimerPair]],
) -> DesignResult:
    pairs = run(target)
    messages = []
    while not pairs and search is not Search.FIXED:
        moved = _step(target, search, len(seq))
        if moved is None:
            messages.append("No primers found: the range cannot be adjusted any further")
            break
        target = moved
        pairs = run(target)
    return DesignResult(pairs, target, seq, messages + _found_message(pairs))
