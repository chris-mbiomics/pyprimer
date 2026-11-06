"""Detection of sequence features: open reading frames and CpG islands."""

from __future__ import annotations

import re
from dataclasses import dataclass

from .sequence import GENETIC_CODE, gc_percent

_ALT_START = re.compile(r"[ATGC]TG")
_MARKED_REGION = re.compile(r"([a-z]*)([A-Z]*)")


@dataclass(frozen=True)
class Region:
    """A half-open region ``[start, end)`` of a sequence (0-based)."""

    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def find_orf(seq: str) -> Region | None:
    """Return the longest open reading frame, or None.

    Translation may initiate from ATG or, failing that, from the rare
    alternative CTG, GTG and TTG start codons; a downstream ATG is preferred
    once found. The returned region ends after the stop codon.
    """
    seq = seq.upper()
    best: list[tuple[int, int, int]] = []  # (codons, start, end) per frame
    for frame in range(3):
        frame_best = (0, -1, -1)
        started = met_started = False
        count = 0
        start = -1
        pos = frame
        while pos < len(seq) - 3:
            codon = seq[pos : pos + 3]
            if GENETIC_CODE.get(codon, "X") != "*":
                if not started:
                    if not _ALT_START.search(codon):
                        pos += 3
                        continue
                    started = True
                    start = pos
                if not met_started and codon == "ATG":
                    met_started = True
                    start = pos
                    count = 0
                count += 1
                pos += 3
                continue
            # stop codon
            if count > frame_best[0]:
                frame_best = (count + 1, start, pos + 3)
            count = 0
            started = met_started = False
            start = -1
            pos += 3
        if count > frame_best[0]:
            frame_best = (count + 1, start, pos + 3)
        best.append(frame_best)

    codons, start, end = sorted(best, key=lambda orf: orf[0], reverse=True)[0]
    if codons <= 1:
        return None
    return Region(start, end)


@dataclass(frozen=True)
class CpGSettings:
    """Parameters for CpG island prediction (Gardiner-Garden & Frommer, 1987).

    Attributes:
        window: Sliding window size (bp).
        min_island: Minimum island length (bp).
        min_gc: Minimum average %GC.
        min_obs_exp: Minimum average observed/expected CpG ratio.
        cpgplot: Emulate EMBOSS cpgplot's (slightly inflated) arithmetic.
    """

    window: int = 200
    min_island: int = 200
    min_gc: float = 50.0
    min_obs_exp: float = 0.6
    cpgplot: bool = False


def find_cpg_islands(seq: str, settings: CpGSettings | None = None) -> list[Region]:
    """Predict CpG islands.

    Window values are averaged over ten consecutive windows. With
    ``settings.cpgplot`` the results reproduce cpgplot, which averages
    eleven windows (dividing by ten) and counts CpGs over window+1 bases.
    """
    settings = settings or CpGSettings()
    seq = seq.lower()
    seq_len = len(seq)
    window = settings.window
    span = 10 if settings.cpgplot else 9

    percent_gc: list[float] = []
    obs_exp: list[float] = []
    for i in range(seq_len - window + 1):
        chunk = seq[i : i + window]
        expected = chunk.count("c") * chunk.count("g") / (len(chunk) - 1)
        if settings.cpgplot:
            chunk = seq[i : i + window + 1]
        observed = chunk.count("cg")
        percent_gc.append(gc_percent(chunk))
        obs_exp.append(observed / expected if expected else 0.0)

    islands: list[Region] = []
    in_island = False
    island_start = 0
    real_pos = 0
    for i in range(seq_len - window - 10 + 1):
        # The reported position is the middle of the averaged windows.
        real_pos = int(i + window / 2)
        if not settings.cpgplot:
            real_pos += span // 2

        av_gc = sum(percent_gc[i : i + span + 1]) / 10
        av_obs_exp = sum(obs_exp[i : i + span + 1]) / 10

        if av_gc > settings.min_gc and av_obs_exp > settings.min_obs_exp:
            if not in_island:
                island_start = real_pos
            in_island = True
            continue
        if in_island:
            in_island = False
            if real_pos - island_start >= settings.min_island:
                end = real_pos - 1 if settings.cpgplot else real_pos + 1
                islands.append(Region(island_start, end))

    # PerlPrimer silently dropped an island still open at the end of the scan.
    if in_island and real_pos + 1 - island_start >= settings.min_island:
        end = real_pos if settings.cpgplot else min(real_pos + 2, seq_len)
        islands.append(Region(island_start, end))
    return islands


def find_marked_regions(seq: str) -> list[Region]:
    """Return upper-case regions embedded in lower-case sequence.

    PerlPrimer optionally lets users mark genes or islands by case: e.g.
    ``acgtATGAAATAGacgt`` marks ``[4, 13)``. Upper-case runs at the very
    start of the sequence are ignored, as in PerlPrimer.
    """
    regions = []
    pos = 0
    for match in _MARKED_REGION.finditer(seq):
        lower, upper = match.group(1), match.group(2)
        if not lower and not upper:
            continue
        if upper and lower:
            regions.append(Region(pos + len(lower), pos + len(lower) + len(upper)))
        pos += len(lower) + len(upper)
    return regions


def highlight_regions(seq: str, regions: list[Region]) -> str:
    """Return the sequence in lower case with the given regions in upper case."""
    out = seq.lower()
    for region in regions:
        out = out[: region.start] + out[region.start : region.end].upper() + out[region.end :]
    return out
