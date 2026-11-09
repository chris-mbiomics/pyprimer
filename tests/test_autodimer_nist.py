"""Replay AutoDimer results recorded from the NIST DNA Analysis web service.

The fixture is produced by ``tools/fetch_autodimer_nist.py``; no network
access is needed here. Every hit must match NIST exactly: the oligos, the
alignment drawing, matches, score and, for primer-dimers, Tm and dG as NIST
formats them.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

from pyprimer.autodimer import (
    AutoDimerSettings,
    Oligo,
    format_value,
    normalise,
    screen_dimers,
    screen_hairpins,
)

FIXTURE = Path(__file__).parent / "data" / "autodimer_nist.json"
CASES: list[dict[str, Any]] = json.loads(FIXTURE.read_text())["cases"]

#: Behaviour of the NIST service that is deliberately not reproduced.
DEVIATIONS = {
    # NIST never pairs lower-case bases; AutoDimer (and pyprimer) read them as upper case.
    "lowercase": "NIST ignores lower-case bases",
    # With R/Y in the 5' part, NIST sizes some hairpin stems past the 3' end.
    "ambiguous_hairpin": "NIST mis-sizes hairpin stems around ambiguity codes",
}

# NIST draws some oligos with their original U and ambiguity codes; pyprimer
# shows them as read (U as T, anything else as N).
_AMBIGUOUS = re.compile(r"[RYKMSWBDHV]")


def _lineup(text: str) -> str:
    return _AMBIGUOUS.sub("N", text.replace("U", "T"))


def _settings(params: dict[str, str]) -> AutoDimerSettings:
    return AutoDimerSettings(
        sodium_m=float(params["monovalentSalt"]),
        strand_um=float(params["strandConcentration"]),
        dg_temp_c=float(params["deltaGTemperature"]),
        min_score=int(params["requiredScore"]),
    )


def _expected(case: dict[str, Any]) -> list[tuple[Any, ...]]:
    rows = []
    for hit in case["response"]:
        if case["type"] == "hairpin":
            labels: tuple[str, ...] = (hit["sequenceLabel"],)
            # Hairpin Tm and dG follow the AutoDimer paper instead (see README).
            values: tuple[str, ...] = ()
        else:
            labels = (hit["sequenceOneLabel"], hit["sequenceTwoLabel"])
            values = (hit["estimatedTemperature"], hit["freeEnergy"])
        common = (hit["matches"], hit["score"], hit["matchSequence"], _lineup(hit["lineup"]))
        rows.append((*labels, *common, *values))
    return rows


def _actual(case: dict[str, Any]) -> list[tuple[Any, ...]]:
    oligos = [Oligo(s["label"], normalise(s["sequence"])) for s in case["sequences"]]
    settings = _settings(case["parameters"])
    if case["type"] == "hairpin":
        return [
            (h.label, h.matches, h.score, h.match_sequence, h.lineup)
            for h in screen_hairpins(oligos, settings).hits
        ]
    return [
        (
            d.label_1,
            d.label_2,
            d.matches,
            d.score,
            d.match_sequence,
            d.lineup,
            format_value(d.tm, "< 0"),
            format_value(d.dg, "> 0"),
        )
        for d in screen_dimers(oligos, settings).hits
    ]


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(c, marks=pytest.mark.xfail(reason=DEVIATIONS[c["name"]], strict=True))
        if c["name"] in DEVIATIONS
        else c
        for c in CASES
    ],
    ids=[c["name"] for c in CASES],
)
def test_matches_nist(case: dict[str, Any]) -> None:
    assert _actual(case) == _expected(case)


def test_fixture_has_numeric_values() -> None:
    """Guard against a fixture where every Tm/dG is "< 0" or "> 0"."""
    values = [
        h["estimatedTemperature"]
        for c in CASES
        if c["type"] == "autodimer"
        for h in c["response"]
        if h["estimatedTemperature"] != "< 0"
    ]
    assert len(values) > 20
