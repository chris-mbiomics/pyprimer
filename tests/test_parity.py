"""Compare the Python port against the original PerlPrimer implementation.

These tests run only when the original source is available: set
``PERLPRIMER_PL`` to the path of ``perlprimer.pl`` (Perl >= 5.14 required).
"""

from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from pyprimer import design as design_module
from pyprimer.design import (
    BisulphiteSettings,
    PCRSettings,
    QPCRSettings,
    SequencingSettings,
    TargetRange,
    design_bisulphite,
    design_pcr,
    design_qpcr,
    design_sequencing,
)
from pyprimer.dimer import primer_dimer
from pyprimer.features import find_cpg_islands, find_orf
from pyprimer.thermo import thermodynamics

from .helpers import cpg_rich, orf_sequence, random_dna

PERL_SOURCE = os.environ.get("PERLPRIMER_PL")
ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(
    not PERL_SOURCE or not shutil.which("perl"),
    reason="set PERLPRIMER_PL to the original perlprimer.pl to run parity tests",
)


class Perl:
    def __init__(self, harness: Path) -> None:
        self.proc = subprocess.Popen(
            ["perl", str(harness)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )

    def __call__(self, **request: Any) -> Any:
        assert self.proc.stdin and self.proc.stdout
        self.proc.stdin.write(json.dumps(request) + "\n")
        self.proc.stdin.flush()
        return json.loads(self.proc.stdout.readline())["result"]

    def close(self) -> None:
        if self.proc.stdin:
            self.proc.stdin.close()
        self.proc.wait(timeout=10)


@pytest.fixture(scope="module")
def perl(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Perl]:
    harness = tmp_path_factory.mktemp("perl") / "harness.pl"
    code = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "perl_harness.py"), str(PERL_SOURCE)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    harness.write_text(code)
    runner = Perl(harness)
    yield runner
    runner.close()


def _pairs(pairs: list[design_module.PrimerPair]) -> list[list[Any]]:
    return [
        [
            p.forward,
            p.forward_pos,
            p.forward_len,
            p.forward_tm,
            p.reverse,
            p.reverse_pos,
            p.reverse_len,
            p.reverse_tm,
            p.reverse_end,
            p.amplicon,
            p.dimer_dg,
            p.full_dimer_dg,
        ]
        for p in pairs
    ]


def _normalise(rows: list[list[Any]]) -> list[list[Any]]:
    return [
        [float(v) if isinstance(v, (int, float, str)) and _num(v) else v for v in row]
        for row in rows
    ]


def _num(v: Any) -> bool:
    try:
        float(v)
    except ValueError:
        return False
    return True


def test_melting_temperature(perl: Perl) -> None:
    rng = random.Random(1)
    for _ in range(200):
        primer = random_dna(rng, rng.randint(12, 35))
        tm, dh, ds = perl(op="tm", seq=primer)
        thermo = thermodynamics(primer)
        assert thermo.tm == pytest.approx(tm, abs=1e-9)
        assert thermo.dh == pytest.approx(dh, abs=1e-9)
        assert thermo.ds == pytest.approx(ds, abs=1e-9)


@pytest.mark.parametrize("full", [False, True])
def test_primer_dimer(perl: Perl, full: bool) -> None:
    rng = random.Random(2)
    for _ in range(300):
        f = random_dna(rng, rng.randint(15, 30))
        # Self-dimers, related primers and unrelated primers.
        r = rng.choice([f, random_dna(rng, rng.randint(15, 30)), f[::-1][:20] + "ACGT"])
        expected = perl(op="dimer", f=f, r=r, full=int(full))
        dimer = primer_dimer(f, r, full=full)
        assert dimer.dg == pytest.approx(expected["dg"], abs=1e-9), (f, r)
        if dimer.found:
            assert dimer.diagram() + "\n" == expected["diagram"], (f, r)


def test_find_orf(perl: Perl) -> None:
    rng = random.Random(3)
    for _ in range(40):
        seq = orf_sequence(rng, rng.randint(0, 200), rng.randint(10, 200), rng.randint(0, 200))
        expected = perl(op="orf", seq=seq)
        orf = find_orf(seq)
        assert ([[orf.start, orf.end]] if orf else []) == expected


def test_find_cpg_islands(perl: Perl) -> None:
    rng = random.Random(4)
    for _ in range(10):
        seq = (
            random_dna(rng, rng.randint(200, 600), gc=0.35)
            + cpg_rich(rng, rng.randint(150, 600))
            + random_dna(rng, rng.randint(400, 600), gc=0.35)
        )
        assert [[r.start, r.end] for r in find_cpg_islands(seq)] == perl(op="cpg", seq=seq)


def _range_vars(t: TargetRange) -> dict[str, int | None]:
    return {
        "max_range_5p": t.outer_start,
        "min_range": t.inner_start,
        "max_range": t.inner_end,
        "max_range_3p": t.outer_end,
    }


@pytest.mark.parametrize(
    "target",
    [
        TargetRange(260, 300, 700, 760),
        TargetRange(200, 300, 700, 700),
        TargetRange(300, 300, 700, 700),
        TargetRange(0, 50, 900, 1000),
    ],
)
def test_design_pcr(perl: Perl, target: TargetRange) -> None:
    rng = random.Random(5)
    seq = orf_sequence(rng, 250, 200, 150)
    settings = PCRSettings()
    expected = perl(
        op="pcr",
        seq=seq,
        vars={
            "min_tm": 57,
            "max_tm": 63,
            "max_diff": 3,
            "pri_win_win": 20,
            "pri_win_max": 24,
            **_range_vars(target),
        },
    )
    result = design_pcr(seq, target, settings)
    assert _normalise(_pairs(result.pairs)) == _normalise(expected)


def test_design_bisulphite(perl: Perl, monkeypatch: pytest.MonkeyPatch) -> None:
    # PerlPrimer counted GpC instead of CpG in reverse primers; use its
    # behaviour here so that everything else can be compared exactly.
    monkeypatch.setattr(design_module, "_count_cpg", lambda window, reverse: window.count("CG"))
    rng = random.Random(6)
    # C-rich, CpG-poor flanks give primers specific for the converted strand.
    flank = "".join(rng.choices("ACGT", weights=[0.15, 0.4, 0.3, 0.15], k=600)).replace("CG", "CA")
    seq = flank[:300] + cpg_rich(rng, 300) + flank[300:]
    target = TargetRange(150, 300, 600, 750)
    expected = perl(
        op="bisulphite",
        seq=seq,
        vars={
            "min_tm": 55,
            "max_tm": 65,
            "max_diff": 5,
            "pri_win_win": 25,
            "pri_win_max": 30,
            **_range_vars(target),
        },
    )
    result = design_bisulphite(seq, target, BisulphiteSettings())
    assert expected, "test sequence should yield primers"
    assert _normalise(_pairs(result.pairs)) == _normalise(expected)


def test_design_qpcr(perl: Perl) -> None:
    rng = random.Random(7)
    seq = random_dna(rng, 900)
    bounds = [180, 330, 520, 700]
    expected = perl(
        op="qpcr",
        seq=seq,
        bounds=bounds,
        vars={
            "min_tm": 58,
            "max_tm": 62,
            "max_diff": 2,
            "pri_win_win": 20,
            "pri_win_max": 24,
            "min_ampsize": 100,
            "max_ampsize": 300,
        },
    )
    result = design_qpcr(seq, bounds, QPCRSettings())
    assert expected, "test sequence should yield primers"
    assert _normalise(_pairs(result.pairs)) == _normalise(expected)


def test_design_sequencing(perl: Perl) -> None:
    rng = random.Random(8)
    seq = random_dna(rng, 4000)
    expected = perl(
        op="sequencing",
        seq=seq,
        vars={
            "min_tm": 56,
            "max_tm": 64,
            "pri_win_win": 20,
            "pri_win_max": 24,
            "min_range": 0,
            "max_range": len(seq),
        },
    )
    result = design_sequencing(seq, settings=SequencingSettings())
    got = [[p.seq, p.pos, p.length, p.tm, p.dimer_dg] for p in result.primers]
    assert len(got) > 3
    assert _normalise(got) == _normalise(expected)
