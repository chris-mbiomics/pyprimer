"""Tests for the AutoDimer port (NIST reference results: test_autodimer_nist.py)."""

from __future__ import annotations

import pytest

from pyprimer.autodimer import (
    AutoDimerSettings,
    Oligo,
    dimer_thermo,
    format_value,
    hairpin_thermo,
    normalise,
    parse_oligos,
    report_text,
    screen_dimers,
    screen_hairpins,
)

PAPER = "GTCAGCTAGTAGCTAGCGCA"  # Figure 3 of Vallone & Butler (2004)
EVERYTHING = AutoDimerSettings(min_score=-100)


class TestInput:
    def test_fasta(self) -> None:
        text = ">F1 forward\nGTCAG CTAG\ntagc\n\n>R1\nACGU\n"
        assert parse_oligos(text) == [Oligo("F1 forward", "GTCAGCTAGTAGC"), Oligo("R1", "ACGT")]

    def test_one_per_line(self) -> None:
        assert parse_oligos("ACGT\n\nGGCC\n") == [Oligo("Seq 1", "ACGT"), Oligo("Seq 2", "GGCC")]

    def test_unknown_bases(self) -> None:
        assert normalise("acgRYn-u") == "ACGNNNNT"

    def test_limit(self) -> None:
        with pytest.raises(ValueError, match="Too many"):
            parse_oligos("ACGT\n" * 101)

    @pytest.mark.parametrize("kwargs", [{"sodium_m": 0}, {"strand_um": -1}, {"dg_temp_c": -300}])
    def test_invalid_settings(self, kwargs: dict[str, float]) -> None:
        with pytest.raises(ValueError):
            AutoDimerSettings(**kwargs)  # type: ignore[arg-type]


class TestDimers:
    def test_paper_example(self) -> None:
        # The paper's self-dimer: 12 matches and 5 mismatches give a score of 7.
        oligos = [Oligo("a", PAPER), Oligo("b", PAPER)]
        hits = screen_dimers(oligos, AutoDimerSettings(min_score=7)).hits
        assert {(h.matches, h.score) for h in hits} == {(12, 7)}

    def test_alignments_per_pair(self) -> None:
        # Equal lengths: every overlap from each end; NIST adds a few extra
        # shifts when the longer oligo is 4+ bases longer, capped here.
        def count(n: int, m: int) -> int:
            screen = screen_dimers([Oligo("a", "A" * n), Oligo("b", "C" * m)], EVERYTHING)
            return len(screen.hits)

        assert count(20, 20) == 40
        assert count(20, 18) == 36
        assert count(15, 9) == 21
        assert count(18, 8) == 23
        assert count(40, 8) == 23

    def test_pairs_and_order(self) -> None:
        oligos = [Oligo(f"S{i}", "GCGCGCGC") for i in range(1, 4)]
        screen = screen_dimers(oligos, AutoDimerSettings(min_score=8))
        assert screen.comparisons == 3
        assert [(h.label_1, h.label_2) for h in screen.hits] == [
            ("S3", "S1"), ("S3", "S1"), ("S3", "S2"), ("S3", "S2"), ("S2", "S1"), ("S2", "S1"),
        ]  # fmt: skip

    def test_no_self_dimers(self) -> None:
        assert screen_dimers([Oligo("a", "GAATTC"), Oligo("b", "AAAAAA")]).hits == []

    def test_threshold_inclusive(self) -> None:
        oligos = [Oligo("a", PAPER), Oligo("b", "TGCGCTAGCTACTAGCTGAC")]
        scores = {h.score for h in screen_dimers(oligos, AutoDimerSettings(min_score=20)).hits}
        assert scores == {20}

    def test_lineup(self) -> None:
        hit = screen_dimers([Oligo("a", "GGCGAAATA"), Oligo("b", "TATTGCGCC")], EVERYTHING).hits[16]
        assert hit.lineup == " 5'-GGCGAAATA-3'\n    ||||x||||\n 3'-CCGCGTTAT-5'"
        assert hit.match_sequence == "GGCGNAATA"
        assert format_value(hit.tm, "< 0") == "4.4" and format_value(hit.dg, "> 0") == "-3.7"


class TestThermodynamics:
    def test_perfect_duplex(self) -> None:
        # Values from the NIST web service.
        tm, dg = dimer_thermo("TTCCTCATGCAATTC", [True] * 15)
        assert (format_value(tm, ""), format_value(dg, "")) == ("45.2", "-12.2")

    def test_self_complementary(self) -> None:
        tm, dg = dimer_thermo("CGCGAATTCGCG", [True] * 12)
        assert (format_value(tm, ""), format_value(dg, "")) == ("51.1", "-12.9")

    def test_unstable(self) -> None:
        assert dimer_thermo("AT", [True, True]) == (None, None)
        assert dimer_thermo("AC", [False, False]) == (None, None)

    def test_hairpin_loop_penalty(self) -> None:
        stem = [True] * 12
        _, dg_loop = hairpin_thermo("GGCCGCAGCGTG", stem)
        _, dg_duplex = dimer_thermo("GGCCGCAGCGTG", stem)
        assert dg_loop is not None and dg_duplex is not None
        assert dg_loop - dg_duplex == pytest.approx(5.6, abs=0.01)

    def test_rounding_half_up(self) -> None:
        assert format_value(0.25, "") == "0.3"
        assert format_value(-0.25, "") == "-0.3"
        assert format_value(None, "> 0") == "> 0"


class TestHairpins:
    def test_folds(self) -> None:
        # Paper figure 2: a 20-mer is folded with 4 and 5 base loops.
        screen = screen_hairpins([Oligo("p", PAPER)], EVERYTHING)
        assert screen.comparisons == 1
        assert len(screen.hits) == 22
        assert {h.loop for h in screen.hits} == {4, 5}
        assert min(h.stem for h in screen.hits) == 2

    def test_paper_hairpin(self) -> None:
        hits = screen_hairpins([Oligo("p", PAPER)], AutoDimerSettings(min_score=3)).hits
        assert [(h.loop, h.stem, h.matches, h.score, h.match_sequence) for h in hits] == [
            (5, 7, 5, 3, "TAGCNNA")
        ]
        assert hits[0].lineup == " TGATCGACTG-5'\nA  ||||XX|\n GCTAGCGCA-3'"

    def test_too_short(self) -> None:
        assert screen_hairpins([Oligo("s", "ACGTACGTA")], EVERYTHING).hits == []


def test_report_text() -> None:
    oligos = [Oligo("F", PAPER), Oligo("R", "TGCGCTAGCTACTAGCTGAC")]
    text = report_text(screen_dimers(oligos, AutoDimerSettings(min_score=20)))
    assert "Total Number of Primer-Primer Comparisons = 1" in text
    assert "R versus F" in text and "Score = 20" in text
    empty = report_text(screen_hairpins(oligos, AutoDimerSettings(min_score=30)))
    assert "NO HAIRPIN HITS" in empty
