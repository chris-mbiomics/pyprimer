"""Unit tests for sequence utilities, thermodynamics and primer-dimers."""

from __future__ import annotations

import pytest

from pyprimer.dimer import primer_dimer
from pyprimer.sequence import (
    bisulphite_convert,
    c_percent,
    clean_sequence,
    complement,
    gc_percent,
    is_degenerate,
    parse_sequence,
    reverse_complement,
    translate,
)
from pyprimer.thermo import Conditions, free_energy_table, melting_temperature, thermodynamics


class TestSequence:
    def test_parse_fasta_uses_first_record(self) -> None:
        name, seq = parse_sequence(">gene one\nACGT acgt\n12 AC-GT\n>second\nTTTT\n")
        assert name == "gene one"
        assert seq == "ACGTacgtACGT"

    def test_parse_raw(self) -> None:
        assert parse_sequence("  ac gt\n") == (None, "acgt")
        assert clean_sequence("1 acgt 5\n") == "acgt"

    def test_complements(self) -> None:
        assert complement("AACGTn") == "TTGCAn"
        assert reverse_complement("AACGt") == "aCGTT"

    def test_bisulphite_converts_upper_case_c(self) -> None:
        assert bisulphite_convert("ACGCc") == "ATGTc"

    def test_composition(self) -> None:
        assert gc_percent("GCat") == 50
        assert c_percent("CCcA") == 50
        assert is_degenerate("ACGN")
        assert not is_degenerate("acgt")

    def test_translate(self) -> None:
        assert translate("ATGTAAxyz", 0, 9) == "M*X"


class TestThermodynamics:
    def test_matches_perlprimer(self) -> None:
        # Reference values from PerlPrimer 1.1.21 with its default conditions.
        t = thermodynamics("ACGTACGTAGCTAGCTAGGA")
        assert t.tm == pytest.approx(59.7261177115384, abs=1e-9)
        assert t.dh == pytest.approx(-153.8)
        assert t.ds == pytest.approx(-428.629807634741, abs=1e-9)

    def test_case_insensitive(self) -> None:
        assert melting_temperature("acgtacgtagct") == melting_temperature("ACGTACGTAGCT")

    def test_more_primer_raises_tm(self) -> None:
        low = melting_temperature("ACGTACGTAGCTAGCTAGGA", Conditions(oligo_nm=50))
        high = melting_temperature("ACGTACGTAGCTAGCTAGGA", Conditions(oligo_nm=1000))
        assert high > low

    def test_magnesium_raises_tm(self) -> None:
        low = melting_temperature("ACGTACGTAGCTAGCTAGGA", Conditions(mg_mm=0))
        high = melting_temperature("ACGTACGTAGCTAGCTAGGA", Conditions(mg_mm=3))
        assert high > low

    def test_excess_dntps_give_no_free_magnesium(self) -> None:
        assert Conditions(mg_mm=0.5, dntp_mm=1).sodium_equivalent == pytest.approx(0.05)

    @pytest.mark.parametrize("bad", ["", "ACGN", "AC GT"])
    def test_rejects_invalid_primers(self, bad: str) -> None:
        with pytest.raises(ValueError):
            thermodynamics(bad)

    def test_rejects_impossible_conditions(self) -> None:
        with pytest.raises(ValueError):
            Conditions(oligo_nm=0)
        with pytest.raises(ValueError):
            Conditions(monovalent_mm=0, mg_mm=0)

    def test_free_energy_table_is_read_only(self) -> None:
        table = free_energy_table()
        assert table["initA"] == 1.03
        with pytest.raises(TypeError):
            table["AATT"] = 0  # type: ignore[index]


class TestPrimerDimer:
    def test_matches_perlprimer(self) -> None:
        # Reference output from PerlPrimer 1.1.21.
        dimer = primer_dimer("GAGAATGTGGGCTCTCAAGC", "GGCGTTAGTTTAAGTCACGG", full=True)
        assert dimer.dg == -0.57
        assert dimer.diagram() == (
            "5' GAGAATGTGGGCTCTCAAGC 3'\n       ..|||...........\n    3' GGCACTGAATTTGATTGCGG 5'\n"
        )

    def test_palindromic_3_prime_end_forms_extensible_dimer(self) -> None:
        dimer = primer_dimer("TTTTTTTTTTTTGAATTC", "TTTTTTTTTTTTGAATTC")
        assert dimer.found
        assert dimer.binding == "111111"  # the six-base GAATTC palindrome pairs
        assert dimer.dg < -2

    def test_no_dimer(self) -> None:
        dimer = primer_dimer("AAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCC")
        assert dimer.dg == 0
        assert not dimer.found
        assert dimer.diagram() == ""

    def test_full_is_at_least_as_stable_as_extensible(self) -> None:
        a, b = "ACGTAGCTAGCATGCAAGTC", "GGATCCATGCTAGCTTTAGC"
        assert primer_dimer(a, b, full=True).dg <= primer_dimer(a, b).dg

    def test_rejects_degenerate(self) -> None:
        with pytest.raises(ValueError):
            primer_dimer("ACGTN", "ACGT")
