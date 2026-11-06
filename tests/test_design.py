"""Tests for primer design, features, cloning, exon mapping and reports."""

from __future__ import annotations

import itertools
import random

import pytest

from pyprimer.cloning import add_cloning_tails, cloning_enzymes, find_enzyme, rebase_enzymes
from pyprimer.design import (
    DesignError,
    PCRSettings,
    QPCRSettings,
    RepeatFilter,
    Search,
    SequencingSettings,
    TargetRange,
    _degenerate_bisulphite,
    design_bisulphite,
    design_pcr,
    design_qpcr,
    design_sequencing,
)
from pyprimer.exons import AlignmentError, map_exons
from pyprimer.features import (
    CpGSettings,
    Region,
    find_cpg_islands,
    find_marked_regions,
    find_orf,
    highlight_regions,
)
from pyprimer.report import MapPrimers, analyse_pair, sequence_map, text_report
from pyprimer.sequence import reverse_complement
from pyprimer.thermo import melting_temperature

from .helpers import cpg_rich, orf_sequence, random_dna


@pytest.fixture(scope="module")
def orf_seq() -> str:
    return orf_sequence(random.Random(11), 250, 200, 250)


class TestFeatures:
    def test_find_orf(self) -> None:
        seq = "ccc" + "ATG" + "GCT" * 20 + "TAA" + "ccc"
        assert find_orf(seq) == Region(3, 69)

    def test_orf_prefers_atg_over_alternative_start(self) -> None:
        # GCC repeats contain no start codons in the other two frames.
        seq = "CTG" + "GCC" * 5 + "ATG" + "GCC" * 10 + "TGA" + "ccc"
        orf = find_orf(seq)
        assert orf is not None and orf.start == 18

    def test_no_orf(self) -> None:
        assert find_orf("TAATAATAATAA") is None

    def test_cpg_island(self) -> None:
        rng = random.Random(4)
        seq = random_dna(rng, 500, gc=0.3) + cpg_rich(rng, 500) + random_dna(rng, 500, gc=0.3)
        islands = find_cpg_islands(seq)
        assert len(islands) == 1
        assert 450 < islands[0].start < 600 and 900 < islands[0].end < 1100

    def test_island_reaching_sequence_end_is_reported(self) -> None:
        rng = random.Random(5)
        seq = random_dna(rng, 400, gc=0.3) + cpg_rich(rng, 600)
        islands = find_cpg_islands(seq)
        assert islands and islands[-1].end <= len(seq)

    def test_cpgplot_mode_runs(self) -> None:
        rng = random.Random(4)
        seq = random_dna(rng, 500, gc=0.3) + cpg_rich(rng, 500) + random_dna(rng, 500, gc=0.3)
        assert find_cpg_islands(seq, CpGSettings(cpgplot=True))

    def test_marked_regions(self) -> None:
        assert find_marked_regions("acgtATGAAATAGacgtGGGcc") == [Region(4, 13), Region(17, 20)]
        assert highlight_regions("ACGTACGT", [Region(2, 4)]) == "acGTacgt"


class TestRepeatFilter:
    def test_runs_and_repeats(self) -> None:
        f = RepeatFilter()
        assert f.run_found("ACGAAAAT")
        assert not f.run_found("ACGAAAT")
        assert f.repeat_found("GCACACACGT") == "CA"
        assert f.rejects("ACACACAC")
        assert not RepeatFilter(enabled=False).rejects("AAAAAAA")

    def test_validation(self) -> None:
        with pytest.raises(ValueError):
            RepeatFilter(repeat=1)


class TestPCR:
    def test_designs_pairs_around_orf(self, orf_seq: str) -> None:
        result = design_pcr(orf_seq, search=Search.OUTWARDS)
        orf = find_orf(orf_seq)
        assert orf and result.pairs
        assert result.target.inner_start == orf.start
        for pair in result.pairs:
            assert pair.forward == orf_seq[pair.forward_pos : pair.forward_pos + pair.forward_len]
            rev_start = pair.reverse_end - pair.reverse_len + 1
            assert pair.reverse == reverse_complement(orf_seq[rev_start : pair.reverse_end + 1])
            assert pair.forward_pos < orf.start and pair.reverse_end > orf.end - 1
            assert 57 < pair.forward_tm < 63 and 57 < pair.reverse_tm < 63
            assert abs(pair.forward_tm - pair.reverse_tm) < 3
            assert pair.amplicon == pair.reverse_end - pair.forward_pos

    def test_pairs_sorted_by_least_stable_dimer(self, orf_seq: str) -> None:
        pairs = design_pcr(orf_seq, search=Search.OUTWARDS).pairs
        dgs = [p.dimer_dg for p in pairs]
        assert dgs == sorted(dgs, reverse=True)

    def test_fixed_search_with_equal_ranges_finds_nothing(self, orf_seq: str) -> None:
        result = design_pcr(orf_seq, search=Search.FIXED)
        assert result.pairs == []

    def test_explicit_target(self, orf_seq: str) -> None:
        target = TargetRange(150, 250, 700, 800)
        result = design_pcr(orf_seq, target, PCRSettings(max_tm_diff=5))
        assert result.pairs
        assert all(150 <= p.forward_pos < 250 and 700 < p.reverse_end <= 800 for p in result.pairs)

    def test_invalid_target(self, orf_seq: str) -> None:
        with pytest.raises(DesignError):
            design_pcr(orf_seq, TargetRange(300, 200, 700, 800))

    def test_no_orf_and_no_target(self) -> None:
        with pytest.raises(DesignError):
            design_pcr("TAA" * 100)

    def test_empty_sequence(self) -> None:
        with pytest.raises(DesignError):
            design_pcr(">only a header\n")

    def test_settings_validation(self) -> None:
        with pytest.raises(ValueError):
            PCRSettings(min_length=25, max_length=20)
        with pytest.raises(ValueError):
            PCRSettings(min_tm=65, max_tm=60)


@pytest.fixture(scope="module")
def bisulphite_seq() -> str:
    """CpG island flanked by C-rich, CpG-free sequence."""
    rng = random.Random(6)
    flank = "".join(rng.choices("ACGT", weights=[0.15, 0.4, 0.3, 0.15], k=600))
    flank = flank.replace("CG", "CA")
    return flank[:300] + cpg_rich(rng, 300) + flank[300:]


class TestBisulphite:
    def test_primers_match_converted_template(self, bisulphite_seq: str) -> None:
        seq = bisulphite_seq
        result = design_bisulphite(seq, TargetRange(150, 300, 600, 750))
        assert result.pairs
        converted = seq.replace("C", "T")
        for p in result.pairs:
            assert p.forward == converted[p.forward_pos : p.forward_pos + p.forward_len]
            rev_start = p.reverse_end - p.reverse_len + 1
            assert p.reverse == reverse_complement(converted[rev_start : p.reverse_end + 1])
            assert p.forward_degenerate and p.reverse_degenerate

    def test_degenerate_bases_mark_cpgs(self) -> None:
        # Forward primers are converted top strand: CpG cytosines become Y.
        assert _degenerate_bisulphite("TTGTTATG", "TCGTCACG", reverse=False) == "TYGTTAYG"
        # Reverse windows are read from the reversed top strand, where a CpG
        # appears as GC; the primer base opposite that C becomes R.
        assert _degenerate_bisulphite("AACAACAA", "TTGCTGCT", reverse=True) == "AACRACRA"

    def test_uses_cpg_island_by_default(self, bisulphite_seq: str) -> None:
        result = design_bisulphite(bisulphite_seq, search=Search.OUTWARDS)
        assert result.target.inner_start is not None


class TestQPCR:
    def test_pairs_span_junctions(self) -> None:
        rng = random.Random(7)
        seq = random_dna(rng, 900)
        bounds = [180, 330, 520, 700]
        settings = QPCRSettings()
        result = design_qpcr(seq, bounds, settings)
        assert result.pairs
        for p in result.pairs:
            assert settings.min_amplicon <= p.amplicon <= settings.max_amplicon
            assert any(p.forward_pos < b < p.reverse_end for b in bounds)

    def test_requires_boundaries(self) -> None:
        with pytest.raises(DesignError):
            design_qpcr(random_dna(random.Random(1), 500), [])
        with pytest.raises(DesignError):
            design_qpcr(random_dna(random.Random(1), 500), [600])

    def test_exon_restriction(self) -> None:
        rng = random.Random(7)
        seq = random_dna(rng, 900)
        bounds = [180, 330, 520, 700]
        result = design_qpcr(seq, bounds, QPCRSettings(first_exon=3, last_exon=4))
        assert all(p.forward_pos >= 330 - 24 for p in result.pairs)


class TestSequencing:
    def test_spacing(self) -> None:
        seq = random_dna(random.Random(8), 4000)
        settings = SequencingSettings()
        result = design_sequencing(seq, settings=settings)
        positions = [p.pos for p in result.primers]
        assert len(positions) > 3
        for a, b in itertools.pairwise(positions):
            assert settings.min_spacing <= b - a <= settings.max_spacing
        for p in result.primers:
            assert p.seq == seq[p.pos : p.pos + p.length]
            assert p.tm == round(melting_temperature(p.seq), 2)

    def test_invalid_range(self) -> None:
        with pytest.raises(DesignError):
            design_sequencing(random_dna(random.Random(1), 500), 400, 100)


class TestCloning:
    def test_rebase_parsed(self) -> None:
        enzymes = rebase_enzymes()
        assert len(enzymes) > 500
        ecori = find_enzyme("ecori")
        assert ecori and ecori.site == "GAATTC" and ecori.pattern == "G'AATT_C"

    def test_cloning_enzymes_exclude_internal_sites(self) -> None:
        names = {e.name for e in cloning_enzymes("AAAGAATTCAAA")}
        assert "EcoRI" not in names
        assert "BamHI" in names
        assert all(len(e.site) == 6 for e in cloning_enzymes())

    def test_tails_and_frame_spacers(self) -> None:
        cloned = add_cloning_tails(
            "ACGTACGTACGTACGTACGT", 101, "TTGCAAGCTTGCAAGCTTGC", 600,
            gene_frame=0,
            forward_tail="GCGCGC_GAATTC",
            reverse_tail="GCGCGC_GGATCC",
            forward_frame=0,
            reverse_frame=0,
        )  # fmt: skip
        # Forward: ((101 - 0) % 3 + 3 + (3 - 6)) % 3 = 2 spacer bases.
        assert cloned.forward == "GCGCGCGAATTCAA_ACGTACGTACGTACGTACGT"
        # Reverse: (2 - 600 % 3 + 0) % 3 = 2 spacer bases.
        assert cloned.reverse == "GCGCGCGGATCCAA_TTGCAAGCTTGCAAGCTTGC"

    def test_tails_without_frame(self) -> None:
        cloned = add_cloning_tails("ACGT", 0, None, 0, forward_tail="GC_GAATTC")
        assert cloned.forward == "GCGAATTC_ACGT"
        assert cloned.reverse is None


class TestExons:
    @staticmethod
    def gene() -> tuple[str, str, list[int]]:
        rng = random.Random(9)
        exons = [random_dna(rng, n) for n in (150, 90, 120, 200)]
        introns = ["GT" + random_dna(rng, n) + "AG" for n in (300, 500, 250)]
        genomic = random_dna(rng, 100) + exons[0]
        for intron, exon in zip(introns, exons[1:], strict=True):
            genomic += intron + exon
        genomic += random_dna(rng, 100)
        return "".join(exons), genomic, [150, 240, 360]

    def test_maps_boundaries(self) -> None:
        mrna, genomic, bounds = self.gene()
        exon_map = map_exons(mrna, genomic)
        assert exon_map.boundaries == bounds
        assert exon_map.unaligned_tail == 0
        for exon in exon_map.exons:
            assert (
                mrna[exon.mrna_start : exon.mrna_end]
                == genomic[exon.genomic_start : exon.genomic_end]
            )

    def test_tolerates_poly_a_tail(self) -> None:
        mrna, genomic, bounds = self.gene()
        exon_map = map_exons(mrna + "A" * 30, genomic)
        assert exon_map.boundaries == bounds
        assert exon_map.unaligned_tail == 30

    def test_unrelated_sequences(self) -> None:
        rng = random.Random(10)
        with pytest.raises(AlignmentError):
            map_exons(random_dna(rng, 300), random_dna(rng, 3000))


class TestReport:
    def test_pair_analysis(self) -> None:
        analysis = analyse_pair("GAGAATGTGGGCTCTCAAGC", "GGCGTTAGTTTAAGTCACGG")
        assert analysis.forward and analysis.forward.tm == 61.23
        labels = [d.label for d in analysis.extensible]
        assert labels == ["Forward vs. Forward"]
        assert "Forward vs. Reverse" in [d.label for d in analysis.non_extensible]
        assert "Forward vs. Forward: -1.26 kcal/mol" in analysis.dimer_text()

    def test_degenerate_primer_reported(self) -> None:
        analysis = analyse_pair("ACGTNACGTACGTACGT", None)
        assert analysis.errors and analysis.forward is None

    def test_sequence_map(self) -> None:
        seq = "cc" + "ATG" + "GCT" * 30 + "TAA" + "c" * 40
        placement = MapPrimers(forward="ATGGCTGCT", forward_pos=2, forward_len=9)
        lines = sequence_map(seq, placement, [Region(2, 98)], width=40).split("\n")
        # Blocks of ruler, primers, sequence and translation, separated by a blank line.
        assert len(seq) == 138 and len(lines) == 4 * 5
        assert lines[0] == "0.........10........20........30........"
        assert lines[5].startswith("40........50")
        assert lines[1] == "5' ATGGCTGCT 3'"
        assert lines[2] == seq[:40]
        assert lines[3].startswith("  M  A  A")

    def test_text_report(self) -> None:
        analysis = analyse_pair("GAGAATGTGGGCTCTCAAGC", "GGCGTTAGTTTAAGTCACGG")
        report = text_report("demo", analysis, "MAP\n", amplicon="500 bases (10 - 510)")
        assert report.startswith("[demo - Report generated at ")
        assert "Forward: 5' GAGAATGTGGGCTCTCAAGC 3'" in report
        assert "Tm: 61.23°C\t\tLength: 20 bases" in report
        assert report.rstrip().endswith("MAP")
