# PyPrimer

Primer design for standard PCR, bisulphite PCR, real-time PCR (qPCR) and
sequencing, with a browser-based interface.

PyPrimer is a Python port of [PerlPrimer](http://perlprimer.sourceforge.net/)
1.1.21 by Owen Marshall. The Perl/Tk desktop GUI has been replaced with a
small web application, and the design algorithms have been ported faithfully.
Tm, ΔH/ΔS, primer-dimers, ORF and CpG prediction and all four design modes
produce identical results to the original (see [Testing](#testing)).

## Features

- **Melting temperature** from SantaLucia (1998) nearest-neighbour parameters,
  with the Mg²⁺/dNTP salt correction of von Ahsen et al. (2001).
- **Primer-dimer prediction** (extensible 3′ dimers and all dimers) using
  mismatch thermodynamics, with alignment diagrams.
- **Standard PCR**: primer pairs flanking a target range (by default, forward
  primers in the first quarter of the sequence and reverse primers in the
  last quarter, each at most 300 bp), with automatic widening or narrowing of the range until pairs
  are found, and restriction-site cloning tails kept in frame (REBASE data).
- **Bisulphite PCR**: primers specific for the converted strand, CpG island
  prediction, and degenerate (Y/R) primers for CpG sites.
- **Real-time PCR**: amplicons spanning exon-exon junctions, from manually
  entered boundaries, a genomic-sequence alignment, or Ensembl annotation.
- **Sequencing**: a tiled set of primers at a chosen spacing.
- **Ensembl** retrieval of cDNA, CDS or genomic sequence by gene name or ID.
- **AutoDimer**: screen a whole multiplex set of oligos for primer-dimers and
  hairpins, a port of NIST AutoDimer (see [AutoDimer](#autodimer)).
- Text reports, TSV export, and one-click BLAST at NCBI.

## Installation

Requires Python 3.10 or later.

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

For a production-grade WSGI server install the optional extra:
`pip install -e '.[server]'` (uses [waitress](https://docs.pylonsproject.org/projects/waitress/)).

## Usage

Start the web interface and open http://localhost:8000:

```bash
.venv/bin/pyprimer serve --open
```

Options: `--host` (default `127.0.0.1`), `--port` (default `8000`), `--debug`.
The server binds to localhost by default. Only use `--host 0.0.0.0` on a
trusted network, because the app has no authentication.

Quick analysis from the command line:

```bash
.venv/bin/pyprimer analyse GAGAATGTGGGCTCTCAAGC GGCGTTAGTTTAAGTCACGG
```

Screen a set of oligos (FASTA-style, `>label` lines) with AutoDimer:

```bash
.venv/bin/pyprimer autodimer primers.txt --score 7
.venv/bin/pyprimer autodimer primers.txt --hairpin --score 3
```

### As a library

```python
from pyprimer.design import design_pcr, Search
from pyprimer.thermo import Conditions, melting_temperature

melting_temperature("ACGTACGTAGCTAGCTAGGA", Conditions(mg_mm=2.0))

result = design_pcr(open("gene.fa").read(), search=Search.OUTWARDS)
for pair in result.pairs[:5]:
    print(pair.forward, pair.reverse, pair.amplicon, pair.dimer_dg)
```

## Project layout

```
src/pyprimer/
  sequence.py    parsing, complementation, composition, genetic code
  thermo.py      nearest-neighbour tables, Tm, dG, reaction conditions
  dimer.py       primer-dimer prediction and diagrams
  features.py    ORF and CpG island detection
  design.py      PCR, bisulphite, qPCR and sequencing primer design
  cloning.py     REBASE parsing and cloning tails
  exons.py       mRNA-to-genome exon mapping (replaces Spidey)
  ensembl.py     Ensembl REST client
  report.py      primer analysis, sequence maps, text reports
  autodimer.py   AutoDimer primer-dimer and hairpin screening of oligo sets
  cli.py         command-line entry point
  web/           Flask app factory, JSON API, template and static assets
  data/gcg.603   REBASE restriction enzyme database
tests/           unit, API, PerlPrimer and NIST AutoDimer parity tests
tools/           Perl harness generator and NIST AutoDimer fixture fetcher
```

## Development

```bash
.venv/bin/pip install -e '.[dev]'
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/mypy
.venv/bin/pytest
```

### Testing

The unit and API tests run anywhere. The parity tests compare this port
against the original Perl code on randomly generated inputs. They extract
PerlPrimer's algorithm subroutines into a headless harness
(`tools/perl_harness.py`), so they need Perl ≥ 5.14 and the original source:

```bash
PERLPRIMER_PL=perlprimer-1.1.21/perlprimer.pl .venv/bin/pytest tests/test_parity.py
```

`tests/test_autodimer_nist.py` replays results recorded from the NIST
AutoDimer web service (`tests/data/autodimer_nist.json`) offline. To record
them again (this sends only the script's synthetic test sequences to NIST):

```bash
.venv/bin/python tools/fetch_autodimer_nist.py > tests/data/autodimer_nist.json
```

## Differences from PerlPrimer

Intentional changes:

- **Bisulphite reverse primers**: PerlPrimer searched reverse primers for
  `CG` on the reversed (not complemented) top strand. That finds GpC sites,
  not CpG sites, so CpG exclusion and the degenerate `R` bases were applied
  to the wrong positions. PyPrimer checks the real CpG sites. The parity test
  disables this fix to confirm everything else matches.
- **CpG islands** that run to the end of the sequence are now reported
  (PerlPrimer silently dropped them).
- **Default target range**: with no range given, primers are placed in the
  first and last quarter of the sequence, each region at most 300 bp.
  PerlPrimer defaulted to the longest ORF (standard PCR) or CpG island
  (bisulphite); the *Find ORF* and *Find CpG islands* buttons still set the
  range that way. Wide ranges produce many pairs and take correspondingly
  longer.
- **Sequencing** uses its own *Exclude %GC* and *GC clamp* options.
  PerlPrimer showed these checkboxes but applied the Standard PCR settings.
- **Exon mapping**: NCBI Spidey is no longer distributed. PyPrimer aligns the
  mRNA to the genomic sequence itself (exact-match exons, junctions placed
  at GT…AG where the sequence is ambiguous), or takes exon structure from
  Ensembl.
- **Ensembl** is queried through its REST API instead of scraping web pages.
- **BLAST** opens NCBI web BLAST with the primers filled in, rather than
  polling and parsing results inside the application.

Not ported: the Contig Viewer TCP link, `.ppr` project files, local BLAST,
and the *Find primers with cloning* range search. Preferences are stored
per browser.

## AutoDimer

The AutoDimer tab is an addition; PerlPrimer had no equivalent. It ports
AutoDimer (Vallone & Butler, *BioTechniques* 37:226–231, 2004), which screens
sets of multiplex primers for primer-dimers and hairpins. The original Visual
Basic program is distributed only as a compiled binary, so the algorithm was
reconstructed from the paper and calibrated against the web version NIST hosts
at [dnaanalysis.nist.gov](https://dnaanalysis.nist.gov/DNAAnalysis/primer-tools).

- Each pair of different oligos is slid against each other, and each oligo is
  folded with 4- and 5-base loops and stems of at least 2 bp. Score = Watson–Crick
  matches − mismatches; alignments scoring at least the threshold are reported.
- Tm and ΔG use SantaLucia (1998) unified nearest-neighbour parameters with
  his empirical salt correction (Na⁺ only). Stacks are counted only between
  adjacent matched pairs, and each mismatch adds 0.438 kcal/mol to ΔH. These
  are rough estimates intended for ranking.

Primer-dimer results (alignments, matches, scores, Tm and ΔG) match the NIST
web service exactly, including some of its quirks:

- Self-dimers are not screened.
- For oligos of different lengths, overlaps are tried from both ends of the
  shorter oligo. If the longer one is 4 or more bases longer, a few extra
  shifts follow. These compare bases offset by one from where they are drawn.
- Unknown bases (N, ambiguity codes) count as mismatches, and the mismatch
  penalty is not applied to the last position of an alignment.

Deliberate differences from the NIST service:

- If the longer oligo is at least 2n + 3 bases (n = shorter length), NIST
  reports an error. PyPrimer caps the extra shifts instead.
- Lower-case bases are read as upper case. NIST never pairs them. U is read as T.
- **Hairpin Tm and ΔG** follow the paper: stem duplex plus a 5.6 kcal/mol
  (entropic) loop penalty, Tm = ΔH/ΔS. The values NIST reports for hairpins
  could not be reproduced. Hairpin folds, matches and scores match NIST, except
  where NIST sizes a stem wrongly next to ambiguity codes.
- The paper's own example (Figure 3) gives Tm 7.5 °C and ΔG −2.42 kcal/mol for
  the VB program, versus 4.0 °C and −1.89 kcal/mol here and on the NIST site.

## Licence and credits

PyPrimer is free software under the GNU General Public License, version 2
or (at your option) any later version. See [LICENSE](LICENSE).

- Original PerlPrimer © 2003-2011 Owen Marshall.
- Restriction enzyme data from [REBASE](http://rebase.neb.com/).
- AutoDimer: P. M. Vallone and J. M. Butler, National Institute of Standards
  and Technology (NIST).
- Marshall OJ. PerlPrimer: cross-platform, graphical primer design for
  standard, bisulphite and real-time PCR. *Bioinformatics* 2004
  20(15):2471-2472.
