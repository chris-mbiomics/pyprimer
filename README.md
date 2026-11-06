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
- **Standard PCR**: primer pairs flanking a target (by default the longest
  ORF), with automatic widening or narrowing of the range until pairs are
  found, and restriction-site cloning tails kept in frame (REBASE data).
- **Bisulphite PCR**: primers specific for the converted strand, CpG island
  prediction, and degenerate (Y/R) primers for CpG sites.
- **Real-time PCR**: amplicons spanning exon-exon junctions, from manually
  entered boundaries, a genomic-sequence alignment, or Ensembl annotation.
- **Sequencing**: a tiled set of primers at a chosen spacing.
- **Ensembl** retrieval of cDNA, CDS or genomic sequence by gene name or ID.
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
  cli.py         command-line entry point
  web/           Flask app factory, JSON API, template and static assets
  data/gcg.603   REBASE restriction enzyme database
tests/           unit, API and PerlPrimer parity tests
tools/           Perl harness generator used by the parity tests
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

## Differences from PerlPrimer

Intentional changes:

- **Bisulphite reverse primers**: PerlPrimer searched reverse primers for
  `CG` on the reversed (not complemented) top strand. That finds GpC sites,
  not CpG sites, so CpG exclusion and the degenerate `R` bases were applied
  to the wrong positions. PyPrimer checks the real CpG sites. The parity test
  disables this fix to confirm everything else matches.
- **CpG islands** that run to the end of the sequence are now reported
  (PerlPrimer silently dropped them).
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

## Licence and credits

PyPrimer is free software under the GNU General Public License, version 2
or (at your option) any later version. See [LICENSE](LICENSE).

- Original PerlPrimer © 2003-2011 Owen Marshall.
- Restriction enzyme data from [REBASE](http://rebase.neb.com/).
- Marshall OJ. PerlPrimer: cross-platform, graphical primer design for
  standard, bisulphite and real-time PCR. *Bioinformatics* 2004
  20(15):2471-2472.
