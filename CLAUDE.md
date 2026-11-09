# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PyPrimer is a Python port of PerlPrimer 1.1.21 (Owen Marshall, GPL-2.0-or-later) — primer design for standard PCR, bisulphite PCR, real-time PCR and sequencing — with a Flask web UI replacing the original Perl/Tk GUI. The original Perl distribution sits in `perlprimer-1.1.21/` (git-ignored, reference only; `perlprimer.pl` is Latin-1 with CRLF line endings, so use `grep -a`).

## Commands

All tooling runs from the project virtualenv (`.venv`, Python 3.12; package installed editable with `pip install -e '.[dev]'`).

```bash
.venv/bin/pyprimer serve --port 8000          # web UI (binds 127.0.0.1 unless --host given; --debug for Flask reloader)
.venv/bin/pyprimer analyse FWD [REV]          # Tm + dimer report on the command line
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/mypy                                 # strict mode, checks src/pyprimer
.venv/bin/pytest                               # parity tests auto-skip without PERLPRIMER_PL
.venv/bin/pytest tests/test_design.py::TestPCR::test_explicit_target   # single test
PERLPRIMER_PL=perlprimer-1.1.21/perlprimer.pl .venv/bin/pytest tests/test_parity.py
```

## Architecture

Library modules under `src/pyprimer/` are pure functions over frozen dataclasses; the web layer only adapts JSON to them.

- `thermo.py` → `dimer.py` → `design.py` is the core dependency chain. `Conditions` (frozen, hashable) is threaded through everything and is part of the `lru_cache` keys for Tm and dimer calculations, which is what keeps pair scoring fast.
- `design.py`: each mode scans candidate windows (`_scan_bounds` → `_windows` → `_pcr_candidates` / `_bisulphite_candidates`), then `_pair_primers` filters and scores pairs and sorts them by full dimer ΔG then extensible ΔG (both descending, stable). `_resolve_target` fills unset target fields from `default_target` (first/last quarter of the sequence, capped at `DEFAULT_REGION` = 300 bp); `_search` implements the "outwards/inwards" 10 bp auto-stepping. `TargetRange` semantics: forward primers start in `[outer_start, inner_start)`, reverse primers end in `(inner_end, outer_end]`; reverse primers carry both `reverse_pos` (on the reverse strand) and `reverse_end` (top-strand coordinate).
- Bisulphite reverse candidates are scanned on the *reversed (not complemented)* top strand, then converted and complemented — keep this in mind when touching CpG/degenerate logic (`_count_cpg`, `_degenerate_bisulphite`).
- `autodimer.py` is a port of NIST AutoDimer (multiplex primer-dimer/hairpin screening). It is independent of `thermo.py`/`Conditions`, because it has its own SantaLucia table and Na⁺-only salt correction. It is calibrated against the NIST web service: `tests/test_autodimer_nist.py` replays `tests/data/autodimer_nist.json` (regenerate with `tools/fetch_autodimer_nist.py`) and requires exact equality. Its deliberate NIST quirks are listed in README.md ("AutoDimer").
- `report.py` builds primer analysis, dimer text and the aligned sequence-map lines (`map_lines`) used both by the text report and by the browser's map renderer.
- `web/api.py`: every endpoint parses its body with `web/schema.build(RequestDataclass, json)`, which validates types strictly against the dataclass type hints and rejects unknown fields; responses use `schema.dump`. Error types map to status codes via blueprint error handlers (`SchemaError`/`DesignError`/`AlignmentError` → 400, `EnsemblError` → 502).
- Front end (`web/templates/index.html`, `web/static/app.js`) is plain Jinja + vanilla JS, no build step. Form inputs declare `data-path="settings.min_tm"` etc.; `readForm` turns them into the nested JSON the API dataclasses expect. Global preferences (conditions, GC/repeat filters, CpG settings) live in localStorage and are merged into each request client-side.

## Faithfulness to PerlPrimer

Algorithms are deliberately literal ports, including Perl quirks (zero treated as "unset" in `_scan_bounds`, amplicon size off by one, ORF scan stopping before the last codon, dimer tie-breaking by last offset). Don't "fix" these without intent — `tests/test_parity.py` runs the original Perl subroutines (extracted headlessly by `tools/perl_harness.py`) on random inputs and requires exact equality. Intentional deviations are listed in README.md ("Differences from PerlPrimer"); if you add one, document it there and, if it breaks parity, monkeypatch the legacy behaviour in the parity test as `test_design_bisulphite` does.

## Data

`src/pyprimer/data/gcg.603` is the REBASE enzyme database (GCG format, read as Latin-1), shipped as package data.
