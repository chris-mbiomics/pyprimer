"""Retrieve gene sequences and exon structures from the Ensembl REST API.

PerlPrimer scraped the Ensembl website, which has long since changed; this
module uses the supported REST service at https://rest.ensembl.org instead.
"""

from __future__ import annotations

import contextlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any

from . import __version__

SERVER = "https://rest.ensembl.org"
TIMEOUT = 30
# Stable IDs such as ENSG00000139618, ENSMUST00000000001.5 or FBgn0000490.
_STABLE_ID = re.compile(r"^[A-Za-z]+\d{5,}(\.\d+)?$")

SEQUENCE_TYPES = ("cdna", "cds", "genomic")


class EnsemblError(RuntimeError):
    """Raised when Ensembl cannot be queried or the gene is not found."""


@dataclass(frozen=True)
class Transcript:
    id: str
    name: str
    biotype: str
    canonical: bool
    exon_lengths: list[int] = field(default_factory=list)

    @property
    def exon_boundaries(self) -> list[int]:
        """cDNA positions at which each exon after the first begins."""
        bounds, pos = [], 0
        for length in self.exon_lengths[:-1]:
            pos += length
            bounds.append(pos)
        return bounds


@dataclass(frozen=True)
class Gene:
    id: str
    name: str
    species: str
    description: str
    location: str
    transcripts: list[Transcript]


def lookup_gene(species: str, query: str) -> Gene:
    """Find a gene by symbol (e.g. ``BRCA2``) or Ensembl stable ID."""
    query = query.strip()
    if not query:
        raise EnsemblError("Enter a gene name or Ensembl ID")
    species = species.strip().lower().replace(" ", "_")
    if _STABLE_ID.match(query):
        path = f"/lookup/id/{urllib.parse.quote(query)}"
    else:
        path = f"/lookup/symbol/{urllib.parse.quote(species)}/{urllib.parse.quote(query)}"
    data = _get_json(path, expand=1)
    if data.get("object_type") == "Transcript":
        data = _get_json(f"/lookup/id/{urllib.parse.quote(data['Parent'])}", expand=1)
    return _parse_gene(data)


def fetch_sequence(stable_id: str, seq_type: str = "cdna") -> str:
    """Fetch a sequence: ``cdna``/``cds`` need a transcript ID, ``genomic`` any ID."""
    if seq_type not in SEQUENCE_TYPES:
        raise EnsemblError(f"Unknown sequence type: {seq_type}")
    data = _get_json(f"/sequence/id/{urllib.parse.quote(stable_id)}", type=seq_type)
    return str(data["seq"])


def _parse_gene(data: dict[str, Any]) -> Gene:
    transcripts = []
    strand = data.get("strand", 1)
    for t in data.get("Transcript", []):
        exons = sorted(t.get("Exon", []), key=lambda e: e["start"], reverse=strand < 0)
        transcripts.append(
            Transcript(
                id=t["id"],
                name=t.get("display_name") or t["id"],
                biotype=t.get("biotype", ""),
                canonical=bool(t.get("is_canonical")),
                exon_lengths=[e["end"] - e["start"] + 1 for e in exons],
            )
        )
    transcripts.sort(key=lambda t: (not t.canonical, t.name))
    return Gene(
        id=data["id"],
        name=data.get("display_name") or data["id"],
        species=data.get("species", ""),
        description=data.get("description") or "",
        location=f"{data.get('seq_region_name')}:{data.get('start')}-{data.get('end')}",
        transcripts=transcripts,
    )


def _get_json(path: str, **params: Any) -> dict[str, Any]:
    url = f"{SERVER}{path}?{urllib.parse.urlencode({**params, 'content-type': 'application/json'})}"
    request = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": f"pyprimer/{__version__}"}
    )
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            result: dict[str, Any] = json.load(response)
            return result
    except urllib.error.HTTPError as exc:
        detail = ""
        with contextlib.suppress(ValueError, AttributeError):
            detail = json.load(exc).get("error", "")
        if exc.code in (400, 404):
            raise EnsemblError(detail or "Not found in Ensembl") from exc
        raise EnsemblError(f"Ensembl request failed ({exc.code}): {detail}") from exc
    except (urllib.error.URLError, TimeoutError) as exc:
        raise EnsemblError(f"Could not reach Ensembl: {exc}") from exc
