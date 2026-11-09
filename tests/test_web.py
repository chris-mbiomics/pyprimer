"""Tests for the web interface and JSON API."""

from __future__ import annotations

import io
import json
import random
from collections.abc import Iterator
from typing import Any

import pytest
from flask import Flask
from flask.testing import FlaskClient

from pyprimer import ensembl
from pyprimer.design import PCRSettings, TargetRange
from pyprimer.web import create_app
from pyprimer.web.schema import SchemaError, build

from .helpers import orf_sequence, random_dna


@pytest.fixture
def app() -> Flask:
    return create_app({"TESTING": True, "MAX_SEQUENCE_LENGTH": 20_000})


@pytest.fixture
def client(app: Flask) -> FlaskClient:
    return app.test_client()


@pytest.fixture(scope="module")
def seq() -> str:
    return orf_sequence(random.Random(11), 250, 200, 250)


def post(client: FlaskClient, path: str, body: Any) -> tuple[int, Any]:
    response = client.post(f"/api{path}", json=body)
    return response.status_code, response.get_json()


def test_index(client: FlaskClient) -> None:
    response = client.get("/")
    assert response.status_code == 200
    assert b"PyPrimer" in response.data
    assert client.get("/static/app.js").status_code == 200


def test_defaults(client: FlaskClient) -> None:
    data = client.get("/api/defaults").get_json()
    assert data["pcr"]["min_tm"] == 57
    assert data["pcr"]["repeats"] == {"enabled": True, "run": 4, "repeat": 3}
    assert data["conditions"]["oligo_nm"] == 200
    assert data["autodimer"]["min_score"] == 7


def test_features(client: FlaskClient, seq: str) -> None:
    status, data = post(client, "/features", {"sequence": seq, "kind": "orf"})
    assert status == 200
    region = data["regions"][0]
    assert data["target"]["inner_start"] == region["start"]
    assert data["highlighted"][region["start"] : region["end"]].isupper()


def test_design_pcr_and_details(client: FlaskClient, seq: str) -> None:
    status, data = post(client, "/design/pcr", {"sequence": seq})
    assert status == 200, data
    assert data["total"] >= len(data["pairs"]) > 0
    assert data["target"] == {
        "outer_start": 0,
        "inner_start": len(seq) // 4,
        "inner_end": len(seq) - len(seq) // 4 - 1,
        "outer_end": len(seq),
    }
    pair = data["pairs"][0]

    status, details = post(
        client,
        "/details",
        {
            "mode": "pcr",
            "sequence": seq,
            "pair": pair,
            "cloning": {"forward_enzyme": "EcoRI", "reverse_enzyme": "BamHI", "forward_frame": 0},
        },
    )
    assert status == 200, details
    assert details["analysis"]["forward"]["seq"] == pair["forward"]
    assert details["cloned"]["forward"].startswith("GCGCGCGAATTC")
    assert details["cloned"]["reverse"].startswith("GCGCGCGGATCC_")
    assert details["map"]["sequence"] == seq.upper()

    response = client.post("/api/report", json={"mode": "pcr", "sequence": seq, "pair": pair})
    assert response.status_code == 200
    assert "attachment" in response.headers["Content-Disposition"]
    assert "SEQUENCE MAP" in response.get_data(as_text=True)


def test_design_pcr_truncates(app: Flask, seq: str) -> None:
    app.config["MAX_PAIRS"] = 2
    status, data = post(app.test_client(), "/design/pcr", {"sequence": seq})
    assert status == 200
    assert len(data["pairs"]) == 2 and data["truncated"]


def test_design_errors(client: FlaskClient, seq: str) -> None:
    assert post(client, "/design/pcr", {"sequence": ""})[0] == 400
    status, data = post(client, "/design/pcr", {"sequence": seq, "bogus": 1})
    assert status == 400 and "bogus" in data["error"]
    status, data = post(client, "/design/pcr", {"sequence": seq, "settings": {"min_tm": "hot"}})
    assert status == 400 and "min_tm" in data["error"]
    status, data = post(client, "/design/pcr", {"sequence": "A" * 30_000})
    assert status == 400 and "too long" in data["error"]
    response = client.post("/api/design/pcr", data="not json", content_type="application/json")
    assert response.status_code == 400


def test_design_qpcr_with_genomic(client: FlaskClient) -> None:
    rng = random.Random(9)
    exons = [random_dna(rng, n) for n in (250, 160, 220)]
    introns = ["GT" + random_dna(rng, n) + "AG" for n in (300, 400)]
    genomic = exons[0] + introns[0] + exons[1] + introns[1] + exons[2]
    status, data = post(client, "/design/qpcr", {"mrna": "".join(exons), "genomic": genomic})
    assert status == 200, data
    assert data["exon_boundaries"] == [250, 410]
    assert len(data["exons"]) == 3


def test_design_bisulphite_and_sequencing(client: FlaskClient) -> None:
    rng = random.Random(8)
    status, data = post(client, "/design/sequencing", {"sequence": random_dna(rng, 3000)})
    assert status == 200 and data["primers"]
    status, data = post(
        client,
        "/design/bisulphite",
        {"sequence": random_dna(rng, 600), "target": {"inner_start": 200, "inner_end": 400}},
    )
    assert status == 200 and "pairs" in data


def test_analyse(client: FlaskClient) -> None:
    status, data = post(client, "/analyse", {"forward": "GAGAATGTGGGCTCTCAAGC"})
    assert status == 200
    assert data["forward"]["tm"] == 61.23 and data["reverse"] is None
    assert post(client, "/analyse", {})[0] == 400


def test_autodimer(client: FlaskClient) -> None:
    oligos = ">F\nGTCAGCTAGTAGCTAGCGCA\n>R\nTGCGCTAGCTACTAGCTGAC\n"
    body = {"sequences": oligos, "settings": {"min_score": 20}}
    status, data = post(client, "/autodimer", body)
    assert status == 200
    assert data["comparisons"] == 1 and data["sequences"] == 2
    assert [(h["label_1"], h["label_2"], h["score"]) for h in data["hits"]] == [("R", "F", 20)] * 2
    assert data["hits"][0]["tm"] == pytest.approx(59.7, abs=0.05)
    assert "Number of Hits = 2" in data["text"]

    status, data = post(client, "/autodimer", {**body, "screen": "hairpin", "settings": {}})
    assert status == 200 and data["comparisons"] == 2
    assert all(h["loop"] in (4, 5) for h in data["hits"])


@pytest.mark.parametrize(
    ("body", "message"),
    [
        ({"sequences": "ACGTACGT"}, "at least two"),
        ({"sequences": ""}, "at least two"),
        ({"sequences": "", "screen": "hairpin"}, "at least one"),
        ({"sequences": "ACGT\nACGT", "screen": "bogus"}, "Unknown screen"),
        ({"sequences": "ACGT\nACGT", "settings": {"sodium_m": 0}}, "greater than zero"),
        ({"sequences": "ACGT\n" * 101}, "Too many"),
    ],
)
def test_autodimer_errors(client: FlaskClient, body: dict[str, object], message: str) -> None:
    status, data = post(client, "/autodimer", body)
    assert status == 400 and message in data["error"]


def test_enzymes(client: FlaskClient) -> None:
    status, data = post(client, "/enzymes", {"region": "GAATTC"})
    assert status == 200
    names = {e["name"] for e in data}
    assert "BamHI" in names and "EcoRI" not in names


def test_exons_endpoint_error(client: FlaskClient) -> None:
    status, data = post(client, "/exons", {"mrna": "ACGT" * 10, "genomic": "TTTT" * 50})
    assert status == 400 and "align" in data["error"]


class FakeResponse(io.BytesIO):
    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


@pytest.fixture
def fake_ensembl(monkeypatch: pytest.MonkeyPatch) -> Iterator[list[str]]:
    gene = {
        "object_type": "Gene", "id": "ENSG01", "display_name": "ABC1", "species": "homo_sapiens",
        "description": "test gene", "seq_region_name": "7", "start": 100, "end": 900, "strand": -1,
        "Transcript": [
            {"id": "ENST02", "display_name": "ABC1-202", "biotype": "protein_coding",
             "is_canonical": 0, "Exon": [{"start": 100, "end": 199}]},
            {"id": "ENST01", "display_name": "ABC1-201", "biotype": "protein_coding",
             "is_canonical": 1, "Exon": [{"start": 100, "end": 149}, {"start": 400, "end": 499},
                                          {"start": 800, "end": 819}]},
        ],
    }  # fmt: skip
    calls: list[str] = []

    def urlopen(request: Any, timeout: float) -> FakeResponse:
        url = request.full_url
        calls.append(url)
        body = {"seq": "ACGT"} if "/sequence/" in url else gene
        return FakeResponse(json.dumps(body).encode())

    monkeypatch.setattr(ensembl.urllib.request, "urlopen", urlopen)
    yield calls


def test_ensembl_lookup(client: FlaskClient, fake_ensembl: list[str]) -> None:
    data = client.get("/api/ensembl/gene?species=Homo sapiens&q=ABC1").get_json()
    assert "/lookup/symbol/homo_sapiens/ABC1" in fake_ensembl[0]
    first = data["transcripts"][0]
    assert first["id"] == "ENST01" and first["canonical"]
    # Minus strand: exons are ordered from the highest coordinate.
    assert first["exon_lengths"] == [20, 100, 50]
    assert first["exon_boundaries"] == [20, 120]

    data = client.get("/api/ensembl/sequence?id=ENST01&type=cdna").get_json()
    assert data["sequence"] == "ACGT"
    assert "type=cdna" in fake_ensembl[-1]


def test_ensembl_lookup_by_id(fake_ensembl: list[str]) -> None:
    ensembl.lookup_gene("homo_sapiens", "ENSG00000139618")
    assert "/lookup/id/ENSG00000139618" in fake_ensembl[0]


def test_ensembl_unavailable(client: FlaskClient, monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(request: Any, timeout: float) -> None:
        raise ensembl.urllib.error.URLError("offline")

    monkeypatch.setattr(ensembl.urllib.request, "urlopen", fail)
    response = client.get("/api/ensembl/gene?q=ABC1")
    assert response.status_code == 502
    assert "Could not reach Ensembl" in response.get_json()["error"]


class TestSchema:
    def test_nested_and_optional(self) -> None:
        settings = build(PCRSettings, {"min_tm": 55, "gc": {"min_gc": 30}, "repeats": None})
        assert settings.min_tm == 55.0 and settings.gc.min_gc == 30.0
        target = build(TargetRange, {"inner_start": 5, "inner_end": None})
        assert target.inner_start == 5 and target.inner_end is None

    @pytest.mark.parametrize(
        "payload",
        [{"min_length": 20.5}, {"exclude_gc": 1}, {"min_tm": True}, {"gc": 5}, {"x": 1}],
    )
    def test_rejects_bad_types(self, payload: dict[str, Any]) -> None:
        with pytest.raises(SchemaError):
            build(PCRSettings, payload)

    def test_validation_errors_become_schema_errors(self) -> None:
        with pytest.raises(SchemaError, match="Tm"):
            build(PCRSettings, {"min_tm": 70})
