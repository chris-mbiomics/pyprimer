"""JSON API used by the web interface."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from flask import Blueprint, Response, current_app, jsonify, request
from werkzeug.exceptions import HTTPException

from .. import ensembl
from ..autodimer import (
    AutoDimerSettings,
    full_report,
    parse_oligos,
    screen_dimers,
    screen_hairpins,
)
from ..cloning import (
    DEFAULT_ANCHOR,
    add_cloning_tails,
    cloning_enzymes,
    cloning_tail,
    find_enzyme,
)
from ..design import (
    BisulphiteSettings,
    DesignError,
    DesignResult,
    PCRSettings,
    PrimerPair,
    QPCRSettings,
    RepeatFilter,
    Search,
    SequencingPrimer,
    SequencingSettings,
    TargetRange,
    design_bisulphite,
    design_pcr,
    design_qpcr,
    design_sequencing,
)
from ..exons import AlignmentError, map_exons
from ..features import CpGSettings, Region, find_cpg_islands, find_orf, highlight_regions
from ..report import (
    MapPrimers,
    PairAnalysis,
    analyse_pair,
    map_lines,
    sequence_map,
    text_report,
)
from ..sequence import parse_sequence
from ..thermo import Conditions
from .schema import SchemaError, build, dump

api = Blueprint("api", __name__)

MODES = ("pcr", "bisulphite", "qpcr", "sequencing")


# --------------------------------------------------------------------------
# Request payloads
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class PCRRequest:
    sequence: str
    target: TargetRange | None = None
    settings: PCRSettings = field(default_factory=PCRSettings)
    conditions: Conditions = field(default_factory=Conditions)
    search: Search = Search.FIXED


@dataclass(frozen=True)
class BisulphiteRequest:
    sequence: str
    target: TargetRange | None = None
    settings: BisulphiteSettings = field(default_factory=BisulphiteSettings)
    conditions: Conditions = field(default_factory=Conditions)
    search: Search = Search.FIXED


@dataclass(frozen=True)
class QPCRRequest:
    mrna: str
    genomic: str = ""
    exon_boundaries: list[int] | None = None
    settings: QPCRSettings = field(default_factory=QPCRSettings)
    conditions: Conditions = field(default_factory=Conditions)


@dataclass(frozen=True)
class SequencingRequest:
    sequence: str
    start: int | None = None
    end: int | None = None
    settings: SequencingSettings = field(default_factory=SequencingSettings)
    conditions: Conditions = field(default_factory=Conditions)


@dataclass(frozen=True)
class FeatureRequest:
    sequence: str
    kind: str = "orf"
    cpg: CpGSettings = field(default_factory=CpGSettings)


@dataclass(frozen=True)
class AnalyseRequest:
    forward: str = ""
    reverse: str = ""
    conditions: Conditions = field(default_factory=Conditions)
    repeats: RepeatFilter = field(default_factory=RepeatFilter)


@dataclass(frozen=True)
class CloningOptions:
    forward_enzyme: str = ""
    reverse_enzyme: str = ""
    anchor: str = DEFAULT_ANCHOR
    forward_frame: int | None = None
    reverse_frame: int | None = None


@dataclass(frozen=True)
class DetailsRequest:
    mode: str
    sequence: str
    pair: PrimerPair | None = None
    primer: SequencingPrimer | None = None
    conditions: Conditions = field(default_factory=Conditions)
    cpg: CpGSettings = field(default_factory=CpGSettings)
    cloning: CloningOptions | None = None
    title: str = ""


@dataclass(frozen=True)
class EnzymeRequest:
    region: str = ""
    exclude_found: bool = True


@dataclass(frozen=True)
class ExonRequest:
    mrna: str
    genomic: str


@dataclass(frozen=True)
class AutoDimerRequest:
    sequences: str = ""
    settings: AutoDimerSettings = field(default_factory=AutoDimerSettings)


# --------------------------------------------------------------------------
# Error handling
# --------------------------------------------------------------------------


class ApiError(Exception):
    def __init__(self, message: str, status: int = 400) -> None:
        super().__init__(message)
        self.message = message
        self.status = status


@api.errorhandler(ApiError)
def _api_error(exc: ApiError) -> tuple[Response, int]:
    return jsonify(error=exc.message), exc.status


@api.errorhandler(SchemaError)
@api.errorhandler(DesignError)
@api.errorhandler(AlignmentError)
def _bad_request(exc: Exception) -> tuple[Response, int]:
    return jsonify(error=str(exc)), 400


@api.errorhandler(ensembl.EnsemblError)
def _ensembl_error(exc: Exception) -> tuple[Response, int]:
    return jsonify(error=str(exc)), 502


@api.errorhandler(HTTPException)
def _http_error(exc: HTTPException) -> tuple[Response, int]:
    return jsonify(error=exc.description), exc.code or 500


def _payload() -> dict[str, Any]:
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        raise ApiError("Expected a JSON object")
    return data


def _sequence(text: str, label: str = "sequence") -> str:
    _, seq = parse_sequence(text)
    if not seq:
        raise ApiError(f"No DNA {label} given")
    limit = current_app.config["MAX_SEQUENCE_LENGTH"]
    if len(seq) > limit:
        raise ApiError(f"The {label} is too long ({len(seq)} bases; the limit is {limit})")
    return seq


# --------------------------------------------------------------------------
# Endpoints
# --------------------------------------------------------------------------


@api.get("/defaults")
def defaults() -> Response:
    """Default settings for every design mode."""
    return jsonify(
        pcr=dump(PCRSettings()),
        bisulphite=dump(BisulphiteSettings()),
        qpcr=dump(QPCRSettings()),
        sequencing=dump(SequencingSettings()),
        conditions=dump(Conditions()),
        cpg=dump(CpGSettings()),
        autodimer=dump(AutoDimerSettings()),
        anchor=DEFAULT_ANCHOR,
    )


@api.post("/features")
def features() -> Response:
    """Find the longest ORF or the CpG islands of a sequence."""
    req = build(FeatureRequest, _payload())
    seq = _sequence(req.sequence)
    regions = _features(seq, req.kind, req.cpg)
    return jsonify(
        length=len(seq),
        regions=dump(regions),
        highlighted=highlight_regions(seq, regions),
        target=dump(_target_from(regions)),
    )


@api.post("/design/pcr")
def pcr() -> Response:
    req = build(PCRRequest, _payload())
    seq = _sequence(req.sequence)
    result = design_pcr(seq, req.target, req.settings, req.conditions, req.search)
    return _design_response(result)


@api.post("/design/bisulphite")
def bisulphite() -> Response:
    req = build(BisulphiteRequest, _payload())
    seq = _sequence(req.sequence)
    result = design_bisulphite(seq, req.target, req.settings, req.conditions, req.search)
    return _design_response(result)


@api.post("/design/qpcr")
def qpcr() -> Response:
    req = build(QPCRRequest, _payload())
    mrna = _sequence(req.mrna, "mRNA sequence")
    exon_map = None
    bounds = req.exon_boundaries
    if bounds is None:
        if not req.genomic.strip():
            raise ApiError("Give either the exon boundaries or a genomic sequence")
        exon_map = map_exons(mrna, _sequence(req.genomic, "genomic sequence"))
        bounds = exon_map.boundaries
    result = design_qpcr(mrna, bounds, req.settings, req.conditions)
    response = _design_payload(result)
    response["exon_boundaries"] = bounds
    response["exons"] = dump(exon_map.exons) if exon_map else None
    return jsonify(response)


@api.post("/design/sequencing")
def sequencing() -> Response:
    req = build(SequencingRequest, _payload())
    seq = _sequence(req.sequence)
    result = design_sequencing(seq, req.start, req.end, req.settings, req.conditions)
    return jsonify(
        primers=dump(result.primers),
        total=len(result.primers),
        truncated=False,
        target=dump(result.target),
        messages=result.messages,
        length=len(result.sequence),
    )


@api.post("/exons")
def exons() -> Response:
    req = build(ExonRequest, _payload())
    exon_map = map_exons(_sequence(req.mrna, "mRNA sequence"), _sequence(req.genomic))
    return jsonify(
        exons=dump(exon_map.exons),
        boundaries=exon_map.boundaries,
        unaligned_tail=exon_map.unaligned_tail,
    )


@api.post("/analyse")
def analyse() -> Response:
    req = build(AnalyseRequest, _payload())
    analysis = analyse_pair(req.forward, req.reverse, req.conditions, req.repeats)
    if not analysis.forward and not analysis.reverse and not analysis.errors:
        raise ApiError("Enter at least one primer sequence")
    return jsonify(_analysis_payload(analysis))


@api.post("/autodimer")
def autodimer() -> Response:
    """AutoDimer screen of a set of oligos for primer-dimers and hairpins."""
    req = build(AutoDimerRequest, _payload())
    try:
        oligos = parse_oligos(req.sequences)
    except ValueError as exc:
        raise ApiError(str(exc)) from exc
    if not oligos:
        raise ApiError("Enter at least one sequence")
    dimers = screen_dimers(oligos, req.settings)
    hairpins = screen_hairpins(oligos, req.settings)
    return jsonify(
        dimers={"hits": dump(dimers.hits), "comparisons": dimers.comparisons},
        hairpins={"hits": dump(hairpins.hits), "comparisons": hairpins.comparisons},
        sequences=len(oligos),
        text=full_report(dimers, hairpins),
    )


@api.post("/details")
def details() -> Response:
    """Analysis, sequence map and cloning primers for a selected result."""
    return jsonify(_details(build(DetailsRequest, _payload()))[0])


@api.post("/report")
def report() -> Response:
    req = build(DetailsRequest, _payload())
    payload, analysis, seq_map = _details(req)
    pair = req.pair
    amplicon = f"{pair.amplicon} bases ({pair.forward_pos} - {pair.reverse_end})" if pair else ""
    cloned = payload["cloned"] or {}
    degenerate = payload["degenerate"] or {}
    text = text_report(
        req.title or "PyPrimer",
        analysis,
        seq_map,
        forward_label=cloned.get("forward") or degenerate.get("forward"),
        reverse_label=cloned.get("reverse") or degenerate.get("reverse"),
        amplicon=amplicon,
    )
    filename = "".join(c if c.isalnum() or c in "-_" else "_" for c in req.title) or "primers"
    return Response(
        text,
        mimetype="text/plain",
        headers={"Content-Disposition": f'attachment; filename="{filename}_report.txt"'},
    )


@api.post("/enzymes")
def enzymes() -> Response:
    """Restriction enzymes suitable for adding cloning sites to an amplicon."""
    req = build(EnzymeRequest, _payload())
    return jsonify(
        [
            {"name": e.name, "pattern": e.pattern, "site": e.site}
            for e in cloning_enzymes(req.region, req.exclude_found)
        ]
    )


@api.get("/ensembl/gene")
def ensembl_gene() -> Response:
    gene = ensembl.lookup_gene(
        request.args.get("species", "homo_sapiens"), request.args.get("q", "")
    )
    return jsonify(
        dump(gene)
        | {
            "transcripts": [
                dump(t) | {"exon_boundaries": t.exon_boundaries} for t in gene.transcripts
            ]
        }
    )


@api.get("/ensembl/sequence")
def ensembl_sequence() -> Response:
    stable_id = request.args.get("id", "").strip()
    if not stable_id:
        raise ApiError("No Ensembl ID given")
    seq_type = request.args.get("type", "cdna")
    return jsonify(
        id=stable_id, type=seq_type, sequence=ensembl.fetch_sequence(stable_id, seq_type)
    )


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def _features(seq: str, kind: str, cpg: CpGSettings) -> list[Region]:
    if kind == "orf":
        orf = find_orf(seq)
        return [orf] if orf else []
    if kind == "cpg":
        return find_cpg_islands(seq, cpg)
    raise ApiError(f"Unknown feature type: {kind}")


def _target_from(regions: list[Region]) -> TargetRange | None:
    if not regions:
        return None
    start, end = regions[0].start, regions[-1].end - 1
    return TargetRange(start, start, end, end)


def _design_payload(result: DesignResult) -> dict[str, Any]:
    limit = current_app.config["MAX_PAIRS"]
    messages = list(result.messages)
    if len(result.pairs) > limit:
        messages.append(f"Showing the first {limit} of {len(result.pairs)} pairs")
    return {
        "pairs": dump(result.pairs[:limit]),
        "total": len(result.pairs),
        "truncated": len(result.pairs) > limit,
        "target": dump(result.target),
        "messages": messages,
        "length": len(result.sequence),
    }


def _design_response(result: DesignResult) -> Response:
    return jsonify(_design_payload(result))


def _analysis_payload(analysis: PairAnalysis) -> dict[str, Any]:
    def dimers(items: list[Any]) -> list[dict[str, Any]]:
        return [{"label": d.label, "dg": d.dimer.dg, "diagram": d.dimer.diagram()} for d in items]

    return {
        "forward": dump(analysis.forward),
        "reverse": dump(analysis.reverse),
        "extensible": dimers(analysis.extensible),
        "non_extensible": dimers(analysis.non_extensible),
        "errors": analysis.errors,
        "text": analysis.dimer_text(),
    }


def _details(req: DetailsRequest) -> tuple[dict[str, Any], PairAnalysis, str]:
    if req.mode not in MODES:
        raise ApiError(f"Unknown mode: {req.mode}")
    seq = _sequence(req.sequence).upper()
    pair, primer = req.pair, req.primer
    if not pair and not primer:
        raise ApiError("Select a primer pair first")

    regions = _features(seq, "cpg" if req.mode == "bisulphite" else "orf", req.cpg)
    translate = req.mode != "bisulphite"

    if pair:
        analysis = analyse_pair(pair.forward, pair.reverse, req.conditions)
        placement = MapPrimers(
            forward=pair.forward,
            forward_pos=pair.forward_pos,
            forward_len=pair.forward_len,
            reverse=pair.reverse,
            reverse_end=pair.reverse_end,
            reverse_len=pair.reverse_len,
        )
    else:
        assert primer is not None
        analysis = analyse_pair(primer.seq, None, req.conditions)
        placement = MapPrimers(
            forward=primer.seq, forward_pos=primer.pos, forward_len=primer.length
        )

    cloned = None
    if pair and req.mode == "pcr" and req.cloning:
        cloned = _cloned(pair, req.cloning, regions)
    degenerate = None
    if pair and req.mode == "bisulphite":
        degenerate = {"forward": pair.forward_degenerate, "reverse": pair.reverse_degenerate}

    seq_map = sequence_map(seq, placement, regions, translate=translate)
    payload = {
        "analysis": _analysis_payload(analysis),
        "map": dump(map_lines(seq, placement, regions, translate=translate)),
        "regions": dump(regions),
        "cloned": cloned,
        "degenerate": degenerate,
    }
    return payload, analysis, seq_map


def _cloned(pair: PrimerPair, options: CloningOptions, orfs: list[Region]) -> dict[str, Any] | None:
    tails = {}
    for side, name in (("forward", options.forward_enzyme), ("reverse", options.reverse_enzyme)):
        if not name:
            continue
        enzyme = find_enzyme(name)
        if not enzyme:
            raise ApiError(f"Unknown restriction enzyme: {name}")
        tails[side] = cloning_tail(enzyme, options.anchor)
    if not tails:
        return None
    result = add_cloning_tails(
        pair.forward,
        pair.forward_pos,
        pair.reverse,
        pair.reverse_end,
        gene_frame=orfs[0].start % 3 if orfs else 0,
        forward_tail=tails.get("forward"),
        reverse_tail=tails.get("reverse"),
        forward_frame=options.forward_frame,
        reverse_frame=options.reverse_frame,
    )
    return {"forward": result.forward, "reverse": result.reverse}
