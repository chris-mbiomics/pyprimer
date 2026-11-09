"""Record reference AutoDimer results from the NIST DNA Analysis web service.

The original AutoDimer (Vallone & Butler 2004) is only distributed as a
compiled Visual Basic program, so ``pyprimer.autodimer`` was reconstructed
from the paper and calibrated against the web version NIST hosts at
https://dnaanalysis.nist.gov/DNAAnalysis/primer-tools. This script sends a
fixed set of synthetic sequences (seeded RNG, nothing user-specific) to that
service and stores the requests and raw responses as a test fixture, which
``tests/test_autodimer_nist.py`` replays offline.

Usage::

    python tools/fetch_autodimer_nist.py > tests/data/autodimer_nist.json
"""

from __future__ import annotations

import datetime
import json
import random
import sys
import urllib.request
from typing import Any

URL = "https://dnaanalysis.nist.gov/DNAAnalysisSvcs/analysis/getresultsfromjson"
PAPER = "GTCAGCTAGTAGCTAGCGCA"  # Figure 3 of the AutoDimer paper
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def params(score: int, na: str = "0.085", ct: str = "1.0", temp: str = "37.0") -> dict[str, str]:
    return {
        "monovalentSalt": na,
        "strandConcentration": ct,
        "deltaGTemperature": temp,
        "requiredScore": str(score),
    }


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(rng: random.Random, seq: str, count: int) -> str:
    bases = list(seq)
    for i in rng.sample(range(len(bases)), count):
        bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def interacting_pair(rng: random.Random, length: int, core: int, mutations: int) -> list[str]:
    """Two oligos whose 3' (or internal) regions are partly complementary."""
    shared = random_dna(rng, core)
    a = random_dna(rng, length - core) + shared
    b = revcomp(mutate(rng, shared, mutations)) + random_dna(rng, length - core)
    return [a, b]


def hairpin_oligo(rng: random.Random, stem: int, loop: int, tail: int) -> str:
    arm = random_dna(rng, stem)
    return random_dna(rng, tail) + arm + random_dna(rng, loop) + revcomp(arm) + random_dna(rng, 2)


def cases() -> list[dict[str, Any]]:
    rng = random.Random(20261109)
    out: list[dict[str, Any]] = []

    def add(name: str, kind: str, seqs: list[str], parameters: dict[str, str]) -> None:
        labels = [f"S{i + 1}" for i in range(len(seqs))]
        out.append(
            {
                "name": name,
                "type": kind,
                "sequences": [
                    {"label": lb, "sequence": s} for lb, s in zip(labels, seqs, strict=True)
                ],
                "parameters": parameters,
            }
        )

    # Exhaustive listings (very low threshold) show every alignment NIST makes.
    add("paper_all", "autodimer", [PAPER, "GCGCGCGCAATTGCGCGC"], params(-100))
    add("equal_length_all", "autodimer", [random_dna(rng, 20), random_dna(rng, 20)], params(-100))
    add("three_seqs_all", "autodimer", [random_dna(rng, n) for n in (12, 15, 9)], params(-100))
    add("paper_hairpin_all", "hairpin", [PAPER, random_dna(rng, 15)], params(-100))
    add("short_hairpin_all", "hairpin", [random_dna(rng, n) for n in (8, 9, 10, 11)], params(-100))
    # Self-dimers: a palindrome next to an unrelated oligo.
    add("palindrome", "autodimer", ["GAATTCCGGAATTC", "AAAAAAAAAAAA"], params(2))
    # Input handling: ambiguity codes, U and lowercase.
    add("ambiguous", "autodimer", ["GCGCGCGCAAUUGCGCGC", "GCGCRYGCAATTGCNCGC"], params(0))
    add("ambiguous_hairpin", "hairpin", ["GCGCRYGCAATTGCNCGC", "GCGCGCAAUUGCGCGC"], params(0))
    add("lowercase", "autodimer", ["gtcagctagtagctagcgca", "TGCGCTAGCTACTAGCTGAC"], params(10))
    # Thresholds and conditions on strongly interacting sets.
    for i, (na, ct, temp) in enumerate(
        [
            ("0.085", "1.0", "37.0"),
            ("0.05", "0.25", "37.0"),
            ("0.2", "5.0", "55.0"),
            ("1.0", "1.0", "25.0"),
            ("0.01", "50.0", "60.0"),
        ]
    ):
        seqs: list[str] = []
        for _ in range(2):
            seqs += interacting_pair(
                rng, rng.randint(18, 26), rng.randint(5, 12), rng.randint(0, 2)
            )
        add(f"dimer_set_{i}", "autodimer", seqs, params(3, na, ct, temp))
        hairpins = [
            hairpin_oligo(rng, rng.randint(3, 8), rng.choice([3, 4, 5, 6]), rng.randint(0, 8))
            for _ in range(4)
        ]
        add(f"hairpin_set_{i}", "hairpin", hairpins, params(1, na, ct, temp))
    for score in (5, 6, 7, 8):
        seqs = [random_dna(rng, rng.randint(18, 28)) for _ in range(5)]
        add(f"threshold_{score}", "autodimer", seqs, params(score))
        add(f"threshold_{score}_hairpin", "hairpin", seqs, params(score - 4))
    # Every alignment for unequal lengths, up to the largest difference NIST
    # handles (it fails once the longer oligo is 2n + 3 bases or more).
    add("lengths_all", "autodimer", [random_dna(rng, n) for n in (8, 11, 15, 18)], params(-100))
    # A self-complementary overlap gets the symmetry correction.
    add("self_complementary", "autodimer", ["CGCGAATTCGCG", "CGCGAATTCGCGA"], params(10))
    # (Nearly) complementary oligos give Tm and dG values to compare.
    for i, (na, ct, temp) in enumerate(
        [("0.085", "1.0", "37.0"), ("0.05", "0.25", "37.0"), ("1.0", "1.0", "25.0")]
    ):
        for j in range(4):
            x = random_dna(rng, rng.randint(8, 16))
            y = revcomp(mutate(rng, x, j % 3) if j % 3 else x)
            add(f"duplex_{i}_{j}", "autodimer", [x, y], params(len(x) - 6, na, ct, temp))
    return out


def fetch(case: dict[str, Any]) -> Any:
    body = {k: case[k] for k in ("sequences", "parameters", "type")}
    request = urllib.request.Request(
        URL,
        json.dumps(body).encode(),
        {"Content-Type": "application/json", "User-Agent": "pyprimer-test-fixtures"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        return {"error": exc.code, "body": exc.read().decode(errors="replace")}


def main() -> None:
    recorded = []
    for case in cases():
        print(f"{case['name']}...", file=sys.stderr)
        recorded.append({**case, "response": fetch(case)})
    json.dump(
        {"source": URL, "fetched": datetime.date.today().isoformat(), "cases": recorded},
        sys.stdout,
        indent=1,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
