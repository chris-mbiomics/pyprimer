"""Nearest-neighbour thermodynamics and melting temperature.

Tm is calculated with SantaLucia's (1998) unified nearest-neighbour
parameters. The salt correction uses the Mg2+/dNTP "sodium equivalent"
formula of von Ahsen et al. (2001) Clin. Chem. 47:1956-1961.

Mismatch free energies (used for primer-dimer scoring) come from:

* Allawi & SantaLucia (1997) Biochemistry 36:10581 (G.T)
* Allawi & SantaLucia (1998) Biochemistry 37:2170 (G.A)
* Allawi & SantaLucia (1998) Biochemistry 37:9435 (A.C)
* Allawi & SantaLucia (1998) Nucleic Acids Res. 26:2694 (C.T)
* Peyret et al. (1999) Biochemistry 38:3468 (A.A, C.C, G.G, T.T)

Duplex keys use PerlPrimer's notation: ``AA`` is 5'-AA-3'/3'-TT-5', and the
four-letter mismatch keys ``AGTG`` mean 5'-AG-3'/3'-TG-5'.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType

from .sequence import is_degenerate

GAS_CONSTANT = 1.987  # cal / (K mol)
KELVIN = 273.15


def _table(spec: str) -> dict[str, float]:
    """Parse a whitespace-separated ``key value`` list (later keys win)."""
    tokens = spec.split()
    return {tokens[i]: float(tokens[i + 1]) for i in range(0, len(tokens), 2)}


# deltaH (kcal/mol) for Watson-Crick nearest neighbours.
DH = _table("""
    AA -7.9 TT -7.9  AT -7.2 TA -7.2  CA -8.5 TG -8.5  GT -8.4 AC -8.4
    CT -7.8 AG -7.8  GA -8.2 TC -8.2  CG -10.6 GC -9.8  GG -8.0 CC -8.0
    initC 0.1 initG 0.1 initA 2.3 initT 2.3
""")

# deltaS (cal/K.mol) for Watson-Crick nearest neighbours.
DS = _table("""
    AA -22.2 TT -22.2  AT -20.4 TA -21.3  CA -22.7 TG -22.7  GT -22.4 AC -22.4
    CT -21.0 AG -21.0  GA -22.2 TC -22.2  CG -27.2 GC -24.4  GG -19.9 CC -19.9
    initC -2.8 initG -2.8 initA 4.1 initT 4.1 sym -1.4
""")

# deltaH (kcal/mol) including internal mismatches.
DH_FULL = _table("""
    AATT -7.9 TTAA -7.9  ATTA -7.2 TAAT -7.2  CAGT -8.5 TGAC -8.5
    GTCA -8.4 ACTG -8.4  CTGA -7.8 AGTC -7.8  GACT -8.2 TCAG -8.2
    CGGC -10.6 GCCG -9.8  GGCC -8.0 CCGG -8.0
    initC 0.1 initG 0.1 initA 2.3 initT 2.3

    AATA 1.2 ATAA 1.2  CAGA -0.9 AGAC -0.9  GACA -2.9 ACAG -2.9  TAAA 4.7 AAAT 4.7
    ACTC 0.0 CTCA 0.0  CCGC -1.5 CGCC -1.5  GCCC 3.6 CCCG 3.6  TCAC 6.1 CACT 6.1
    AGTG -3.1 GTGA -3.1  CGGG -4.9 GGGC -4.9  GGCG -6.0 GCGG -6.0  TGAG 1.6 GAGT 1.6
    ATTT -2.7 TTTA -2.7  CTGT -5.0 TGTC -5.0  GTCT -2.2 TCTG -2.2  TTAT 0.2 TATT 0.2

    AGTT 1.0 TTGA 1.0  ATTG -2.5 GTTA -2.5  CGGT -4.1 TGGC -4.1  CTGG -2.8 GGTC -2.8
    GGCT 3.3 TCGG 3.3  GGTT 5.8 TTGG 5.8  GTCG -4.4 GCTG -4.4  GTTG 4.1 GTTG 4.1
    TGAT -0.1 TAGT -0.1  TGGT -1.4 TGGT -1.4  TTAG -1.3 GATT -1.3

    AATG -0.6 GTAA -0.6  AGTA -0.7 ATGA -0.7  CAGG -0.7 GGAC -0.7  CGGA -4.0 AGGC -4.0
    GACG -0.6 GCAG -0.6  GGCA 0.5 ACGG 0.5  TAAG 0.7 GAAT 0.7  TGAA 3.0 AAGT 3.0

    ACTT 0.7 TTCA 0.7  ATTC -1.2 CTTA -1.2  CCGT -0.8 TGCC -0.8  CTGC -1.5 CGTC -1.5
    GCCT 2.3 TCCG 2.3  GTCC 5.2 CCTG 5.2  TCAT 1.2 TACT 1.2  TTAC 1.0 CATT 1.0

    AATC 2.3 CTAA 2.3  ACTA 5.3 ATCA 5.3  CAGC 1.9 CGAC 1.9  CCGA 0.6 AGCC 0.6
    GACC 5.2 CCAG 5.2  GCCA -0.7 ACCG -0.7  TAAC 3.4 CAAT 3.4  TCAA 7.6 AACT 7.6
""")

# deltaS (cal/K.mol) including internal mismatches.
DS_FULL = _table("""
    AATT -22.2 TTAA -22.2  ATTA -20.4 TAAT -21.3  CAGT -22.7 TGAC -22.7
    GTCA -22.4 ACTG -22.4  CTGA -21.0 AGTC -21.0  GACT -22.2 TCAG -22.2
    CGGC -27.2 GCCG -24.4  GGCC -19.9 CCGG -19.9
    initC -2.8 initG -2.8 initA 4.1 initT 4.1 sym -1.4

    AATA 1.7 ATAA 1.7  CAGA -4.2 AGAC -4.2  GACA -9.8 ACAG -9.8  TAAA 12.9 AAAT 12.9
    ACTC -4.4 CTCA -4.4  CCGC -7.2 CGCC -7.2  GCCC 8.9 CCCG 8.9  TCAC 16.4 CACT 16.4
    AGTG -9.5 GTGA -9.5  CGGG -15.3 GGGC -15.3  GGCG -15.8 GCGG -15.8  TGAG 3.6 GAGT 3.6
    ATTT -10.8 TTTA -10.8  CTGT -15.8 TGTC -15.8  GTCT -8.4 TCTG -8.4  TTAT -1.5 TATT -1.5

    AGTT 0.9 TTGA 0.9  ATTG -8.3 GTTA -8.3  CGGT -11.7 TGGC -11.7  CTGG -8.0 GGTC -8.0
    GGCT 10.4 TCGG 10.4  GGTT 16.3 TTGG 16.3  GTCG -12.3 GCTG -12.3  GTTG 9.5 GTTG 9.5
    TGAT -1.7 TAGT -1.7  TGGT -6.2 TGGT -6.2  TTAG -5.3 GATT -5.3

    AATG -2.3 GTAA -2.3  AGTA -2.3 ATGA -2.3  CAGG -2.3 GGAC -2.3  CGGA -13.2 AGGC -13.2
    GACG -1.0 GCAG -1.0  GGCA 3.2 ACGG 3.2  TAAG 0.7 GAAT 0.7  TGAA 7.4 AAGT 7.4

    ACTT 0.2 TTCA 0.2  ATTC -6.2 CTTA -6.2  CCGT -4.5 TGCC -4.5  CTGC -6.1 CGTC -6.1
    GCCT 5.4 TCCG 5.4  GTCC 13.5 CCTG 13.5  TCAT 0.7 TACT 0.7  TTAC 0.7 CATT 0.7

    AATC 4.6 CTAA 4.6  ACTA 14.6 ATCA 14.6  CAGC 3.7 CGAC 3.7  CCGA -0.6 AGCC -0.6
    GACC 14.2 CCAG 14.2  GCCA -3.8 ACCG -3.8  TAAC 8.0 CAAT 8.0  TCAA 20.2 AACT 20.2
""")

# Salt-independent initiation free energies (kcal/mol).
DG_INIT = {"initC": 0.98, "initG": 0.98, "initA": 1.03, "initT": 1.03}


@dataclass(frozen=True)
class Conditions:
    """PCR reaction conditions.

    Attributes:
        oligo_nm: Primer concentration (nM).
        mg_mm: Mg2+ concentration (mM).
        monovalent_mm: Monovalent cation (Na+/K+) concentration (mM).
        dntp_mm: Total dNTP concentration (mM).
        dimer_temp_c: Temperature at which primer-dimer dG is evaluated (C).
    """

    oligo_nm: float = 200.0
    mg_mm: float = 1.5
    monovalent_mm: float = 50.0
    dntp_mm: float = 0.2
    dimer_temp_c: float = 37.0

    def __post_init__(self) -> None:
        if self.oligo_nm <= 0:
            raise ValueError("Primer concentration must be positive")
        if min(self.mg_mm, self.monovalent_mm, self.dntp_mm) < 0:
            raise ValueError("Salt and dNTP concentrations cannot be negative")
        if self.sodium_equivalent <= 0:
            raise ValueError("Monovalent cation or free Mg2+ concentration must be positive")

    @property
    def sodium_equivalent(self) -> float:
        """Sodium-equivalent concentration (M) after von Ahsen et al. (2001)."""
        # Free Mg2+ is taken as [Mg2+] - [dNTP]; dNTPs in excess leave none.
        free_mg = self.mg_mm - self.dntp_mm
        correction = math.sqrt(free_mg) if free_mg > 0 else 0.0
        return (self.monovalent_mm + 120 * correction) / 1000


DEFAULT_CONDITIONS = Conditions()


@dataclass(frozen=True)
class Thermo:
    """Melting temperature and duplex thermodynamics of a primer."""

    tm: float  # C
    dh: float  # kcal/mol
    ds: float  # cal/(K mol), salt-corrected
    dg: float  # kcal/mol at Conditions.dimer_temp_c


def melting_temperature(primer: str, conditions: Conditions = DEFAULT_CONDITIONS) -> float:
    """Return the salt-corrected Tm (C) of a primer."""
    return thermodynamics(primer, conditions).tm


def thermodynamics(primer: str, conditions: Conditions = DEFAULT_CONDITIONS) -> Thermo:
    """Calculate Tm, dH, dS and dG of a primer hybridised to its complement.

    Raises:
        ValueError: if the primer is empty or contains degenerate bases.
    """
    primer = primer.upper()
    if not primer:
        raise ValueError("Empty primer sequence")
    if is_degenerate(primer):
        raise ValueError(f"Cannot calculate the Tm of a degenerate sequence: {primer}")
    return _thermodynamics(primer, conditions)


@lru_cache(maxsize=65536)
def _thermodynamics(primer: str, conditions: Conditions) -> Thermo:
    dh = 0.0
    ds = 0.0
    for i in range(len(primer) - 1):
        pair = primer[i : i + 2]
        dh += DH[pair]
        ds += DS[pair]
    for terminal in (primer[0], primer[-1]):
        dh += DH["init" + terminal]
        ds += DS["init" + terminal]

    ds += 0.368 * (len(primer) - 1) * math.log(conditions.sodium_equivalent)

    # For PCR, [target] << [primer], so the strand concentration term is [primer]/4.
    oligo_molar = conditions.oligo_nm / 1e9
    tm = (dh * 1000) / (ds + GAS_CONSTANT * math.log(oligo_molar / 4)) - KELVIN
    dg = dh - (KELVIN + conditions.dimer_temp_c) * (ds / 1000)
    return Thermo(tm=tm, dh=dh, ds=ds, dg=dg)


def free_energy_table(conditions: Conditions = DEFAULT_CONDITIONS) -> Mapping[str, float]:
    """Return nearest-neighbour dG values (kcal/mol) for the given conditions.

    Includes mismatch pairs and the salt-independent initiation terms.
    """
    return _free_energy_table(conditions.sodium_equivalent, conditions.dimer_temp_c)


@lru_cache(maxsize=32)
def _free_energy_table(sodium_equivalent: float, temperature: float) -> Mapping[str, float]:
    entropy_adjust = 0.368 * math.log(sodium_equivalent)
    table = dict(DG_INIT)
    for key, dh in DH_FULL.items():
        if key.startswith("init"):
            continue
        ds = DS_FULL[key] + entropy_adjust
        table[key] = dh - (KELVIN + temperature) * (ds / 1000)
    return MappingProxyType(table)
