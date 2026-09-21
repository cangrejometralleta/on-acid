"""
dissociation.py — acid-base equilibrium without solving a quadratic by hand.

The original idea: for a weak monoprotic acid the quadratic

    x^2 + Ka*x - Ka*C = 0        (x = [H+])

ALWAYS has the same coefficients (a=1, b=Ka, c=-Ka*C), so it is solved
once, symbolically, and what remains is a direct operation — exact, with
no C - x ~= C approximation:

    [H+] = (Ka/2) * (sqrt(1 + 4*C/Ka) - 1)

For several constants (polyprotic acids) this generalizes to a sum over
the species: every proton released contributes j * alpha_j * C, and the
proton balance is closed with water included.

For bases the same algebra runs on [OH-] and closes with pH = pKw - pOH.

No dependencies: standard library only.

Background
----------
- Acid dissociation constant:
  https://en.wikipedia.org/wiki/Acid_dissociation_constant
- Self-ionization of water (the Kw term in the balance):
  https://en.wikipedia.org/wiki/Self-ionization_of_water
- Loss of significance, which the rationalized form avoids:
  https://en.wikipedia.org/wiki/Loss_of_significance
- Bisection method, used to close the polyprotic balance:
  https://en.wikipedia.org/wiki/Bisection_method
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Sequence

__all__ = [
    "KW",
    "monoprotic_h",
    "polyprotic_h",
    "acid_ph",
    "base_ph",
    "ph",
    "henderson_hasselbalch",
    "buffer_h",
    "buffer_ph",
    "equivalence_volume",
    "titration_ph",
    "titration_curve",
    "fractions",
    "pk",
    "k_from_pk",
    "Species",
    "ACIDS",
    "BASES",
    "find_species",
    "demo",
]

KW = 1.0e-14  # ion product of water at 25 C


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def pk(k: float) -> float:
    """pKa/pKb/pH from a constant or a concentration."""
    if k <= 0:
        raise ValueError("the constant must be > 0")
    return -math.log10(k)


def k_from_pk(value: float) -> float:
    """Constant from its pK."""
    return 10.0 ** (-value)


def _pkw(kw: float) -> float:
    return -math.log10(kw)


# --------------------------------------------------------------------------
# monoprotic case: the direct formula
# --------------------------------------------------------------------------

def monoprotic_h(ka: float, c: float) -> float:
    """
    Exact [H+] for a weak monoprotic acid, ignoring water.

    Equivalent to (Ka/2)*(sqrt(1 + 4C/Ka) - 1), written in rationalized
    form: when Ka << C, subtracting two nearly equal numbers loses digits
    in floating point, and this form does not.
    """
    if ka <= 0 or c < 0:
        raise ValueError("Ka must be > 0 and C must be >= 0")
    if c == 0:
        return 0.0
    return 2.0 * ka * c / (ka + math.sqrt(ka * ka + 4.0 * ka * c))


def _classic_monoprotic_h(ka: float, c: float) -> float:
    """The same formula straight out of the quadratic, kept for comparison."""
    return (ka / 2.0) * (math.sqrt(1.0 + 4.0 * c / ka) - 1.0)


# --------------------------------------------------------------------------
# general case: a sum over the species
# --------------------------------------------------------------------------

def _betas(kas: Sequence[float]) -> list[float]:
    """Cumulative constants beta_j = Ka1*Ka2*...*Kaj, with beta_0 = 1."""
    products = [1.0]
    for ka in kas:
        if ka <= 0:
            raise ValueError("every Ka must be > 0")
        products.append(products[-1] * ka)
    return products


def fractions(kas: Sequence[float], h: float) -> list[float]:
    """
    Fraction alpha_j of each species at a given [H+].

    alpha_0 is the fully protonated form (H_nA), alpha_n is A^n-.
    They sum to 1.
    """
    return _fractions_from_betas(_betas(kas), h)


def _fractions_from_betas(betas: Sequence[float], h: float) -> list[float]:
    """Evaluate the fractions from cumulative constants already prepared."""
    n = len(betas) - 1
    terms = [betas[j] * h ** (n - j) for j in range(n + 1)]
    total = sum(terms)
    return [t / total for t in terms]


def _released_protons(betas: Sequence[float], c: float, h: float) -> float:
    """The sum: sum_j j * alpha_j * C — protons the acid contributes."""
    alphas = _fractions_from_betas(betas, h)
    return c * sum(j * a for j, a in enumerate(alphas))


def _solve_charge_balance(
    betas: Sequence[float],
    c: float,
    kw: float,
    spectator: float = 0.0,
) -> float:
    """
    Close the charge balance for h, by bisection on a logarithmic scale.

        h + spectator = sum_j j*alpha_j*C + Kw/h

    `spectator` is the net charge of ions that do not react: positive for
    the cation a strong base brings in, negative for the anion of a strong
    acid. The residual grows with h, so bisection always converges.
    """
    def residual(h: float) -> float:
        return h + spectator - _released_protons(betas, c, h) - kw / h

    lo, hi = 1e-20, max(1.0, c * (len(betas) - 1) + abs(spectator) + 1.0)
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if residual(mid) > 0:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def polyprotic_h(kas: Sequence[float], c: float, kw: float = KW) -> float:
    """
    Exact [H+] for an acid with one or several constants, water included.
    The proton balance

        h = sum_j j*alpha_j*C + Kw/h

    is closed by bisection on h, which is monotone and always converges.
    """
    kas = list(kas)
    if not kas:
        raise ValueError("at least one Ka is required")
    if c < 0:
        raise ValueError("C must be >= 0")
    return _solve_charge_balance(_betas(kas), c, kw)


# --------------------------------------------------------------------------
# pH
# --------------------------------------------------------------------------

def acid_ph(
    ka: float | Sequence[float],
    c: float,
    kw: float = KW,
    exact: bool = False,
) -> float:
    """
    pH of a weak acid. `ka` is one constant or a sequence of them.

    With a single Ka and `exact` off, the direct formula answers and water
    is ignored, which is correct except in very dilute solutions or with an
    acid weak enough that water competes. With `exact=True`, or with several
    constants, the full balance is solved instead.
    """
    if isinstance(ka, (int, float)):
        if not exact:
            return pk(monoprotic_h(float(ka), c))
        kas: Sequence[float] = [float(ka)]
    else:
        kas = list(ka)
    return pk(polyprotic_h(kas, c, kw))


def base_ph(
    kb: float | Sequence[float],
    c: float,
    kw: float = KW,
    exact: bool = False,
) -> float:
    """
    pH of a weak base: [OH-] comes from the same formula and the result
    closes with pH = pKw - pOH — the familiar "subtract from 14", with pKw
    left as a parameter for temperatures other than 25 C.
    """
    poh = acid_ph(kb, c, kw=kw, exact=exact)  # same algebra, on Kb and [OH-]
    return _pkw(kw) - poh


def ph(species: "Species | str", c: float, kw: float = KW, exact: bool = False) -> float:
    """pH of a tabulated species, given as an object or as a name/formula."""
    if isinstance(species, str):
        species = find_species(species)
    if species.is_base:
        return base_ph(species.k, c, kw=kw, exact=exact)
    return acid_ph(species.k, c, kw=kw, exact=exact)


# --------------------------------------------------------------------------
# buffers
#
# https://en.wikipedia.org/wiki/Henderson%E2%80%93Hasselbalch_equation
# --------------------------------------------------------------------------

def henderson_hasselbalch(ka: float, acid_c: float, base_c: float) -> float:
    """
    pH = pKa + log10(base/acid) — the classroom approximation.

    It assumes both members of the pair survive mixing untouched, which
    fails once either one is dilute enough to compete with water or close
    enough to Ka to be consumed. Use `buffer_ph` when that matters.
    """
    if acid_c <= 0 or base_c <= 0:
        raise ValueError("a buffer needs both members of the pair")
    return pk(ka) + math.log10(base_c / acid_c)


def buffer_h(ka: float, acid_c: float, base_c: float, kw: float = KW) -> float:
    """
    Exact [H+] for a weak acid and its conjugate base in the same solution.

    The salt contributes a spectator cation equal to `base_c`, and the two
    members share one analytical concentration, so the charge balance is
    the same one the free acid solves.
    """
    if ka <= 0:
        raise ValueError("Ka must be > 0")
    if acid_c < 0 or base_c < 0:
        raise ValueError("concentrations must be >= 0")
    return _solve_charge_balance(_betas([ka]), acid_c + base_c, kw, spectator=base_c)


def buffer_ph(ka: float, acid_c: float, base_c: float, kw: float = KW) -> float:
    """pH of a buffer, solved exactly rather than approximated."""
    return pk(buffer_h(ka, acid_c, base_c, kw))


# --------------------------------------------------------------------------
# titration with a strong titrant
#
# https://en.wikipedia.org/wiki/Titration_curve
# --------------------------------------------------------------------------

def equivalence_volume(
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    protons: int = 1,
) -> float:
    """
    Titrant volume that neutralizes `protons` equivalents of the analyte.

    Volumes carry whatever unit you pass in, and the result carries it back.
    """
    if titrant_c <= 0:
        raise ValueError("the titrant must have a concentration > 0")
    return protons * analyte_c * analyte_volume / titrant_c


def titration_ph(
    kas: Sequence[float],
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    titrant_volume: float,
    kw: float = KW,
) -> float:
    """
    pH at one point of a weak acid titrated with a strong base.

    Dilution is applied to both the analyte and the titrant, and the strong
    base enters only as its spectator cation. Nothing here is piecewise, so
    the buffer region, the equivalence point and the excess-titrant tail all
    come out of the same balance.
    """
    total_volume = analyte_volume + titrant_volume
    if total_volume <= 0:
        raise ValueError("the total volume must be > 0")
    diluted_analyte = analyte_c * analyte_volume / total_volume
    spectator = titrant_c * titrant_volume / total_volume
    return pk(_solve_charge_balance(_betas(kas), diluted_analyte, kw, spectator))


def titration_curve(
    kas: Sequence[float],
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    up_to: float | None = None,
    points: int = 101,
    kw: float = KW,
) -> list[tuple[float, float]]:
    """
    The whole curve, as (titrant volume, pH) pairs.

    `up_to` defaults to twice the final equivalence volume, which shows every
    plateau and leaves room for the tail. `points` counts the samples,
    endpoints included.
    """
    if points < 2:
        raise ValueError("a curve needs at least two points")
    if up_to is None:
        up_to = 2.0 * equivalence_volume(analyte_c, analyte_volume, titrant_c, len(kas))
    step = up_to / (points - 1)
    return [
        (i * step, titration_ph(kas, analyte_c, analyte_volume, titrant_c, i * step, kw))
        for i in range(points)
    ]


# --------------------------------------------------------------------------
# common constants at 25 C
#
# Textbook values, rounded to two significant figures; see the CRC Handbook
# of Chemistry and Physics for the primary tables:
# https://en.wikipedia.org/wiki/CRC_Handbook_of_Chemistry_and_Physics
#
# Spanish names are kept as aliases so existing lookups keep resolving.
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Species:
    name: str
    formula: str
    k: tuple[float, ...]
    is_base: bool = False
    alias: tuple[str, ...] = field(default=())

    @property
    def pk(self) -> tuple[float, ...]:
        return tuple(pk(x) for x in self.k)

    def __repr__(self) -> str:
        label = "Kb" if self.is_base else "Ka"
        values = ", ".join(f"{x:.3g}" for x in self.k)
        return f"<{self.name} ({self.formula}) {label}=[{values}]>"


def _make_acid(name: str, formula: str, *k: float, alias: Iterable[str] = ()) -> Species:
    return Species(name, formula, tuple(k), False, tuple(alias))


def _make_base(name: str, formula: str, *k: float, alias: Iterable[str] = ()) -> Species:
    return Species(name, formula, tuple(k), True, tuple(alias))


ACIDS: dict[str, Species] = {s.formula: s for s in (
    _make_acid("acetic acid",       "CH3COOH",  1.8e-5,
               alias=("ácido acético", "acetico", "vinegar", "vinagre")),
    _make_acid("formic acid",       "HCOOH",    1.8e-4,
               alias=("ácido fórmico", "formico")),
    _make_acid("lactic acid",       "C3H6O3",   1.4e-4,
               alias=("ácido láctico", "lactico")),
    _make_acid("benzoic acid",      "C6H5COOH", 6.3e-5,
               alias=("ácido benzoico", "benzoico")),
    _make_acid("hydrofluoric acid", "HF",       6.8e-4,
               alias=("ácido fluorhídrico", "fluorhidrico")),
    _make_acid("nitrous acid",      "HNO2",     4.5e-4,
               alias=("ácido nitroso", "nitroso")),
    _make_acid("hydrocyanic acid",  "HCN",      6.2e-10,
               alias=("ácido cianhídrico", "cianhidrico")),
    _make_acid("hypochlorous acid", "HClO",     3.0e-8,
               alias=("ácido hipocloroso", "hipocloroso")),
    _make_acid("carbonic acid",     "H2CO3",    4.3e-7, 4.7e-11,
               alias=("ácido carbónico", "carbonico")),
    _make_acid("hydrosulfuric acid", "H2S",     8.9e-8, 1.0e-19,
               alias=("ácido sulfhídrico", "sulfhidrico")),
    _make_acid("sulfurous acid",    "H2SO3",    1.5e-2, 1.0e-7,
               alias=("ácido sulfuroso", "sulfuroso")),
    _make_acid("oxalic acid",       "H2C2O4",   5.9e-2, 6.4e-5,
               alias=("ácido oxálico", "oxalico")),
    _make_acid("ascorbic acid",     "C6H8O6",   7.9e-5, 1.6e-12,
               alias=("ácido ascórbico", "ascorbico", "vitamin c", "vitamina c")),
    _make_acid("phosphoric acid",   "H3PO4",    7.5e-3, 6.2e-8, 4.2e-13,
               alias=("ácido fosfórico", "fosforico")),
    _make_acid("citric acid",       "C6H8O7",   7.4e-4, 1.7e-5, 4.0e-7,
               alias=("ácido cítrico", "citrico")),
)}

BASES: dict[str, Species] = {s.formula: s for s in (
    _make_base("ammonia",      "NH3",     1.8e-5,
               alias=("amoníaco", "amoniaco", "amonio")),
    _make_base("methylamine",  "CH3NH2",  4.4e-4, alias=("metilamina",)),
    _make_base("ethylamine",   "C2H5NH2", 4.5e-4, alias=("etilamina",)),
    _make_base("trimethylamine", "(CH3)3N", 6.3e-5, alias=("trimetilamina",)),
    _make_base("pyridine",     "C5H5N",   1.7e-9, alias=("piridina",)),
    _make_base("aniline",      "C6H5NH2", 4.3e-10, alias=("anilina",)),
    _make_base("hydrazine",    "N2H4",    1.7e-6, alias=("hidrazina",)),
    _make_base("hydroxylamine", "NH2OH",  1.1e-8, alias=("hidroxilamina",)),
)}


def find_species(key: str) -> Species:
    """Find a species by formula, name or alias, ignoring accents and case."""
    def normalize(s: str) -> str:
        table = str.maketrans("áéíóúÁÉÍÓÚ", "aeiouAEIOU")
        return s.translate(table).strip().lower()

    target = normalize(key)
    for table in (ACIDS, BASES):
        for s in table.values():
            candidates = {normalize(s.formula), normalize(s.name), *(normalize(a) for a in s.alias)}
            if target in candidates:
                return s
    raise KeyError(f"no constants available for {key!r}")


# --------------------------------------------------------------------------

def _demo_species() -> None:
    """Prints four pH values and one speciation vector."""
    print("acetic acid 0.1 M        -> pH", round(ph("CH3COOH", 0.1), 3))
    print("ammonia 0.1 M            -> pH", round(ph("NH3", 0.1), 3))
    print("phosphoric acid 0.1 M    -> pH", round(ph("H3PO4", 0.1), 3))
    print("HCN 1e-6 M (exact)       -> pH", round(ph("HCN", 1e-6, exact=True), 3))
    print("H3PO4 speciation at pH 7 ->",
          [round(x, 4) for x in fractions(find_species("H3PO4").k, 1e-7)])


def _demo_buffer() -> None:
    """Prints the acetate pair where the approximation holds and where it drifts."""
    acetic = find_species("CH3COOH").k[0]
    print("acetate buffer 1:1       -> pH", round(buffer_ph(acetic, 0.1, 0.1), 3),
          "(Henderson-Hasselbalch:", round(henderson_hasselbalch(acetic, 0.1, 0.1), 3), ")")
    print("the same pair at 1e-5 M  -> pH", round(buffer_ph(acetic, 1e-5, 1e-5), 3),
          "(Henderson-Hasselbalch:", round(henderson_hasselbalch(acetic, 1e-5, 1e-5), 3), ")")


def _demo_titration() -> None:
    """Prints five points of the acetic acid curve against strong base."""
    acetic = find_species("CH3COOH").k[0]
    print("0.1 M acetic acid, 50 mL, titrated with 0.1 M NaOH:")
    for volume in (0.0, 25.0, 49.9, 50.0, 60.0):
        print(f"  {volume:5.1f} mL -> pH",
              round(titration_ph([acetic], 0.1, 50.0, 0.1, volume), 3))


def demo() -> None:
    """Runs the worked example. This is the console entry point."""
    _demo_species()
    _demo_buffer()
    _demo_titration()


if __name__ == "__main__":
    demo()
