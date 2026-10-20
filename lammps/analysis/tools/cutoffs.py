"""Bond cutoffs (A) for the OVITO scripts (not meant to be run).

TABLE: about 1.25 x the equilibrium bond length, from compute_bonds.py
(E. Pedretti, diamond/Si3N4 work). Reference distances:
  C-C 1.544 (diamond)   C-Si 1.888 (3C-SiC)   Si-Si 2.352 (c-Si)
  Si-N 1.70-1.77 (beta-Si3N4)   Si-O 1.61 (alpha-quartz)
  C-N 1.47, C-O 1.43, N-O 1.40   C-H 1.09, N-H 1.02, O-H 0.97, Si-H 1.48
ALTERNATIVE: the values of ovito_cho.py (detached atoms) where they
differ; they changed some earlier counts.
Other pairs: 1.25 x (sum of the covalent radii of ASE), printed when
used. Every cutoff must sit in the valley after the first peak of the
partial RDF of your system: check with structure.py.
"""
from ase.data import atomic_numbers, covalent_radii

TABLE = {
    ("C", "C"): 1.95, ("C", "H"): 1.45, ("C", "O"): 1.85, ("C", "N"): 1.90,
    ("C", "Si"): 2.35, ("H", "H"): 1.20, ("H", "O"): 1.30, ("H", "N"): 1.35,
    ("H", "Si"): 1.85, ("O", "O"): 1.80, ("O", "N"): 1.75, ("O", "Si"): 2.05,
    ("N", "N"): 1.80, ("N", "Si"): 2.20, ("Si", "Si"): 2.85,
}
ALTERNATIVE = {("C", "C"): 1.90, ("C", "H"): 1.40, ("C", "O"): 1.90, ("H", "H"): 1.00}
MULT = 1.25


def cutoff(a, b, overrides=None, alternative=False):
    """(cutoff, source) for the pair a-b."""
    for table, src in ((overrides or {}, "given"),
                       (ALTERNATIVE if alternative else {}, "ovito_cho"), (TABLE, "table")):
        for key in ((a, b), (b, a)):
            if key in table:
                return table[key], src
    r = MULT * (covalent_radii[atomic_numbers[a]] + covalent_radii[atomic_numbers[b]])
    return round(r, 2), f"{MULT} x covalent radii"


def parse(items):
    """['Fe-Fe:2.9', 'C-H:1.4'] -> {('Fe', 'Fe'): 2.9, ('C', 'H'): 1.4}"""
    out = {}
    for it in items or []:
        pair, r = it.split(":")
        a, b = pair.split("-")
        out[(a, b)] = float(r)
    return out


def table_for(elements, overrides=None, alternative=False, verbose=True):
    """{(a, b): r} for all pairs of the elements present; prints them."""
    els = sorted(set(elements), key=lambda e: atomic_numbers[e])
    out = {}
    for i, a in enumerate(els):
        for b in els[i:]:
            r, src = cutoff(a, b, overrides, alternative)
            out[(a, b)] = r
            if verbose:
                print(f"  bond cutoff {a}-{b}: {r:.2f} A ({src})")
    return out
