"""
core/connectivity.py
--------------------
Bond detection from atomic positions using the covalent-radius criterion.

Two atoms A and B are considered bonded when:
    dist(A, B)  <  (r_cov(A) + r_cov(B)) × tolerance

Reference radii: Alvarez (2008), DOI 10.1039/b801115j
Default tolerance: 1.20  (standard in VESTA, Mercury, ASE)

The module uses scipy.spatial.cKDTree for O(N log N) neighbour search,
making it practical for unit cells with tens of thousands of atoms.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field

import numpy as np

# ── Covalent radii (Å) — Alvarez 2008 ────────────────────────────────────────
COVALENT_RADII: dict[str, float] = {
    "H":  0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B":  0.84, "C":  0.76, "N":  0.71,
    "O":  0.66, "F":  0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P":  1.07,
    "S":  1.05, "Cl": 1.02, "Ar": 1.06,
    "K":  2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V":  1.53,
    "Cr": 1.39, "Mn": 1.61, "Fe": 1.52, "Co": 1.50, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y":  1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39,
    "Te": 1.38, "I":  1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03,
    "Hf": 1.75, "Ta": 1.70, "W":  1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45,
    "Pb": 1.46, "Bi": 1.48,
}
_R_DEFAULT = 0.90   # fallback for unknown elements


def get_radius(symbol: str) -> float:
    return COVALENT_RADII.get(symbol, _R_DEFAULT)


# ─────────────────────────────────────────────────────────────────────────────
# Bond settings (per-pair cutoffs, editable by the user)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BondSettings:
    """
    Persistent bond-detection parameters, editable per species pair.

    Default radii source
    --------------------
    Alvarez, S. (2008) "Dalton Transactions", 2832–2838.
    DOI: 10.1039/b801115j
    Default tolerance factor: 1.20 (same as VESTA / Mercury / ASE).

    Custom cutoffs
    --------------
    ``custom_max``  : frozenset({sym_A, sym_B}) → max distance (Å)
    ``custom_min``  : frozenset({sym_A, sym_B}) → min distance (Å)  [optional]
    When a pair is absent from these dicts the default formula is used.
    """

    tolerance: float = 1.20   # global scale factor on (r_A + r_B)
    min_dist:  float = 0.40   # global minimum (avoids self-bonds)

    # Per-pair overrides
    custom_max: dict = dc_field(default_factory=dict)
    custom_min: dict = dc_field(default_factory=dict)

    # ── helpers ──────────────────────────────────────────────────────────────

    def default_max(self, sym_a: str, sym_b: str) -> float:
        """Default max distance using Alvarez 2008 + tolerance factor."""
        return (get_radius(sym_a) + get_radius(sym_b)) * self.tolerance

    def max_dist(self, sym_a: str, sym_b: str) -> float:
        """Effective max distance (custom override if set, otherwise default)."""
        return self.custom_max.get(
            frozenset([sym_a, sym_b]),
            self.default_max(sym_a, sym_b),
        )

    def min_dist_pair(self, sym_a: str, sym_b: str) -> float:
        """Effective min distance for this pair."""
        return self.custom_min.get(
            frozenset([sym_a, sym_b]),
            self.min_dist,
        )

    def reset(self):
        """Clear all custom overrides (restore factory defaults)."""
        self.custom_max.clear()
        self.custom_min.clear()

    def pairs_for(self, species: list[str]) -> list[tuple[str, str]]:
        """
        Return all unique (A, B) pairs (including A==A) for a list of species,
        sorted for display.
        """
        unique = sorted(set(species))
        pairs  = []
        for i, a in enumerate(unique):
            for b in unique[i:]:
                pairs.append((a, b))
        return pairs


# ─────────────────────────────────────────────────────────────────────────────
# Bond detection
# ─────────────────────────────────────────────────────────────────────────────

def compute_bonds(
    coords:    np.ndarray,
    symbols:   list[str],
    tolerance: float = 1.20,
    min_dist:  float = 0.40,
    settings:  "BondSettings | None" = None,
) -> list[tuple[int, int]]:
    """
    Find all bonded atom pairs in a structure.

    Parameters
    ----------
    coords    : (N, 3) Cartesian coordinates (Å)
    symbols   : element symbols, length N
    tolerance : global scale factor on (r_i + r_j); default 1.20 (VESTA)
    min_dist  : global minimum bond length (Å); avoids self-bonds
    settings  : optional BondSettings for per-pair overrides

    Returns
    -------
    List of (i, j) index pairs with i < j.
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError as e:
        raise ImportError(
            "scipy is required for bond detection.\n"
            "Install with: pip install scipy"
        ) from e

    if settings is not None:
        tolerance = settings.tolerance
        min_dist  = settings.min_dist

    radii   = np.array([get_radius(s) for s in symbols])
    max_cut = float((radii.max() * 2) * tolerance)   # global search radius

    tree  = cKDTree(coords)
    pairs = tree.query_pairs(max_cut, output_type="ndarray")
    if len(pairs) == 0:
        return []

    if settings is None:
        # Vectorised fast path (no per-pair overrides): compute every
        # candidate distance and cutoff in one NumPy pass.  A Python loop
        # over query_pairs is O(pairs) with a per-pair np.linalg.norm and
        # becomes the dominant cost when the pair count is large — e.g. the
        # big chiral tubes the polar-map spurious-bond scan builds for a
        # buckled lattice (MoS₂ …).  Results are identical to the loop.
        i_idx = pairs[:, 0]
        j_idx = pairs[:, 1]
        d  = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
        hi = (radii[i_idx] + radii[j_idx]) * tolerance
        mask = (d > min_dist) & (d < hi)
        return [(int(i), int(j)) for i, j in pairs[mask]]

    # Per-pair override path (user-edited cutoffs): keep the exact per-pair
    # logic, which cannot be expressed as a single vectorised comparison.
    bonds = []
    for i, j in pairs:
        d   = float(np.linalg.norm(coords[i] - coords[j]))
        si, sj = symbols[i], symbols[j]
        lo = settings.min_dist_pair(si, sj)
        hi = settings.max_dist(si, sj)
        if lo < d < hi:
            bonds.append((int(i), int(j)))

    return bonds


# ─────────────────────────────────────────────────────────────────────────────
# Build rendering arrays for lines (fast, any system size)
# ─────────────────────────────────────────────────────────────────────────────

def bond_line_arrays(
    coords:  np.ndarray,
    symbols: list[str],
    bonds:   list[tuple[int, int]],
    cpk:     dict[str, tuple],
    default_color: tuple = (0.7, 0.7, 0.7, 1.0),
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build vertex and colour arrays for half-bond line rendering.

    Each bond i→j is split at its midpoint into two segments:
      - i   → mid  coloured with CPK(i)
      - mid → j    coloured with CPK(j)

    This gives the standard VESTA "bicolour bond" appearance.

    Returns
    -------
    pts    : (2 * 2 * B, 3)  float32   — line endpoints (pairs)
    colors : (2 * 2 * B, 4)  float32   — RGBA per vertex
    """
    if not bonds:
        return np.zeros((0, 3), dtype=np.float32), np.zeros((0, 4), dtype=np.float32)

    pts_list = []
    col_list = []

    for i, j in bonds:
        ci  = np.array(cpk.get(symbols[i],  default_color), dtype=np.float32)
        cj  = np.array(cpk.get(symbols[j],  default_color), dtype=np.float32)
        pi  = coords[i]
        pj  = coords[j]
        mid = (pi + pj) * 0.5

        # Half-bond i → midpoint
        pts_list.extend([pi,  mid])
        col_list.extend([ci,  ci])
        # Half-bond midpoint → j
        pts_list.extend([mid, pj])
        col_list.extend([cj,  cj])

    pts    = np.array(pts_list, dtype=np.float32)
    colors = np.array(col_list, dtype=np.float32)
    return pts, colors
