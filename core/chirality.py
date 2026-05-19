"""
core/chirality.py
-----------------
Chirality engine for nanotube generation from arbitrary 2D lattices.

Key concepts
------------
For a 2D lattice with vectors a1, a2 (angle γ between them):

  Chiral vector:      Ch = n·a1 + m·a2
  Translational vec:  T  = t1·a1 + t2·a2   such that  Ch · T = 0

The perpendicularity condition expands to:
  n·t1·|a1|²  +  (n·t2 + m·t1)·(a1·a2)  +  m·t2·|a2|²  =  0

For hexagonal lattices (|a1|=|a2|, γ=60°):
  An exact integer solution always exists (by symmetry).

For rectangular / oblique lattices:
  An exact solution generally does not exist. This module finds the
  best integer approximation via a bounded search, and computes a
  *strain* metric — the fractional angular residual — that quantifies
  the periodicity error. This is a key scientific contribution of the tool.

References
----------
Dresselhaus et al., Phys. Rep. 1995 (hexagonal theory)
Frey & Doren, TubeGen 3.4, 2011 (hexagonal implementation)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from math import gcd as _gcd

import numpy as np

from .io import LatticeStructure


# ─────────────────────────────────────────────────────────────────────────────
# Result dataclass
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ChiralityResult:
    """All geometric information for a (n, m) nanotube."""

    # Chiral indices
    n: int
    m: int

    # Lattice info
    a1: np.ndarray
    a2: np.ndarray

    # Derived vectors (Å)
    Ch_vec:  np.ndarray = field(repr=False)  # chiral vector
    T_vec:   np.ndarray = field(repr=False)  # translational vector
    t1: int = 0
    t2: int = 0

    # Scalar properties
    diameter:      float = 0.0   # Å
    Ch_norm:       float = 0.0   # |Ch|, Å
    T_norm:        float = 0.0   # |T|,  Å
    theta_deg:     float = 0.0   # chiral angle, degrees
    n_atoms:       int   = 0     # atoms per nanotube unit cell
    n_atoms_cell:  int   = 2     # atoms in the 2D unit cell (set by compute_chirality)
    strain:        float = 0.0   # periodicity strain, %

    def __post_init__(self):
        self.Ch_norm  = float(np.linalg.norm(self.Ch_vec))
        self.T_norm   = float(np.linalg.norm(self.T_vec))
        self.diameter = self.Ch_norm / math.pi
        th = math.atan2(
            np.dot(self.Ch_vec, _perp(self.a1)),
            np.dot(self.Ch_vec, self.a1)
        )
        self.theta_deg = math.degrees(th)
        # General formula: n_atoms_cell × |n·t2 − m·t1|
        self.n_atoms = self.n_atoms_cell * abs(self.n * self.t2 - self.m * self.t1)

    def __repr__(self) -> str:
        return (
            f"ChiralityResult(n={self.n}, m={self.m} | "
            f"D={self.diameter:.4f} Å | θ={self.theta_deg:.2f}° | "
            f"atoms={self.n_atoms} | strain={self.strain:.4f}%)"
        )


def _perp(v: np.ndarray) -> np.ndarray:
    """90° counter-clockwise rotation of a 2D vector."""
    return np.array([-v[1], v[0]])


# ─────────────────────────────────────────────────────────────────────────────
# Core algorithm: find best translational vector T
# ─────────────────────────────────────────────────────────────────────────────

def _exact_hexagonal_T(n: int, m: int) -> tuple[int, int]:
    """
    For hexagonal lattices, the exact (t1, t2) is given analytically.
    Uses the standard Dresselhaus formula.
    """
    from math import gcd
    d_R = gcd(2 * m + n, 2 * n + m)
    t1  =  (2 * m + n) // d_R
    t2  = -(2 * n + m) // d_R
    return t1, t2


def _search_T(
    n: int, m: int, a1: np.ndarray, a2: np.ndarray,
    limit: int = 300,
) -> tuple[int, int, float]:
    """
    Find the best integer (t1, t2) minimising |Ch · T| / (|Ch| · |T|).

    Works for all lattice types (hexagonal γ=60° and γ=120°, rectangular,
    oblique) by using the actual dot products — no index-based shortcuts.

    Returns (t1, t2, strain_fraction).
    strain_fraction = 0 means perfectly periodic (e.g. hexagonal or rectangular
    zigzag/armchair).
    """
    Ch      = n * a1 + m * a2
    Ch_norm = np.linalg.norm(Ch)

    if Ch_norm < 1e-12:
        return 0, 1, 0.0  # degenerate (n=m=0), handled upstream

    # General perpendicularity condition:
    #   (Ch · a1)·t1  +  (Ch · a2)·t2  =  0
    dot_Ch_a1 = float(np.dot(Ch, a1))
    dot_Ch_a2 = float(np.dot(Ch, a2))

    # Geometry-based shortcuts: if Ch is already perpendicular to a lattice
    # vector, that vector *is* the exact T (works for rectangular zigzag/armchair
    # and as a degenerate case of the general formula).
    if abs(dot_Ch_a1) < 1e-8:
        return 1, 0, 0.0   # T = a1 is exactly perpendicular to Ch
    if abs(dot_Ch_a2) < 1e-8:
        return 0, 1, 0.0   # T = a2 is exactly perpendicular to Ch

    target = -dot_Ch_a2 / dot_Ch_a1   # ideal t1/t2 ratio

    best_err = float("inf")
    best_t   = (1, 1)
    found    = False

    for t2 in range(1, limit + 1):
        t1 = round(target * t2)
        if t1 == 0:
            # T = t2·a2 only works if dot_Ch_a2 ≈ 0 (handled above).
            # Here it means rounding collapsed to 0; try ±1 instead.
            candidates = [(-1, t2), (1, t2)]
        else:
            T      = t1 * a1 + t2 * a2
            T_norm = np.linalg.norm(T)
            if T_norm > 1e-12:
                err = abs(np.dot(Ch, T)) / (Ch_norm * T_norm)
                if err < best_err:
                    best_err = err
                    best_t   = (t1, t2)
                if err < 1e-9:
                    found = True
                    break
            # Also explore t1 ± 1 (catches rounding boundary cases)
            candidates = [(t1 - 1, t2), (t1 + 1, t2)]

        for t1_alt, t2_alt in candidates:
            if t1_alt == 0:
                continue
            T_alt      = t1_alt * a1 + t2_alt * a2
            T_norm_alt = np.linalg.norm(T_alt)
            if T_norm_alt < 1e-12:
                continue
            err_alt = abs(np.dot(Ch, T_alt)) / (Ch_norm * T_norm_alt)
            if err_alt < best_err:
                best_err = err_alt
                best_t   = (t1_alt, t2_alt)
            if err_alt < 1e-9:
                found = True
                break   # exact solution found via alt

        if found:
            break

    # Reduce (t1, t2) by their GCD to get the primitive translation vector.
    # Without this, a later t2 iteration (e.g. t2=31 instead of 1) can give
    # the same angular error via floating-point luck on a multiple of the true
    # T, yielding a tube cell many times longer than necessary.
    t1_out, t2_out = best_t
    g = _gcd(abs(t1_out), abs(t2_out))
    if g > 1:
        t1_out //= g
        t2_out //= g
        # Recompute error for the reduced vector (should be same or smaller)
        T_red  = t1_out * a1 + t2_out * a2
        T_norm = np.linalg.norm(T_red)
        if T_norm > 1e-12:
            best_err = abs(np.dot(Ch, T_red)) / (Ch_norm * T_norm)

    return t1_out, t2_out, best_err


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def compute_chirality(
    n: int,
    m: int,
    structure: LatticeStructure,
    search_limit: int = 300,
) -> ChiralityResult | None:
    """
    Compute the full chirality description for index pair (n, m).

    Parameters
    ----------
    n, m         : chiral indices
    structure    : LatticeStructure from core.io
    search_limit : max iterations for T-vector search (higher = more accurate)

    Returns None for the degenerate (0, 0) case.
    """
    if n == 0 and m == 0:
        return None

    a1, a2 = structure.a1, structure.a2

    # Use the general search for ALL lattice types.
    # _exact_hexagonal_T (Dresselhaus formula) only works for γ=60° convention;
    # _search_T handles γ=60°, γ=120°, rectangular, and oblique correctly,
    # and returns strain≈0 for lattices with an exact perpendicular T.
    t1, t2, strain = _search_T(n, m, a1, a2, limit=search_limit)

    Ch_vec = n * a1 + m * a2
    T_vec  = t1 * a1 + t2 * a2

    return ChiralityResult(
        n=n, m=m,
        a1=a1, a2=a2,
        Ch_vec=Ch_vec, T_vec=T_vec,
        t1=t1, t2=t2,
        strain=strain * 100,          # convert to %
        n_atoms_cell=len(structure.atoms),
    )


def unique_sector_deg(structure: "LatticeStructure") -> float:
    """
    Return the maximum chiral angle (°) of the symmetry-unique sector
    for the given lattice.

    Hexagonal (a=b, γ=60°/120°)  →  30°  (zigzag 0° ↔ armchair 30°)
    Square    (a=b, γ=90°)       →  45°  (zigzag 0° ↔ armchair 45°)
    Rectangular (a≠b, γ=90°)    →  90°  (full quadrant)
    Oblique                      →  γ    (sector from a₁ direction to a₂ direction)

    For oblique lattices the two boundary directions are exactly a₁ (θ=0°) and
    a₂ (θ=γ).  All (n≥0, m≥0) integer pairs map to directions in [0°, γ], so
    the sector opening angle equals γ.  For the common case γ > 90° (e.g. the
    pza-C10 structure with γ = 102.4°) some points have negative x in Cartesian
    polar coordinates; the panel accommodates this.
    """
    lt = structure.lattice_type
    if lt == "hexagonal":
        return 30.0
    if lt == "rectangular":
        if abs(structure.a - structure.b) < 1e-3:
            return 45.0
        return 90.0
    # oblique: sector is exactly [0°, γ]
    return structure.gamma_deg


def scan_chirality(
    structure:    "LatticeStructure",
    n_max:        int = 30,
    m_max:        int | None = None,
    max_diameter: float = 25.0,
    max_atoms:    int | None = None,
    max_T_norm:   float | None = None,
    search_limit: int = 50,
    unique_only:  bool = True,
) -> list[ChiralityResult]:
    """
    Scan all (n, m) pairs and return a sorted list of ChiralityResults.

    Parameters
    ----------
    structure    : input 2D structure
    n_max, m_max : index range (m_max defaults to n_max)
    max_diameter : filter by diameter (Å)
    max_atoms    : optional upper limit on atoms/cell
    max_T_norm   : optional upper limit on translational vector length (Å);
                   useful for oblique lattices where T can be very long
    search_limit : accuracy of T-vector brute-force search.
                   Default is 50, which is sufficient for polar-map scanning
                   (keeps T vectors short and tube cells manageable).
                   Use 300 for high-accuracy single-tube computation.
    unique_only  : if True (default), skip (n,m) pairs that are equivalent
                   to (m,n) by the lattice symmetry (applies to hexagonal
                   and square lattices where |a₁|=|a₂|).
    """
    if m_max is None:
        m_max = n_max

    # Determine the unique-sector angular cutoff.
    # For symmetric lattices (hexagonal, square) we filter by the chiral angle
    # of Ch = n·a₁ + m·a₂ rather than by an index heuristic.
    # This is convention-agnostic: it works for both γ=60° and γ=120°.
    #   γ=60°  hexagonal : armchair at (n,n)  → θ=30° (same n→2n rule)
    #   γ=120° hexagonal : armchair at (2n,n) → θ=30° (m>n gives θ>30°, excluded)
    #   square (γ=90°)   : armchair at (n,n)  → θ=45°
    _lt        = structure.lattice_type
    _symmetric = unique_only and (
        _lt == "hexagonal"
        or (_lt == "rectangular" and abs(structure.a - structure.b) < 1e-3)
    )
    a1, a2 = structure.a1, structure.a2
    if _symmetric:
        # Compute the actual armchair boundary angle from the lattice vectors
        # (a1 + a2 is always the armchair direction for hexagonal/square).
        # This avoids false rejections when γ deviates slightly from the ideal
        # value (e.g. 60.0001° instead of 60.0000°).
        arm_vec   = a1 + a2
        _theta_max = math.degrees(math.atan2(arm_vec[1], arm_vec[0]))
    else:
        _theta_max = 90.0

    results = []
    for n in range(n_max + 1):
        for m in range(m_max + 1):
            # Skip pairs outside the symmetry-unique angular sector.
            # We compute the actual chiral angle from the lattice vectors so
            # the cutoff is correct regardless of γ convention.
            if _symmetric:
                Ch_test = n * a1 + m * a2
                theta_test = math.degrees(math.atan2(Ch_test[1], Ch_test[0]))
                if theta_test > _theta_max + 1e-6 or theta_test < -1e-6:
                    continue
            res = compute_chirality(n, m, structure,
                                    search_limit=search_limit)
            if res is None:
                continue
            if res.diameter > max_diameter:
                continue
            if max_atoms is not None and res.n_atoms > max_atoms:
                continue
            if max_T_norm is not None and res.T_norm > max_T_norm:
                continue
            results.append(res)

    results.sort(key=lambda r: (r.diameter, r.theta_deg))
    return results
