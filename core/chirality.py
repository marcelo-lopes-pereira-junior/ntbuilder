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
  best integer approximation and computes a *strain* metric — the
  fractional angular residual — that quantifies the periodicity error.
  This is a key scientific contribution of the tool.

  Perpendicularity fixes t1 once t2 is chosen (the best integer is
  round(x·t2) with x = -(Ch·a2)/(Ch·a1)), so choosing T is choosing a
  DENOMINATOR — the classical problem of approximating x by a rational of
  bounded denominator. Its record holders are the continued-fraction
  convergents of x and their largest fitting semiconvergents, so the search
  visits O(log limit) candidates instead of sweeping all of them. Searching
  (t1, t2) as a pair adds nothing: on biphenylene (4,1) at limit 2000 the
  exhaustive 2-D search over 1 678 000 pairs returns the same t = (-320,
  1543) as the 6000 candidates of the 1-D form, 269 times slower.

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
        # Folded into [0, 180).  Rolling makes Ch and -Ch the same tube, so the
        # group that acts on the indices always contains -I and the wedge of
        # distinct chiral directions is 180 degrees wide -- but on a lattice
        # too poor in symmetry to bring every orbit into the first quadrant the
        # representative comes back with m < 0, and its raw angle is negative.
        # Left negative, 41 % of the points of an oblique polar map (measured
        # on the reduced AgBr3 cell: 90 of 221, down to -78.8 degrees) fell
        # below the axis and were clipped, while the wedge drawn from 0 to 180
        # sat empty above gamma.  The fold reports the direction of (-n, -m),
        # which is the same tube.
        self.theta_deg = math.degrees(th) % 180.0
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


def _search_T_scan(
    n: int, m: int, a1: np.ndarray, a2: np.ndarray,
    limit: int = 300,
) -> tuple[int, int, float]:
    """
    Reference implementation: sweep every t2 from 1 to *limit*.

    Kept because it is obviously correct and because the fast path is checked
    against it.  It is O(limit) and visits ~15 useful points in 45 369 steps
    on biphenylene (4,1); ``_search_T`` reaches the same answer in O(log
    limit).  Not called in normal operation.

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
# Fast path: the candidates a bounded search can possibly win at
# ─────────────────────────────────────────────────────────────────────────────

def _bounded_denominators(x: float, limit: int) -> list[int]:
    """Every t2 <= limit at which the approximation of x can set a new record.

    Perpendicularity fixes t1 for a given t2 -- the best integer is
    round(x*t2) -- so choosing T is choosing a denominator, and since |T| grows
    with t2 the criterion |Ch.T|/(|Ch||T|) is, to leading order, |x - t1/t2|.
    Its record holders are the convergents of the continued fraction of x plus,
    at each level, the semiconvergent with the largest multiplier that still
    fits under the bound; nothing else can win.  A sweep over every t2 spends
    O(limit) steps to visit O(log limit) useful ones: on biphenylene (4,1) the
    sweep's 15 record holders are exactly [1, 4, 5, 24, 29, 53, 82, 135, 352,
    1191, 1543, 5820, 13183, 32186, 45369], and this recursion lists them in
    15 steps.

    a0 = floor(x) belongs to the NUMERATOR: the denominator recursion starts at
    a1, from the reciprocal of the fractional part.  Started at a0 the list
    came out wrong and lost to the sweep in 984 of 1440 test cases.

    The recursion runs on a float and is stopped when the remainder falls to
    float noise -- past that point the digits are not x's any more.  It cannot
    claim rationality, which is not decidable in floating point, only that a
    denominator was reached within the limit.
    """
    out: list[int] = [1]                      # q_0 = 1 is always a candidate
    x = abs(float(x))
    frac = x - math.floor(x)
    if frac <= 1e-13:
        return [1]
    x = 1.0 / frac
    q_prev, q = 0, 1                          # q_{-1}, q_0
    for _ in range(64):
        a = math.floor(x)
        if a < 1:
            break
        # The largest semiconvergent of this level that fits.  Between two
        # convergents these mediants are the only other possible records.
        room = (limit - q_prev) // q if q else 0
        j = int(min(a, room))
        if j >= 1:
            out.append(int(j * q + q_prev))
        q_next = int(a) * q + q_prev
        if q_next > limit:
            break
        out.append(int(q_next))
        frac = x - a
        if frac <= 1e-13:                     # exact hit, or float noise
            break
        q_prev, q = q, q_next
        x = 1.0 / frac
    return sorted({t for t in out if 1 <= t <= limit})


def T_options(
    n: int, m: int, a1: np.ndarray, a2: np.ndarray,
    limit: int = 300,
) -> list[dict]:
    """The trade-off front between cell length and periodicity residual.

    Every candidate that can hold the record, cheapest cell first.  Each entry
    is ``{t1, t2, T_norm, strain}`` with strain in per cent, and the list is a
    Pareto front: |T| grows down it and the residual falls.

    This is not an extra computation -- it is the search itself.  Only the
    continued-fraction convergents of the ideal ratio can improve on their
    predecessors, so the candidate list IS the front, with O(log limit)
    entries instead of thousands.  On the AgBr3 monolayer's (4, 1) it runs
    from 112 atoms at 64 % down to 12 208 atoms at 0.0004 %, and a user who
    wants a cell small enough for DFT picks a row rather than hunting for a
    search_limit that happens to produce it.
    """
    Ch = n * a1 + m * a2
    Ch_norm = float(np.linalg.norm(Ch))
    if Ch_norm < 1e-12:
        return [{"t1": 0, "t2": 1, "T_norm": float(np.linalg.norm(a2)),
                 "strain": 0.0}]
    dot_Ch_a1 = float(np.dot(Ch, a1))
    dot_Ch_a2 = float(np.dot(Ch, a2))
    if abs(dot_Ch_a1) < 1e-8:
        return [{"t1": 1, "t2": 0, "T_norm": float(np.linalg.norm(a1)),
                 "strain": 0.0}]
    if abs(dot_Ch_a2) < 1e-8:
        return [{"t1": 0, "t2": 1, "T_norm": float(np.linalg.norm(a2)),
                 "strain": 0.0}]

    target = -dot_Ch_a2 / dot_Ch_a1
    out: list[dict] = []
    best = float("inf")
    for t2 in _bounded_denominators(target, limit):
        t1_ideal = round(target * t2)
        for t1 in (t1_ideal - 1, t1_ideal, t1_ideal + 1):
            if t1 == 0:
                continue
            g = _gcd(abs(t1), abs(t2))
            t1r, t2r = (t1 // g, t2 // g) if g > 1 else (t1, t2)
            T = t1r * a1 + t2r * a2
            T_norm = float(np.linalg.norm(T))
            if T_norm < 1e-12:
                continue
            err = abs(float(np.dot(Ch, T))) / (Ch_norm * T_norm)
            if err < best - 1e-15:
                best = err
                out.append({"t1": int(t1r), "t2": int(t2r),
                            "T_norm": T_norm, "strain": err * 100.0})
        if best < 1e-12:
            break
    return out


def _search_T(
    n: int, m: int, a1: np.ndarray, a2: np.ndarray,
    limit: int = 300,
    max_strain: float | None = None,
    max_T_norm: float | None = None,
) -> tuple[int, int, float]:
    """Best integer (t1, t2) minimising |Ch . T| / (|Ch| . |T|).

    With ``max_strain`` (per cent) the choice inverts: the SHORTEST cell whose
    residual is within tolerance, rather than the smallest residual whatever
    the length.  ``max_T_norm`` caps the length in angstroms instead.  Both
    read the front from :func:`T_options`; neither is on by default, so the
    unconstrained answer is unchanged.

    Same criterion and same tie-breaking as :func:`_search_T_scan`, evaluated
    only where a record is possible.  Searching (t1, t2) as a pair is not
    smarter: measured on biphenylene (4,1) with limit 2000, the exhaustive
    2-D search over 1 678 000 pairs returns the very same t = (-320, 1543) as
    the 6000 candidates of the 1-D form, 269 times slower, because for each
    t2 the optimal t1 is forced to round(x*t2) and only +-1 around it can
    compete once |T| enters the denominator.
    """
    Ch      = n * a1 + m * a2
    Ch_norm = float(np.linalg.norm(Ch))

    if Ch_norm < 1e-12:
        return 0, 1, 0.0

    dot_Ch_a1 = float(np.dot(Ch, a1))
    dot_Ch_a2 = float(np.dot(Ch, a2))

    # Ch already perpendicular to a lattice vector: that vector IS T.  This is
    # the axial case -- (n,0) and (0,m) on any rectangular cell -- and it costs
    # one primitive cell, whatever the cell's ratio.
    if abs(dot_Ch_a1) < 1e-8:
        return 1, 0, 0.0
    if abs(dot_Ch_a2) < 1e-8:
        return 0, 1, 0.0

    if max_strain is not None or max_T_norm is not None:
        front = T_options(n, m, a1, a2, limit)
        pick = None
        if max_strain is not None:
            # Cheapest cell that is good enough.  The front is ordered by
            # length, so the first hit is the shortest.
            pick = next((c for c in front if c["strain"] <= max_strain), None)
        if pick is None and max_T_norm is not None:
            fits = [c for c in front if c["T_norm"] <= max_T_norm]
            pick = fits[-1] if fits else None
        if pick is None:
            pick = front[-1]          # nothing satisfies it: the best there is
        return pick["t1"], pick["t2"], pick["strain"] / 100.0

    target = -dot_Ch_a2 / dot_Ch_a1

    best_err = float("inf")
    best_t = (1, 1)
    for t2 in _bounded_denominators(target, limit):
        t1_ideal = round(target * t2)
        for t1 in (t1_ideal - 1, t1_ideal, t1_ideal + 1):
            if t1 == 0:
                continue
            T = t1 * a1 + t2 * a2
            T_norm = float(np.linalg.norm(T))
            if T_norm < 1e-12:
                continue
            err = abs(float(np.dot(Ch, T))) / (Ch_norm * T_norm)
            if err < best_err:
                best_err = err
                best_t = (t1, t2)
        if best_err < 1e-12:
            break

    # Reduce by the GCD: a multiple of the true T scores the same angle and
    # would hand back a cell several periods long.
    t1_out, t2_out = best_t
    g = _gcd(abs(t1_out), abs(t2_out))
    if g > 1:
        t1_out //= g
        t2_out //= g
        T_red = t1_out * a1 + t2_out * a2
        T_norm = float(np.linalg.norm(T_red))
        if T_norm > 1e-12:
            best_err = abs(float(np.dot(Ch, T_red))) / (Ch_norm * T_norm)

    return t1_out, t2_out, best_err


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def compute_chirality(
    n: int,
    m: int,
    structure: LatticeStructure,
    search_limit: int = 300,
    max_strain: float | None = None,
    max_T_norm: float | None = None,
) -> ChiralityResult | None:
    """
    Compute the full chirality description for index pair (n, m).

    Parameters
    ----------
    n, m         : chiral indices
    structure    : LatticeStructure from core.io
    max_strain   : per cent. Trade exactness for a shorter cell: return the
                   SHORTEST T whose periodicity residual is within this
                   tolerance instead of the most exact one. On AgBr3 (4,1),
                   0.1 %% buys 432 atoms in a 53.7 A cell where the exact
                   answer needs 12 208 atoms in 1518 A.
    max_T_norm   : angstroms. Cap the cell length instead, taking the best
                   residual that fits.
    search_limit : bound on |t2|, the denominator of the ratio t1/t2 that
                   approximates perpendicularity. Not an iteration count: the
                   search evaluates only the O(log search_limit) candidates
                   where the residual can improve. An exact T is reached when
                   this bound reaches the reduced denominator of the ratio —
                   45 369 for biphenylene (4,1).

    Returns None for the degenerate (0, 0) case.
    """
    if n == 0 and m == 0:
        return None

    a1, a2 = structure.a1, structure.a2

    # Use the general search for ALL lattice types.
    # _exact_hexagonal_T (Dresselhaus formula) only works for γ=60° convention;
    # _search_T handles γ=60°, γ=120°, rectangular, and oblique correctly,
    # and returns strain≈0 for lattices with an exact perpendicular T.
    t1, t2, strain = _search_T(n, m, a1, a2, limit=search_limit,
                               max_strain=max_strain, max_T_norm=max_T_norm)

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
    # Delegated to core.planegroup, which derives the wedge from the order of
    # the group that actually acts on the indices instead of tabulating one
    # answer per lattice name.  The old table had no row for the centred
    # rectangular lattice and returned gamma for it, the oblique answer.
    from core.planegroup import sector_deg
    return sector_deg(structure)


def basis_swap_invariant(structure: "LatticeStructure", tol: float = 2e-3) -> bool:
    """
    Is the atomic basis invariant under exchanging the two lattice axes?

    Folding the map by identifying (n, m) with (m, n) is only legitimate when
    the *whole crystal* — lattice **and** basis — is symmetric under the mirror
    that swaps a₁ and a₂.  Equal lattice constants are necessary but not
    sufficient: the basis can break the mirror even on a perfectly square or
    hexagonal lattice.

    Penta-graphene is the canonical counter-example.  Its lattice is square
    (a = b = 3.63 Å, γ = 90°), but the sp³ carbons at z = ±0.687 Å sit at
    fractional positions such as (0.634, 0.134) whose mirror image (0.134,
    0.634) is occupied by an atom at the *opposite* height.  The (5,0) and
    (0,5) nanotubes are therefore genuinely inequivalent structures, and
    folding the map would hide half of them from the user.

    An atom at fractional (u, v) with out-of-plane offset z must be matched by
    an atom of the same species at (v, u) with the same z, modulo one lattice
    translation.

    Returns True when the basis is invariant (folding is safe).
    """
    from core.planegroup import structure_point_group
    swap = np.array([[0, 1], [1, 0]], dtype=np.int64)
    return any(np.array_equal(U, swap)
               for U in structure_point_group(structure, tol))


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

    # Which (n, m) are distinct is decided by the orbits of the group that acts
    # on the indices -- the structure's point group plus inversion -- rather
    # than by an angular cutoff read off the lattice name.  That is what makes
    # every one of the 17 plane groups come out right, and it also widens the
    # scan from the first quadrant to the half plane: on a lattice whose only
    # automorphisms are +-I the quadrant spans chiral directions [0, gamma]
    # and hides every tube with m < 0.  Where the symmetry does justify the
    # quadrant, orbit reduction returns it unchanged.
    from core.planegroup import unique_indices
    a1, a2 = structure.a1, structure.a2
    if unique_only:
        pairs = unique_indices(structure, n_max, m_max)
    else:
        pairs = [(n, m) for n in range(n_max + 1) for m in range(m_max + 1)
                 if n or m]

    results = []
    for n, m in pairs:
        # Cheap diameter pre-filter.  |Ch|/π is *exactly* the diameter that
        # ChiralityResult computes, so discarding out-of-range pairs here is
        # identical to the ``res.diameter > max_diameter`` check below — but it
        # skips the expensive T-vector search for the vast majority of pairs
        # when n_max is large (the polar map at n_max=100 scans ~10 k pairs and
        # keeps only a few hundred).
        D_cheap = float(np.linalg.norm(n * a1 + m * a2)) / math.pi
        if D_cheap > max_diameter:
            continue
        res = compute_chirality(n, m, structure, search_limit=search_limit)
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
