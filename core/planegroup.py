"""
core/planegroup.py
------------------
The symmetry of a 2D crystal, computed rather than classified by cases.

Everything the nanotube construction needs from symmetry splits in two:

* **What the tube is** — Ch, T, the diameter, whether an exact perpendicular
  translation exists — depends only on the metric (|a1|, |a2|, gamma), that is
  on the Bravais class.  The basis never enters.
* **How many distinct tubes there are** — the fundamental domain of the (n, m)
  map — depends on the point group of the *structure* in its setting, so it
  ranges over all 17 plane groups.

This module computes the second one directly instead of testing for the four
lattice names one at a time.  The route is standard crystallography:

1. Lagrange-reduce the basis, so every lattice automorphism has small integer
   entries and a bounded search is exhaustive.
2. Enumerate the integer matrices ``U`` with ``det U = +-1`` and
   ``U G U^T = G`` (``G`` the Gram matrix).  Those are the automorphisms of the
   *lattice* — its holohedry.
3. Keep the ``U`` for which the *basis* also maps onto itself, allowing an
   arbitrary fractional translation.  Allowing that translation is what makes
   a glide count: pm and pg differ by half a lattice vector, and half a lattice
   vector is only a shift of the origin of the sheet before it is rolled, so
   the two enumerate identically.

Index vectors and fractional coordinates transform by the *same* integer
matrix, which is what makes the whole thing cheap.  With rows ``a1, a2`` in
``A`` and a row vector ``v = (n, m)``,

    Ch = v A,     v -> v U,     f -> f U.

Why this replaces the old heuristics
------------------------------------
``lattice_type`` knows hexagonal, rectangular and oblique, so a rhombic cell
(|a1| = |a2|, gamma not 60/90/120) — a centred rectangular lattice, the fifth
Bravais class — was called oblique and its mirror was lost.  And the old
``basis_swap_invariant`` matched each atom's image modulo a *lattice*
translation, so it saw a mirror through the origin and missed a glide.  Both
errors were conservative: the map came out too large, never wrong.  A third
error was not: for a lattice whose only automorphism is +-I, scanning
``n, m >= 0`` covers directions [0, gamma] only, and every tube with m < 0 was
missing.  Orbit reduction over the computed group fixes all three at once.
"""
from __future__ import annotations

import numpy as np

from core.io import LatticeStructure

# Entries of an automorphism of a Lagrange-reduced 2D basis lie in {-1, 0, 1}:
# a1 must map to a vector of its own length, and in a reduced basis those are
# among +-a1, +-a2, +-(a1 +- a2).  The search runs to 2 anyway -- it costs 625
# matrix products -- so a basis that reduction leaves slightly off ideal, as a
# noisy CIF does, cannot silently lose an operation.
_SEARCH = 2

# Fractional-coordinate tolerance.  Loose enough for experimental CIFs, tight
# enough that a genuine 1/3 is not mistaken for a 1/2.
_TOL = 2e-3


# ─────────────────────────────────────────────────────────────────────────────
# Lattice
# ─────────────────────────────────────────────────────────────────────────────

def gram(structure: LatticeStructure) -> np.ndarray:
    """Metric tensor G = A A^T of the 2D cell."""
    A = np.array([structure.a1, structure.a2], dtype=float)
    return A @ A.T


def reduce_basis(G: np.ndarray) -> np.ndarray:
    """
    Lagrange--Gauss reduction.

    Returns the integer matrix ``R`` (det = +-1) for which ``R G R^T`` is
    reduced, i.e. ``|a1| <= |a2|`` and ``|a1.a2| <= |a1|^2 / 2``.

    Reduction is what makes the automorphism search finite.  Cells in the wild
    are not reduced -- one of the structures shipped with NTBuilder has
    gamma = 164.13 degrees -- and in such a basis an automorphism can have
    large entries.
    """
    R = np.eye(2, dtype=np.int64)
    for _ in range(100):                      # terminates in a few steps
        g = R @ G @ R.T
        if g[0, 0] > g[1, 1]:
            R = R[::-1].copy()
            g = R @ G @ R.T
        if g[0, 0] <= 1e-12:
            break
        mu = int(round(g[0, 1] / g[0, 0]))
        if mu == 0:
            break
        R[1] -= mu * R[0]
    return R


def lattice_point_group(structure: LatticeStructure,
                        tol: float = 1e-6) -> list[np.ndarray]:
    """
    The holohedry: every integer ``U`` with ``det U = +-1`` and
    ``U G U^T = G``, expressed in the structure's own basis.

    Its order is 2, 4, 8 or 12 for the oblique, rectangular (primitive or
    centred), square and hexagonal lattices respectively.
    """
    G = gram(structure)
    R = reduce_basis(G)
    Gr = R @ G @ R.T
    scale = float(np.trace(Gr)) or 1.0
    Rinv = np.round(np.linalg.inv(R.astype(float))).astype(np.int64)

    ops = []
    rng = range(-_SEARCH, _SEARCH + 1)
    for a in rng:
        for b in rng:
            for c in rng:
                for d in rng:
                    U = np.array([[a, b], [c, d]], dtype=np.int64)
                    if abs(a * d - b * c) != 1:
                        continue
                    if np.max(np.abs(U @ Gr @ U.T - Gr)) > tol * scale:
                        continue
                    ops.append(Rinv @ U @ R)
    return ops


# ─────────────────────────────────────────────────────────────────────────────
# Structure
# ─────────────────────────────────────────────────────────────────────────────

def _fractional(structure: LatticeStructure) -> tuple[np.ndarray, list, list]:
    A = np.array([structure.a1, structure.a2], dtype=float)
    pos = np.array([a["pos"] for a in structure.atoms], dtype=float)
    frac = np.linalg.solve(A.T, pos.T).T % 1.0
    syms = [a["symbol"] for a in structure.atoms]
    zs = [float(a.get("z", 0.0)) for a in structure.atoms]
    return frac, syms, zs


def _maps_onto_itself(U: np.ndarray, frac: np.ndarray, syms: list,
                      zs: list, tol: float) -> bool:
    """
    Is the basis invariant under ``U`` followed by *some* translation?

    The translation is a free parameter, not a lattice vector, which is exactly
    what lets a glide (and a mirror that does not pass through the chosen
    origin) be recognised.  Candidates are the offsets that map atom 0 onto a
    compatible atom, so at most N of them.

    The out-of-plane offset must be preserved, not merely matched: an operation
    that maps the sheet onto itself only after flipping z relates the tube to
    its radially inverted partner, and for a Janus sheet like MoSSe those are
    different structures.
    """
    img = frac @ U
    for j in range(len(frac)):
        if syms[j] != syms[0] or abs(zs[j] - zs[0]) > 1e-3:
            continue
        t = frac[j] - img[0]
        shifted = (img + t) % 1.0
        ok = True
        for i in range(len(frac)):
            d = np.abs((shifted[i] - frac + 0.5) % 1.0 - 0.5).max(axis=1)
            hit = [k for k in np.where(d < tol)[0]
                   if syms[k] == syms[i] and abs(zs[k] - zs[i]) <= 1e-3]
            if not hit:
                ok = False
                break
        if ok:
            return True
    return False


def structure_point_group(structure: LatticeStructure,
                          tol: float = _TOL) -> list[np.ndarray]:
    """
    The point group of the crystal: the subgroup of the holohedry that the
    atomic basis also respects, up to a translation.

    Always a subgroup of :func:`lattice_point_group`, which is the
    crystallographic statement that a structure can only lose symmetry
    relative to its lattice, never gain it.  Penta-graphene is the standard
    example: a square lattice (holohedry 4mm) carrying an sp3 basis that
    breaks every mirror, so the crystal is 4 and (5,0) is not (0,5).
    """
    frac, syms, zs = _fractional(structure)
    if len(frac) == 0:
        return lattice_point_group(structure)
    return [U for U in lattice_point_group(structure)
            if _maps_onto_itself(U, frac, syms, zs, tol)]


# ─────────────────────────────────────────────────────────────────────────────
# The group that acts on the chiral indices
# ─────────────────────────────────────────────────────────────────────────────

def chirality_group(structure: LatticeStructure,
                    tol: float = _TOL) -> list[np.ndarray]:
    """
    The group by which two (n, m) pairs may be identified.

    It is the structure's point group together with ``-I``.  Inversion is
    included whether or not the crystal has it, because it is a property of
    rolling and not of the sheet: Ch and -Ch wrap the same atoms around the
    same cylinder, so they can never denote different tubes.
    """
    ops = structure_point_group(structure, tol)
    minus = -np.eye(2, dtype=np.int64)
    out = list(ops)
    for U in ops:
        V = minus @ U
        if not any(np.array_equal(V, W) for W in out):
            out.append(V)
    return out


def unique_indices(structure: LatticeStructure, n_max: int,
                   m_max: int | None = None,
                   tol: float = _TOL) -> list[tuple[int, int]]:
    """
    One representative of every symmetry-distinct (n, m) within the box.

    The domain scanned is the half plane ``n > 0`` plus ``(0, m > 0)``, not the
    first quadrant.  On a lattice whose only automorphisms are +-I -- any truly
    oblique one -- the quadrant spans chiral directions [0, gamma] and leaves
    everything between gamma and 180 degrees unreachable, so tubes with m < 0
    simply do not exist for the user.  Reducing the half plane by the computed
    group gives the quadrant back automatically whenever the symmetry is there
    to justify it.
    """
    if m_max is None:
        m_max = n_max
    ops = chirality_group(structure, tol)
    in_box = lambda v: (0 <= v[0] <= n_max and abs(v[1]) <= m_max
                        and not (v[0] == 0 and v[1] <= 0))

    seen: set[tuple[int, int]] = set()
    out: list[tuple[int, int]] = []
    for n in range(0, n_max + 1):
        for m in range(-m_max, m_max + 1):
            if not in_box((n, m)) or (n, m) in seen:
                continue
            orbit = {(int(o[0]), int(o[1]))
                     for o in (np.array((n, m), dtype=np.int64) @ U for U in ops)}
            seen |= orbit
            # Report the conventional member of the orbit.  The literature
            # names tubes in the first quadrant with n >= m, so prefer m >= 0
            # first and the largest n after that; (5, 0) rather than (0, 5),
            # and (1, 2) rather than the (2, -1) a 4-fold rotation also allows.
            # Only an orbit with no first-quadrant member at all -- which
            # happens exactly when the lattice is too poor in symmetry to
            # bring it there -- falls back to a negative m.
            cands = [v for v in orbit if in_box(v)] or [(n, m)]
            out.append(max(cands, key=lambda v: (v[1] >= 0, v[0], v[1])))
    return sorted(out)


def sector_deg(structure: LatticeStructure, tol: float = _TOL) -> float:
    """
    Opening angle, in degrees, of the symmetry-unique wedge of chiral angles.

    A finite group acting on the circle of directions has a fundamental domain
    of 360 / |G| degrees, so the wedge is read off the group order rather than
    tabulated per lattice.  Graphene gives 30 (order 12, 6mm), MoS2 also 30
    (3m, order 6, plus inversion from the rolling), penta-graphene and every
    primitive rectangular sheet 90, and a bare oblique one 180.

    This is an opening angle only.  Where the wedge begins depends on how the
    mirrors sit relative to a1, which is what separates p3m1 from p31m; use
    :func:`unique_indices` when the actual set of tubes is wanted.
    """
    return 360.0 / len(chirality_group(structure, tol))
