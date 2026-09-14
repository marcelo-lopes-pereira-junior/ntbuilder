"""
core/builder.py
---------------
Constructs the 3D nanotube from a LatticeStructure and a ChiralityResult.

Algorithm
---------
1. Enumerate the lattice sites that fall in one nanotube unit cell: the
   parallelogram spanned by Ch and T in the sheet.
2. Give each site its fractional coordinates (s, w) in the (Ch, T) basis and
   keep those with 0 ≤ s < 1 and 0 ≤ w < 1; set u = s·|Ch| and v = w·|T|.
   For an exact T (perpendicular to Ch) these are the plain projections onto
   the two axes; for an approximate T they are not, and projections would cut
   the wrong region (see _iter_atoms).
4. Roll: map the circumferential coordinate u → (x, y) on a cylinder.
   • For flat structures: radius R = |Ch| / 2π (all atoms at same r).
   • For buckled / multi-layer structures: each atom's radial distance is
         r_atom = R + sign * z_offset
     where z_offset is the out-of-plane displacement stored in the structure,
     and sign = +1 (outward) or −1 (inward) depending on roll_inward.
   The angle is always θ = u / R (circumferential coordinate based on Ch).
5. The axial coordinate v becomes z.

Rolling direction for buckled structures
-----------------------------------------
Buckled 2D materials (pentagraphene, silicene, …) and layered materials
(MoS₂, MoSSe, …) have atoms at different z-offsets in the 2D plane.
When rolling:
  roll_inward=False  (default, "outward"):
      atoms with z>0 are placed at larger r (outer wall)
  roll_inward=True ("inward"):
      atoms with z>0 are placed at smaller r (inner wall)
For Janus materials (MoSSe), the two options produce chemically
distinct nanotubes (S-outer vs. Se-outer).

Bond validation
---------------
check_spurious_bonds(structure, nt) compares the set of bonded species
pairs in the flat 2D structure with those in the 3D nanotube.  Any pair
that appears only in the 3D structure signals an unphysical bond caused
by the curvature (typically for small-diameter tubes of multi-layer materials).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterator

import numpy as np

from .io import LatticeStructure
from .chirality import ChiralityResult


# ─────────────────────────────────────────────────────────────────────────────
# Result container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class NanotubeStructure:
    """3D nanotube ready for export."""

    chirality: ChiralityResult
    symbols:   list[str]
    coords:    np.ndarray   # shape (N, 3), Cartesian Å
    box:       np.ndarray   # shape (3,)  — [Lx, Ly, Lz]
    vacuum:    float = 10.0

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)

    @property
    def diameter(self) -> float:
        return self.chirality.diameter

    @property
    def length(self) -> float:
        return float(self.box[2])

    def __repr__(self) -> str:
        r = self.chirality
        return (
            f"NanotubeStructure(({r.n},{r.m}) | "
            f"D={self.diameter:.4f} Å | L={self.length:.4f} Å | "
            f"atoms={self.n_atoms})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Builder
# ─────────────────────────────────────────────────────────────────────────────

_TOL = 1e-4   # Å.  A site this close below a cell boundary counts as ON it.
              # The SAME shift at both ends of the half-open window keeps
              # exactly one image of every site.  The previous asymmetric pair
              # (1e-4 below, 5e-4 above) dropped any site that landed between
              # the two, together with its image.


def _index_range(k: int, c: float, lo: float, hi: float,
                 fallback: tuple[int, int]) -> tuple[int, int]:
    """Integers i with lo <= (k*i + c) < hi as a padded half-open range.

    ``k`` is an integer coefficient; when it is 0 the condition does not
    depend on i and the row is either empty or given by ``fallback``.  The
    range is padded by one on each side -- the caller checks every i exactly.
    """
    if k == 0:
        return fallback if lo <= c < hi else (0, 0)
    a, b = (lo - c) / k, (hi - c) / k
    if k < 0:
        a, b = b, a
    return math.floor(a) - 1, math.ceil(b) + 2


def _iter_atoms(
    structure: LatticeStructure,
    chirality: ChiralityResult,
) -> Iterator[tuple[str, float, float, float]]:
    """Yield (symbol, u, v, z_offset) for every atom of one nanotube unit cell."""
    Ch_norm, T_norm = chirality.Ch_norm, chirality.T_norm
    for _k, sym, s, w, _x, _y, z_off in _iter_sites(structure, chirality):
        yield sym, s * Ch_norm, w * T_norm, z_off


def _iter_sites(
    structure: LatticeStructure,
    chirality: ChiralityResult,
) -> Iterator[tuple[int, str, float, float, float, float, float]]:
    """
    Yield (k, symbol, s, w, x, y, z_offset) for every atom of one nanotube
    unit cell: k is the atom of the 2D cell it comes from and (x, y) its
    lattice-index coordinates in the sheet, so each tube atom keeps its
    identity in the flat layer (what check_curvature_bonds compares).

    The cell is the parallelogram spanned by Ch and T in the sheet.  A site
    belongs to it when its fractional coordinates (s, w) in the (Ch, T) basis
    lie in [0, 1), and u = s*|Ch|, v = w*|T| are returned, so the roll sends
    Ch onto the whole circumference and T onto the axial period exactly.

    For an exact T, perpendicular to Ch, this is the rectangle 0 <= u < |Ch|,
    0 <= v < |T| of plain projections, which is how the cell used to be cut.
    For an approximate T it is not: projections onto two oblique axes cut a
    region 1/sin^2(Ch, T) times the cell, so sites were repeated -- 184 atoms
    instead of 6, down to 0.13 A apart, on biphenylene (5,6) with t = (1,1) --
    and the tube had a seam.  Fractional coordinates hold exactly
    n_atoms_cell * |n*t2 - m*t1| atoms, each once, with no seam; the cost of
    the approximation becomes a uniform shear of the sheet, the residual that
    is reported as strain.

    With x = i + f1 and y = j + f2 the lattice-index coordinates of a site,
        s = (x*t2 - y*t1) / D,   w = (n*y - m*x) / D,   D = n*t2 - m*t1,
    both linear in i, so each row j gets analytic bounds on i and the
    supercell is never materialised.
    """
    a1, a2 = structure.a1, structure.a2
    n, m = chirality.n, chirality.m
    t1, t2 = chirality.t1, chirality.t2
    D = n * t2 - m * t1
    if D == 0:
        raise ValueError(
            f"T = ({t1},{t2}) is parallel to Ch = ({n},{m}): no nanotube cell.")
    Ch_norm = chirality.Ch_norm
    T_norm  = chirality.T_norm
    lo_s, hi_s = -_TOL / Ch_norm, 1.0 - _TOL / Ch_norm
    lo_w, hi_w = -_TOL / T_norm,  1.0 - _TOL / T_norm

    # The cell's corners in index space bound the rows and columns to visit.
    xs = (0, n, t1, n + t1)
    ys = (0, m, t2, m + t2)
    x_range = (min(xs) - 2, max(xs) + 3)
    j_lo, j_hi = min(ys) - 2, max(ys) + 2

    basis = np.column_stack([a1, a2])

    for k, atom in enumerate(structure.atoms):
        sym   = atom["symbol"]
        z_off = atom.get("z", 0.0)
        f1, f2 = (float(f) for f in np.linalg.lstsq(
            basis, np.asarray(atom["pos"], dtype=float), rcond=None)[0])

        for j in range(j_lo, j_hi + 1):
            y = j + f2
            # s = (t2*i + cs)/D and w = (-m*i + cw)/D
            cs = f1 * t2 - y * t1
            cw = n * y - m * f1
            if D > 0:
                r_s = _index_range(t2, cs, lo_s * D, hi_s * D, x_range)
                r_w = _index_range(-m, cw, lo_w * D, hi_w * D, x_range)
            else:
                # Dividing by a negative D flips both inequalities.
                r_s = _index_range(-t2, -cs, lo_s * -D, hi_s * -D, x_range)
                r_w = _index_range(m, -cw, lo_w * -D, hi_w * -D, x_range)

            for i in range(max(r_s[0], r_w[0], x_range[0]),
                           min(r_s[1], r_w[1], x_range[1])):
                x = i + f1
                s = (x * t2 - y * t1) / D
                w = (n * y - m * x) / D
                if lo_s <= s < hi_s and lo_w <= w < hi_w:
                    yield k, sym, s, w, x, y, z_off


def build_nanotube(
    structure: LatticeStructure,
    chirality: ChiralityResult,
    vacuum: float = 10.0,
    roll_inward: bool = False,
) -> NanotubeStructure:
    """
    Build a 3D nanotube from a 2D structure and its chirality description.

    Parameters
    ----------
    structure   : 2D unit cell (LatticeStructure)
    chirality   : output of core.chirality.compute_chirality()
    vacuum      : vacuum padding around the nanotube in the xy plane (Å)
    roll_inward : if True, atoms with positive z-offset are placed at
                  smaller radius (inner wall); default False places them
                  at larger radius (outer wall).

    Returns
    -------
    NanotubeStructure with Cartesian coordinates centred in the simulation box.
    """
    Ch_norm   = chirality.Ch_norm
    T_norm    = chirality.T_norm
    radius    = Ch_norm / (2.0 * np.pi)   # mean rolling radius

    # +1 → z>0 atoms go outward; -1 → z>0 atoms go inward
    roll_sign = -1.0 if roll_inward else +1.0

    sym_out: list[str]   = []
    x_out:   list[float] = []
    y_out:   list[float] = []
    z_out:   list[float] = []

    for sym, u, v, z_off in _iter_atoms(structure, chirality):
        r_atom = radius + roll_sign * z_off
        angle  = u / radius          # circumferential → azimuthal angle
        sym_out.append(sym)
        x_out.append(r_atom * math.cos(angle))
        y_out.append(r_atom * math.sin(angle))
        z_out.append(v)

    if not sym_out:
        raise RuntimeError(
            f"No atoms selected for ({chirality.n},{chirality.m}). "
            "Check that the structure and chirality indices are compatible."
        )

    coords = np.column_stack([
        np.array(x_out, dtype=float),
        np.array(y_out, dtype=float),
        np.array(z_out, dtype=float),
    ])

    # Centre nanotube in the xy box
    # Use the outermost radius for box sizing
    max_r  = radius + abs(structure.max_z_offset)
    box_xy = max_r * 2.0 + vacuum
    box_z  = T_norm
    shift  = box_xy / 2.0
    coords[:, 0] += shift
    coords[:, 1] += shift

    box = np.array([box_xy, box_xy, box_z])

    return NanotubeStructure(
        chirality=chirality,
        symbols=sym_out,
        coords=coords,
        box=box,
        vacuum=vacuum,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Bond validation
# ─────────────────────────────────────────────────────────────────────────────

SPURIOUS_MIN_SHORTENING = 0.10


def check_spurious_bonds(
    structure:  LatticeStructure,
    nt:         NanotubeStructure,
    tolerance:  float = 1.20,
    settings:   "BondSettings | None" = None,
    min_shortening: float = SPURIOUS_MIN_SHORTENING,
) -> set[frozenset]:
    """
    Detect bonds in the 3D nanotube that do not exist in the flat 2D structure.

    Returns a set of frozensets, each frozenset being an unordered pair of
    element symbols, e.g. {frozenset({'S', 'S'})} for a spurious S–S bond.
    An empty set means no spurious bonds were found.

    A species pair is spurious when, in the tube, it comes within the bond
    cutoff (``tolerance`` times the sum of covalent radii) although it is not
    bonded in the flat sheet, AND its shortest distance in the tube is at least
    ``min_shortening`` (a fraction) below its shortest distance in the flat
    sheet.  The second condition is a hysteresis on the cutoff.  Without it a
    pair that already sits just outside the cutoff in the sheet is flagged by
    any curvature at all: in C2DB's Fe2Mo2F2O8 the Fe-Mo pair across the
    Fe-O-Mo bridge is 3.71 A in the sheet against a 3.67 A cutoff, and every
    one of its 125 tubes was flagged, the (17,2) tube of D = 36.6 A for a
    3.5 % squeeze.  With 10 % the gentle squeezes are released (115 of those
    125 tubes, and 10.4 % of a random sample of 250 flagged catalogue tubes)
    while real contacts stay flagged: the median flagged pair is 26 % shorter
    than in the sheet, and the (3,1) tube of Fe2Mo2F2O8 is at -21 %.
    ``min_shortening=0`` restores the plain cutoff test.

    Parameters
    ----------
    structure : original flat 2D unit cell
    nt        : built 3D nanotube
    tolerance : bond-detection tolerance factor (used when settings=None)
    settings  : optional BondSettings for per-pair cutoffs (takes precedence)
    min_shortening : fraction the pair must shorten relative to the sheet
    """
    from .connectivity import compute_bonds, BondSettings as _BS

    # ── 2D bonds: use the flat atom positions including z ────────────────────
    flat_coords = np.array([
        [a["pos"][0], a["pos"][1], a.get("z", 0.0)]
        for a in structure.atoms
    ], dtype=float)
    flat_syms = [a["symbol"] for a in structure.atoms]

    # Build a 3×3 supercell of the 2D cell to capture all intra-cell bonds.
    a1, a2 = structure.a1, structure.a2
    supercell_coords: list = []
    supercell_syms:   list = []
    for di in range(-1, 2):
        for dj in range(-1, 2):
            shift = di * np.array([a1[0], a1[1], 0.0]) + \
                    dj * np.array([a2[0], a2[1], 0.0])
            supercell_coords.append(flat_coords + shift)
            supercell_syms.extend(flat_syms)
    sc_coords = np.vstack(supercell_coords)

    flat_bonds = compute_bonds(sc_coords, supercell_syms,
                               tolerance=tolerance, settings=settings)
    flat_species_pairs: set[frozenset] = set()
    for i, j in flat_bonds:
        flat_species_pairs.add(frozenset([supercell_syms[i], supercell_syms[j]]))

    def flat_min_distance(pair: frozenset) -> float:
        """Shortest distance of this species pair in the sheet (periodic in-plane)."""
        from scipy.spatial import cKDTree
        a_, b_ = (tuple(pair) * 2)[:2]
        sc_syms = np.array(supercell_syms)
        ia = np.where(np.array(flat_syms) == a_)[0]
        ib = np.where(sc_syms == b_)[0]
        if len(ia) == 0 or len(ib) == 0:
            return float("inf")
        # A wider supercell for the reference distance than for the bonds: a
        # long cell edge can put the nearest image beyond the 3 x 3 block.
        big, big_syms = [], []
        for di in range(-2, 3):
            for dj in range(-2, 3):
                shift = di * np.array([a1[0], a1[1], 0.0]) + dj * np.array([a2[0], a2[1], 0.0])
                big.append(flat_coords + shift)
                big_syms.extend(flat_syms)
        big = np.vstack(big)
        jb = np.where(np.array(big_syms) == b_)[0]
        d, _ = cKDTree(big[jb]).query(flat_coords[ia], k=min(2, len(jb)))
        d = np.atleast_2d(d)
        d = d[d > 0.4]
        return float(d.min()) if d.size else float("inf")

    # ── 3D bonds, periodic along the tube axis ───────────────────────────────
    # The tube repeats along z with period box[2].  Searching one cell on its
    # own misses every close pair that straddles the cell boundary, and where
    # that boundary falls depends on how the sheet was cut: in the nanotube
    # catalogue, 2 727 tubes that are identical in both rolling senses got
    # different spurious-bond verdicts for the two cuts, and 58 of 60 of those
    # agree once the neighbouring cells are searched too.  Enough images are
    # added to cover the longest bond the cutoffs allow, since a short cell
    # (T of 2.5 A against a 6 A cutoff) needs more than one on each side.
    from .connectivity import get_radius
    X = np.asarray(nt.coords, dtype=float)
    syms = list(nt.symbols)
    N = len(syms)
    Lz = float(nt.box[2])
    if settings is not None:
        tol_scale = settings.tolerance
    else:
        tol_scale = tolerance
    reach = 2.0 * max(get_radius(sy) for sy in set(syms)) * tol_scale
    if settings is not None and getattr(settings, "pair_cutoffs", None):
        reach = max([reach] + [float(v) for v in settings.pair_cutoffs.values()])
    k_img = max(1, int(np.ceil(reach / Lz))) if Lz > 1e-9 else 0
    # The original cell first, so that index < N means "an atom of this cell".
    order = [0] + [k for k in range(-k_img, k_img + 1) if k != 0]
    blocks = [X + np.array([0.0, 0.0, k * Lz]) for k in order]
    aug = np.vstack(blocks)
    nt_bonds = compute_bonds(aug, syms * len(blocks), tolerance=tolerance, settings=settings)
    shortest: dict[frozenset, float] = {}
    for i, j in nt_bonds:
        if i < N or j < N:                    # at least one atom of the cell itself
            pair = frozenset([syms[i % N], syms[j % N]])
            d = float(np.linalg.norm(aug[i] - aug[j]))
            shortest[pair] = min(d, shortest.get(pair, float("inf")))

    # ── Spurious = bonded in 3D, not in 2D, and genuinely shortened ──────────
    spurious: set[frozenset] = set()
    for pair, d_tube in shortest.items():
        if pair in flat_species_pairs:
            continue
        if min_shortening > 0:
            d_flat = flat_min_distance(pair)
            if d_tube > (1.0 - min_shortening) * d_flat:
                continue
        spurious.add(pair)
    return spurious


CURVATURE_MIN_CHANGE = 0.10


def tube_sites(structure: LatticeStructure, chirality: ChiralityResult) -> tuple[np.ndarray, ...]:
    """(k, s, w, x, y) arrays of one tube cell, in build_nanotube's atom order."""
    rows = [(k, s, w, x, y) for k, _sym, s, w, x, y, _z in _iter_sites(structure, chirality)]
    if not rows:
        raise RuntimeError(f"No atoms selected for ({chirality.n},{chirality.m}).")
    k, s, w, x, y = (np.array(c) for c in zip(*rows))
    return k.astype(np.int64), s, w, x, y


def _roll_sites(structure, chirality, s, w, z_off, roll_inward):
    """Cartesian positions of sheet points (s, w, z) after rolling, axis at the origin."""
    R = chirality.Ch_norm / (2.0 * np.pi)
    r = R + (-1.0 if roll_inward else 1.0) * z_off
    ang = 2.0 * np.pi * s
    return np.column_stack([r * np.cos(ang), r * np.sin(ang), w * chirality.T_norm])


def check_curvature_bonds(
    structure:   LatticeStructure,
    chirality:   ChiralityResult,
    roll_inward: bool = False,
    tolerance:   float = 1.20,
    settings:    "BondSettings | None" = None,
    min_change:  float = CURVATURE_MIN_CHANGE,
    sites:       "tuple[np.ndarray, ...] | None" = None,
) -> tuple[set[frozenset], set[frozenset]]:
    """
    Bonds that rolling forms and breaks, compared atom by atom with the sheet.

    Returns (formed, broken), each a set of species pairs.  Every tube atom
    keeps its identity in the flat layer (the 2D-cell atom and the lattice
    cell it comes from), so every pair of atoms has a distance in the sheet,
    d_flat, and one in the tube, d_tube, with the same bond cutoff of
    compute_bonds (``tolerance`` times the sum of covalent radii):

    * formed: d_tube is below the cutoff, d_flat is not, and d_tube is at least
      ``min_change`` shorter than d_flat;
    * broken: d_flat is below the cutoff, d_tube is not, and d_tube is at least
      ``min_change`` longer than d_flat.

    The margin is the same hysteresis as check_spurious_bonds, now per pair
    of atoms and in both directions.  check_spurious_bonds sees neither a
    broken bond nor a new bond of a species pair that is bonded elsewhere in
    the sheet; both are curvature artefacts all the same.  A tube pair that
    comes closer than the minimum bond length counts as formed, not as
    "no bond".

    The comparison is analytic in the sheet coordinates: a flat bond from
    atom k to atom k' of lattice offset (dx, dy) sends a tube atom at (s, w)
    to its partner at (s + ds, w + dw), with ds = (dx t2 - dy t1)/D and
    dw = (n dy - m dx)/D, so no pairs are searched on the seam.  Formed bonds
    come from a neighbour search in the tube with its axial images, where
    an image is the lattice shift T and a turn around the tube the shift Ch.
    ``sites`` (from tube_sites) saves enumerating the cell twice when both
    senses are checked.
    """
    from scipy.spatial import cKDTree
    from .connectivity import get_radius

    syms = [a["symbol"] for a in structure.atoms]
    nb = len(syms)
    zb = np.array([a.get("z", 0.0) for a in structure.atoms], dtype=float)
    if settings is None:
        hi = np.array([[(get_radius(a) + get_radius(b)) * tolerance for b in syms] for a in syms])
        lo = np.full((nb, nb), 0.40)
    else:
        hi = np.array([[settings.max_dist(a, b) for b in syms] for a in syms])
        lo = np.array([[settings.min_dist_pair(a, b) for b in syms] for a in syms])
    reach = float(hi.max())

    if sites is None:
        sites = tube_sites(structure, chirality)
    k, s, w, x, y = sites
    N = len(k)
    n, m, t1, t2 = chirality.n, chirality.m, chirality.t1, chirality.t2
    D = n * t2 - m * t1
    a1 = np.asarray(structure.a1, dtype=float)[:2]
    a2 = np.asarray(structure.a2, dtype=float)[:2]
    P = _roll_sites(structure, chirality, s, w, zb[k], roll_inward)

    def pair_name(i, j):
        return frozenset([syms[i], syms[j]])

    # ── Broken: every bond of the sheet, followed onto the tube ──────────────
    basis = np.column_stack([a1, a2])
    F = np.array([np.linalg.lstsq(basis, np.asarray(a["pos"], dtype=float)[:2], rcond=None)[0]
                  for a in structure.atoms])
    area = abs(a1[0] * a2[1] - a1[1] * a2[0])
    ni = int(np.ceil(reach / (area / np.linalg.norm(a2)))) + 1
    nj = int(np.ceil(reach / (area / np.linalg.norm(a1)))) + 1
    di, dj = np.meshgrid(np.arange(-ni, ni + 1), np.arange(-nj, nj + 1), indexing="ij")
    di, dj = di.ravel(), dj.ravel()
    broken: set[frozenset] = set()
    by_k = [np.nonzero(k == kk)[0] for kk in range(nb)]
    for kk in range(nb):
        idx = by_k[kk]
        if len(idx) == 0:
            continue
        for kp in range(nb):
            dx = di + F[kp, 0] - F[kk, 0]
            dy = dj + F[kp, 1] - F[kk, 1]
            vec = np.outer(dx, a1) + np.outer(dy, a2)
            d_flat = np.sqrt((vec ** 2).sum(axis=1) + (zb[kp] - zb[kk]) ** 2)
            sel = (d_flat > lo[kk, kp]) & (d_flat < hi[kk, kp])
            if not sel.any():
                continue
            if pair_name(kk, kp) in broken:
                continue
            for bx, by, bd in zip(dx[sel], dy[sel], d_flat[sel]):
                ds = (bx * t2 - by * t1) / D
                dw = (n * by - m * bx) / D
                Q = _roll_sites(structure, chirality, s[idx] + ds, w[idx] + dw,
                                np.full(len(idx), zb[kp]), roll_inward)
                d_tube = np.linalg.norm(P[idx] - Q, axis=1)
                if np.any((d_tube >= hi[kk, kp]) & (d_tube >= (1.0 + min_change) * bd)):
                    broken.add(pair_name(kk, kp))
                    break

    # ── Formed: every close pair of the tube, traced back to the sheet ───────
    Lz = float(chirality.T_norm)
    k_img = max(1, int(np.ceil(reach / Lz))) if Lz > 1e-9 else 0
    order = np.array([0] + [b for b in range(-k_img, k_img + 1) if b != 0])
    aug = np.vstack([P + np.array([0.0, 0.0, b * Lz]) for b in order])
    pairs = cKDTree(aug).query_pairs(reach, output_type="ndarray")
    formed: set[frozenset] = set()
    if len(pairs):
        pairs = pairs[(pairs[:, 0] < N) | (pairs[:, 1] < N)]
        i, j = pairs[:, 0], pairs[:, 1]
        ii, jj = i % N, j % N
        bi, bj = order[i // N], order[j // N]
        ki, kj = k[ii], k[jj]
        dx = (x[jj] + bj * t1) - (x[ii] + bi * t1)
        dy = (y[jj] + bj * t2) - (y[ii] + bi * t2)
        turns = np.round((dx * t2 - dy * t1) / D)
        dx = dx - turns * n
        dy = dy - turns * m
        vec = np.outer(dx, a1) + np.outer(dy, a2)
        d_flat = np.sqrt((vec ** 2).sum(axis=1) + (zb[kj] - zb[ki]) ** 2)
        d_tube = np.linalg.norm(aug[i] - aug[j], axis=1)
        h, l = hi[ki, kj], lo[ki, kj]
        new = (d_tube < h) & ~((d_flat > l) & (d_flat < h)) & (d_tube <= (1.0 - min_change) * d_flat)
        for a, b in set(zip(ki[new].tolist(), kj[new].tolist())):
            formed.add(pair_name(a, b))
    return formed, broken


def curvature_tokens(formed: set[frozenset], broken: set[frozenset]) -> list[str]:
    """Catalogue spelling: "+A-B" for a bond rolling forms, "-A-B" for one it breaks."""
    def name(p):
        a = sorted(p)
        return f"{a[0]}-{a[-1]}"
    return sorted("+" + name(p) for p in formed) + sorted("-" + name(p) for p in broken)
