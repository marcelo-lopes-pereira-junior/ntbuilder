"""
core/builder.py
---------------
Constructs the 3D nanotube from a LatticeStructure and a ChiralityResult.

Algorithm
---------
1. Tile the 2D unit cell over a supercell large enough to cover at least
   one full (Ch × T) unit cell of the nanotube.
2. Project each atom onto the Ch and T axes using dot products.
3. Keep only atoms whose (u, v) coordinates satisfy:
       0 ≤ u < |Ch|  and  0 ≤ v < |T|
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

_TOL_LO = 1e-4   # lower-boundary inclusion tolerance (v ≥ -_TOL_LO)
_TOL_HI = 5e-4   # upper-boundary exclusion tolerance (v < ceil - _TOL_HI)
                 # Larger than _TOL_LO to catch floating-point images that land
                 # just below the upper boundary due to rounding.


def _iter_atoms(
    structure: LatticeStructure,
    chirality: ChiralityResult,
    margin: int,
) -> Iterator[tuple[str, float, float, float]]:
    """
    Yield (symbol, u_coord, v_coord, z_offset) for every atom in the nanotube
    unit cell.

    Uses analytic i-bounds for each (j, atom) pair so the supercell is never
    materialised as a 2D array.  Memory complexity is O(n_atoms_output) instead
    of O(margin²), making arbitrarily large unit cells feasible.
    """
    a1, a2  = structure.a1, structure.a2
    Ch_norm = chirality.Ch_norm
    T_norm  = chirality.T_norm
    u_hat   = chirality.Ch_vec / Ch_norm
    v_hat   = chirality.T_vec  / T_norm

    # Scalar projections of lattice vectors onto the nanotube axes
    a1_u = float(a1 @ u_hat)
    a1_v = float(a1 @ v_hat)
    a2_u = float(a2 @ u_hat)
    a2_v = float(a2 @ v_hat)

    for atom in structure.atoms:
        sym   = atom["symbol"]
        pos   = atom["pos"]
        z_off = atom.get("z", 0.0)
        p_u   = float(pos @ u_hat)
        p_v   = float(pos @ v_hat)

        for j in range(-margin, margin + 1):
            # Contributions from j and the atom offset
            ju = j * a2_u + p_u
            jv = j * a2_v + p_v

            # ── Analytic i-bounds from the u-constraint ───────────────────────
            # -_TOL_LO ≤ i*a1_u + ju < Ch_norm - _TOL_HI
            if abs(a1_u) > 1e-12:
                lo_u = (-_TOL_LO - ju) / a1_u
                hi_u = (Ch_norm - _TOL_HI - ju) / a1_u
                if a1_u > 0:
                    i_lo_u = math.ceil(lo_u  - 1e-9)
                    i_hi_u = math.floor(hi_u + 1e-9) + 1   # exclusive
                else:
                    i_lo_u = math.ceil(hi_u  - 1e-9)
                    i_hi_u = math.floor(lo_u + 1e-9) + 1
            else:
                # a1_u ≈ 0: u independent of i — skip row if out of range
                if not (-_TOL_LO <= ju < Ch_norm - _TOL_HI):
                    continue
                i_lo_u, i_hi_u = -margin, margin + 1

            # ── Analytic i-bounds from the v-constraint ───────────────────────
            # -_TOL_LO ≤ i*a1_v + jv < T_norm - _TOL_HI
            if abs(a1_v) > 1e-12:
                lo_v = (-_TOL_LO - jv) / a1_v
                hi_v = (T_norm - _TOL_HI - jv) / a1_v
                if a1_v > 0:
                    i_lo_v = math.ceil(lo_v  - 1e-9)
                    i_hi_v = math.floor(hi_v + 1e-9) + 1
                else:
                    i_lo_v = math.ceil(hi_v  - 1e-9)
                    i_hi_v = math.floor(lo_v + 1e-9) + 1
            else:
                if not (-_TOL_LO <= jv < T_norm - _TOL_HI):
                    continue
                i_lo_v, i_hi_v = -margin, margin + 1

            # Intersect the two ranges with the overall margin guard
            i_lo = max(i_lo_u, i_lo_v, -margin)
            i_hi = min(i_hi_u, i_hi_v, margin + 1)

            for i in range(i_lo, i_hi):
                u = i * a1_u + ju
                v = i * a1_v + jv
                # Final exact check using asymmetric tolerances
                if (-_TOL_LO <= u < Ch_norm - _TOL_HI and
                        -_TOL_LO <= v < T_norm  - _TOL_HI):
                    yield sym, u, v, z_off


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

    # Supercell margin — guarantees the iterator covers at least one full
    # (Ch × T) unit cell regardless of the lattice geometry.
    margin = max(
        abs(chirality.n) + abs(chirality.t1),
        abs(chirality.m) + abs(chirality.t2),
    ) + 3

    sym_out: list[str]   = []
    x_out:   list[float] = []
    y_out:   list[float] = []
    z_out:   list[float] = []

    for sym, u, v, z_off in _iter_atoms(structure, chirality, margin):
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

def check_spurious_bonds(
    structure:  LatticeStructure,
    nt:         NanotubeStructure,
    tolerance:  float = 1.20,
    settings:   "BondSettings | None" = None,
) -> set[frozenset]:
    """
    Detect bonds in the 3D nanotube that do not exist in the flat 2D structure.

    Returns a set of frozensets, each frozenset being an unordered pair of
    element symbols, e.g. {frozenset({'S', 'S'})} for a spurious S–S bond.
    An empty set means no spurious bonds were found.

    Parameters
    ----------
    structure : original flat 2D unit cell
    nt        : built 3D nanotube
    tolerance : bond-detection tolerance factor (used when settings=None)
    settings  : optional BondSettings for per-pair cutoffs (takes precedence)
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

    # ── 3D bonds ─────────────────────────────────────────────────────────────
    nt_coords_centered = nt.coords - nt.coords.mean(axis=0)
    nt_bonds = compute_bonds(nt_coords_centered, list(nt.symbols),
                             tolerance=tolerance, settings=settings)
    nt_species_pairs: set[frozenset] = set()
    for i, j in nt_bonds:
        nt_species_pairs.add(frozenset([nt.symbols[i], nt.symbols[j]]))

    # ── Spurious = in 3D but not in 2D ───────────────────────────────────────
    return nt_species_pairs - flat_species_pairs
