"""
core/deformations.py
--------------------
Post-build coordinate transformations.

All functions take a NanotubeStructure and return a *new* NanotubeStructure
(immutable-style: originals are never mutated).

Transformations
---------------
apply_axial_strain(nt, strain)
    Stretch or compress the nanotube along its axis (Z by convention).
    strain > 0 → tensile; strain < 0 → compressive.

apply_torsion(nt, twist_rate)
    Apply a uniform twist φ(z) = twist_rate · z around the tube axis.
    twist_rate in degrees per Å.

torsion_closes(nt, twist_rate, n_rep)
    True when the total twist over the cell is a rotation symmetry of the
    untwisted structure, so the twisted cell stays periodic along Z.

apply_radial_strain(nt, strain)  [experimental]
    Uniform radial scaling of (x,y) coordinates — models hydrostatic
    in-plane pressure. Use with caution: does not relax bonds.
"""

from __future__ import annotations

import math
import numpy as np

from .builder import NanotubeStructure


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _copy(nt: NanotubeStructure, new_coords: np.ndarray,
          new_box: np.ndarray | None = None) -> NanotubeStructure:
    """Return a new NanotubeStructure sharing metadata but with new coords/box."""
    return NanotubeStructure(
        chirality=nt.chirality,
        symbols=list(nt.symbols),
        coords=new_coords.copy(),
        box=new_box.copy() if new_box is not None else nt.box.copy(),
        vacuum=nt.vacuum,
    )


def _tube_centre(nt: NanotubeStructure) -> tuple[float, float]:
    """XY centre of the simulation box (= tube axis position)."""
    return float(nt.box[0]) / 2.0, float(nt.box[1]) / 2.0


def _xy_centroid(nt: NanotubeStructure) -> tuple[float, float]:
    """Compute the atomic XY centroid (centre of mass with unit masses).

    Used as the torsion axis so the operation behaves correctly for
    bundles or off-centred tubes, where the geometric box centre and
    the actual centre of the atomic distribution differ.
    """
    if nt.coords.shape[0] == 0:
        return _tube_centre(nt)
    cx = float(np.mean(nt.coords[:, 0]))
    cy = float(np.mean(nt.coords[:, 1]))
    return cx, cy


def replicate_z(nt: NanotubeStructure, n_rep: int) -> NanotubeStructure:
    """Tile *nt* ``n_rep`` times along the Z axis.

    The simulation box length along Z is scaled accordingly.  Used by
    torsion to apply the twist to the *displayed* (multi-cell) structure
    rather than only the single unit cell.
    """
    if n_rep <= 1:
        return nt
    Lz   = float(nt.box[2])
    base = nt.coords.copy()
    tiles = [base + np.array([0.0, 0.0, k * Lz]) for k in range(n_rep)]
    new_coords  = np.vstack(tiles)
    new_symbols = list(nt.symbols) * n_rep
    new_box     = nt.box.copy()
    new_box[2] *= n_rep
    return NanotubeStructure(
        chirality = nt.chirality,
        symbols   = new_symbols,
        coords    = new_coords,
        box       = new_box,
        vacuum    = nt.vacuum,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Axial strain / compression
# ─────────────────────────────────────────────────────────────────────────────

def apply_axial_strain(nt: NanotubeStructure, strain: float) -> NanotubeStructure:
    """
    Scale the nanotube uniformly along the Z axis by factor (1 + strain).

    Parameters
    ----------
    nt     : source NanotubeStructure
    strain : fractional strain, e.g. +0.05 = 5 % tension, -0.05 = 5 % compression.
             Allowed range: −0.5 to +2.0.

    Returns
    -------
    New NanotubeStructure with scaled Z coordinates and updated box[2].
    """
    if not (-0.5 <= strain <= 2.0):
        raise ValueError(f"strain must be in [-0.5, 2.0], got {strain:.3f}")

    factor    = 1.0 + strain
    new_coords = nt.coords.copy()
    new_coords[:, 2] *= factor

    new_box      = nt.box.copy()
    new_box[2]  *= factor

    return _copy(nt, new_coords, new_box)


# ─────────────────────────────────────────────────────────────────────────────
# Torsion (twist)
# ─────────────────────────────────────────────────────────────────────────────

TORSION_CLOSE_TOL = 0.02   # Å, atom-matching tolerance of torsion_closes


def torsion_closes(
    nt:         NanotubeStructure,
    twist_rate: float,
    n_rep:      int   = 1,
    tol:        float = TORSION_CLOSE_TOL,
) -> bool:
    """
    Return True when a uniform twist keeps the structure periodic along Z.

    ``apply_torsion`` rotates the atom at height z by φ(z) = twist_rate · z
    about the atomic XY centroid.  The image of that atom one period higher,
    at z + L, is rotated by φ(z) + Φ with Φ = twist_rate · L.  The twisted
    cell therefore joins itself across the Z boundary exactly when the
    rotation R(Φ) about the same axis maps the *untwisted* periodic structure
    onto itself: Φ must be a rotation symmetry of the structure (a multiple
    of 360°/g for a single tube of rotational order g, of lcm(60°, 360°/g)
    for a hexagonal bundle, and so on).

    The test is numerical: every atom rotated by Φ, wrapped into the cell
    along Z, must have an atom of the same species within ``tol`` Å
    (cKDTree, periodic images across the Z boundary included).  An isometry
    that sends each atom within ``tol`` of a distinct partner is a bijection
    whenever ``tol`` is below half the shortest interatomic distance.

    Parameters
    ----------
    nt          : untwisted structure, periodic along Z with period box[2].
    twist_rate  : rotation rate in degrees per Å.
    n_rep       : number of unit cells the twist spans (the same ``n_rep``
                  passed to ``apply_torsion``), so L = n_rep · box[2].
    tol         : matching tolerance in Å.
    """
    from scipy.spatial import cKDTree

    n_atoms = nt.coords.shape[0]
    if abs(twist_rate) < 1e-12 or n_atoms == 0:
        return True

    Lz_cell = float(nt.box[2])
    if Lz_cell <= 0.0:
        return False
    phi = twist_rate * Lz_cell * max(1, int(n_rep))      # degrees

    cx, cy = _xy_centroid(nt)
    dx = nt.coords[:, 0] - cx
    dy = nt.coords[:, 1] - cy
    r_max = float(np.sqrt(dx * dx + dy * dy).max())

    # Φ ≡ 0 (mod 360°): the identity, periodic for any structure.  The
    # residual angle moves the outermost atom by at most r_max · |Δφ|.
    residual = (phi + 180.0) % 360.0 - 180.0
    if abs(math.radians(residual)) * r_max <= tol:
        return True

    ref = nt.coords.copy()
    ref[:, 2] %= Lz_cell
    symbols = np.asarray(nt.symbols)

    # Periodic images of the atoms that sit within tol of either Z face.
    lo = ref[:, 2] < tol
    hi = ref[:, 2] > Lz_cell - tol
    img_idx = np.concatenate([np.arange(n_atoms), np.flatnonzero(lo),
                              np.flatnonzero(hi)])
    img = np.vstack([ref, ref[lo] + [0.0, 0.0, Lz_cell],
                     ref[hi] - [0.0, 0.0, Lz_cell]])

    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    rot = np.empty_like(ref)
    rot[:, 0] = cx + c * dx - s * dy
    rot[:, 1] = cy + s * dx + c * dy
    rot[:, 2] = ref[:, 2]

    dist, j = cKDTree(img).query(rot, distance_upper_bound=tol)
    if not np.all(np.isfinite(dist)):
        return False
    return bool(np.all(symbols[img_idx[j]] == symbols))


def apply_torsion(
    nt:        NanotubeStructure,
    twist_rate: float,
    z_vacuum:   float | None = None,
    n_rep:      int = 1,
    closes:     bool | None = None,
) -> NanotubeStructure:
    """
    Apply a uniform helical twist φ(z) = twist_rate × z to the nanotube.

    Each atom at axial position z is rotated by φ(z) degrees around the
    **atomic XY centroid** (centre of mass with unit weights) — this
    matches the geometric centre for a single SWNT and also handles
    bundles correctly, where the box centre and the actual centre of the
    atomic distribution differ.

    Atomic z-coordinates are unchanged.  The total angle over the twisted
    length L is Φ = twist_rate × L.  When Φ is a rotation symmetry of the
    untwisted structure about the twist axis (see ``torsion_closes``) the
    twisted cell joins itself across the Z boundary and stays periodic, so
    the Z box is kept.  Otherwise the **axial (Z) periodicity is broken** —
    periodic boundary conditions would join an end at angle 0 to an end at
    angle Φ — and a vacuum slab is added along Z.

    Parameters
    ----------
    nt          : source NanotubeStructure.
    twist_rate  : rotation rate in degrees per Å.  Positive → right-hand
                  (conventional) twist.
    z_vacuum    : Z padding added to each end of the simulation box, in Å.
                  When ``None`` the function keeps the Z box if the twist
                  closes the cell and otherwise uses ``nt.vacuum`` (the
                  same lateral vacuum used during the radial build).  An
                  explicit value is always honoured; ``0`` keeps the
                  original Z box length.
    n_rep       : tile the input ``n_rep`` times along Z **before**
                  applying the twist.  Useful when the user has set
                  ``Reps > 1`` in the viewer and expects the torsion to
                  act on the supercell that they see, not on a single
                  unit cell.  Default 1 = no replication.
    closes      : result of ``torsion_closes(nt, twist_rate, n_rep)`` when
                  the caller has already computed it; ``None`` computes it
                  if needed (only when ``z_vacuum`` is ``None``).

    Returns
    -------
    New NanotubeStructure with twisted XY coordinates and a Z box
    extended by ``2·z_vacuum``.
    """
    if z_vacuum is None:
        if closes is None:
            closes = torsion_closes(nt, twist_rate, n_rep=n_rep)
        z_vacuum = 0.0 if closes else float(nt.vacuum)

    # ── Optional pre-replication along Z ────────────────────────────────────
    if n_rep > 1:
        nt = replicate_z(nt, n_rep)

    cx, cy = _xy_centroid(nt)

    coords    = nt.coords.copy()
    dx        = coords[:, 0] - cx
    dy        = coords[:, 1] - cy
    z         = coords[:, 2]

    phi_rad   = np.deg2rad(twist_rate * z)    # per-atom rotation angle
    cos_phi   = np.cos(phi_rad)
    sin_phi   = np.sin(phi_rad)

    coords[:, 0] = cx + cos_phi * dx - sin_phi * dy
    coords[:, 1] = cy + sin_phi * dx + cos_phi * dy

    new_box = nt.box.copy()
    if z_vacuum > 0.0:
        # Pad symmetrically along Z and re-centre the atoms in the new box.
        new_box[2] = float(nt.box[2]) + 2.0 * float(z_vacuum)
        shift_z    = (new_box[2] - float(nt.box[2])) / 2.0
        coords[:, 2] += shift_z

    return _copy(nt, coords, new_box)


def torsion_warning(twist_rate: float, z_vacuum: float,
                    n_rep: int = 1, closes: bool = False,
                    total_angle: float | None = None) -> str | None:
    """
    Return a short user-facing message describing the consequence of
    applying torsion on a periodic structure, or ``None`` when no torsion
    is applied.

    With ``closes=True`` (see ``torsion_closes``) the message says that the
    twisted cell stays periodic along Z; ``total_angle`` (degrees), when
    given, is quoted.  Otherwise it warns about the loss of axial
    periodicity and the vacuum slab.

    Mentions ``n_rep`` when the torsion was applied to a multi-cell
    supercell so the user is aware that the displayed Reps were absorbed
    into the geometry.  For a non-closing twist the resulting structure is
    *no longer* a periodic unit cell and the Reps spinbox should be hidden
    until the torsion is undone.
    """
    if abs(twist_rate) < 1e-9:
        return None
    if closes:
        angle = (f" turns the cell by {total_angle:.2f}° along its length, a "
                 f"rotation symmetry of the untwisted structure, so it"
                 if total_angle is not None else "")
        parts = [
            f"Torsion of {twist_rate:+.4f} °/Å{angle} closes the cell: the "
            f"twisted cell stays periodic along Z."
        ]
        if n_rep > 1:
            parts.append(
                f"The twist was applied to the supercell of {n_rep} unit "
                f"cells that you had on screen, which is now the periodic "
                f"cell; Reps replicates it."
            )
        return "  ".join(parts)
    parts = [
        f"Torsion of {twist_rate:+.4f} °/Å breaks the axial periodicity of "
        f"the nanotube — the structure is no longer commensurate along Z."
    ]
    if n_rep > 1:
        parts.append(
            f"The twist was applied to the supercell of {n_rep} unit cells "
            f"that you had on screen; the resulting geometry is a single "
            f"finite segment and the Reps control has been disabled."
        )
    parts.append(
        f"A vacuum slab of {z_vacuum:.2f} Å was added to each end of the "
        f"simulation box.  Treat the result as a finite-length segment in "
        f"PBC-aware codes."
    )
    return "  ".join(parts)


# ─────────────────────────────────────────────────────────────────────────────
# Radial strain (experimental)
# ─────────────────────────────────────────────────────────────────────────────

def apply_radial_strain(nt: NanotubeStructure, strain: float) -> NanotubeStructure:
    """
    Scale (x,y) coordinates uniformly by (1 + strain) around the tube axis.

    This models isotropic in-plane compression or expansion without atomic
    relaxation. Useful for generating a series of structures for DFT,
    but the resulting bond lengths will be unphysical until relaxed.

    Parameters
    ----------
    nt     : source NanotubeStructure
    strain : fractional radial strain (−0.3 to +0.3 recommended).
    """
    if not (-0.5 <= strain <= 1.0):
        raise ValueError(f"radial strain must be in [-0.5, 1.0], got {strain:.3f}")

    factor = 1.0 + strain
    cx, cy = _tube_centre(nt)

    coords        = nt.coords.copy()
    coords[:, 0]  = cx + (coords[:, 0] - cx) * factor
    coords[:, 1]  = cy + (coords[:, 1] - cy) * factor

    # Adjust box XY (the tube is larger/smaller).
    # Builder convention: box_xy = diameter + vacuum (vacuum = total padding, not per-side).
    new_box        = nt.box.copy()
    old_diam       = nt.diameter
    new_box[0]     = old_diam * factor + nt.vacuum
    new_box[1]     = new_box[0]

    # Re-centre atoms in the new box
    new_cx = new_box[0] / 2.0
    new_cy = new_box[1] / 2.0
    coords[:, 0] += new_cx - cx
    coords[:, 1] += new_cy - cy

    return _copy(nt, coords, new_box)


# ─────────────────────────────────────────────────────────────────────────────
# Deformation descriptor (for UI / methods text)
# ─────────────────────────────────────────────────────────────────────────────

def deformation_description(
    axial_strain:  float = 0.0,
    twist_rate:    float = 0.0,
    radial_strain: float = 0.0,
) -> str:
    """Return a one-line human-readable description of applied deformations."""
    parts = []
    if abs(axial_strain) > 1e-6:
        sign = "tensile" if axial_strain > 0 else "compressive"
        parts.append(f"axial {sign} strain {abs(axial_strain) * 100:.2f}%")
    if abs(twist_rate) > 1e-6:
        parts.append(f"torsion {twist_rate:+.4f} °/Å")
    if abs(radial_strain) > 1e-6:
        sign = "expansive" if radial_strain > 0 else "compressive"
        parts.append(f"radial {sign} strain {abs(radial_strain) * 100:.2f}%")
    return ", ".join(parts) if parts else "none"
