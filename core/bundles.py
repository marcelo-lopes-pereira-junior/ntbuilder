"""
core/bundles.py
---------------
Nanotube bundle builder.

Given a single NanotubeStructure (already built), replicate it on a 2D
lattice to produce a periodic bundle supercell.  The result is a new
NanotubeStructure whose XY box vectors correspond to the bundle lattice
vectors — ready for periodic DFT/MD calculations.

Supported geometries
--------------------
  "linear"      : 1-D chain along X  (1 × N_y tubes, default N_y=2)
  "triangle"    : equilateral triangle (3 tubes)
  "hexagonal7"  : 1 central + 6 surrounding tubes (7 total)
  "square4"     : 2×2 grid (4 tubes)
  "grid"        : arbitrary N_x × N_y rectangular grid

In all cases the bundle lattice pitch is:
    p = diameter + spacing     (centre-to-centre)
where *spacing* is the surface-to-surface gap (default 3.4 Å, vdW).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .builder import NanotubeStructure


# ─────────────────────────────────────────────────────────────────────────────
# Data container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BundleResult:
    """Result of build_bundle()."""
    nanotube:   NanotubeStructure   # merged supercell, export-ready
    geometry:   str
    n_tubes:    int
    pitch:      float               # Å — centre-to-centre distance
    spacing:    float               # Å — surface-to-surface gap


# ─────────────────────────────────────────────────────────────────────────────
# Geometry definitions  (offsets in units of *pitch*)
# ─────────────────────────────────────────────────────────────────────────────

def _grid_offsets(nx: int, ny: int) -> list[tuple[float, float]]:
    """Rectangular N×M grid, centred at origin."""
    offs = []
    for iy in range(ny):
        for ix in range(nx):
            offs.append((ix - (nx - 1) / 2.0, iy - (ny - 1) / 2.0))
    return offs


_GEOMETRY_OFFSETS: dict[str, list[tuple[float, float]]] = {
    # fractional units of (pitch_x, pitch_y)
    "linear":     [(0.0, 0.0), (1.0, 0.0)],
    "triangle":   [(0.0, 0.0), (1.0, 0.0), (0.5, math.sqrt(3) / 2.0)],
    "hexagonal7": [
        (0.0, 0.0),                                 # centre
        (1.0, 0.0), (-1.0, 0.0),                    # ±x
        (0.5,  math.sqrt(3) / 2.0),
        (-0.5, math.sqrt(3) / 2.0),
        (0.5, -math.sqrt(3) / 2.0),
        (-0.5,-math.sqrt(3) / 2.0),
    ],
    "square4":    _grid_offsets(2, 2),
}


# ─────────────────────────────────────────────────────────────────────────────
# Main builder
# ─────────────────────────────────────────────────────────────────────────────

def build_bundle(
    nt:       NanotubeStructure,
    geometry: str   = "hexagonal7",
    spacing:  float = 3.4,
    vacuum:   float = 10.0,
    nx:       int   = 2,
    ny:       int   = 2,
) -> BundleResult:
    """
    Replicate *nt* on a 2D lattice to form a nanotube bundle supercell.

    Parameters
    ----------
    nt        : single-nanotube structure (axis along Z).
    geometry  : one of ``"linear"``, ``"triangle"``, ``"hexagonal7"``,
                ``"square4"``, ``"grid"``.
    spacing   : surface-to-surface intertube gap in Å (default 3.4 Å).
                Sets the bundle's internal pitch.
    vacuum    : lateral vacuum padding around the **outermost** tube
                surfaces in Å (default 10 Å).  The simulation box is sized
                so that the surface of each external tube has *vacuum* Å
                of free space to the box boundary.  Set to ``0`` to obtain
                a tightly packed periodic bundle (touching the box at the
                intertube spacing).
    nx, ny    : grid dimensions (used only when geometry="grid").

    Returns
    -------
    BundleResult with merged NanotubeStructure ready for export.

    Notes
    -----
    The tube axis (Z) is unchanged.  The XY box is rebuilt from scratch
    so that the *vacuum* parameter is correctly applied to the bundle as
    a whole — not inherited from the single-tube ``nt`` box, which was
    sized for an isolated tube.
    """
    geometry = geometry.lower().strip()
    if geometry not in _GEOMETRY_OFFSETS and geometry != "grid":
        raise ValueError(
            f"Unknown geometry '{geometry}'. "
            f"Choose from: {', '.join(list(_GEOMETRY_OFFSETS) + ['grid'])}."
        )

    # ── Effective outer diameter from the atoms' centre of mass ──────────────
    # ``nt.diameter`` returns the diameter of the chirality *of the innermost
    # wall* — for an MWNT or an already-built bundle this is far smaller
    # than the actual radial extent of the atoms, so using it as the pitch
    # produces overlapping replicas.  We also pivot on the *atomic centroid*
    # in XY (not the box centre) so that off-centred input structures —
    # bundles fed to /api/bundle a second time, MWNTs whose merging path
    # left them slightly off, etc. — are still tiled symmetrically around
    # their geometric centre.
    if nt.coords.shape[0] > 0:
        cx_in = float(np.mean(nt.coords[:, 0]))
        cy_in = float(np.mean(nt.coords[:, 1]))
    else:
        cx_in = float(nt.box[0]) / 2.0
        cy_in = float(nt.box[1]) / 2.0
    radii = np.hypot(nt.coords[:, 0] - cx_in, nt.coords[:, 1] - cy_in)
    d_outer = 2.0 * float(np.max(radii)) if radii.size else float(nt.diameter)
    # Fall back to chirality.diameter if for some reason the radial scan
    # yields a smaller value than the analytic single-tube diameter
    # (numerical pathology).
    d     = max(d_outer, float(nt.diameter))
    pitch = d + spacing         # centre-to-centre distance

    # ── Raw offsets in fractional pitch units ────────────────────────────────
    if geometry == "grid":
        offsets = _grid_offsets(nx, ny)
    else:
        offsets = _GEOMETRY_OFFSETS[geometry]

    n_tubes = len(offsets)

    # ── Convert offsets to Å ─────────────────────────────────────────────────
    abs_offsets = [(ox * pitch, oy * pitch) for ox, oy in offsets]

    # ── Re-centre offsets on the mid-span of their bounding box ──────────────
    # The geometry tables intentionally store some configurations with
    # asymmetric origins — ``linear`` lives in ``[(0,0), (1,0)]`` and
    # ``triangle`` in ``[(0,0), (1,0), (0.5, sqrt(3)/2)]`` — to keep the
    # tabulated coordinates readable.  Without this re-centring step the
    # bundle's atomic centroid would land at the mid-span of those offsets
    # (e.g. ``(pitch/2, 0)`` for ``linear``) rather than at the box centre,
    # producing a bundle visibly shifted toward one corner of the
    # simulation box.  Geometries that are already centred on the origin
    # (``hexagonal7`` and the ``square4``/``grid`` outputs of
    # :func:`_grid_offsets`) are unaffected because ``mid_x`` and ``mid_y``
    # evaluate to zero for them.
    xs = [o[0] for o in abs_offsets]
    ys = [o[1] for o in abs_offsets]
    mid_x = (max(xs) + min(xs)) / 2.0
    mid_y = (max(ys) + min(ys)) / 2.0
    if abs(mid_x) > 1e-12 or abs(mid_y) > 1e-12:
        abs_offsets = [(ox - mid_x, oy - mid_y) for ox, oy in abs_offsets]
        xs = [o[0] for o in abs_offsets]
        ys = [o[1] for o in abs_offsets]

    # ── Bounding box: span of centres + tube diameter + 2·vacuum ─────────────
    # The "+ d" accounts for the half-diameter on each side; the "+ 2·vacuum"
    # adds the user-requested padding measured from the outer tube *surface*.
    x_span = (max(xs) - min(xs)) + d + 2.0 * vacuum
    y_span = (max(ys) - min(ys)) + d + 2.0 * vacuum

    box_x = x_span
    box_y = y_span
    box_z = float(nt.box[2])

    cx_bundle = box_x / 2.0
    cy_bundle = box_y / 2.0

    # ── Place each tube relative to the new bundle centre ────────────────────
    # Use the same atomic centroid as above (rather than the box centre)
    # so that off-centred inputs still produce a properly centred bundle.
    cx_nt = cx_in
    cy_nt = cy_in

    all_symbols: list[str] = []
    all_coords:  list[np.ndarray] = []

    for ox, oy in abs_offsets:
        coords_k = nt.coords.copy()
        coords_k[:, 0] += (cx_bundle + ox) - cx_nt
        coords_k[:, 1] += (cy_bundle + oy) - cy_nt
        all_symbols.extend(nt.symbols)
        all_coords.append(coords_k)

    merged_coords = np.vstack(all_coords)
    merged_box    = np.array([box_x, box_y, box_z])

    merged = NanotubeStructure(
        chirality=nt.chirality,
        symbols=all_symbols,
        coords=merged_coords,
        box=merged_box,
        vacuum=vacuum,
    )

    return BundleResult(
        nanotube=merged,
        geometry=geometry,
        n_tubes=n_tubes,
        pitch=pitch,
        spacing=spacing,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Available geometry names (for UI dropdowns)
# ─────────────────────────────────────────────────────────────────────────────

GEOMETRIES = list(_GEOMETRY_OFFSETS.keys()) + ["grid"]

GEOMETRY_LABELS = {
    "linear":     "Linear (2 tubes)",
    "triangle":   "Triangle (3 tubes)",
    "square4":    "Square 2×2 (4 tubes)",
    "hexagonal7": "Hexagonal (7 tubes)",
    "grid":       "Custom grid (N×M)",
}
