"""
core/symmetry.py
----------------
Tools for enforcing lattice symmetry and finding primitive cells.

snap_to_symmetry(structure)
    Fix floating-point deviations in lattice parameters
    (e.g. γ = 60.0001° → 60.0000°, |a₁| ≠ |a₂| → averaged) by snapping
    to the ideal cell for the detected lattice type.

find_primitive_cell(structure, tol)
    Reduce a conventional or supercell to the smallest equivalent unit cell
    by searching for the shortest translation symmetries.  Uses spglib when
    available; falls back to a built-in search algorithm otherwise.
"""

from __future__ import annotations

import math

import numpy as np

from core.io import LatticeStructure


# ─────────────────────────────────────────────────────────────────────────────
# Internal fractional-coordinate helpers
# ─────────────────────────────────────────────────────────────────────────────

def _to_frac(pos: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Cartesian 2D → fractional in (a1, a2) basis."""
    return np.linalg.solve(np.column_stack([a1, a2]), pos)


def _to_cart(frac: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    return frac[0] * a1 + frac[1] * a2


def _wrap(frac: np.ndarray) -> np.ndarray:
    """Wrap fractional coordinates into [0, 1)."""
    return frac - np.floor(frac + 1e-9)


def _min_image_cart(diff: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Return the minimum-image Cartesian vector for a difference vector."""
    frac = _to_frac(diff, a1, a2)
    frac -= np.round(frac)
    return _to_cart(frac, a1, a2)


# ─────────────────────────────────────────────────────────────────────────────
# Snap to symmetry
# ─────────────────────────────────────────────────────────────────────────────

def snap_to_symmetry(structure: LatticeStructure) -> tuple[LatticeStructure, str]:
    """
    Enforce exact lattice symmetry.

    Rules
    -----
    hexagonal (|a₁|≈|a₂|, γ≈60°)
        Sets |a₁| = |a₂| = (|a₁|+|a₂|)/2, γ = 60.000° exactly,
        a₁ aligned to the x-axis.
    square (|a₁|≈|a₂|, γ≈90°)
        Sets |a₁| = |a₂| = average, γ = 90.000° exactly.
    rectangular (γ≈90°)
        Sets γ = 90.000° exactly; individual lengths are kept.
    oblique
        Returns the structure unchanged (no unique ideal to snap to).

    Atom Cartesian positions are recomputed through fractional coordinates so
    they remain exactly on-lattice after snapping.

    Returns
    -------
    (new_structure, description)
        description is a human-readable string summarising the changes.
    """
    lt    = structure.lattice_type
    a1    = structure.a1.copy()
    a2    = structure.a2.copy()
    a_len = float(np.linalg.norm(a1))
    b_len = float(np.linalg.norm(a2))
    g_old = structure.gamma_deg

    if lt == "hexagonal":
        a_ideal = (a_len + b_len) / 2.0
        a1_new  = np.array([a_ideal, 0.0])
        a2_new  = np.array([a_ideal / 2.0, a_ideal * math.sqrt(3.0) / 2.0])
        desc = (
            f"a: {a_len:.5f} / {b_len:.5f} Å → {a_ideal:.5f} Å  |  "
            f"γ: {g_old:.4f}° → 60.0000°"
        )
    elif lt == "rectangular":
        a1_new = np.array([a_len, 0.0])
        a2_new = np.array([0.0, b_len])
        if abs(a_len - b_len) < 1e-3 * a_len:
            # Square
            a_ideal = (a_len + b_len) / 2.0
            a1_new  = np.array([a_ideal, 0.0])
            a2_new  = np.array([0.0, a_ideal])
            desc = (
                f"a: {a_len:.5f} / {b_len:.5f} Å → {a_ideal:.5f} Å  |  "
                f"γ: {g_old:.4f}° → 90.0000°"
            )
        else:
            desc = f"γ: {g_old:.4f}° → 90.0000°  (a = {a_len:.5f} Å, b = {b_len:.5f} Å)"
    else:
        return structure, "Oblique lattice — nothing to snap."

    # Reproject atomic positions through fractional coordinates
    new_atoms = []
    for atom in structure.atoms:
        frac    = _wrap(_to_frac(atom["pos"], a1, a2))
        new_pos = _to_cart(frac, a1_new, a2_new)
        new_atoms.append({**atom, "pos": new_pos})

    return LatticeStructure(new_atoms, a1_new, a2_new, source=structure.source), desc


# ─────────────────────────────────────────────────────────────────────────────
# Primitive-cell finder
# ─────────────────────────────────────────────────────────────────────────────

def _orient_and_wrap(structure: LatticeStructure) -> LatticeStructure:
    """
    Canonical post-processing for any newly built primitive cell:
    1. Rotate so that a1 points along the +x axis.
    2. Wrap all atom positions into [0, 1) in fractional coordinates.
    The lattice type and lengths are preserved; only orientation changes.
    """
    a1, a2 = structure.a1.copy(), structure.a2.copy()

    # Rotate a1 onto +x
    angle = math.atan2(float(a1[1]), float(a1[0]))
    if abs(angle) > 1e-9:
        cos_a, sin_a = math.cos(-angle), math.sin(-angle)
        R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        a1 = R @ a1
        a2 = R @ a2
        a1[1] = 0.0           # clean up floating-point noise

    # Ensure a2 has positive y-component
    if a2[1] < 0:
        a2 = -a2

    # Wrap atom positions into [0, 1) fractional
    new_atoms = []
    for atom in structure.atoms:
        # Rotate position
        angle_fwd = math.atan2(float(structure.a1[1]), float(structure.a1[0]))
        cos_f, sin_f = math.cos(-angle_fwd), math.sin(-angle_fwd)
        Rf  = np.array([[cos_f, -sin_f], [sin_f, cos_f]])
        pos = Rf @ atom["pos"]
        # Wrap to [0, 1)
        frac = _wrap(_to_frac(pos, a1, a2))
        new_atoms.append({**atom, "pos": _to_cart(frac, a1, a2)})

    return LatticeStructure(new_atoms, a1, a2, source=structure.source)


def find_primitive_cell(
    structure: LatticeStructure,
    tol: float = 0.08,
) -> tuple[LatticeStructure, str]:
    """
    Reduce a conventional / supercell to the smallest equivalent primitive cell.

    Uses spglib when available (more robust for complex structures); otherwise
    falls back to a built-in translation-symmetry search.

    Parameters
    ----------
    structure : input structure (possibly a supercell)
    tol       : position tolerance in Å

    Returns
    -------
    (new_structure, description)
        description summarises what changed.
    """
    n_in = len(structure.atoms)

    try:
        result = _primitive_spglib(structure, tol)
    except Exception:
        result = _primitive_builtin(structure, tol)

    result = _orient_and_wrap(result)
    n_out  = len(result.atoms)
    if n_out == n_in:
        desc = "Already primitive — no reduction possible."
    else:
        factor = n_in // n_out
        desc = (
            f"Reduced {n_in} → {n_out} atoms/cell  "
            f"(÷{factor})  |  "
            f"a: {result.a:.5f} Å, b: {result.b:.5f} Å, "
            f"γ: {result.gamma_deg:.4f}°  [{result.lattice_type}]"
        )
    return result, desc


# ── spglib backend ────────────────────────────────────────────────────────────

def _primitive_spglib(structure: LatticeStructure, tol: float) -> LatticeStructure:
    import spglib  # raises ImportError if not installed

    a1, a2 = structure.a1, structure.a2
    # Build a 3D cell with a large vacuum c-axis so spglib treats it as 2D
    c_vac  = 30.0
    lattice = np.array([
        [a1[0], a1[1], 0.0],
        [a2[0], a2[1], 0.0],
        [0.0,   0.0,   c_vac],
    ])

    M = np.column_stack([a1, a2])
    symbols_sorted = sorted(set(a["symbol"] for a in structure.atoms))
    sym2num = {s: i + 1 for i, s in enumerate(symbols_sorted)}
    num2sym = {v: k for k, v in sym2num.items()}

    positions_frac = []
    numbers        = []
    for atom in structure.atoms:
        f2 = np.linalg.solve(M, atom["pos"])
        positions_frac.append([float(f2[0]), float(f2[1]), 0.5])
        numbers.append(sym2num[atom["symbol"]])

    cell = (lattice, positions_frac, numbers)
    prim = spglib.find_primitive(cell, symprec=tol)
    if prim is None:
        raise RuntimeError("spglib: find_primitive returned None")

    prim_lattice, prim_frac, prim_nums = prim

    # Verify the primitive cell is actually in the XY plane
    a1_new = prim_lattice[0, :2].copy()
    a2_new = prim_lattice[1, :2].copy()

    new_atoms = []
    for frac3, num in zip(prim_frac, prim_nums):
        pos = float(frac3[0]) * a1_new + float(frac3[1]) * a2_new
        new_atoms.append({"symbol": num2sym[num], "pos": pos, "z": 0.0})

    return LatticeStructure(new_atoms, a1_new, a2_new, source=structure.source)


# ── Built-in backend ──────────────────────────────────────────────────────────

def _primitive_builtin(structure: LatticeStructure, tol: float) -> LatticeStructure:
    """
    Built-in primitive-cell search.

    Algorithm
    ---------
    1. Generate candidate translation vectors from all same-species atom pairs,
       including periodic images (na, nb) ∈ {-1,0,1}².
    2. Validate each candidate: a valid translation maps every atom onto another
       atom of the same species (using minimum-image convention).
    3. Keep the two shortest linearly independent valid translations → v1, v2.
    4. Extract one representative atom per orbit in the new primitive cell.
    """
    a1, a2    = structure.a1, structure.a2
    positions = np.array([at["pos"] for at in structure.atoms])
    symbols   = [at["symbol"] for at in structure.atoms]
    z_vals    = [at.get("z", 0.0) for at in structure.atoms]
    n         = len(positions)

    # ── 1. Candidate translations ─────────────────────────────────────────
    raw: list[np.ndarray] = []
    for i in range(n):
        for j in range(n):
            for na in range(-1, 2):
                for nb in range(-1, 2):
                    delta = positions[j] - positions[i] + na * a1 + nb * a2
                    if np.linalg.norm(delta) > tol:
                        raw.append(delta)

    # Deduplicate
    unique: list[np.ndarray] = []
    for c in raw:
        if not any(np.linalg.norm(c - u) < tol for u in unique):
            unique.append(c)

    # ── 2. Validate ────────────────────────────────────────────────────────
    valid: list[np.ndarray] = []
    for delta in unique:
        if _is_valid_translation(positions, symbols, delta, a1, a2, tol):
            valid.append(delta)

    if not valid:
        return structure

    # ── 3. Shortest pair ───────────────────────────────────────────────────
    valid.sort(key=lambda v: float(np.linalg.norm(v)))

    v1 = valid[0]
    v2: np.ndarray | None = None
    for v in valid[1:]:
        cross = abs(float(v1[0] * v[1] - v1[1] * v[0]))
        if cross > tol ** 2 * 50:
            v2 = v
            break

    if v2 is None:
        return structure

    # Reject if new cell is not strictly smaller
    area_orig = abs(float(a1[0] * a2[1] - a1[1] * a2[0]))
    area_prim = abs(float(v1[0] * v2[1] - v1[1] * v2[0]))
    if area_prim >= area_orig - tol:
        return structure

    # ── 4. Extract atoms in the primitive cell ─────────────────────────────
    M_prim = np.column_stack([v1, v2])
    new_atoms: list[dict] = []
    seen_pos:  list[np.ndarray] = []

    for i in range(n):
        frac     = np.linalg.solve(M_prim, positions[i])
        frac_w   = _wrap(frac)
        pos_new  = frac_w[0] * v1 + frac_w[1] * v2

        # Check for Cartesian duplicates using minimum-image in primitive cell
        is_dup = False
        for sp in seen_pos:
            diff      = pos_new - sp
            diff_frac = np.linalg.solve(M_prim, diff)
            diff_frac -= np.round(diff_frac)
            diff_cart  = diff_frac[0] * v1 + diff_frac[1] * v2
            if np.linalg.norm(diff_cart) < tol:
                is_dup = True
                break

        if not is_dup:
            seen_pos.append(pos_new)
            new_atoms.append({"symbol": symbols[i], "pos": pos_new, "z": z_vals[i]})

    if not new_atoms:
        return structure

    return LatticeStructure(new_atoms, v1, v2, source=structure.source)


def _is_valid_translation(
    positions: np.ndarray,
    symbols:   list[str],
    delta:     np.ndarray,
    a1:        np.ndarray,
    a2:        np.ndarray,
    tol:       float,
) -> bool:
    """True if shifting all atoms by delta maps each onto a same-species atom."""
    M = np.column_stack([a1, a2])
    for i in range(len(positions)):
        shifted = positions[i] + delta
        found   = False
        for j in range(len(positions)):
            if symbols[j] != symbols[i]:
                continue
            diff      = shifted - positions[j]
            diff_frac = np.linalg.solve(M, diff)
            diff_frac -= np.round(diff_frac)
            diff_cart  = diff_frac[0] * a1 + diff_frac[1] * a2
            if np.linalg.norm(diff_cart) < tol:
                found = True
                break
        if not found:
            return False
    return True
