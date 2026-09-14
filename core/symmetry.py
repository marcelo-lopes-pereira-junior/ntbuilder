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
import os

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

# Tolerance ladder for the LATTICE, tight to loose: (relative length, degrees).
# A relaxation never lands exactly on the ideal metric -- C2DB's 1AgBrHfIO-1
# has a and b 0.3 % apart and gamma 119.70 deg, which is p3m1 with the
# distortion of the calculation -- so the search asks the same question at
# several tolerances and keeps the highest symmetry it finds, reporting which
# tolerance it needed.  Measured on 3000 catalogue systems this never exceeds
# 0.23 A (median 0.001 A).  An angstrom ladder up to 0.5 A applied to the
# lattice made biphenylene square and promoted 321 of those 3000 systems.
SYM_LADDER = ((1e-4, 0.02), (1e-3, 0.1), (3e-3, 0.3), (1e-2, 0.6), (2e-2, 1.2))

# Tolerance ladder for the ATOMS, in angstroms, as in Materials Studio's Find
# Symmetry.  The one it replaces was a multiple of the lattice tolerance in
# fractional coordinates, which let atoms be matched up to 1.18 A apart (2.2 %
# of the catalogue beyond 0.5 A).  BASIS_MAX_A caps the climb.
BASIS_LADDER_A = (0.001, 0.01, 0.1, 0.5)
BASIS_MAX_A = float(os.environ.get("NTB_BASIS_MAX_A", 0.5))


def symmetry_search(structure: LatticeStructure,
                    ladder=SYM_LADDER) -> tuple[str, float, float]:
    """Highest-symmetry Bravais class within the ladder, and the tolerance for it.

    Returns (class, relative length tolerance, angle tolerance).  Ties go to
    the tightest tolerance that already gives that class, so a cell that is
    exactly hexagonal is reported as hexagonal at 1e-4, not at 2e-2.
    """
    return _symmetry_search_basis(structure, ladder)[:3]


def _symmetry_search_basis(structure: LatticeStructure, ladder=SYM_LADDER):
    """:func:`symmetry_search` plus the basis change that shows the class.

    Every rung is asked in every equivalent basis of the lattice, so the class
    no longer depends on which basis the file used.  Ties keep the basis as
    given and the tightest rung.
    """
    from .io import lattice_type_at, SYM_RANK, equivalent_bases, _metric
    bases = equivalent_bases(structure.a1, structure.a2)
    best = None
    for rel, ang in ladder:
        for k, (M, B) in enumerate(bases):
            a, b, g = _metric(B)
            lt = lattice_type_at(a, b, g, len_tol=rel * 0.5 * (a + b), ang_tol=ang)
            key = (SYM_RANK[lt], k == 0)
            if best is None or key > best[0]:
                best = (key, lt, rel, ang, M)
    return best[1], best[2], best[3], best[4]


def basis_tolerance(structure: LatticeStructure, start_A: float = 0.0,
                    max_A: float | None = None,
                    ladder=BASIS_LADDER_A) -> float:
    """The step of the atom ladder at which the basis shows the most symmetry.

    Climbs the angstrom ladder up to ``max_A`` and keeps the tightest step that
    gives the largest point group, never below the step that covers
    ``start_A`` (how far the lattice snap moved things: a basis test stricter
    than that would reject the symmetry just snapped in).  Stored on the
    structure in fractional units of the longest edge, which is what
    :func:`core.planegroup.structure_point_group` compares; returned in A.
    """
    from .planegroup import structure_point_group
    cap = BASIS_MAX_A if max_A is None else max_A
    steps = [t for t in ladder if t <= cap + 1e-12] or [ladder[0]]
    floor = next((t for t in steps if t >= start_A - 1e-12), steps[-1])
    edge = max(structure.a, structure.b)
    best_tol, best_order = steps[-1], -1
    for tol in steps:
        if tol < floor:
            continue
        structure.sym_tol = tol / edge
        order = len(structure_point_group(structure))
        if order > best_order:
            best_tol, best_order = tol, order
    structure.sym_tol = best_tol / edge
    return best_tol


def snap_to_symmetry(structure: LatticeStructure,
                     ladder=SYM_LADDER) -> tuple[LatticeStructure, str]:
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
    lt, rel_tol, ang_tol, M = _symmetry_search_basis(structure, ladder)
    # The snap is done in the basis that shows the class and carried back to
    # the basis the structure uses, so indices (n, m) keep their meaning.  For
    # M = I this is exactly the previous behaviour.
    A_orig = np.array([structure.a1, structure.a2], dtype=float)
    moved_basis = not np.array_equal(M, np.eye(2, dtype=int))
    B = M @ A_orig
    a1    = B[0].copy()
    a2    = B[1].copy()
    a_len = float(np.linalg.norm(a1))
    b_len = float(np.linalg.norm(a2))
    g_old = math.degrees(math.acos(max(-1.0, min(1.0, float(a1 @ a2) / (a_len * b_len)))))

    if lt == "hexagonal":
        a_ideal = (a_len + b_len) / 2.0
        a1_new  = np.array([a_ideal, 0.0])
        a2_new  = np.array([a_ideal / 2.0, a_ideal * math.sqrt(3.0) / 2.0])
        desc = (
            f"a: {a_len:.5f} / {b_len:.5f} Å → {a_ideal:.5f} Å  |  "
            f"γ: {g_old:.4f}° → 60.0000°"
        )
    elif lt == "square":
        a_ideal = (a_len + b_len) / 2.0
        a1_new  = np.array([a_ideal, 0.0])
        a2_new  = np.array([0.0, a_ideal])
        desc = (
            f"a: {a_len:.5f} / {b_len:.5f} Å → {a_ideal:.5f} Å  |  "
            f"γ: {g_old:.4f}° → 90.0000°"
        )
    elif lt == "rectangular":
        a1_new = np.array([a_len, 0.0])
        a2_new = np.array([0.0, b_len])
        desc = f"γ: {g_old:.4f}° → 90.0000°  (a = {a_len:.5f} Å, b = {b_len:.5f} Å)"
    elif lt == "centred rectangular":
        # a = b is the symmetry; gamma is free and stays.  This branch did not
        # exist -- the rhombic cell fell through to "oblique" and was returned
        # untouched -- and it is the one that matters most for an exact T: the
        # mirror direction (1,1) closes exactly only while a = b holds, and
        # C2DB's WI3 ships a and b agreeing to the 13th decimal and no
        # further, which formally destroys it.
        a_ideal = (a_len + b_len) / 2.0
        g = math.radians(g_old)
        a1_new = np.array([a_ideal, 0.0])
        a2_new = np.array([a_ideal * math.cos(g), a_ideal * math.sin(g)])
        desc = (
            f"a: {a_len:.5f} / {b_len:.5f} Å → {a_ideal:.5f} Å  |  "
            f"γ: {g_old:.4f}° kept (rhombic)"
        )
    else:
        basis_tolerance(structure)
        return structure, "Oblique lattice — nothing to snap."

    if moved_basis:
        Bn = np.array([a1_new, a2_new], dtype=float)
        if np.linalg.det(B) < 0:                 # keep the handedness of the cell
            Bn[:, 1] *= -1.0
        Minv = np.round(np.linalg.inv(M.astype(float))).astype(int)
        An = Minv @ Bn
        th = math.atan2(An[0, 1], An[0, 0])      # a1 back on the x axis
        c, s_ = math.cos(th), math.sin(th)
        An = An @ np.array([[c, -s_], [s_, c]])
        a1_new, a2_new = An[0], An[1]
        a1, a2 = A_orig[0].copy(), A_orig[1].copy()
        desc = f"{desc}  (base equivalente {M.tolist()})"

    # Reproject atomic positions through fractional coordinates
    new_atoms = []
    for atom in structure.atoms:
        frac    = _wrap(_to_frac(atom["pos"], a1, a2))
        new_pos = _to_cart(frac, a1_new, a2_new)
        new_atoms.append({**atom, "pos": new_pos})

    out = LatticeStructure(new_atoms, a1_new, a2_new, source=structure.source)
    moved = rel_tol * max(a_len, b_len)
    tol_atoms = basis_tolerance(out, start_A=moved)
    desc = f"{desc}  [tol {rel_tol * 100:.2f} % / {ang_tol:.2f}°, átomos {tol_atoms:g} Å]"
    return out, desc


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

    def area_per_atom(x) -> float:
        b1 = np.asarray(x.a1, float)[:2]
        b2 = np.asarray(x.a2, float)[:2]
        return abs(float(b1[0] * b2[1] - b1[1] * b2[0])) / max(len(x.atoms), 1)

    # A reduction removes atoms AND area in the same proportion, so the area
    # per atom is conserved exactly.  It is the cheapest possible check and it
    # catches every way this can go wrong: a wrong basis vector, an invented
    # symmetry, a dropped sublattice.  Without it a broken cell reached the
    # builder silently -- the AgBr3 monolayer came back with four times the
    # area per atom and built a single strand instead of a tube.
    ref = area_per_atom(structure)
    candidates = []
    try:
        candidates.append(_primitive_spglib(structure, tol))
    except Exception:
        pass
    try:
        candidates.append(_primitive_builtin(structure, tol))
    except Exception:
        pass

    result = None
    for cand in candidates:
        cand = _orient_and_wrap(cand)
        if abs(area_per_atom(cand) - ref) <= 0.02 * ref:
            result = cand
            break
    if result is None:
        # Nothing trustworthy: hand back the input rather than a wrong cell.
        return structure, ("Primitive-cell search rejected: no candidate "
                           "conserved the area per atom "
                           f"({ref:.3f} A^2/atom). Structure unchanged.")

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
    """Primitive cell via spglib, with the slab's third dimension respected.

    Three things here are not optional, and each was a live bug.

    The out-of-plane offsets go IN.  Forcing every atom to z = 0.5 hands
    spglib a flattened structure, which has symmetry the slab does not: on the
    AgBr3 monolayer the flat projection reduces 16 atoms to 4, while the real
    puckered slab reduces to 8.

    The vacuum axis is identified, not assumed.  spglib is free to return the
    basis in any order, and it does: for that same AgBr3 cell the 30 A vacuum
    vector comes back as row 2, with the two in-plane vectors in rows 1 and 3.
    Taking rows 1 and 2 as the plane adopted the vacuum as a lattice vector --
    a 30 A period with nothing in it -- and the tube built on it came out as a
    single strand.  ``no_idealize=True`` keeps the input Cartesian frame, so
    the vacuum row is the one along the original z and the other two carry
    z = 0 exactly.

    The offsets come back OUT, recovered from the fractional coordinate along
    the vacuum axis.  Returning z = 0 for every atom flattened MoSSe and
    penta-graphene, which silently turns roll_inward into a no-op -- the Janus
    non-equivalence is exactly what those systems are in the test set for.
    """
    import spglib  # raises ImportError if not installed

    a1, a2 = structure.a1, structure.a2
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
        zf = 0.5 + float(atom.get("z", 0.0)) / c_vac
        positions_frac.append([float(f2[0]), float(f2[1]), zf])
        numbers.append(sym2num[atom["symbol"]])

    cell = (lattice, positions_frac, numbers)
    prim = spglib.standardize_cell(cell, symprec=tol,
                                   to_primitive=True, no_idealize=True)
    if prim is None:
        prim = spglib.find_primitive(cell, symprec=tol)
    if prim is None:
        raise RuntimeError("spglib: no primitive cell returned")

    prim_lattice, prim_frac, prim_nums = prim
    prim_lattice = np.asarray(prim_lattice, float)
    prim_frac    = np.asarray(prim_frac, float)

    # Which row is the vacuum?  The one most nearly along the original z.
    zhat = np.array([0.0, 0.0, 1.0])
    k_vac = int(np.argmax(np.abs(prim_lattice @ zhat)))
    plane = [i for i in range(3) if i != k_vac]
    if any(abs(prim_lattice[i, 2]) > 1e-6 * c_vac for i in plane):
        raise RuntimeError("spglib: primitive cell is not a slab in the XY plane")

    i, j = plane
    v1 = prim_lattice[i, :2].copy()
    v2 = prim_lattice[j, :2].copy()
    # Keep the pair right-handed, or the cell comes back mirrored.
    if v1[0] * v2[1] - v1[1] * v2[0] < 0.0:
        i, j = j, i
        v1, v2 = v2, v1

    new_atoms = []
    for frac3, num in zip(prim_frac, prim_nums):
        pos = float(frac3[i]) * v1 + float(frac3[j]) * v2
        z   = (float(frac3[k_vac]) - 0.5) * c_vac
        new_atoms.append({"symbol": num2sym[num], "pos": pos, "z": z})

    return LatticeStructure(new_atoms, v1, v2, source=structure.source)


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
