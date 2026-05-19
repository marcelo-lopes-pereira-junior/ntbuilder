"""
core/io.py
----------
Readers for 2D crystal structures.

Supported formats:
  - XYZ  (.xyz)  — simple, no cell info; cell inferred from data or passed explicitly
  - PDB  (.pdb)  — reads CRYST1 for cell, ATOM/HETATM for positions
  - CIF  (.cif)  — requires `gemmi` (pip install gemmi)

All readers return a LatticeStructure object with:
  - atoms : list of {'symbol': str, 'pos': np.ndarray([x, y]), 'z': float}
  - a1    : np.ndarray([ax, ay])  — first lattice vector (in-plane)
  - a2    : np.ndarray([bx, by])  — second lattice vector (in-plane)
  - source: str  — original file path

Auto-orientation
----------------
If the periodic 2D layer is not parallel to the XY plane (e.g. input file uses
the XZ or YZ plane), the loader automatically permutes coordinates so the layer
normal aligns with Z.  This is detected from the 3D lattice vectors (a1 × a2)
for XYZ files and is a no-op for CIF/PDB (which follow standard conventions).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Data container
# ─────────────────────────────────────────────────────────────────────────────

class LatticeStructure:
    """
    Holds the 2D unit cell of a layered material.

    Each atom is a dict:
        {'symbol': str, 'pos': np.array([x, y]), 'z': float}

    The 'z' key is the out-of-plane (buckling) offset in Å.
    It is 0.0 for flat structures and non-zero for buckled ones
    (e.g. silicene, pentagraphene) or multi-layer materials (e.g. MoS₂).
    When rolling into a nanotube, 'z' is added to (or subtracted from)
    the cylinder radius, placing buckled atoms at different radii.
    """

    _BUCKLING_THRESHOLD = 0.05   # Å — below this is considered flat

    def __init__(
        self,
        atoms: list[dict],
        a1: np.ndarray,
        a2: np.ndarray,
        source: str = "",
    ):
        # Ensure every atom has a 'z' key (default 0.0)
        for atom in atoms:
            atom.setdefault("z", 0.0)
        self.atoms  = atoms   # [{'symbol': str, 'pos': np.array([x, y]), 'z': float}]
        self.source = source

        # ── Normalise hexagonal γ=120° → γ=60° (Dresselhaus convention) ──────
        # Input files from VASP/QE/VESTA often use the γ=120° primitive cell for
        # hexagonal layers.  In that frame, (n,n) falls at 60° from a₁ and is
        # geometrically equivalent to the zigzag direction — NOT armchair.
        # Replacing a₂ → a₁+a₂ converts to γ=60° where:
        #   (n,0) = zigzag (θ=0°),  (n,n) = armchair (θ=30°)  ← standard convention
        # Atom Cartesian positions are stored absolutely and are unaffected.
        a1_ = np.asarray(a1, dtype=float)
        a2_ = np.asarray(a2, dtype=float)
        len1 = float(np.linalg.norm(a1_))
        len2 = float(np.linalg.norm(a2_))
        if len1 > 1e-10 and len2 > 1e-10:
            cos_g   = float(np.dot(a1_, a2_) / (len1 * len2))
            gamma   = math.degrees(math.acos(max(-1.0, min(1.0, cos_g))))
            if abs(len1 - len2) < 1e-3 * len1 and abs(gamma - 120.0) < 0.5:
                a2_ = a1_ + a2_   # γ=120° → γ=60°
        self.a1 = a1_
        self.a2 = a2_

    @property
    def has_buckling(self) -> bool:
        """True if any atom has a significant out-of-plane offset."""
        return any(abs(a["z"]) > self._BUCKLING_THRESHOLD for a in self.atoms)

    @property
    def max_z_offset(self) -> float:
        """Maximum absolute z-offset (Å); 0.0 for flat structures."""
        return max((abs(a["z"]) for a in self.atoms), default=0.0)

    @property
    def d_min(self) -> float:
        """
        Minimum physically meaningful nanotube diameter (Å).

        For buckled/layered structures the innermost atom shell must lie at
        positive radius:  R_inner = R − max_z_offset > 0
        ⟹ D > 2 × max_z_offset.

        Returns 0.0 for flat (non-buckled) structures.
        """
        dz = self.max_z_offset
        return 2.0 * dz if dz > self._BUCKLING_THRESHOLD else 0.0

    # Derived geometry
    @property
    def a(self) -> float:
        return float(np.linalg.norm(self.a1))

    @property
    def b(self) -> float:
        return float(np.linalg.norm(self.a2))

    @property
    def gamma_deg(self) -> float:
        cos_g = np.dot(self.a1, self.a2) / (self.a * self.b)
        return math.degrees(math.acos(np.clip(cos_g, -1, 1)))

    @property
    def lattice_type(self) -> str:
        """
        Crystal system of the 2D lattice.

        'hexagonal'   |a₁|=|a₂| and γ=60° *or* γ=120° (both conventions used
                       in XYZ/CIF files for the hexagonal lattice).
        'rectangular' γ=90° (includes the square, a=b case).
        'oblique'     everything else.
        """
        g = self.gamma_deg
        if abs(self.a - self.b) < 1e-3 and (abs(g - 60) < 0.5 or abs(g - 120) < 0.5):
            return "hexagonal"
        if abs(g - 90) < 0.5:
            return "rectangular"
        return "oblique"

    @property
    def is_square(self) -> bool:
        """True for square lattices (a=b, γ=90°)."""
        return self.lattice_type == "rectangular" and abs(self.a - self.b) < 1e-3

    def __repr__(self) -> str:
        return (
            f"LatticeStructure({len(self.atoms)} atoms | "
            f"a={self.a:.4f} Å, b={self.b:.4f} Å, γ={self.gamma_deg:.2f}° "
            f"[{self.lattice_type}] | {Path(self.source).name})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# XYZ reader
# ─────────────────────────────────────────────────────────────────────────────

def read_xyz(
    path: str | Path,
    a1: np.ndarray | None = None,
    a2: np.ndarray | None = None,
) -> LatticeStructure:
    """
    Read an XYZ file. Lattice vectors must be provided explicitly unless
    a comment line contains 'Lattice="ax ay az bx by bz ..."' (extended XYZ).

    The reader automatically detects the out-of-plane axis from the lattice
    normal (a1 × a2) and permutes coordinates to the XY plane if needed.

    Parameters
    ----------
    path : path to .xyz file
    a1, a2 : lattice vectors (Å, 2D or 3D) — required if not embedded in file
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1]

    # Parse full 3D lattice from extended XYZ (preferred)
    a1_file_3d, a2_file_3d = _parse_extxyz_lattice_3d(comment)

    # Resolve a1
    if a1 is not None:
        a_arr = np.asarray(a1, dtype=float)
        a1_3d = np.array([a_arr[0], a_arr[1], 0.0]) if len(a_arr) == 2 else a_arr[:3]
    elif a1_file_3d is not None:
        a1_3d = a1_file_3d
    else:
        raise ValueError(
            "Lattice vector a1 not found in XYZ file. "
            "Pass a1= and a2= explicitly, or use extended XYZ format."
        )

    # Resolve a2
    if a2 is not None:
        a_arr = np.asarray(a2, dtype=float)
        a2_3d = np.array([a_arr[0], a_arr[1], 0.0]) if len(a_arr) == 2 else a_arr[:3]
    elif a2_file_3d is not None:
        a2_3d = a2_file_3d
    else:
        raise ValueError(
            "Lattice vector a2 not found in XYZ file. "
            "Pass a1= and a2= explicitly, or use extended XYZ format."
        )

    # Read all 3D atom positions
    symbols      = []
    positions_3d = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 3:
            continue
        symbols.append(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3]) if len(parts) >= 4 else 0.0
        positions_3d.append([x, y, z])

    if not symbols:
        raise ValueError("No atoms found in XYZ file.")

    positions_3d = np.array(positions_3d, dtype=float)

    # Auto-orient so the layer lies in the XY plane
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a1_3d, a2_3d, symbols)

    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


def _parse_extxyz_lattice_3d(comment: str):
    """
    Extract full 3D lattice vectors from an extended XYZ comment line.

    Expected format:  Lattice="ax ay az  bx by bz  cx cy cz"

    Returns (a1_3d, a2_3d) as numpy arrays of shape (3,), or (None, None).
    """
    import re
    m = re.search(r'Lattice\s*=\s*"([^"]+)"', comment, re.IGNORECASE)
    if not m:
        return None, None
    vals = list(map(float, m.group(1).split()))
    if len(vals) >= 9:
        return np.array(vals[0:3]), np.array(vals[3:6])
    if len(vals) >= 6:
        return np.array([vals[0], vals[1], 0.0]), np.array([vals[3], vals[4], 0.0])
    return None, None


def _parse_extxyz_lattice(comment: str):
    """Compatibility wrapper — returns 2D (xy) projections of the lattice."""
    a1_3d, a2_3d = _parse_extxyz_lattice_3d(comment)
    if a1_3d is None:
        return None, None
    return a1_3d[:2], a2_3d[:2]


# ─────────────────────────────────────────────────────────────────────────────
# 2D plane auto-orientation
# ─────────────────────────────────────────────────────────────────────────────

def _auto_orient_to_xy(
    positions_3d: np.ndarray,
    a1_3d:        np.ndarray,
    a2_3d:        np.ndarray,
    symbols:      list[str],
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    """
    Detect the out-of-plane axis from the lattice normal (a1 × a2) and
    permute coordinates so that axis becomes Z.

    Parameters
    ----------
    positions_3d : (N, 3) array of Cartesian coordinates
    a1_3d, a2_3d : full 3D lattice vectors
    symbols      : element symbols, length N

    Returns
    -------
    atoms  : list of {'symbol', 'pos': array([x,y]), 'z': float}
    a1_2d  : (2,) in-plane projection of a1
    a2_2d  : (2,) in-plane projection of a2
    """
    # Determine out-of-plane axis from the lattice normal
    normal = np.cross(a1_3d, a2_3d)
    norm   = np.linalg.norm(normal)
    if norm > 1e-10:
        # Axis most aligned with the layer normal
        out_axis = int(np.argmax(np.abs(normal / norm)))
    else:
        # Degenerate lattice → fall back to minimum-spread axis
        spread   = positions_3d.std(axis=0) if len(positions_3d) > 1 else np.zeros(3)
        out_axis = int(np.argmin(spread))

    # Permutation: out_axis → index 2 (Z)
    if out_axis == 2:
        perm = [0, 1, 2]          # already XY — nothing to do
    elif out_axis == 1:
        perm = [0, 2, 1]          # Y out-of-plane: (x,y,z) → (x,z,y)
    else:                          # out_axis == 0
        perm = [1, 2, 0]          # X out-of-plane: (x,y,z) → (y,z,x)

    if out_axis != 2:
        axis_names = ["X", "Y", "Z"]
        print(
            f"[NTBuilder] Auto-rotating structure: "
            f"{axis_names[out_axis]}-axis is out-of-plane. "
            f"Permuting to XY plane."
        )

    pos  = positions_3d[:, perm]  # (N, 3) — new z is out-of-plane
    a1p  = a1_3d[perm]
    a2p  = a2_3d[perm]

    # Centre z so the layer sits at z=0
    z_mean = pos[:, 2].mean()
    z_vals = pos[:, 2] - z_mean

    atoms = [
        {"symbol": sym, "pos": np.array(pos[i, :2], dtype=float), "z": float(z_vals[i])}
        for i, sym in enumerate(symbols)
    ]

    # ── Normalise hexagonal convention ──────────────────────────────────────
    # Some DFT codes (VASP, QE) output hexagonal cells with γ=120° instead of
    # the Dresselhaus γ=60° standard.  In the γ=120° frame, the "armchair"
    # direction (n=m) falls at 60° from a₁ — NOT at 30° — so it is
    # *geometrically equivalent to the zigzag direction*.
    #
    # Fix: if |a₁|≈|a₂| and γ≈120°, replace a₂ → a₁+a₂.  This leaves all
    # Cartesian atom positions unchanged (they are stored absolutely) while
    # converting to the γ=60° basis where:
    #   (n,0)  → zigzag  (θ=0°)
    #   (n,n)  → armchair (θ=30°)   ← matches Dresselhaus convention
    a1_2d = a1p[:2]
    a2_2d = a2p[:2]
    len_a1 = float(np.linalg.norm(a1_2d))
    len_a2 = float(np.linalg.norm(a2_2d))
    if len_a1 > 1e-10 and len_a2 > 1e-10:
        cos_gamma = float(np.dot(a1_2d, a2_2d) / (len_a1 * len_a2))
        gamma_deg = math.degrees(math.acos(max(-1.0, min(1.0, cos_gamma))))
        if abs(len_a1 - len_a2) < 1e-3 * len_a1 and abs(gamma_deg - 120.0) < 0.5:
            a2_2d = a1_2d + a2_2d   # γ=120° → γ=60°

    return atoms, a1_2d.copy(), a2_2d.copy()


# ─────────────────────────────────────────────────────────────────────────────
# PDB reader
# ─────────────────────────────────────────────────────────────────────────────

def read_pdb(path: str | Path) -> LatticeStructure:
    """
    Read a PDB file.
    Lattice vectors are taken from the CRYST1 record (a, b, gamma).
    Only ATOM and HETATM records are parsed; the z coordinate is ignored.
    """
    path = Path(path)
    a1 = a2 = None
    atoms = []

    for line in path.read_text(encoding="utf-8").splitlines():
        tag = line[:6].strip()

        if tag == "CRYST1":
            # CRYST1  a  b  c  alpha  beta  gamma  sGroup  Z
            # PDB fixed-width columns (0-indexed): a[6:15] b[15:24] c[24:33]
            # alpha[33:40] beta[40:47] gamma[47:54] sGroup[55:66]
            a_len   = float(line[6:15])
            b_len   = float(line[15:24])
            gamma   = float(line[47:54])
            gamma_r = math.radians(gamma)
            a1 = np.array([a_len, 0.0])
            a2 = np.array([b_len * math.cos(gamma_r),
                           b_len * math.sin(gamma_r)])

        elif tag in ("ATOM", "HETATM"):
            sym = line[76:78].strip() or line[12:16].strip().lstrip("0123456789")
            x = float(line[30:38])
            y = float(line[38:46])
            atoms.append({"symbol": sym, "pos": np.array([x, y])})

    if a1 is None:
        raise ValueError("No CRYST1 record found in PDB file.")
    if not atoms:
        raise ValueError("No ATOM/HETATM records found in PDB file.")

    return LatticeStructure(atoms, a1, a2, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# CIF symmetry-expansion helpers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_sym_ops(text: str) -> list[str]:
    """
    Extract symmetry operation strings from CIF text.
    Looks for '_symmetry_equiv_pos_as_xyz' or '_space_group_symop_operation_xyz'.
    Returns a list of strings like ['-y,x-y,z', 'x,x-y,-z', ...].
    Falls back to ['x,y,z'] (P1) if none are found.
    """
    import re
    for key in (
        r"_symmetry_equiv_pos_as_xyz",
        r"_space_group_symop_operation_xyz",
    ):
        m = re.search(key + r"\s+((?:.|\n)*?)(?=\n_|\nloop_|\Z)", text)
        if not m:
            continue
        ops = []
        for line in m.group(1).splitlines():
            line = line.strip().strip("'\"")
            if "," in line and not line.startswith("_"):
                ops.append(line)
        if ops:
            return ops
    return ["x,y,z"]


def _expand_asu(
    sym_ops:  list[str],
    asu_fracs: list[tuple[float, float, float]],
    asu_syms:  list[str],
    tol: float = 1e-3,
) -> tuple[list[str], list[tuple[float, float, float]]]:
    """
    Apply every symmetry operation in *sym_ops* to every ASU atom and
    return deduplicated (symbols, fractional-coords) covering the full
    primitive cell.

    Fractional coordinates are wrapped into [0, 1).
    The z-coordinate is filtered to the dominant z-layer (the layer with
    the most atoms after rounding to tol) so that only one 2D layer is kept
    even when the CIF represents a 3D bulk structure.
    """
    import re

    def _eval_op(op_str: str, x: float, y: float, z: float) -> tuple[float, ...]:
        result = []
        for part in op_str.lower().split(","):
            part = part.strip()
            # Evaluate safely: only x, y, z and arithmetic are allowed
            try:
                val = eval(          # noqa: S307  (trusted CIF data)
                    part,
                    {"__builtins__": {}},
                    {"x": x, "y": y, "z": z},
                )
                result.append(float(val) % 1.0)
            except Exception:
                result.append(0.0)
        # Pad to 3 components if a 2D op omitted z
        while len(result) < 3:
            result.append(z % 1.0)
        return tuple(result)

    all_syms:  list[str]                     = []
    all_fracs: list[tuple[float, float, float]] = []

    for sym, (fx, fy, fz) in zip(asu_syms, asu_fracs):
        for op in sym_ops:
            nx, ny, nz = _eval_op(op, fx, fy, fz)
            # Deduplicate
            dup = any(
                abs(nx - ex) < tol and abs(ny - ey) < tol and abs(nz - ez) < tol
                for (ex, ey, ez) in all_fracs
            )
            if not dup:
                all_syms.append(sym)
                all_fracs.append((nx, ny, nz))

    if not all_fracs:
        return asu_syms, list(asu_fracs)

    # Delegate layer selection to the shared helper
    return _select_2d_layer(all_syms, all_fracs)


# ─────────────────────────────────────────────────────────────────────────────
# Shared 2D-layer selection (used by both CIF backends)
# ─────────────────────────────────────────────────────────────────────────────

def _select_2d_layer(
    syms:  list[str],
    fracs: list[tuple[float, float, float]],
) -> tuple[list[str], list[tuple[float, float, float]]]:
    """
    From a fully-expanded set of fractional coordinates, keep only the atoms
    belonging to a single 2D slab.

    Heuristic
    ---------
    • Compute the occupied z-span accounting for periodic wrapping.
    • If span < 0.45 of the c-axis → single slab: keep ALL atoms (their
      different z values encode physical buckling / multi-layer offsets).
    • If span ≥ 0.45 → bulk-like: isolate the densest z-cluster by placing
      the largest gap at the periodic boundary and keeping the remainder.
    """
    if not fracs:
        return syms, fracs

    z_arr = sorted(frac[2] for frac in fracs)

    if len(z_arr) > 1:
        gaps    = [(z_arr[(i + 1) % len(z_arr)] - z_arr[i]) % 1.0
                   for i in range(len(z_arr))]
        max_gap = max(gaps)
    else:
        max_gap = 0.0
    z_spread = 1.0 - max_gap

    if z_spread < 0.45:
        return syms, fracs   # single slab — keep everything

    # Bulk: keep atoms on the dense side of the largest gap
    cut_idx = gaps.index(max_gap)
    z_lo    = z_arr[cut_idx]
    z_hi    = z_arr[(cut_idx + 1) % len(z_arr)]

    cluster_syms:  list[str]   = []
    cluster_fracs: list[tuple] = []
    for sym, frac in zip(syms, fracs):
        z = frac[2]
        if z_lo < z_hi:
            in_gap = z_lo < z <= z_hi
        else:
            in_gap = not (z_hi < z <= z_lo)
        if not in_gap:
            cluster_syms.append(sym)
            cluster_fracs.append(frac)

    return (cluster_syms, cluster_fracs) if cluster_syms else (syms, fracs)


# ─────────────────────────────────────────────────────────────────────────────
# VASP POSCAR / CONTCAR reader
# ─────────────────────────────────────────────────────────────────────────────

def read_poscar(path: str | Path) -> LatticeStructure:
    """
    Read a VASP 5 POSCAR or CONTCAR file.

    Format
    ------
    Line 1  : comment
    Line 2  : universal scale factor
    Lines 3-5: lattice vectors a, b, c (Å × scale)
    Line 6  : element symbols (VASP 5)
    Line 7  : atom counts per species
    Line 8  : 'Direct' or 'Cartesian'
    Lines 9+: fractional or Cartesian coordinates
    """
    path  = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    scale = float(lines[1].split()[0])

    def _vec(line):
        return np.array([float(x) * scale for x in line.split()[:3]])

    a_vec = _vec(lines[2])
    b_vec = _vec(lines[3])
    c_vec = _vec(lines[4])

    # VASP 5: element names on line 6
    species_line = lines[5].split()
    count_line   = lines[6].split()
    species = species_line
    counts  = [int(c) for c in count_line]
    symbols = []
    for sym, cnt in zip(species, counts):
        symbols.extend([sym] * cnt)
    coord_line_idx = 7

    mode = lines[coord_line_idx].strip()[0].lower()   # 'd' or 'c'
    n_atoms = sum(counts)

    positions = []
    for i in range(n_atoms):
        vals = [float(x) for x in lines[coord_line_idx + 1 + i].split()[:3]]
        if mode == 'd':   # fractional → Cartesian
            cart = vals[0] * a_vec + vals[1] * b_vec + vals[2] * c_vec
        else:
            cart = np.array(vals)
        positions.append(cart)

    positions_3d = np.array(positions)
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a_vec, b_vec, symbols)
    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# XSF (XCrysDen) reader
# ─────────────────────────────────────────────────────────────────────────────

def read_xsf(path: str | Path) -> LatticeStructure:
    """
    Read an XSF (XCrysDen Structure File) for a periodic crystal.

    Handles CRYSTAL / PRIMVEC / PRIMCOORD blocks.
    Atomic numbers are converted to element symbols automatically.
    """
    import re

    _Z_TO_SYM = {
        1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O",
        9:"F", 10:"Ne", 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P",
        16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca", 22:"Ti", 23:"V",
        24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 30:"Zn",
        31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 42:"Mo", 46:"Pd",
        47:"Ag", 48:"Cd", 49:"In", 50:"Sn", 74:"W", 78:"Pt", 79:"Au",
        82:"Pb",
    }

    path  = Path(path)
    text  = path.read_text(encoding="utf-8")
    lines = [l.split("#")[0].strip() for l in text.splitlines()]
    lines = [l for l in lines if l]

    def _next_keyword(keyword):
        for i, l in enumerate(lines):
            if l.upper().startswith(keyword.upper()):
                return i
        return None

    # Lattice vectors
    pv = _next_keyword("PRIMVEC")
    if pv is None:
        raise ValueError("No PRIMVEC block found in XSF file.")
    a_vec = np.array([float(x) for x in lines[pv + 1].split()[:3]])
    b_vec = np.array([float(x) for x in lines[pv + 2].split()[:3]])
    c_vec = np.array([float(x) for x in lines[pv + 3].split()[:3]])

    # Atoms
    pc = _next_keyword("PRIMCOORD")
    if pc is None:
        raise ValueError("No PRIMCOORD block found in XSF file.")
    n_atoms = int(lines[pc + 1].split()[0])

    symbols, positions = [], []
    for i in range(n_atoms):
        parts = lines[pc + 2 + i].split()
        # First token: atomic number (int) or element symbol (str)
        tok = parts[0]
        try:
            sym = _Z_TO_SYM.get(int(tok), f"X{tok}")
        except ValueError:
            sym = tok
        symbols.append(sym)
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

    positions_3d = np.array(positions)
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a_vec, b_vec, symbols)
    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS data file reader
# ─────────────────────────────────────────────────────────────────────────────

def read_lammps(path: str | Path) -> LatticeStructure:
    """
    Read a LAMMPS data file (atom_style full or atomic).

    Element names are inferred from inline comments in the Masses section
    (e.g. '# C' or '# Carbon') or from the atomic mass itself.
    """
    _MASS_TO_SYM = {
        1.008: "H",  4.003: "He", 6.941: "Li", 9.012: "Be",
        10.811: "B", 12.011: "C", 14.007: "N", 15.999: "O",
        18.998: "F", 22.990: "Na", 24.305: "Mg", 26.982: "Al",
        28.086: "Si", 30.974: "P", 32.06: "S", 35.453: "Cl",
        39.948: "Ar", 39.098: "K", 40.078: "Ca", 47.867: "Ti",
        50.942: "V",  51.996: "Cr", 54.938: "Mn", 55.845: "Fe",
        58.933: "Co", 58.693: "Ni", 63.546: "Cu", 65.38: "Zn",
        69.723: "Ga", 72.63: "Ge", 74.922: "As", 78.971: "Se",
        79.904: "Br", 95.96: "Mo", 106.42: "Pd", 107.87: "Ag",
        112.41: "Cd", 114.82: "In", 118.71: "Sn", 183.84: "W",
        195.08: "Pt", 196.97: "Au", 207.2: "Pb",
    }

    path  = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    # ── Parse box ──────────────────────────────────────────────────────────────
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0
    for line in lines:
        clean = line.split("#")[0].strip()
        if "xlo xhi" in clean:
            xlo, xhi = float(clean.split()[0]), float(clean.split()[1])
        elif "ylo yhi" in clean:
            ylo, yhi = float(clean.split()[0]), float(clean.split()[1])
        elif "zlo zhi" in clean:
            zlo, zhi = float(clean.split()[0]), float(clean.split()[1])

    Lx = xhi - xlo;  Ly = yhi - ylo;  Lz = zhi - zlo

    # ── Parse Masses section → type_id → element symbol ───────────────────────
    type_sym: dict[int, str] = {}
    in_masses = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower() == "masses":
            in_masses = True;  continue
        if in_masses:
            if stripped == "" or (stripped and stripped[0].isalpha() and stripped != ""):
                if stripped and not stripped[0].isdigit():
                    in_masses = False
                    continue
            if stripped == "":
                continue
            parts   = stripped.split()
            if not parts or not parts[0].isdigit():
                in_masses = False;  continue
            tid  = int(parts[0])
            mass = float(parts[1])
            # Check inline comment for symbol first
            comment_sym = ""
            if "#" in line:
                comment = line.split("#", 1)[1].strip().split()[0] if line.split("#", 1)[1].strip() else ""
                if comment and comment.isalpha() and len(comment) <= 2:
                    comment_sym = comment.capitalize()
            if comment_sym:
                type_sym[tid] = comment_sym
            else:
                # Match by mass (nearest within 0.1 amu)
                best = min(_MASS_TO_SYM.items(), key=lambda kv: abs(kv[0] - mass))
                type_sym[tid] = best[1] if abs(best[0] - mass) < 0.1 else f"T{tid}"

    # ── Parse Atoms section ────────────────────────────────────────────────────
    # atom_style atomic:  atom_id  type  x  y  z
    # atom_style full:    atom_id  mol   type  charge  x  y  z
    in_atoms = False
    raw_atoms: list[tuple[int, float, float, float]] = []  # (type, x, y, z)
    for line in lines:
        stripped = line.split("#")[0].strip()
        if stripped.lower().startswith("atoms") and stripped.lower() != "atom types":
            in_atoms = True;  continue
        if in_atoms:
            if not stripped:
                continue
            parts = stripped.split()
            if not parts[0].isdigit():
                in_atoms = False;  continue
            # Detect style by column count
            if len(parts) >= 7:   # full style: id mol type charge x y z
                tid = int(parts[2])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
            else:                  # atomic style: id type x y z
                tid = int(parts[1])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            raw_atoms.append((tid, x, y, z))

    if not raw_atoms:
        raise ValueError("No Atoms section found in LAMMPS data file.")

    symbols      = [type_sym.get(tid, f"T{tid}") for tid, *_ in raw_atoms]
    positions_3d = np.array([[x - xlo, y - ylo, z - zlo]
                              for _, x, y, z in raw_atoms])

    a_vec = np.array([Lx, 0.0, 0.0])
    b_vec = np.array([0.0, Ly, 0.0])
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a_vec, b_vec, symbols)
    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# Quantum ESPRESSO pw.x input reader
# ─────────────────────────────────────────────────────────────────────────────

def read_qe(path: str | Path) -> LatticeStructure:
    """
    Read a Quantum ESPRESSO pw.x input file (.in / .pwi).

    Parses CELL_PARAMETERS (angstrom or bohr), ATOMIC_SPECIES (for element
    labels), and ATOMIC_POSITIONS (angstrom, bohr, or crystal).
    """
    BOHR = 0.52917721067   # Å per Bohr

    path  = Path(path)
    text  = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    # ── alat from celldm(1) or A= ─────────────────────────────────────────────
    import re
    alat = None
    m = re.search(r"celldm\s*\(\s*1\s*\)\s*=\s*([\d.eEdD+-]+)", text, re.I)
    if m:
        alat = float(m.group(1).replace("d", "e").replace("D", "e")) * BOHR
    else:
        m = re.search(r"\bA\s*=\s*([\d.eEdD+-]+)", text, re.I)
        if m:
            alat = float(m.group(1).replace("d", "e").replace("D", "e"))

    # ── CELL_PARAMETERS ───────────────────────────────────────────────────────
    cp_match = re.search(
        r"CELL_PARAMETERS\s*\{?\s*(\w+)\s*\}?(.*?)(?=ATOMIC_|K_POINTS|\Z)",
        text, re.I | re.S
    )
    if cp_match is None:
        raise ValueError("No CELL_PARAMETERS block found in QE file.")
    cp_unit = cp_match.group(1).lower()
    cp_body = cp_match.group(2).strip().splitlines()
    cp_rows = []
    for line in cp_body:
        nums = line.split()
        if len(nums) >= 3:
            try:
                cp_rows.append([float(x) for x in nums[:3]])
            except ValueError:
                pass
        if len(cp_rows) == 3:
            break
    if len(cp_rows) < 3:
        raise ValueError("Could not parse 3 cell vectors from CELL_PARAMETERS.")
    cp_mat = np.array(cp_rows)

    if cp_unit in ("bohr",):
        cp_mat *= BOHR
    elif cp_unit in ("alat",) and alat is not None:
        cp_mat *= alat
    # else: angstrom — no conversion

    a_vec, b_vec, c_vec = cp_mat[0], cp_mat[1], cp_mat[2]

    # ── ATOMIC_SPECIES → label alias map ──────────────────────────────────────
    # Maps QE labels (e.g. 'C1', 'Mo_sv') → element symbol
    as_match = re.search(
        r"ATOMIC_SPECIES(.*?)(?=ATOMIC_POSITIONS|CELL_PARAMETERS|K_POINTS|\Z)",
        text, re.I | re.S
    )
    label_to_sym: dict[str, str] = {}
    if as_match:
        for line in as_match.group(1).splitlines():
            parts = line.split()
            if len(parts) >= 2:
                label = parts[0]
                # Strip trailing digits/underscores to get element
                sym = re.sub(r"[_\d].*$", "", label)
                label_to_sym[label] = sym.capitalize()

    # ── ATOMIC_POSITIONS ──────────────────────────────────────────────────────
    ap_match = re.search(
        r"ATOMIC_POSITIONS\s*\{?\s*(\w+)\s*\}?(.*?)(?=CELL_PARAMETERS|K_POINTS|\Z)",
        text, re.I | re.S
    )
    if ap_match is None:
        raise ValueError("No ATOMIC_POSITIONS block found in QE file.")
    ap_unit = ap_match.group(1).lower()
    symbols, positions = [], []
    for line in ap_match.group(2).splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            label = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        sym = label_to_sym.get(label, re.sub(r"[_\d].*$", "", label).capitalize())
        symbols.append(sym)

        if ap_unit == "crystal":
            cart = x * a_vec + y * b_vec + z * c_vec
        elif ap_unit == "bohr":
            cart = np.array([x, y, z]) * BOHR
        elif ap_unit in ("alat",) and alat is not None:
            cart = np.array([x, y, z]) * alat
        else:
            cart = np.array([x, y, z])   # angstrom
        positions.append(cart)

    if not symbols:
        raise ValueError("No atomic positions found in QE file.")

    positions_3d = np.array(positions)
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a_vec, b_vec, symbols)
    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# CIF reader  (gemmi preferred, built-in fallback)
# ─────────────────────────────────────────────────────────────────────────────

def read_cif(path: str | Path) -> LatticeStructure:
    """
    Read a CIF file.

    Backend selection
    -----------------
    gemmi (preferred)
        Uses gemmi's own SpaceGroup database to look up every symmetry
        operation for the declared space group, then applies them via
        ``gemmi.Op.apply_to_xyz``.  This handles all 230 space groups
        (including triclinic, oblique, and complex symmetries) correctly.
    Built-in fallback (no gemmi)
        Parses ``_symmetry_equiv_pos_as_xyz`` by regex and evaluates each
        operation with a restricted ``eval``.  Adequate for P1 and simple
        space groups.
    """
    path = Path(path)
    try:
        return _read_cif_gemmi(path)
    except ImportError:
        return _read_cif_builtin(path)


def _read_cif_gemmi(path: Path) -> LatticeStructure:
    """
    CIF reader using gemmi for robust symmetry expansion.

    Key improvement over the previous version
    ------------------------------------------
    Symmetry is now expanded via ``gemmi.make_small_structure_from_block``
    + ``SpaceGroup.operations()``, which uses gemmi's built-in International
    Tables database.  The old custom ``_parse_sym_ops`` / ``_expand_asu``
    regex path is no longer used here, so space groups with non-trivial
    generators (oblique, trigonal, monoclinic with non-standard settings,
    etc.) are handled correctly.
    """
    import re as _re
    try:
        import gemmi
    except ImportError as e:
        raise ImportError(
            "gemmi not available. Install with: pip install gemmi"
        ) from e

    doc   = gemmi.cif.read(str(path))
    block = doc.sole_block()

    # ── Cell parameters ────────────────────────────────────────────────────
    def _cell(key: str) -> float:
        val = block.find_value(key)
        if val is None:
            raise ValueError(f"CIF key '{key}' not found in {path.name}")
        return float(val.strip("()").split("(")[0])   # strip e.s.d. like 2.46(1)

    a_len = _cell("_cell_length_a")
    b_len = _cell("_cell_length_b")
    c_len = _cell("_cell_length_c")
    alpha = math.radians(_cell("_cell_angle_alpha"))
    beta  = math.radians(_cell("_cell_angle_beta"))
    gamma = math.radians(_cell("_cell_angle_gamma"))

    cos_a, cos_b, cos_g = math.cos(alpha), math.cos(beta), math.cos(gamma)
    sin_g = math.sin(gamma)
    cx    = c_len * cos_b
    cy    = c_len * (cos_a - cos_b * cos_g) / (sin_g if abs(sin_g) > 1e-10 else 1.0)
    cz    = math.sqrt(max(c_len**2 - cx**2 - cy**2, 0.0))

    a1_3d = np.array([a_len, 0.0,         0.0])
    a2_3d = np.array([b_len * cos_g, b_len * sin_g, 0.0])

    # ── ASU atoms via gemmi SmallStructure ────────────────────────────────
    small = gemmi.make_small_structure_from_block(block)

    # ── Space-group symmetry operations ───────────────────────────────────
    # gemmi API changed between versions:
    #   < 0.6  : small.find_spacegroup() returns SpaceGroup or None
    #   >= 0.6 : SmallStructure has no find_spacegroup(); use
    #            gemmi.find_spacegroup_by_name(small.spacegroup_hm) instead
    sg = None
    try:
        # old API (gemmi < 0.6)
        sg = small.find_spacegroup()
    except AttributeError:
        # new API (gemmi >= 0.6)
        hm = getattr(small, "spacegroup_hm", None) or getattr(small, "spacegroup", None)
        if hm:
            try:
                sg = gemmi.find_spacegroup_by_name(str(hm))
            except Exception:
                pass
        if sg is None:
            # last resort: scan the CIF block for the H-M symbol directly
            for key in ("_symmetry_space_group_name_H-M",
                        "_space_group_name_H-M_alt",
                        "_symmetry_Int_Tables_number"):
                raw = block.find_value(key)
                if raw:
                    raw = raw.strip().strip("'\"")
                    try:
                        sg = gemmi.find_spacegroup_by_name(raw)
                    except Exception:
                        pass
                    if sg is None:
                        try:
                            sg = gemmi.find_spacegroup_by_number(int(raw))
                        except Exception:
                            pass
                    if sg is not None:
                        break

    if sg is not None:
        ops = list(sg.operations())
        sg_name = sg.hm
    else:
        # P1 fallback: identity only
        ops     = [gemmi.Op()]
        sg_name = "P1 (unknown)"
    print(f"[NTBuilder] CIF space group: {sg_name}  ({len(ops)} operations)")

    # ── Expand ASU using gemmi operations ─────────────────────────────────
    tol = 1e-3

    all_syms:  list[str]   = []
    all_fracs: list[tuple] = []

    for site in small.sites:
        # Prefer type_symbol (clean element); fall back to label stripping
        elem = site.type_symbol.strip()
        if not elem:
            elem = _re.sub(r"[^A-Za-z]", "", site.label)[:2]
        if not elem:
            continue

        fxyz = [site.fract.x, site.fract.y, site.fract.z]

        for op in ops:
            nxyz = op.apply_to_xyz(fxyz)
            nx   = float(nxyz[0]) % 1.0
            ny   = float(nxyz[1]) % 1.0
            nz   = float(nxyz[2]) % 1.0

            # Deduplication with minimum-image distance in fractional space
            dup = False
            for (ex, ey, ez) in all_fracs:
                dx = abs(nx - ex); dx = min(dx, 1.0 - dx)
                dy = abs(ny - ey); dy = min(dy, 1.0 - dy)
                dz = abs(nz - ez); dz = min(dz, 1.0 - dz)
                if dx < tol and dy < tol and dz < tol:
                    dup = True
                    break
            if not dup:
                all_syms.append(elem)
                all_fracs.append((nx, ny, nz))

    if not all_syms:
        raise ValueError(f"No atoms found in CIF file {path.name}.")

    # ── Keep only one 2D slab ─────────────────────────────────────────────
    all_syms, all_fracs = _select_2d_layer(all_syms, all_fracs)

    # ── Fractional → Cartesian ────────────────────────────────────────────
    symbols:      list[str]  = []
    positions_3d: list[list] = []
    for sym, (fx, fy, fz) in zip(all_syms, all_fracs):
        cart_x = fx * a_len + fy * b_len * cos_g + fz * cx
        cart_y =              fy * b_len * sin_g  + fz * cy
        cart_z =                                    fz * cz
        symbols.append(sym)
        positions_3d.append([cart_x, cart_y, cart_z])

    positions_3d_arr = np.array(positions_3d, dtype=float)

    # Auto-orient so the 2D layer lies in XY
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d_arr, a1_3d, a2_3d, symbols)

    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


def _read_cif_builtin(path: Path) -> LatticeStructure:
    """
    Minimal CIF parser (no gemmi required).
    Reads _cell_length_*, _cell_angle_*, _atom_site_fract_* columns, and
    symmetry operations from _symmetry_equiv_pos_as_xyz (or equivalent) so
    that non-P1 space groups (e.g. Materials Studio P6/MMM) are handled.
    """
    import re
    text = path.read_text(encoding="utf-8")

    def _val(key: str) -> float:
        # Accept values that may start with a minus sign
        m = re.search(rf"{re.escape(key)}\s+(-?[\d.]+)", text)
        if not m:
            raise ValueError(f"CIF key '{key}' not found in {path.name}")
        return float(m.group(1))

    a_len   = _val("_cell_length_a")
    b_len   = _val("_cell_length_b")
    gamma   = _val("_cell_angle_gamma")
    gamma_r = math.radians(gamma)

    # c axis — used for fractional-z → Cartesian conversion and the full lattice
    c_len_raw = re.search(r"_cell_length_c\s+(-?[\d.]+)", text)
    c_len     = float(c_len_raw.group(1)) if c_len_raw else 0.0
    alpha_raw = re.search(r"_cell_angle_alpha\s+(-?[\d.]+)", text)
    beta_raw  = re.search(r"_cell_angle_beta\s+(-?[\d.]+)", text)
    alpha_r   = math.radians(float(alpha_raw.group(1))) if alpha_raw else math.pi / 2
    beta_r    = math.radians(float(beta_raw.group(1)))  if beta_raw  else math.pi / 2

    # Full 3D c vector
    cos_a, cos_b, cos_g = math.cos(alpha_r), math.cos(beta_r), math.cos(gamma_r)
    sin_g = math.sin(gamma_r)
    cx = c_len * cos_b
    cy = c_len * (cos_a - cos_b * cos_g) / sin_g if abs(sin_g) > 1e-10 else 0.0
    cz = math.sqrt(max(c_len**2 - cx**2 - cy**2, 0.0))

    a1_3d = np.array([a_len, 0.0, 0.0])
    a2_3d = np.array([b_len * math.cos(gamma_r), b_len * math.sin(gamma_r), 0.0])

    # ── Read ASU atoms from the _atom_site loop ────────────────────────────────
    # Find the _atom_site loop — handle multiple 'loop_' blocks by finding the
    # one that contains _atom_site_fract_x
    loop_pattern = re.compile(
        r"loop_\s*((?:_atom_site_\S+\s*)+)((?:(?!loop_|_\w).*(?:\n|$))*)",
        re.MULTILINE,
    )
    header_block = data_block = None
    for m in loop_pattern.finditer(text):
        if "_atom_site_fract_x" in m.group(1):
            header_block = m.group(1)
            data_block   = m.group(2)
            break

    if header_block is None:
        raise ValueError(f"No _atom_site loop with fractional coords in {path.name}")

    columns = re.findall(r"_atom_site_(\S+)", header_block)

    def _col(name: str) -> int | None:
        try:
            return columns.index(name)
        except ValueError:
            return None

    # Prefer type_symbol over label for element symbols
    sym_col = _col("type_symbol")
    if sym_col is None:
        sym_col = _col("label") or 0
    fx_col  = _col("fract_x")
    fy_col  = _col("fract_y")
    fz_col  = _col("fract_z")

    if fx_col is None or fy_col is None:
        raise ValueError(
            f"CIF file {path.name} does not contain fractional coordinates."
        )

    asu_syms:  list[str] = []
    asu_fracs: list[tuple[float, float, float]] = []

    for line in data_block.splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("_") or parts[0].startswith("#"):
            continue
        if len(parts) <= max(sym_col, fx_col, fy_col):
            continue
        raw_sym = parts[sym_col]
        # Strip digits, underscores, parentheses from label (e.g. "C1" → "C")
        sym = re.sub(r"[0-9_()\[\]]", "", raw_sym)
        if not sym:
            continue
        try:
            fx = float(parts[fx_col])
            fy = float(parts[fy_col])
            fz = float(parts[fz_col]) if fz_col is not None and len(parts) > fz_col else 0.0
        except ValueError:
            continue
        asu_syms.append(sym)
        asu_fracs.append((fx, fy, fz))

    if not asu_syms:
        raise ValueError(f"No atoms parsed from CIF file {path.name}.")

    # ── Expand ASU using symmetry operations ──────────────────────────────────
    sym_ops = _parse_sym_ops(text)
    symbols_exp, fracs_exp = _expand_asu(sym_ops, asu_fracs, asu_syms)

    # ── Convert fractional → Cartesian ────────────────────────────────────────
    symbols      = []
    positions_3d = []
    for sym, (fx, fy, fz) in zip(symbols_exp, fracs_exp):
        cart_x = fx * a_len + fy * b_len * math.cos(gamma_r) + fz * cx
        cart_y =              fy * b_len * math.sin(gamma_r)  + fz * cy
        cart_z =                                                fz * cz
        symbols.append(sym)
        positions_3d.append([cart_x, cart_y, cart_z])

    positions_3d = np.array(positions_3d, dtype=float)

    # Auto-orient so the 2D layer lies in XY (usually a no-op for CIF)
    atoms, a1_2d, a2_2d = _auto_orient_to_xy(positions_3d, a1_3d, a2_3d, symbols)

    return LatticeStructure(atoms, a1_2d, a2_2d, source=str(path))


# ─────────────────────────────────────────────────────────────────────────────
# Unified loader
# ─────────────────────────────────────────────────────────────────────────────

def load_structure(
    path: str | Path,
    a1: np.ndarray | None = None,
    a2: np.ndarray | None = None,
) -> LatticeStructure:
    """
    Auto-detect format and load a 2D crystal structure.

    Parameters
    ----------
    path : file path
    a1, a2 : lattice vectors (only needed for plain XYZ without embedded cell)

    Supported formats
    -----------------
    .xyz                 — plain or extended XYZ
    .pdb                 — Protein Data Bank (CRYST1 record)
    .cif                 — Crystallographic Information File
    .poscar / .contcar / .vasp — VASP POSCAR / CONTCAR
    .xsf                 — XCrysDen Structure File
    .lammps / .data      — LAMMPS data file
    .in / .pwi           — Quantum ESPRESSO pw.x input
    """
    path = Path(path)
    ext  = path.suffix.lower()
    name = path.name.lower()

    if ext == ".xyz":
        return read_xyz(path, a1=a1, a2=a2)
    elif ext == ".pdb":
        return read_pdb(path)
    elif ext == ".cif":
        return read_cif(path)
    elif ext in (".poscar", ".contcar", ".vasp") or name in ("poscar", "contcar"):
        return read_poscar(path)
    elif ext == ".xsf":
        return read_xsf(path)
    elif ext in (".lammps", ".data"):
        return read_lammps(path)
    elif ext in (".in", ".pwi"):
        return read_qe(path)
    else:
        raise ValueError(
            f"Unsupported file format: '{ext}'. "
            "Supported: .xyz, .pdb, .cif, .poscar/.vasp, .xsf, .lammps/.data, .in/.pwi"
        )
