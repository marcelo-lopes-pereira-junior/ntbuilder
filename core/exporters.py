"""
core/exporters.py
-----------------
Export a NanotubeStructure to common simulation formats.

Supported formats
-----------------
  .pdb    — Protein Data Bank (CRYST1 + ATOM records)
  .xyz    — plain XYZ / extended XYZ
  POSCAR  — VASP 5 format
  .lammps — LAMMPS data file (full atom style)
  .pwi    — Quantum ESPRESSO pw.x input skeleton
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .builder import NanotubeStructure


# ─────────────────────────────────────────────────────────────────────────────
# PDB
# ─────────────────────────────────────────────────────────────────────────────

def write_pdb(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write a PDB file with CRYST1 box and ATOM records."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"REMARK  BPN nanotube ({r.n},{r.m}) | "
            f"D={nt.diameter:.4f} Ang | T={nt.length:.4f} Ang | "
            f"atoms={nt.n_atoms} | strain={r.strain:.4f}%\n"
        )
        f.write(
            f"CRYST1{Lx:9.3f}{Ly:9.3f}{Lz:9.3f}"
            f"  90.00  90.00  90.00 P 1           1\n"
        )
        for idx, (sym, (x, y, z)) in enumerate(
            zip(nt.symbols, nt.coords), start=1
        ):
            f.write(
                f"ATOM  {idx:5d}  {sym:<3s} MOL A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {sym:>2s}\n"
            )
        f.write("END\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# XYZ
# ─────────────────────────────────────────────────────────────────────────────

def write_xyz(nt: NanotubeStructure, path: str | Path,
              extended: bool = True) -> Path:
    """
    Write an XYZ file.

    Performance: the body of the file is assembled into a single string
    (atom rows joined with "\\n") and flushed with one ``f.write`` call.
    For nanotubes with tens of thousands of atoms this is ~10× faster
    than the per-row loop, which matters on syncing folders (OneDrive,
    Dropbox) where every write triggers a sync event.

    Parameters
    ----------
    extended : if True, embed lattice in the comment line (extended XYZ).
    """
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    if extended:
        lattice = (
            f'{Lx:.6f} 0.000000 0.000000 '
            f'0.000000 {Ly:.6f} 0.000000 '
            f'0.000000 0.000000 {Lz:.6f}'
        )
        header = (
            f"{nt.n_atoms}\n"
            f'Lattice="{lattice}" '
            f'Properties=species:S:1:pos:R:3 '
            f'nanotube=({r.n},{r.m}) diameter={nt.diameter:.4f} '
            f'strain={r.strain:.4f}\n'
        )
    else:
        header = (
            f"{nt.n_atoms}\n"
            f"({r.n},{r.m}) nanotube | "
            f"D={nt.diameter:.4f} Ang | strain={r.strain:.4f}%\n"
        )

    # Build the atom block as a single string.  Pre-allocate a Python list
    # of size N (avoids list.append re-allocations) and join at the end.
    n = nt.n_atoms
    lines = [None] * n
    coords = nt.coords
    syms   = nt.symbols
    for i in range(n):
        s = syms[i]
        x, y, z = coords[i, 0], coords[i, 1], coords[i, 2]
        lines[i] = f"{s:<4s}  {x:14.8f}  {y:14.8f}  {z:14.8f}"

    body = "\n".join(lines) + "\n"
    path.write_text(header + body, encoding="utf-8")
    return path


# ─────────────────────────────────────────────────────────────────────────────
# VASP POSCAR
# ─────────────────────────────────────────────────────────────────────────────

def write_poscar(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write a VASP 5 POSCAR file."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    # Collect species in order of first appearance
    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)
    counts = [nt.symbols.count(s) for s in species_order]

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"BPN nanotube ({r.n},{r.m}) D={nt.diameter:.4f} Ang "
            f"strain={r.strain:.4f}%\n"
        )
        f.write("1.0\n")
        f.write(f"  {Lx:16.10f}   0.0000000000   0.0000000000\n")
        f.write(f"   0.0000000000  {Ly:16.10f}   0.0000000000\n")
        f.write(f"   0.0000000000   0.0000000000  {Lz:16.10f}\n")
        f.write("  " + "  ".join(species_order) + "\n")
        f.write("  " + "  ".join(str(c) for c in counts) + "\n")
        f.write("Cartesian\n")

        # Write atoms grouped by species
        coords_by_species = {s: [] for s in species_order}
        for sym, coord in zip(nt.symbols, nt.coords):
            coords_by_species[sym].append(coord)
        for s in species_order:
            for x, y, z in coords_by_species[s]:
                f.write(f"  {x:16.10f}  {y:16.10f}  {z:16.10f}\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS data file
# ─────────────────────────────────────────────────────────────────────────────

def write_lammps(nt: NanotubeStructure, path: str | Path) -> Path:
    """
    Write a LAMMPS data file (full atom style).
    Atom types are assigned in order of first appearance.
    """
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)
    type_map = {s: i + 1 for i, s in enumerate(species_order)}

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"# LAMMPS data — BPN nanotube ({r.n},{r.m}) "
            f"D={nt.diameter:.4f} Ang strain={r.strain:.4f}%\n\n"
        )
        f.write(f"{nt.n_atoms} atoms\n")
        f.write(f"{len(species_order)} atom types\n\n")
        f.write(f"0.0 {Lx:.6f} xlo xhi\n")
        f.write(f"0.0 {Ly:.6f} ylo yhi\n")
        f.write(f"0.0 {Lz:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        masses = {"C": 12.011, "N": 14.007, "B": 10.811, "H": 1.008,
                  "O": 15.999, "S": 32.06}
        for s in species_order:
            m = masses.get(s, 1.0)
            f.write(f"  {type_map[s]}  {m:.3f}  # {s}\n")
        f.write("\nAtoms  # full\n\n")
        for idx, (sym, (x, y, z)) in enumerate(
            zip(nt.symbols, nt.coords), start=1
        ):
            # atom_id molecule_id atom_type charge x y z
            f.write(
                f"  {idx}  1  {type_map[sym]}  0.0  "
                f"{x:.8f}  {y:.8f}  {z:.8f}\n"
            )

    return path


# ─────────────────────────────────────────────────────────────────────────────
# Quantum ESPRESSO pw.x skeleton
# ─────────────────────────────────────────────────────────────────────────────

def write_qe(nt: NanotubeStructure, path: str | Path,
             ecutwfc: float = 60.0, ecutrho: float = 480.0) -> Path:
    """Write a Quantum ESPRESSO pw.x input skeleton."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)

    # Bohr conversion
    BOHR = 0.52917721067
    Lx_b, Ly_b, Lz_b = Lx / BOHR, Ly / BOHR, Lz / BOHR

    with path.open("w", encoding="utf-8") as f:
        f.write(f"! BPN nanotube ({r.n},{r.m})\n")
        f.write(f"! D={nt.diameter:.4f} Ang | strain={r.strain:.4f}%\n\n")
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write("  prefix      = 'nanotube'\n")
        f.write("/\n\n")
        f.write("&SYSTEM\n")
        f.write(f"  ibrav      = 0\n")
        f.write(f"  nat        = {nt.n_atoms}\n")
        f.write(f"  ntyp       = {len(species_order)}\n")
        f.write(f"  ecutwfc    = {ecutwfc}\n")
        f.write(f"  ecutrho    = {ecutrho}\n")
        f.write("/\n\n")
        f.write("&ELECTRONS\n")
        f.write("  conv_thr   = 1.0d-8\n")
        f.write("/\n\n")
        f.write("CELL_PARAMETERS angstrom\n")
        f.write(f"  {Lx:.10f}   0.0000000000   0.0000000000\n")
        f.write(f"   0.0000000000  {Ly:.10f}   0.0000000000\n")
        f.write(f"   0.0000000000   0.0000000000  {Lz:.10f}\n\n")
        f.write("ATOMIC_SPECIES\n")
        masses = {"C": 12.011, "N": 14.007, "B": 10.811, "H": 1.008,
                  "O": 15.999, "S": 32.06}
        pseudo = {"C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                  "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF",
                  "B": "B.pbe-n-kjpaw_psl.1.0.0.UPF"}
        for s in species_order:
            m = masses.get(s, 1.0)
            pp = pseudo.get(s, f"{s}.pseudo.UPF")
            f.write(f"  {s}  {m:.3f}  {pp}\n")
        f.write("\nATOMIC_POSITIONS angstrom\n")
        for sym, (x, y, z) in zip(nt.symbols, nt.coords):
            f.write(f"  {sym}  {x:.10f}  {y:.10f}  {z:.10f}\n")
        f.write("\nK_POINTS gamma\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# XSF (XCrysDen Structure File)
# ─────────────────────────────────────────────────────────────────────────────

_SYM_TO_Z = {
    "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9,
    "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17,
    "Ar":18, "K":19, "Ca":20, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26,
    "Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga":31, "Ge":32, "As":33, "Se":34,
    "Br":35, "Mo":42, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50,
    "W":74, "Pt":78, "Au":79, "Pb":82,
}

def write_xsf(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write an XSF file (XCrysDen / VESTA compatible)."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    with path.open("w", encoding="utf-8") as f:
        f.write(f"# NTBuilder — nanotube ({r.n},{r.m})\n")
        f.write(f"# D={nt.diameter:.4f} Å  strain={r.strain:.4f}%\n\n")
        f.write("CRYSTAL\n\n")
        f.write("PRIMVEC\n")
        f.write(f"  {Lx:.10f}  0.0000000000  0.0000000000\n")
        f.write(f"  0.0000000000  {Ly:.10f}  0.0000000000\n")
        f.write(f"  0.0000000000  0.0000000000  {Lz:.10f}\n\n")
        f.write("PRIMCOORD\n")
        f.write(f"  {nt.n_atoms}  1\n")
        for sym, (x, y, z) in zip(nt.symbols, nt.coords):
            Z = _SYM_TO_Z.get(sym, 6)   # default to C if unknown
            f.write(f"  {Z:3d}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# CP2K input skeleton
# ─────────────────────────────────────────────────────────────────────────────

def write_cp2k(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write a CP2K input skeleton (&FORCE_EVAL / &SUBSYS section)."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)

    masses = {"C":12.011, "N":14.007, "B":10.811, "H":1.008,
              "O":15.999, "S":32.06,  "Si":28.086,"Mo":95.96}

    with path.open("w", encoding="utf-8") as f:
        f.write(f"# NTBuilder — nanotube ({r.n},{r.m})\n")
        f.write(f"# D={nt.diameter:.4f} Å  strain={r.strain:.4f}%\n\n")
        f.write("&GLOBAL\n")
        f.write("  PROJECT  nanotube\n")
        f.write("  RUN_TYPE  ENERGY\n")
        f.write("&END GLOBAL\n\n")
        f.write("&FORCE_EVAL\n")
        f.write("  METHOD  Quickstep\n")
        f.write("  &SUBSYS\n")
        f.write("    &CELL\n")
        f.write(f"      A  {Lx:.10f}  0.0  0.0\n")
        f.write(f"      B  0.0  {Ly:.10f}  0.0\n")
        f.write(f"      C  0.0  0.0  {Lz:.10f}\n")
        f.write("      PERIODIC  XYZ\n")
        f.write("    &END CELL\n")
        f.write("    &COORD\n")
        for sym, (x, y, z) in zip(nt.symbols, nt.coords):
            f.write(f"      {sym:<4s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")
        f.write("    &END COORD\n")
        for sym in species_order:
            f.write(f"    &KIND {sym}\n")
            f.write(f"      ELEMENT  {sym}\n")
            f.write(f"      MASS  {masses.get(sym, 1.0):.3f}\n")
            f.write( "      BASIS_SET  DZVP-MOLOPT-GTH\n")
            f.write( "      POTENTIAL  GTH-PBE\n")
            f.write( "    &END KIND\n")
        f.write("  &END SUBSYS\n")
        f.write("&END FORCE_EVAL\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# SIESTA .fdf skeleton
# ─────────────────────────────────────────────────────────────────────────────

def write_siesta(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write a SIESTA .fdf input skeleton."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = nt.box

    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)
    type_map = {s: i + 1 for i, s in enumerate(species_order)}

    # Atomic numbers for ChemicalSpeciesLabel
    _Z = {"H":1,"C":6,"N":7,"O":8,"B":5,"Si":14,"S":16,"Mo":42,"P":15,"F":9}

    with path.open("w", encoding="utf-8") as f:
        f.write(f"# NTBuilder — nanotube ({r.n},{r.m})\n")
        f.write(f"# D={nt.diameter:.4f} Å  strain={r.strain:.4f}%\n\n")
        f.write(f"SystemLabel  nanotube_{r.n}_{r.m}\n")
        f.write(f"NumberOfAtoms  {nt.n_atoms}\n")
        f.write(f"NumberOfSpecies  {len(species_order)}\n\n")
        f.write("%block LatticeVectors\n")
        f.write(f"  {Lx:.10f}  0.000000  0.000000\n")
        f.write(f"  0.000000  {Ly:.10f}  0.000000\n")
        f.write(f"  0.000000  0.000000  {Lz:.10f}\n")
        f.write("%endblock LatticeVectors\n\n")
        f.write("AtomicCoordinatesFormat  Ang\n\n")
        f.write("%block ChemicalSpeciesLabel\n")
        for sym in species_order:
            f.write(f"  {type_map[sym]}  {_Z.get(sym, 6)}  {sym}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\n")
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for sym, (x, y, z) in zip(nt.symbols, nt.coords):
            f.write(f"  {x:14.8f}  {y:14.8f}  {z:14.8f}  {type_map[sym]}\n")
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# CIF writer
# ─────────────────────────────────────────────────────────────────────────────

def write_cif(nt: NanotubeStructure, path: str | Path) -> Path:
    """Write a CIF file (space group P1, Cartesian → fractional)."""
    path = Path(path)
    r    = nt.chirality
    Lx, Ly, Lz = [float(x) for x in nt.box]

    species_order: list[str] = []
    for s in nt.symbols:
        if s not in species_order:
            species_order.append(s)

    with path.open("w", encoding="utf-8") as f:
        f.write(f"data_nanotube_{r.n}_{r.m}\n\n")
        f.write(f"_cell_length_a   {Lx:.6f}\n")
        f.write(f"_cell_length_b   {Ly:.6f}\n")
        f.write(f"_cell_length_c   {Lz:.6f}\n")
        f.write( "_cell_angle_alpha   90.00000\n")
        f.write( "_cell_angle_beta    90.00000\n")
        f.write( "_cell_angle_gamma   90.00000\n\n")
        f.write( "_symmetry_space_group_name_H-M   'P 1'\n")
        f.write( "_symmetry_Int_Tables_number        1\n\n")
        f.write( "loop_\n")
        f.write( "  _symmetry_equiv_pos_as_xyz\n")
        f.write( "  'x,y,z'\n\n")
        f.write( "loop_\n")
        f.write( "  _atom_site_label\n")
        f.write( "  _atom_site_type_symbol\n")
        f.write( "  _atom_site_fract_x\n")
        f.write( "  _atom_site_fract_y\n")
        f.write( "  _atom_site_fract_z\n")
        counts: dict[str, int] = {}
        for sym, (x, y, z) in zip(nt.symbols, nt.coords):
            counts[sym] = counts.get(sym, 0) + 1
            label = f"{sym}{counts[sym]}"
            fx, fy, fz = x / Lx, y / Ly, z / Lz
            f.write(f"  {label:<8s}  {sym:<4s}  {fx:.8f}  {fy:.8f}  {fz:.8f}\n")

    return path


# ─────────────────────────────────────────────────────────────────────────────
# Unified export
# ─────────────────────────────────────────────────────────────────────────────

_WRITERS = {
    ".pdb":    write_pdb,
    ".xyz":    write_xyz,
    ".poscar": write_poscar,
    "poscar":  write_poscar,
    ".vasp":   write_poscar,
    "vasp":    write_poscar,
    ".lammps": write_lammps,
    ".data":   write_lammps,
    ".pwi":    write_qe,
    ".in":     write_qe,
    ".xsf":    write_xsf,
    "xsf":     write_xsf,
    ".inp":    write_cp2k,
    "cp2k":    write_cp2k,
    ".fdf":    write_siesta,
    "siesta":  write_siesta,
    ".cif":    write_cif,
    "cif":     write_cif,
}


def export(
    nt: NanotubeStructure,
    path: str | Path,
    fmt: str | None = None,
) -> Path:
    """
    Export a NanotubeStructure to the specified file.

    Parameters
    ----------
    nt   : nanotube structure
    path : output file path
    fmt  : format override (e.g. 'poscar', '.pdb'); inferred from extension if None
    """
    path = Path(path)
    key  = (fmt or path.suffix).lower()

    writer = _WRITERS.get(key)
    if writer is None:
        raise ValueError(
            f"Unknown format '{key}'. "
            f"Supported: {sorted(_WRITERS.keys())}"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    return writer(nt, path)
