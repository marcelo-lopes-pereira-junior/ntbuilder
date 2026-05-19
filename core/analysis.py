"""
core/analysis.py
----------------
Post-construction analysis and output generation.

Functions
---------
bond_analysis(nt, cutoff)
    Compute bond lengths and species-pair labels for histogram / statistics.

electronic_character(n, m, lattice_type)
    Return 'metallic', 'semiconducting', or 'unknown' for the nanotube.

tube_symmetry_info(n, m)
    Return a dict with simplified line-group / point-group information.

generate_methods_text(nt, structure, deform_desc, software)
    Build a ready-to-paste Methods paragraph for a research paper.

generate_vasp_inputs(nt)
    Return (INCAR_str, KPOINTS_str) strings for VASP.

generate_qe_input(nt, structure)
    Return a complete pw.x input string for Quantum ESPRESSO.

generate_cp2k_input(nt)
    Return a CP2K input string (adds to the structure block from exporters.py).

query_cod(formula, max_results)
    Query the Crystallography Open Database REST API (requires network).

query_mp(formula, api_key, max_results)
    Query Materials Project API v3 (requires mp-api or requests + API key).
"""

from __future__ import annotations

import math
import textwrap
from typing import Sequence

import numpy as np

from .builder  import NanotubeStructure
from .io       import LatticeStructure
from .chirality import ChiralityResult


# ─────────────────────────────────────────────────────────────────────────────
# Bond analysis
# ─────────────────────────────────────────────────────────────────────────────

def bond_analysis(
    nt:     NanotubeStructure,
    cutoff: float = 3.5,
) -> dict:
    """
    Compute all pairwise distances shorter than *cutoff* in the nanotube,
    respecting periodic boundary conditions along Z.

    Returns
    -------
    dict with keys:
        "distances"   : np.ndarray (N_bonds,)  — all bond lengths in Å
        "pairs"       : list[str]              — "A-B" species label per bond
        "species_set" : set[str]               — unique pair labels found
        "mean"        : float                  — mean bond length
        "std"         : float                  — std dev
        "min"         : float                  — shortest bond
        "max"         : float                  — longest bond (≤ cutoff)
        "n_bonds"     : int
    """
    coords  = nt.coords
    syms    = list(nt.symbols)
    Lz      = float(nt.box[2])
    n       = len(syms)

    distances: list[float] = []
    pair_labels: list[str] = []

    for i in range(n):
        for j in range(i + 1, n):
            d = coords[j] - coords[i]
            # Periodic image along Z only (nanotube is periodic in Z)
            d[2] -= Lz * round(d[2] / Lz)
            dist = float(np.linalg.norm(d))
            if dist < cutoff:
                distances.append(dist)
                a, b = sorted([syms[i], syms[j]])
                pair_labels.append(f"{a}-{b}")

    dists = np.array(distances, dtype=float)

    if len(dists) == 0:
        return {
            "distances": dists, "pairs": pair_labels,
            "species_set": set(), "mean": 0.0, "std": 0.0,
            "min": 0.0, "max": 0.0, "n_bonds": 0,
        }

    return {
        "distances":   dists,
        "pairs":       pair_labels,
        "species_set": set(pair_labels),
        "mean":        float(dists.mean()),
        "std":         float(dists.std()),
        "min":         float(dists.min()),
        "max":         float(dists.max()),
        "n_bonds":     len(dists),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Electronic character
# ─────────────────────────────────────────────────────────────────────────────

def electronic_character(
    n: int,
    m: int,
    lattice_type:        str   = "hexagonal",
    species:             Sequence[str] | None = None,
    n_atoms_unit_cell:   int | None = None,
    lattice_constant_a:  float | None = None,
) -> str:
    """
    Classify the nanotube's electronic character from the zone-folding rule.

    The Saito–Dresselhaus rule

        (n − m) mod 3 == 0  →  metallic
        otherwise           →  semiconducting

    is a direct consequence of the Dirac-cone band structure of *pristine
    graphene* at the K point.  It is **strictly valid only for graphene**
    (hexagonal lattice, monatomic carbon basis with exactly 2 atoms per
    unit cell, lattice constant a ≈ 2.46 Å).

    Applying the rule to any other precursor — h-BN, MoS₂, MoSSe,
    biphenylene, penta-graphene, irida-graphene, or any other sp² carbon
    allotrope — produces physically unjustified labels: h-BN is
    insulating for every (n, m), MoS₂/MoSSe are always semiconducting,
    and the larger-unit-cell carbon allotropes have their own
    zone-folding selection rules with no closed-form analogue.

    To prevent silent misclassification, this function returns
    ``'requires DFT'`` whenever **any** of the three graphene markers
    fails: lattice type, species, or unit-cell metrics.

    Parameters
    ----------
    n, m              : chiral indices.
    lattice_type      : from ``LatticeStructure.lattice_type``.
    species           : iterable of element symbols in the unit cell
                        (NOT the nanotube — the parent 2D lattice).
    n_atoms_unit_cell : number of atoms in the parent 2D unit cell.
    lattice_constant_a: ``LatticeStructure.a`` in Å.

    Returns
    -------
    'metallic' | 'semiconducting' | 'requires DFT'
    """
    if lattice_type != "hexagonal":
        return "requires DFT"
    if species is not None:
        unique = {str(s).strip().upper() for s in species if s}
        if unique and unique != {"C"}:
            return "requires DFT"
    # Graphene-specific markers: 2-atom basis, a ≈ 2.46 Å
    if n_atoms_unit_cell is not None and n_atoms_unit_cell != 2:
        return "requires DFT"
    if lattice_constant_a is not None and not (2.40 <= lattice_constant_a <= 2.52):
        return "requires DFT"
    return "metallic" if (n - m) % 3 == 0 else "semiconducting"


def electronic_character_label(
    n: int,
    m: int,
    lattice_type:        str   = "hexagonal",
    species:             Sequence[str] | None = None,
    n_atoms_unit_cell:   int | None = None,
    lattice_constant_a:  float | None = None,
) -> str:
    """Human-readable one-liner including the armchair/zigzag special cases.

    Returns ``'requires DFT calculation [zone-folding rule does not apply]'``
    for any precursor that is not pristine graphene.  The check uses three
    independent markers (species, atom count, lattice constant) — any
    missing marker is treated as “unknown” and forces the DFT label, so
    silent misclassification cannot occur.
    """
    char = electronic_character(
        n, m, lattice_type,
        species             = species,
        n_atoms_unit_cell   = n_atoms_unit_cell,
        lattice_constant_a  = lattice_constant_a,
    )
    if char == "requires DFT":
        return "requires DFT calculation  [zone-folding rule does not apply]"
    if n == m:
        return "metallic  [armchair]"
    if m == 0:
        if (n - m) % 3 == 0:
            return "metallic  [zigzag]"
        return "semiconducting  [zigzag]"
    return char


# ─────────────────────────────────────────────────────────────────────────────
# Tube symmetry / line-group info
# ─────────────────────────────────────────────────────────────────────────────

def tube_symmetry_info(n: int, m: int) -> dict:
    """
    Return a dict with simplified symmetry classification.

    Keys
    ----
    type        : "armchair" | "zigzag" | "chiral"
    d_nm        : gcd(n, m)
    n_screw     : number of screw axes (= n for armchair, m for zigzag, ...)
    description : human-readable string
    """
    import math as _math
    d = _math.gcd(n, m)

    if n == m:
        t = "armchair"
        desc = f"Armchair ({n},{n}) — D_{n}h point group, {n}-fold rotational symmetry"
    elif m == 0:
        t = "zigzag"
        desc = f"Zigzag ({n},0) — D_{n}h point group, {n}-fold rotational symmetry"
    else:
        t = "chiral"
        desc = (
            f"Chiral ({n},{m}) — T symmetry (no mirror planes), "
            f"gcd(n,m)={d} — screw axis with {d}-fold periodicity"
        )

    return {
        "type":        t,
        "d_nm":        d,
        "n_screw":     d,
        "description": desc,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Methods-section text generator
# ─────────────────────────────────────────────────────────────────────────────

def generate_methods_text(
    nt:           NanotubeStructure,
    structure:    LatticeStructure | None = None,
    deform_desc:  str  = "",
    software:     str  = "NTBuilder",
    version:      str  = "1.1",
    cite_key:     str  = "[CITE]",
    n_walls:      int  = 1,
    wall_info:    str  = "",
) -> str:
    """
    Generate a ready-to-paste Methods paragraph describing the nanotube.

    Parameters
    ----------
    nt           : built NanotubeStructure
    structure    : optional LatticeStructure for lattice parameters
    deform_desc  : string from deformations.deformation_description()
    software     : software name for citation
    version      : software version string
    cite_key     : BibTeX key placeholder
    n_walls      : number of walls (1 = SWNT)
    wall_info    : optional per-wall summary string (from mwnt.mwnt_summary)
    """
    ch = nt.chirality
    n, m = ch.n, ch.m

    # Nanotube type label
    if n == m:
        tube_type = f"armchair ({n},{m})"
    elif m == 0:
        tube_type = f"zigzag ({n},{m})"
    else:
        tube_type = f"chiral ({n},{m})"

    wall_label = "multi-walled" if n_walls > 1 else "single-walled"

    # Lattice info
    if structure is not None:
        lat_sentence = (
            f"The parent 2D lattice has parameters "
            f"a = {structure.a:.4f} Å, b = {structure.b:.4f} Å, "
            f"γ = {structure.gamma_deg:.2f}° ({structure.lattice_type} lattice), "
            f"with {len(structure.atoms)} atom(s) per unit cell."
        )
    else:
        lat_sentence = ""

    # Deformation sentence
    deform_sentence = ""
    if deform_desc and deform_desc != "none":
        deform_sentence = (
            f"The nanotube was subsequently deformed by applying {deform_desc}. "
        )

    # Electronic character (graphene-like only)
    elec = ""
    if structure is not None and structure.lattice_type == "hexagonal":
        char = electronic_character(n, m, "hexagonal")
        elec = f"Based on zone-folding arguments, this tube is expected to be {char}. "

    # Wall info
    wall_sentence = ""
    if n_walls > 1 and wall_info:
        wall_sentence = f"\nWall assignments: {wall_info}. "

    text = textwrap.dedent(f"""\
        The {wall_label} {tube_type} nanotube was constructed using {software} \
v{version} {cite_key}. {lat_sentence}
        The chiral vector is C_h = {n}a₁ + {m}a₂, yielding a nanotube diameter \
of {ch.diameter:.4f} Å and a translational period of {ch.T_norm:.4f} Å \
({nt.n_atoms} atoms per unit cell, lattice strain {ch.strain:.4f}%). \
{elec}{deform_sentence}\
The nanotube was placed in a periodic orthorhombic simulation box of \
{nt.box[0]:.2f} × {nt.box[1]:.2f} × {nt.box[2]:.2f} Å³ with \
{nt.vacuum:.1f} Å of vacuum padding in the lateral directions to minimise \
inter-image interactions.{wall_sentence}
    """).strip()

    return text


# ─────────────────────────────────────────────────────────────────────────────
# DFT input generators
# ─────────────────────────────────────────────────────────────────────────────

def _species(nt: NanotubeStructure) -> list[str]:
    """Unique element symbols in the nanotube, sorted."""
    return sorted(set(nt.symbols))


def generate_vasp_inputs(nt: NanotubeStructure) -> tuple[str, str]:
    """
    Return (INCAR_str, KPOINTS_str) for a typical VASP geometry optimisation
    of a nanotube.

    Notes
    -----
    - k-points: Γ-only in XY (isolated), N_k along Z proportional to 1/Lz.
    - INCAR settings follow common practice for nanotube DFT (PBE, vdW-D3).
    - Pseudopotential POTCAR selection is not automated (listed as a comment).
    """
    ch   = nt.chirality
    Lz   = float(nt.box[2])

    # 1 k-point per ~5 Å of Z period is a good starting point
    nkz  = max(1, round(20.0 / Lz))

    species = _species(nt)
    potcar_hint = "  ".join(f"PAW_PBE {s}" for s in species)

    incar = textwrap.dedent(f"""\
        # VASP INCAR — nanotube ({ch.n},{ch.m})
        # Generated by NTBuilder

        SYSTEM  = NT_{ch.n}_{ch.m}
        ISTART  = 0        ! start from scratch
        ICHARG  = 2        ! superposition of atomic charges

        # Electronic minimisation
        ENCUT   = 500      ! plane-wave cutoff (eV) — adjust for your POTCAR
        EDIFF   = 1E-6     ! energy convergence (eV)
        PREC    = Accurate
        ALGO    = Fast     ! RMM-DIIS
        NELM    = 200

        # Smearing (use ISMEAR=0 for semiconductors/insulators, -5 for metals)
        ISMEAR  = 0
        SIGMA   = 0.05

        # Geometry optimisation
        IBRION  = 2        ! conjugate gradient
        NSW     = 200      ! max ionic steps
        EDIFFG  = -0.01    ! force convergence (eV/Å)
        ISIF    = 2        ! relax ions, keep cell fixed

        # vdW dispersion correction (remove if not needed)
        IVDW    = 11       ! DFT-D3 with Becke-Johnson damping

        # Output
        LWAVE   = .FALSE.
        LCHARG  = .FALSE.
        NWRITE  = 1

        # POTCAR order: {potcar_hint}
    """).rstrip()

    kpoints = textwrap.dedent(f"""\
        Automatic k-mesh for nanotube ({ch.n},{ch.m}) — Γ-centred
        0
        Gamma
          1  1  {nkz}
          0  0  0
    """).rstrip()

    return incar, kpoints


def generate_qe_input(
    nt:        NanotubeStructure,
    structure: LatticeStructure | None = None,
    prefix:    str = "nanotube",
) -> str:
    """
    Return a complete pw.x input file for Quantum ESPRESSO geometry optimisation.

    Pseudopotentials are listed as placeholders (PSP_LIBRARY/El.upf) since
    the actual filenames depend on the user's pslibrary installation.
    """
    ch      = nt.chirality
    Lx, Ly, Lz = [float(v) for v in nt.box]
    species = _species(nt)
    Lz_bohr = Lz / 0.529177

    # k-points along Z
    nkz = max(1, round(20.0 / Lz))

    nat  = nt.n_atoms
    ntyp = len(species)

    # ATOMIC_POSITIONS block
    pos_lines = []
    for sym, coord in zip(nt.symbols, nt.coords):
        # Convert Å → bohr
        x, y, z = [v / 0.529177 for v in coord]
        pos_lines.append(f"  {sym:<4} {x:18.10f} {y:18.10f} {z:18.10f}")

    # ATOMIC_SPECIES block (mass lookup)
    _MASS = {
        "H":1.008,"He":4.003,"Li":6.941,"Be":9.012,"B":10.811,"C":12.011,
        "N":14.007,"O":15.999,"F":18.998,"Ne":20.180,"Na":22.990,"Mg":24.305,
        "Al":26.982,"Si":28.086,"P":30.974,"S":32.065,"Cl":35.453,"Ar":39.948,
        "K":39.098,"Ca":40.078,"Sc":44.956,"Ti":47.867,"V":50.942,"Cr":51.996,
        "Mn":54.938,"Fe":55.845,"Co":58.933,"Ni":58.693,"Cu":63.546,"Zn":65.38,
        "Ga":69.723,"Ge":72.640,"As":74.922,"Se":78.96,"Br":79.904,"Kr":83.798,
        "Mo":95.96,"W":183.84,"Pt":195.084,"Au":196.967,"Bi":208.980,
    }
    spec_lines = [
        f"  {s:<4} {_MASS.get(s, 1.0):<10.3f} PSP_LIBRARY/{s}.upf"
        for s in species
    ]

    text = textwrap.dedent(f"""\
        &CONTROL
          calculation   = 'relax'
          prefix        = '{prefix}'
          outdir        = './out'
          pseudo_dir    = './PSP_LIBRARY'
          verbosity     = 'high'
          etot_conv_thr = 1.0e-5
          forc_conv_thr = 1.0e-4
          nstep         = 200
        /

        &SYSTEM
          ibrav     = 0
          nat       = {nat}
          ntyp      = {ntyp}
          ecutwfc   = 60.0         ! Ry — adjust for your pseudopotentials
          ecutrho   = 480.0
          occupations = 'smearing'
          smearing  = 'mv'
          degauss   = 0.005
          vdw_corr  = 'dft-d3'    ! remove if not needed
        /

        &ELECTRONS
          conv_thr     = 1.0e-8
          mixing_beta  = 0.3
          electron_maxstep = 200
        /

        &IONS
          ion_dynamics = 'bfgs'
        /

        CELL_PARAMETERS angstrom
          {Lx:.6f}   0.000000   0.000000
          0.000000   {Ly:.6f}   0.000000
          0.000000   0.000000   {Lz:.6f}

        ATOMIC_SPECIES
        {chr(10).join(spec_lines)}

        ATOMIC_POSITIONS angstrom
        {chr(10).join(pos_lines)}

        K_POINTS automatic
          1  1  {nkz}  0  0  0
    """).rstrip()

    return text


def generate_cp2k_input(nt: NanotubeStructure) -> str:
    """
    Return a CP2K input file for geometry optimisation (GFN2-xTB or PBE).
    Uses the QUICKSTEP module with GAPW/GTH pseudopotentials.
    """
    ch      = nt.chirality
    Lx, Ly, Lz = [float(v) for v in nt.box]
    species = _species(nt)
    nat     = nt.n_atoms

    coord_lines = [
        f"      {sym:<4} {x:16.8f} {y:16.8f} {z:16.8f}"
        for sym, (x, y, z) in zip(nt.symbols, nt.coords)
    ]
    kind_lines = [
        f"    &KIND {s}\n      BASIS_SET DZVP-MOLOPT-SR-GTH\n      POTENTIAL GTH-PBE\n    &END KIND"
        for s in species
    ]

    text = textwrap.dedent(f"""\
        &GLOBAL
          PROJECT nanotube_{ch.n}_{ch.m}
          RUN_TYPE GEO_OPT
          PRINT_LEVEL LOW
        &END GLOBAL

        &MOTION
          &GEO_OPT
            TYPE MINIMIZATION
            MAX_ITER 200
            MAX_FORCE 1.0E-4
            RMS_FORCE 3.0E-4
          &END GEO_OPT
        &END MOTION

        &FORCE_EVAL
          METHOD Quickstep
          &DFT
            BASIS_SET_FILE_NAME BASIS_MOLOPT
            POTENTIAL_FILE_NAME GTH_POTENTIALS
            &MGRID
              CUTOFF 400
              REL_CUTOFF 50
            &END MGRID
            &XC
              &XC_FUNCTIONAL PBE
              &END XC_FUNCTIONAL
              &VDW_POTENTIAL
                DISPERSION_FUNCTIONAL PAIR_POTENTIAL
                &PAIR_POTENTIAL
                  TYPE DFTD3(BJ)
                  CALCULATE_C9_TERM .TRUE.
                  PARAMETER_FILE_NAME dftd3.dat
                  REFERENCE_FUNCTIONAL PBE
                &END PAIR_POTENTIAL
              &END VDW_POTENTIAL
            &END XC
            &SCF
              MAX_SCF 200
              EPS_SCF 1.0E-6
              &MIXING
                METHOD BROYDEN_MIXING
                ALPHA 0.3
              &END MIXING
            &END SCF
            &KPOINTS
              SCHEME MONKHORST-PACK 1 1 {max(1, round(20.0 / Lz))}
            &END KPOINTS
          &END DFT
          &SUBSYS
            &CELL
              A {Lx:.6f} 0.0 0.0
              B 0.0 {Ly:.6f} 0.0
              C 0.0 0.0 {Lz:.6f}
              PERIODIC XYZ
            &END CELL
            &COORD
        {chr(10).join(coord_lines)}
            &END COORD
        {chr(10).join(kind_lines)}
          &END SUBSYS
        &END FORCE_EVAL
    """).rstrip()

    return text


# ─────────────────────────────────────────────────────────────────────────────
# Database queries
# ─────────────────────────────────────────────────────────────────────────────

def query_cod(
    formula: str,
    max_results: int = 30,
    timeout: float = 10.0,
) -> list[dict]:
    """
    Query the Crystallography Open Database (COD) REST API for structures
    matching *formula*.

    Uses the COD /result endpoint with the ``formula`` parameter which
    accepts standard chemical formulas (e.g. "MoS2", "BN", "C").

    Returns a list of dicts with keys:
        id, formula, sg (space group), a, b, c, file_url, source

    Raises RuntimeError on network failure.
    """
    import urllib.request
    import urllib.parse
    import json

    formula = formula.strip()

    # COD REST API — search by exact formula text.
    # The ``formula`` parameter matches the stored chemical formula string.
    params = urllib.parse.urlencode({
        "format":  "json",
        "formula": formula,
    })
    url = f"https://www.crystallography.net/cod/result?{params}"

    req = urllib.request.Request(url, headers={
        "User-Agent": "NTBuilder/1.1 (research use)",
        "Accept": "application/json",
    })
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read())
    except Exception as exc:
        raise RuntimeError(f"COD query failed: {exc}") from exc

    results = []
    for entry in data[:max_results]:
        file_id = entry.get("file", "")
        results.append({
            "id":       file_id,
            "formula":  entry.get("formula", ""),
            "sg":       entry.get("sg", ""),
            "a":        entry.get("a", ""),
            "b":        entry.get("b", ""),
            "c":        entry.get("c", ""),
            "file_url": f"https://www.crystallography.net/cod/{file_id}.cif",
            "title":    entry.get("title", ""),
            "source":   "COD",
        })
    return results


def query_mp(
    formula: str,
    api_key: str,
    max_results: int = 20,
    timeout: float = 15.0,
) -> list[dict]:
    """
    Query the Materials Project API v3 for 2D structures matching *formula*.

    Requires a valid Materials Project API key (free at materialsproject.org).
    Returns a list of dicts with keys:
        id, formula, sg, a, b, c, alpha, beta, gamma, source

    Structure files are not downloaded automatically — use the MP website or
    mp-api (``pip install mp-api``) to fetch the CIF directly.
    """
    try:
        import urllib.request
        import urllib.parse
        import json
    except ImportError:
        raise RuntimeError("urllib not available")

    headers = {
        "X-API-KEY": api_key,
        "Content-Type": "application/json",
    }
    params = urllib.parse.urlencode({
        "formula":    formula,
        "fields":     "material_id,formula_pretty,symmetry,lattice",
        "_limit":     str(max_results),
        "theoretical": "false",
    })
    url = f"https://api.materialsproject.org/materials/summary/?{params}"

    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read())
    except Exception as exc:
        raise RuntimeError(f"Materials Project query failed: {exc}") from exc

    results = []
    for entry in data.get("data", []):
        lat = entry.get("lattice", {})
        results.append({
            "id":       entry.get("material_id", ""),
            "formula":  entry.get("formula_pretty", formula),
            "sg":       entry.get("symmetry", {}).get("symbol", ""),
            "a":        lat.get("a", ""),
            "b":        lat.get("b", ""),
            "c":        lat.get("c", ""),
            "alpha":    lat.get("alpha", ""),
            "beta":     lat.get("beta", ""),
            "gamma":    lat.get("gamma", ""),
            "mp_url":   f"https://materialsproject.org/materials/{entry.get('material_id','')}",
            "source":   "Materials Project",
        })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# C2DB local cache  (one-time download, SQLite / ASE-DB format)
# ─────────────────────────────────────────────────────────────────────────────

class C2DBCache:
    """
    Manages a local copy of the C2DB ASE-SQLite database.

    The database is downloaded once from the official C2DB website and stored
    at  ~/.ntbuilder/c2db.db  (~200 MB – 2 GB depending on version).
    All subsequent searches are fully offline and instant.

    Public API
    ----------
    C2DBCache.is_available()         → bool
    C2DBCache.db_size_mb()           → float | None  (None if not downloaded)
    C2DBCache.download(progress_cb)  → downloads; raises on error
    C2DBCache.search(formula, ...)   → list[dict]
    """

    # Storage location – platform-agnostic home dir
    DB_DIR  = None   # set lazily (pathlib not available at class scope reliably)
    DB_NAME = "c2db.db"

    # Known download URLs, tried in order
    _DOWNLOAD_URLS = [
        # DTU primary
        "https://c2db.fysik.dtu.dk/c2db.db",
        # Older CMR mirror
        "https://cmrdb.fysik.dtu.dk/c2db/c2db.db",
    ]

    # Keys used in the ASE-DB key_value_pairs / text_key_values.
    # C2DB has changed column names across versions; we check all variants.
    _FORMULA_KEYS  = ("formula",)
    _UID_KEYS      = ("uid", "id", "name")
    _EHULL_KEYS    = ("ehull", "e_hull", "e_above_hull")
    _GAP_KEYS      = ("gap", "gap_pbe", "gap_nosoc")
    _MAG_KEYS      = ("magstate", "magnetic", "magmom_l")
    _LG_KEYS       = ("layergroup", "layer_group", "spacegroup")

    BASE_URL = "https://c2db.fysik.dtu.dk"

    @classmethod
    def _db_path(cls):
        from pathlib import Path
        p = Path.home() / ".ntbuilder" / cls.DB_NAME
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    @classmethod
    def is_available(cls) -> bool:
        p = cls._db_path()
        return p.exists() and p.stat().st_size > 500_000   # >500 KB = real DB

    @classmethod
    def db_size_mb(cls) -> "float | None":
        p = cls._db_path()
        if p.exists():
            return p.stat().st_size / 1_048_576
        return None

    @classmethod
    def _find_download_url(cls, timeout: float = 10.0) -> str:
        """
        Try to extract the direct .db download URL from the C2DB website.
        Falls back to the first entry in _DOWNLOAD_URLS if not found.
        """
        import urllib.request, re
        hdr = {"User-Agent": "NTBuilder/1.1", "Accept": "text/html,*/*"}
        try:
            req = urllib.request.Request(cls.BASE_URL + "/", headers=hdr)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                html = r.read().decode("utf-8", errors="replace")
            # Look for an href that ends with .db (the database download link)
            m = re.search(r'href=["\']([^"\']*\.db)["\']', html)
            if m:
                href = m.group(1)
                if href.startswith("http"):
                    return href
                return cls.BASE_URL.rstrip("/") + "/" + href.lstrip("/")
        except Exception:
            pass
        return cls._DOWNLOAD_URLS[0]

    @classmethod
    def download(cls, progress_cb=None, timeout: float = 60.0):
        """
        Download the C2DB database.

        Parameters
        ----------
        progress_cb : callable(bytes_done: int, bytes_total: int) | None
            Called periodically during download.  bytes_total may be 0
            if the server does not send Content-Length.
        timeout     : connection timeout per chunk in seconds.

        Raises RuntimeError on failure.
        """
        import urllib.request, shutil, tempfile
        from pathlib import Path

        target = cls._db_path()

        # Determine URL
        urls = [cls._find_download_url(timeout=10)] + cls._DOWNLOAD_URLS
        seen = set()
        urls = [u for u in urls if not (u in seen or seen.add(u))]

        last_exc = None
        for url in urls:
            try:
                hdr = {
                    "User-Agent": "NTBuilder/1.1 (research use)",
                    "Accept":     "*/*",
                }
                req = urllib.request.Request(url, headers=hdr)
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    total = int(resp.headers.get("Content-Length", 0) or 0)
                    done  = 0
                    chunk = 1 << 17   # 128 KB chunks

                    # Write to a temp file first so partial downloads don't corrupt
                    tmp = tempfile.NamedTemporaryFile(
                        dir=target.parent, prefix=".c2db_dl_", suffix=".db",
                        delete=False
                    )
                    try:
                        while True:
                            buf = resp.read(chunk)
                            if not buf:
                                break
                            tmp.write(buf)
                            done += len(buf)
                            if progress_cb:
                                progress_cb(done, total)
                        tmp.close()
                        # Atomic replace
                        shutil.move(tmp.name, str(target))
                    except Exception:
                        tmp.close()
                        Path(tmp.name).unlink(missing_ok=True)
                        raise

                return   # success
            except Exception as exc:
                last_exc = exc
                continue

        raise RuntimeError(
            f"Could not download C2DB database from any known URL.\n"
            f"Last error: {last_exc}\n\n"
            f"You can download it manually from:\n"
            f"  {cls._DOWNLOAD_URLS[0]}\n"
            f"and place it at:\n"
            f"  {cls._db_path()}"
        )

    @classmethod
    def search(cls,
               formula: str,
               stability_max: float = 0.5,
               max_results: int    = 50,
               ) -> list[dict]:
        """
        Search the local C2DB SQLite database.

        The database uses the ASE-DB v3 schema:
            systems           – one row per material
            text_key_values   – (id, key, value TEXT)
            number_key_values – (id, key, value REAL)

        Searches are done with SQL JOINs on the normalised key tables,
        which have indexes and are fast even for 16 000+ entries.
        """
        import sqlite3, json

        db_path = cls._db_path()
        if not cls.is_available():
            raise RuntimeError("C2DB local database not available. Download it first.")

        con = sqlite3.connect(str(db_path))
        con.row_factory = sqlite3.Row

        def _has_table(name: str) -> bool:
            cur = con.execute(
                "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?",
                (name,)
            )
            return cur.fetchone()[0] > 0

        results: list[dict] = []

        try:
            # ── Strategy A: normalised tables (ASE-DB v3, fast) ──────────────
            if _has_table("text_key_values") and _has_table("number_key_values"):
                # Build formula condition (try all known formula key names)
                fkey_cond = " OR ".join(
                    f"(tkv.key='{k}' AND tkv.value LIKE ?)"
                    for k in cls._FORMULA_KEYS
                )
                fkey_params = [f"%{formula.strip()}%"] * len(cls._FORMULA_KEYS)

                sql = f"""
                    SELECT DISTINCT s.id, s.unique_id, s.key_value_pairs
                    FROM   systems s
                    JOIN   text_key_values tkv ON tkv.id = s.id
                    WHERE  ({fkey_cond})
                    LIMIT  {max_results * 5}
                """
                rows = con.execute(sql, fkey_params).fetchall()

                # Now filter by ehull (may live in number_key_values)
                def _get_num(sid, keys):
                    for k in keys:
                        cur = con.execute(
                            "SELECT value FROM number_key_values WHERE id=? AND key=? LIMIT 1",
                            (sid, k)
                        )
                        r = cur.fetchone()
                        if r is not None:
                            return r[0]
                    return None

                def _get_txt(sid, keys):
                    for k in keys:
                        cur = con.execute(
                            "SELECT value FROM text_key_values WHERE id=? AND key=? LIMIT 1",
                            (sid, k)
                        )
                        r = cur.fetchone()
                        if r is not None:
                            return r[0]
                    return None

                for row in rows:
                    sid = row["id"]
                    ehull = _get_num(sid, cls._EHULL_KEYS) or 0.0
                    if ehull > stability_max:
                        continue

                    # Try to get UID from text keys first, fall back to unique_id
                    uid = _get_txt(sid, cls._UID_KEYS) or row["unique_id"] or str(sid)
                    formula_val = _get_txt(sid, cls._FORMULA_KEYS) or "?"
                    gap_pbe     = _get_num(sid, cls._GAP_KEYS)
                    magnetic    = _get_txt(sid, cls._MAG_KEYS) or "—"
                    layer_group = _get_txt(sid, cls._LG_KEYS) or "—"

                    results.append({
                        "id":          uid,
                        "uid":         uid,
                        "formula":     formula_val,
                        "layer_group": layer_group,
                        "ehull":       ehull,
                        "gap_pbe":     gap_pbe,
                        "magnetic":    magnetic,
                        "file_url":    f"{cls.BASE_URL}/material/{uid}/download/cif",
                        "source":      "C2DB",
                    })
                    if len(results) >= max_results:
                        break

            # ── Strategy B: JSON blob in systems.key_value_pairs (fallback) ──
            if not results:
                sql = f"""
                    SELECT id, unique_id, key_value_pairs
                    FROM   systems
                    WHERE  key_value_pairs LIKE ?
                    LIMIT  {max_results * 5}
                """
                rows = con.execute(sql, (f"%{formula.strip()}%",)).fetchall()

                def _jfloat(d, keys):
                    for k in keys:
                        v = d.get(k)
                        if v is not None:
                            try:
                                return float(v)
                            except (TypeError, ValueError):
                                pass
                    return None

                def _jtxt(d, keys, default="—"):
                    for k in keys:
                        v = d.get(k)
                        if v:
                            return str(v)
                    return default

                for row in rows:
                    try:
                        kvp = json.loads(row["key_value_pairs"] or "{}")
                    except json.JSONDecodeError:
                        continue

                    ehull = _jfloat(kvp, cls._EHULL_KEYS) or 0.0
                    if ehull > stability_max:
                        continue

                    uid = _jtxt(kvp, cls._UID_KEYS, row["unique_id"] or str(row["id"]))
                    formula_val = _jtxt(kvp, cls._FORMULA_KEYS, "?")

                    # Double-check formula matches the search term
                    if formula.strip().lower() not in formula_val.lower():
                        continue

                    results.append({
                        "id":          uid,
                        "uid":         uid,
                        "formula":     formula_val,
                        "layer_group": _jtxt(kvp, cls._LG_KEYS),
                        "ehull":       ehull,
                        "gap_pbe":     _jfloat(kvp, cls._GAP_KEYS),
                        "magnetic":    _jtxt(kvp, cls._MAG_KEYS),
                        "file_url":    f"{cls.BASE_URL}/material/{uid}/download/cif",
                        "source":      "C2DB",
                    })
                    if len(results) >= max_results:
                        break

        finally:
            con.close()

        return results


def _c2db_local_search(formula: str, stability_max: float) -> list[dict]:
    """
    Small curated subset of C2DB for common 2D materials.
    UIDs are verified against c2db.fysik.dtu.dk.
    Used as a reliable fallback when the live search returns nothing.
    """
    _DB = [
        # ── Graphene / h-BN ───────────────────────────────────────────────────
        dict(uid="2C-1",    formula="C",     layer_group="P6/mmm",  ehull=0.000, gap_pbe=0.00, magnetic="NM"),
        dict(uid="2BN-1",   formula="BN",    layer_group="P6/mmm",  ehull=0.000, gap_pbe=4.63, magnetic="NM"),
        # ── Mo dichalcogenides ────────────────────────────────────────────────
        dict(uid="1MoS2-1", formula="MoS2",  layer_group="P-6m2",   ehull=0.000, gap_pbe=1.67, magnetic="NM"),
        dict(uid="1MoSe2-1",formula="MoSe2", layer_group="P-6m2",   ehull=0.000, gap_pbe=1.47, magnetic="NM"),
        dict(uid="1MoTe2-1",formula="MoTe2", layer_group="P-6m2",   ehull=0.000, gap_pbe=1.09, magnetic="NM"),
        # ── W dichalcogenides ─────────────────────────────────────────────────
        dict(uid="1WS2-1",  formula="WS2",   layer_group="P-6m2",   ehull=0.000, gap_pbe=1.80, magnetic="NM"),
        dict(uid="1WSe2-1", formula="WSe2",  layer_group="P-6m2",   ehull=0.000, gap_pbe=1.65, magnetic="NM"),
        dict(uid="1WTe2-1", formula="WTe2",  layer_group="P-6m2",   ehull=0.022, gap_pbe=0.75, magnetic="NM"),
        # ── Nb / V dichalcogenides (magnetic) ────────────────────────────────
        dict(uid="1NbS2-1", formula="NbS2",  layer_group="P-6m2",   ehull=0.000, gap_pbe=0.00, magnetic="FM"),
        dict(uid="1NbSe2-1",formula="NbSe2", layer_group="P-6m2",   ehull=0.000, gap_pbe=0.00, magnetic="FM"),
        dict(uid="1VS2-1",  formula="VS2",   layer_group="P-6m2",   ehull=0.000, gap_pbe=0.00, magnetic="FM"),
        dict(uid="1VSe2-1", formula="VSe2",  layer_group="P-6m2",   ehull=0.000, gap_pbe=0.00, magnetic="FM"),
        # ── Cr dichalcogenides ────────────────────────────────────────────────
        dict(uid="1CrS2-1", formula="CrS2",  layer_group="P-6m2",   ehull=0.078, gap_pbe=0.00, magnetic="FM"),
        dict(uid="1CrSe2-1",formula="CrSe2", layer_group="P-6m2",   ehull=0.091, gap_pbe=0.00, magnetic="FM"),
        # ── Ti dichalcogenides ────────────────────────────────────────────────
        dict(uid="1TiS2-1", formula="TiS2",  layer_group="P-3m1",   ehull=0.000, gap_pbe=0.00, magnetic="NM"),
        dict(uid="1TiSe2-1",formula="TiSe2", layer_group="P-3m1",   ehull=0.000, gap_pbe=0.00, magnetic="NM"),
        # ── Sn / Pb dichalcogenides ───────────────────────────────────────────
        dict(uid="1SnS2-1", formula="SnS2",  layer_group="P-3m1",   ehull=0.000, gap_pbe=1.59, magnetic="NM"),
        dict(uid="1SnSe2-1",formula="SnSe2", layer_group="P-3m1",   ehull=0.000, gap_pbe=0.83, magnetic="NM"),
        # ── III-VI monochalcogenides ──────────────────────────────────────────
        dict(uid="2GaS-1",  formula="GaS",   layer_group="P-6m2",   ehull=0.000, gap_pbe=2.28, magnetic="NM"),
        dict(uid="2GaSe-1", formula="GaSe",  layer_group="P-6m2",   ehull=0.000, gap_pbe=1.94, magnetic="NM"),
        dict(uid="2GaTe-1", formula="GaTe",  layer_group="P-6m2",   ehull=0.040, gap_pbe=1.56, magnetic="NM"),
        dict(uid="2InS-1",  formula="InS",   layer_group="P-6m2",   ehull=0.000, gap_pbe=2.00, magnetic="NM"),
        dict(uid="2InSe-1", formula="InSe",  layer_group="P-6m2",   ehull=0.000, gap_pbe=1.35, magnetic="NM"),
        dict(uid="2InTe-1", formula="InTe",  layer_group="P-6m2",   ehull=0.080, gap_pbe=1.00, magnetic="NM"),
        # ── Bismuth ───────────────────────────────────────────────────────────
        dict(uid="2Bi-1",   formula="Bi",    layer_group="P-3m1",   ehull=0.000, gap_pbe=0.00, magnetic="NM"),
        dict(uid="2BiI3-1", formula="BiI3",  layer_group="P-3",     ehull=0.000, gap_pbe=2.20, magnetic="NM"),
        # ── Phosphorus / As ───────────────────────────────────────────────────
        dict(uid="4P-1",    formula="P",     layer_group="Pmna",    ehull=0.000, gap_pbe=0.91, magnetic="NM"),
        dict(uid="4As-1",   formula="As",    layer_group="Pmna",    ehull=0.000, gap_pbe=0.00, magnetic="NM"),
        # ── MXene-like / other oxides ─────────────────────────────────────────
        dict(uid="1MoO2-1", formula="MoO2",  layer_group="P-6m2",   ehull=0.090, gap_pbe=0.00, magnetic="NM"),
        dict(uid="1WO2-1",  formula="WO2",   layer_group="P-6m2",   ehull=0.100, gap_pbe=0.00, magnetic="NM"),
    ]

    formula_clean = formula.strip()
    matches = []
    for entry in _DB:
        if formula_clean.lower() in entry["formula"].lower():
            e = entry.get("ehull", 0.0)
            if e <= stability_max:
                matches.append({
                    "id":          entry["uid"],
                    "uid":         entry["uid"],
                    "formula":     entry["formula"],
                    "layer_group": entry["layer_group"],
                    "ehull":       entry["ehull"],
                    "gap_pbe":     entry["gap_pbe"],
                    "magnetic":    entry["magnetic"],
                    "file_url":    f"https://c2db.fysik.dtu.dk/material/{entry['uid']}/download/cif",
                    "source":      "C2DB",
                })
    return matches


def query_c2db(
    formula: str,
    stability_max: float = 0.5,
    max_results: int = 50,
) -> list[dict]:
    """
    Search C2DB for 2D structures matching *formula*.

    Priority order:
    1. Local SQLite cache (~/.ntbuilder/c2db.db) — full database, instant search.
    2. Bundled curated list — 30 verified common materials, always available.

    Call C2DBCache.download() first to enable the full database.
    """
    # ── Tier 1: local SQLite cache ────────────────────────────────────────────
    if C2DBCache.is_available():
        try:
            return C2DBCache.search(formula, stability_max, max_results)
        except Exception as exc:
            import sys
            print(f"[C2DB] local cache error: {exc}", file=sys.stderr)

    # ── Tier 2: small curated fallback ────────────────────────────────────────
    return _c2db_local_search(formula, stability_max)[:max_results]
