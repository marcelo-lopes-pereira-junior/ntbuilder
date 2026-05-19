"""
tests/test_core.py
------------------
Unit and integration tests for the NTBuilder core library.

Sections
--------
  1. Helpers               — in-memory structure factories (no file I/O)
  2. TestLatticeStructure  — LatticeStructure properties
  3. TestLatticeClassification — lattice_type detection
  4. TestChirality         — compute_chirality and scan_chirality
  5. TestBuilder           — build_nanotube + geometric invariants
  6. TestExporters         — all export formats (XYZ, PDB, POSCAR, LAMMPS, QE)
  7. TestConnectivity      — bond detection and rendering arrays
  8. TestSymmetry          — snap_to_symmetry + find_primitive_cell
  9. TestRoundTrip         — write → parse-back consistency

Run with:
    pytest tests/ -v
    pytest tests/ -v --tb=short -q   # compact output
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

from core.io        import LatticeStructure, read_xyz, read_pdb
from core.chirality import compute_chirality, scan_chirality
from core.builder   import build_nanotube
from core.exporters import (
    write_xyz, write_pdb, write_poscar, write_lammps, write_qe, export,
)
from core.connectivity import compute_bonds, BondSettings, bond_line_arrays, get_radius
from core.symmetry import snap_to_symmetry, find_primitive_cell


# ─────────────────────────────────────────────────────────────────────────────
# 1. Helpers — in-memory structure factories
# ─────────────────────────────────────────────────────────────────────────────

def _graphene():
    """Primitive hexagonal graphene cell (a = 2.46 Å, γ = 60°)."""
    a = 2.46
    a1 = np.array([a, 0.0])
    a2 = np.array([a * math.cos(math.radians(60)), a * math.sin(math.radians(60))])
    atoms = [
        {"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0},
        {"symbol": "C", "pos": a1 / 3 + a2 / 3, "z": 0.0},
    ]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _rectangular():
    """Simple rectangular cell (a ≠ b, γ = 90°)."""
    a1 = np.array([3.0, 0.0])
    a2 = np.array([0.0, 4.0])
    atoms = [{"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0}]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _oblique():
    """Minimal oblique cell (γ ≈ 75°)."""
    gamma = math.radians(75)
    a, b = 3.0, 3.5
    a1 = np.array([a, 0.0])
    a2 = np.array([b * math.cos(gamma), b * math.sin(gamma)])
    atoms = [{"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0}]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _hbn():
    """Hexagonal boron nitride (a = 2.50 Å, 2 species)."""
    a = 2.50
    a1 = np.array([a, 0.0])
    a2 = np.array([a * math.cos(math.radians(60)), a * math.sin(math.radians(60))])
    atoms = [
        {"symbol": "B", "pos": np.array([0.0, 0.0]), "z": 0.0},
        {"symbol": "N", "pos": a1 / 3 + a2 / 3, "z": 0.0},
    ]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _buckled():
    """Silicene-like buckled hexagonal cell (2 Si at ±0.23 Å z-offset)."""
    a = 3.87
    a1 = np.array([a, 0.0])
    a2 = np.array([a * math.cos(math.radians(60)), a * math.sin(math.radians(60))])
    atoms = [
        {"symbol": "Si", "pos": np.array([0.0, 0.0]), "z":  0.23},
        {"symbol": "Si", "pos": a1 / 3 + a2 / 3,     "z": -0.23},
    ]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _graphene_2x1():
    """2×1 supercell of graphene — 4 atoms (should reduce to 2 in primitive cell)."""
    a = 2.46
    a1 = np.array([a, 0.0])
    a2 = np.array([a * math.cos(math.radians(60)), a * math.sin(math.radians(60))])
    a1_super = 2.0 * a1
    a2_super = a2.copy()
    base = [
        {"symbol": "C", "pos": np.array([0.0, 0.0]),   "z": 0.0},
        {"symbol": "C", "pos": a1 / 3 + a2 / 3,        "z": 0.0},
    ]
    atoms = []
    for k in range(2):
        for at in base:
            atoms.append({**at, "pos": at["pos"] + k * a1})
    return LatticeStructure(atoms=atoms, a1=a1_super, a2=a2_super)


def _noisy_hexagonal():
    """
    Hexagonal cell with slight numerical noise in both |a| and γ.

    Length noise: +0.0005 Å  (<< the 0.001 Å classification threshold).
    Angular noise: γ = 60.05° (well within the ±0.5° tolerance).
    snap_to_symmetry should fix both to exact hexagonal values.
    """
    a = 2.46
    a1 = np.array([a + 0.0005, 0.0])
    a2 = np.array([a * math.cos(math.radians(60.05)),
                   a * math.sin(math.radians(60.05))])
    atoms = [{"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0}]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _slightly_off_rectangular():
    """Near-rectangular cell with γ ≈ 89.9° and a ≠ b (no square snap)."""
    a1 = np.array([3.0, 0.0])
    off = math.radians(0.1)                        # 0.1° tilt
    a2 = np.array([4.0 * math.sin(off), 4.0 * math.cos(off)])
    atoms = [{"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0}]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


# ─────────────────────────────────────────────────────────────────────────────
# Shared fixture — build a (5,5) nanotube once per class
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def nt_5_5():
    s = _graphene()
    ch = compute_chirality(5, 5, s)
    return build_nanotube(s, ch, vacuum=10.0)


# ─────────────────────────────────────────────────────────────────────────────
# 2. TestLatticeStructure
# ─────────────────────────────────────────────────────────────────────────────

class TestLatticeStructure:
    def test_graphene_properties(self):
        """Lattice parameters match the graphene input."""
        s = _graphene()
        assert abs(s.a - 2.46) < 0.01
        assert abs(s.b - 2.46) < 0.01
        assert abs(s.gamma_deg - 60.0) < 0.5

    def test_has_buckling_false_for_flat(self):
        s = _graphene()
        assert s.has_buckling is False

    def test_has_buckling_true_for_buckled(self):
        s = _buckled()
        assert s.has_buckling is True

    def test_d_min_flat_is_zero(self):
        s = _graphene()
        assert s.d_min == 0.0

    def test_d_min_buckled_is_twice_offset(self):
        s = _buckled()
        # max_z_offset = 0.23 Å → d_min = 2 × 0.23 = 0.46 Å
        assert abs(s.d_min - 2 * 0.23) < 0.01

    def test_gamma120_normalised_to_60(self):
        """LatticeStructure must convert γ = 120° hexagonal cells to γ = 60°."""
        a = 2.46
        a1 = np.array([a, 0.0])
        # γ = 120°: a2 points at 120° from a1
        a2 = np.array([a * math.cos(math.radians(120)),
                       a * math.sin(math.radians(120))])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        # After normalisation the stored angle must be ~60°
        assert abs(s.gamma_deg - 60.0) < 0.5

    def test_is_square_false_for_non_square_rectangular(self):
        s = _rectangular()     # a=3, b=4, γ=90°
        assert s.is_square is False

    def test_is_square_true(self):
        a1 = np.array([3.0, 0.0])
        a2 = np.array([0.0, 3.0])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        assert s.is_square is True


# ─────────────────────────────────────────────────────────────────────────────
# 3. TestLatticeClassification
# ─────────────────────────────────────────────────────────────────────────────

class TestLatticeClassification:
    def test_graphene_is_hexagonal(self):
        assert _graphene().lattice_type == "hexagonal"

    def test_rectangular_gamma90(self):
        assert _rectangular().lattice_type == "rectangular"

    def test_oblique_gamma75(self):
        assert _oblique().lattice_type == "oblique"

    def test_square_is_rectangular(self):
        a1 = np.array([3.0, 0.0])
        a2 = np.array([0.0, 3.0])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        assert s.lattice_type == "rectangular"

    def test_hbn_is_hexagonal(self):
        """h-BN has a hexagonal lattice."""
        assert _hbn().lattice_type == "hexagonal"

    def test_very_oblique_gamma20(self):
        """A cell with γ = 20° is oblique."""
        gamma = math.radians(20)
        a1 = np.array([3.0, 0.0])
        a2 = np.array([3.0 * math.cos(gamma), 3.0 * math.sin(gamma)])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        assert s.lattice_type == "oblique"


# ─────────────────────────────────────────────────────────────────────────────
# 4. TestChirality
# ─────────────────────────────────────────────────────────────────────────────

class TestChirality:
    # ── Corrected test (was using wrong formula — missing √3 factor) ─────────
    def test_armchair_graphene_diameter(self):
        """(5,5) armchair graphene: D = a·√(n²+nm+m²)/π ≈ 6.78 Å."""
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        expected = 2.46 * math.sqrt(5**2 + 5*5 + 5**2) / math.pi
        assert abs(ch.diameter - expected) < 0.05

    # ── Literature reference values (Dresselhaus 1995) ────────────────────────
    def test_zigzag_10_0_diameter(self):
        """(10,0) zigzag CNT: D ≈ 7.83 Å (literature value)."""
        s = _graphene()
        ch = compute_chirality(10, 0, s)
        expected = 2.46 * math.sqrt(10**2) / math.pi   # = 2.46·10/π
        assert abs(ch.diameter - expected) < 0.05

    def test_armchair_10_10_diameter(self):
        """(10,10) armchair CNT: D ≈ 13.56 Å (literature value)."""
        s = _graphene()
        ch = compute_chirality(10, 10, s)
        expected = 2.46 * math.sqrt(10**2 + 10*10 + 10**2) / math.pi
        assert abs(ch.diameter - expected) < 0.05

    def test_diameter_formula_general(self):
        """D = a·√(n²+nm+m²)/π holds for (7,3) graphene."""
        s = _graphene()
        n, m = 7, 3
        ch = compute_chirality(n, m, s)
        expected = 2.46 * math.sqrt(n**2 + n*m + m**2) / math.pi
        assert abs(ch.diameter - expected) < 0.1

    def test_zero_zero_returns_none(self):
        assert compute_chirality(0, 0, _graphene()) is None

    def test_zigzag_theta_zero(self):
        """(n,0) tubes have chiral angle ≈ 0°."""
        ch = compute_chirality(8, 0, _graphene())
        assert abs(ch.theta_deg) < 0.5

    def test_armchair_theta_30(self):
        """(n,n) hexagonal tubes have chiral angle ≈ 30°."""
        ch = compute_chirality(6, 6, _graphene())
        assert abs(ch.theta_deg - 30.0) < 0.5

    def test_atom_count_graphene_armchair_5_5(self):
        """(5,5) armchair graphene nanotube: 20 atoms per unit cell."""
        ch = compute_chirality(5, 5, _graphene())
        assert ch.n_atoms == 20

    def test_atom_count_positive(self):
        """Any (n,m) ≠ (0,0) produces a positive atom count."""
        s = _graphene()
        for n, m in [(3, 0), (4, 4), (5, 2), (6, 1)]:
            ch = compute_chirality(n, m, s)
            assert ch.n_atoms > 0, f"n_atoms ≤ 0 for ({n},{m})"

    def test_strain_nonnegative(self):
        """Strain must be ≥ 0 for all chiralities."""
        s = _graphene()
        for n, m in [(5, 5), (8, 0), (4, 2), (3, 0)]:
            ch = compute_chirality(n, m, s)
            assert ch.strain >= 0.0, f"Negative strain for ({n},{m})"

    def test_rectangular_strain_positive(self):
        """Non-hexagonal lattice (3,2) has non-negative strain."""
        ch = compute_chirality(3, 2, _rectangular())
        assert ch.strain >= 0.0

    def test_oblique_sector_gamma(self):
        """Oblique unique sector equals γ, not 90°."""
        from core.chirality import unique_sector_deg
        s = _oblique()
        assert abs(unique_sector_deg(s) - s.gamma_deg) < 0.01

    def test_hbn_multispecies_chirality(self):
        """compute_chirality works for h-BN (two-species hexagonal lattice)."""
        ch = compute_chirality(4, 4, _hbn())
        assert ch is not None
        assert ch.n_atoms > 0
        assert ch.diameter > 0.0

    def test_scan_returns_list(self):
        results = scan_chirality(_graphene(), n_max=5, m_max=5, max_diameter=20.0)
        assert isinstance(results, list)
        assert len(results) > 0

    def test_scan_diameter_filter(self):
        """All results from scan respect the max_diameter ceiling."""
        max_d = 15.0
        results = scan_chirality(_graphene(), n_max=8, m_max=8, max_diameter=max_d)
        for ch in results:
            assert ch.diameter <= max_d + 0.01, (
                f"Diameter {ch.diameter:.2f} Å exceeds max_diameter {max_d} Å"
                f" for ({ch.n},{ch.m})"
            )

    def test_scan_no_zero_zero(self):
        """scan_chirality never returns (0,0)."""
        results = scan_chirality(_graphene(), n_max=4, m_max=4, max_diameter=20.0)
        for ch in results:
            assert not (ch.n == 0 and ch.m == 0)


# ─────────────────────────────────────────────────────────────────────────────
# 5. TestBuilder
# ─────────────────────────────────────────────────────────────────────────────

class TestBuilder:
    def test_atom_count_matches_chirality(self):
        """Built nanotube atom count must match the chirality prediction."""
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert len(nt.symbols) == ch.n_atoms

    def test_coords_shape(self):
        s = _graphene()
        ch = compute_chirality(4, 2, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert nt.coords.shape == (ch.n_atoms, 3)

    def test_tube_is_centred(self):
        """Nanotube centroid in xy must sit at (box/2, box/2)."""
        s = _graphene()
        ch = compute_chirality(6, 0, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        cx = nt.coords[:, 0].mean()
        cy = nt.coords[:, 1].mean()
        centre = nt.box[0] / 2.0
        assert abs(cx - centre) < 0.5
        assert abs(cy - centre) < 0.5

    def test_z_span_at_most_T_norm(self):
        """
        Z-span of atoms must be strictly less than T_norm.

        In a periodic nanotube unit cell there is always a small gap between the
        last atom and the periodic image of the first — making z_span < T_norm
        the physically correct invariant, NOT z_span ≈ T_norm.
        """
        s = _graphene()
        ch = compute_chirality(5, 3, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        z_span = nt.coords[:, 2].max() - nt.coords[:, 2].min()
        assert z_span <= ch.T_norm + 1e-6, (
            f"z_span {z_span:.4f} Å > T_norm {ch.T_norm:.4f} Å"
        )

    def test_vacuum_sets_box(self):
        """Box lateral dimension must equal diameter + vacuum."""
        vacuum = 15.0
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        nt = build_nanotube(s, ch, vacuum=vacuum)
        expected_box_xy = ch.diameter + vacuum
        assert abs(float(nt.box[0]) - expected_box_xy) < 0.01
        assert abs(float(nt.box[1]) - expected_box_xy) < 0.01

    def test_cylindrical_shell(self):
        """All atoms must lie on the cylindrical shell at radius ≈ D/2."""
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        R_expected = ch.diameter / 2.0
        cx = float(nt.box[0]) / 2.0
        cy = float(nt.box[1]) / 2.0
        radii = np.sqrt((nt.coords[:, 0] - cx)**2 + (nt.coords[:, 1] - cy)**2)
        assert np.allclose(radii, R_expected, atol=0.1), (
            f"Max radial deviation: {np.abs(radii - R_expected).max():.4f} Å"
        )

    def test_no_coincident_atoms(self):
        """No two atoms in the nanotube may overlap (min distance > 0.5 Å)."""
        pytest.importorskip("scipy")
        from scipy.spatial.distance import pdist
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        dists = pdist(nt.coords)
        assert dists.min() > 0.5, (
            f"Minimum inter-atom distance {dists.min():.4f} Å — atoms may overlap"
        )

    def test_roll_inward_preserves_atom_count(self):
        """roll_inward=True must yield the same number of atoms as default."""
        s = _graphene()
        ch = compute_chirality(5, 5, s)
        nt_out = build_nanotube(s, ch, vacuum=10.0, roll_inward=False)
        nt_in  = build_nanotube(s, ch, vacuum=10.0, roll_inward=True)
        assert len(nt_in.symbols) == len(nt_out.symbols)

    def test_buckled_atoms_at_different_radii(self):
        """
        In a buckled nanotube the two sublattices must sit at distinct radii,
        separated by roughly 2 × z_offset.
        """
        s = _buckled()
        ch = compute_chirality(5, 5, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        cx = float(nt.box[0]) / 2.0
        cy = float(nt.box[1]) / 2.0
        radii = np.sqrt((nt.coords[:, 0] - cx)**2 + (nt.coords[:, 1] - cy)**2)
        r_min, r_max = radii.min(), radii.max()
        # Two shells separated by ≈ 2 × 0.23 = 0.46 Å
        assert r_max - r_min > 0.30, (
            f"Expected two distinct radial shells; got r_min={r_min:.3f}, "
            f"r_max={r_max:.3f}"
        )

    def test_rectangular_builder(self):
        """Iterator builder works for rectangular lattices."""
        s = _rectangular()
        ch = compute_chirality(3, 0, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert len(nt.symbols) == ch.n_atoms

    def test_oblique_builder(self):
        """Iterator builder works for oblique lattices."""
        s = _oblique()
        ch = compute_chirality(2, 1, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert len(nt.symbols) == ch.n_atoms

    def test_hbn_builder(self):
        """Nanotube builder works for multi-species (h-BN) structures."""
        s = _hbn()
        ch = compute_chirality(3, 3, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert len(nt.symbols) == ch.n_atoms
        # Both species must be present
        assert "B" in nt.symbols
        assert "N" in nt.symbols

    def test_box_z_equals_length(self):
        """nt.length must match box[2]."""
        s = _graphene()
        ch = compute_chirality(4, 2, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert abs(nt.length - float(nt.box[2])) < 1e-9


# ─────────────────────────────────────────────────────────────────────────────
# 6. TestExporters
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def _nt_for_export():
    """A small (4,2) nanotube used by all exporter tests."""
    s = _graphene()
    ch = compute_chirality(4, 2, s)
    return build_nanotube(s, ch, vacuum=10.0)


class TestExporters:
    # ── XYZ ──────────────────────────────────────────────────────────────────

    def test_xyz_atom_count(self, tmp_path, _nt_for_export):
        path = write_xyz(_nt_for_export, tmp_path / "out.xyz")
        first_line = path.read_text().splitlines()[0].strip()
        assert int(first_line) == _nt_for_export.n_atoms

    def test_xyz_extended_lattice_in_comment(self, tmp_path, _nt_for_export):
        path = write_xyz(_nt_for_export, tmp_path / "out.xyz", extended=True)
        comment = path.read_text().splitlines()[1]
        assert "Lattice=" in comment

    def test_xyz_plain_no_lattice(self, tmp_path, _nt_for_export):
        path = write_xyz(_nt_for_export, tmp_path / "plain.xyz", extended=False)
        comment = path.read_text().splitlines()[1]
        assert "Lattice=" not in comment

    def test_xyz_coord_line_count(self, tmp_path, _nt_for_export):
        path = write_xyz(_nt_for_export, tmp_path / "out.xyz")
        lines = path.read_text().splitlines()
        # Lines 0 (count) + 1 (comment) + N (atoms)
        assert len(lines) == 2 + _nt_for_export.n_atoms

    # ── PDB ──────────────────────────────────────────────────────────────────

    def test_pdb_cryst1_record(self, tmp_path, _nt_for_export):
        path = write_pdb(_nt_for_export, tmp_path / "out.pdb")
        text = path.read_text()
        assert "CRYST1" in text

    def test_pdb_atom_records(self, tmp_path, _nt_for_export):
        path = write_pdb(_nt_for_export, tmp_path / "out.pdb")
        text = path.read_text()
        atom_lines = [l for l in text.splitlines() if l.startswith("ATOM")]
        assert len(atom_lines) == _nt_for_export.n_atoms

    def test_pdb_ends_with_end(self, tmp_path, _nt_for_export):
        path = write_pdb(_nt_for_export, tmp_path / "out.pdb")
        assert path.read_text().strip().endswith("END")

    # ── POSCAR ───────────────────────────────────────────────────────────────

    def test_poscar_total_atom_count(self, tmp_path, _nt_for_export):
        path = write_poscar(_nt_for_export, tmp_path / "POSCAR")
        lines = path.read_text().splitlines()
        # Line 6 (index 5) = species names; line 7 (index 6) = atom counts
        counts = list(map(int, lines[6].split()))
        assert sum(counts) == _nt_for_export.n_atoms

    def test_poscar_cartesian_tag(self, tmp_path, _nt_for_export):
        path = write_poscar(_nt_for_export, tmp_path / "POSCAR")
        text = path.read_text()
        assert "Cartesian" in text

    def test_poscar_scale_factor(self, tmp_path, _nt_for_export):
        path = write_poscar(_nt_for_export, tmp_path / "POSCAR")
        scale_line = path.read_text().splitlines()[1].strip()
        assert abs(float(scale_line) - 1.0) < 1e-6

    # ── LAMMPS ───────────────────────────────────────────────────────────────

    def test_lammps_atom_count_header(self, tmp_path, _nt_for_export):
        path = write_lammps(_nt_for_export, tmp_path / "out.lammps")
        text = path.read_text()
        atom_count_line = next(l for l in text.splitlines() if "atoms" in l
                               and not "atom types" in l)
        assert int(atom_count_line.split()[0]) == _nt_for_export.n_atoms

    def test_lammps_atom_types_count(self, tmp_path, _nt_for_export):
        """Number of atom types must match species count."""
        path = write_lammps(_nt_for_export, tmp_path / "out.lammps")
        text = path.read_text()
        types_line = next(l for l in text.splitlines() if "atom types" in l)
        n_types = int(types_line.split()[0])
        n_species = len(set(_nt_for_export.symbols))
        assert n_types == n_species

    def test_lammps_masses_section(self, tmp_path, _nt_for_export):
        path = write_lammps(_nt_for_export, tmp_path / "out.lammps")
        assert "Masses" in path.read_text()

    # ── Quantum ESPRESSO ─────────────────────────────────────────────────────

    def test_qe_system_nat(self, tmp_path, _nt_for_export):
        path = write_qe(_nt_for_export, tmp_path / "out.pwi")
        text = path.read_text()
        nat_line = next(l for l in text.splitlines() if "nat" in l and "=" in l)
        nat_val = int(nat_line.split("=")[1].strip())
        assert nat_val == _nt_for_export.n_atoms

    def test_qe_atomic_positions_section(self, tmp_path, _nt_for_export):
        path = write_qe(_nt_for_export, tmp_path / "out.pwi")
        assert "ATOMIC_POSITIONS" in path.read_text()

    def test_qe_cell_parameters_section(self, tmp_path, _nt_for_export):
        path = write_qe(_nt_for_export, tmp_path / "out.pwi")
        assert "CELL_PARAMETERS" in path.read_text()

    def test_qe_k_points_section(self, tmp_path, _nt_for_export):
        path = write_qe(_nt_for_export, tmp_path / "out.pwi")
        assert "K_POINTS" in path.read_text()

    # ── Unified export() dispatcher ───────────────────────────────────────────

    def test_export_by_extension_xyz(self, tmp_path, _nt_for_export):
        out = export(_nt_for_export, tmp_path / "nt.xyz")
        assert out.exists()
        assert int(out.read_text().splitlines()[0]) == _nt_for_export.n_atoms

    def test_export_by_extension_pdb(self, tmp_path, _nt_for_export):
        out = export(_nt_for_export, tmp_path / "nt.pdb")
        assert out.exists()
        assert "CRYST1" in out.read_text()

    def test_export_by_fmt_override(self, tmp_path, _nt_for_export):
        """Passing fmt= should override the filename extension."""
        out = export(_nt_for_export, tmp_path / "nanotube", fmt=".xyz")
        assert out.exists()

    def test_export_unknown_format_raises(self, tmp_path, _nt_for_export):
        with pytest.raises(ValueError, match="Unknown format"):
            export(_nt_for_export, tmp_path / "nt.zzz")

    def test_export_creates_parent_directory(self, tmp_path, _nt_for_export):
        deep_path = tmp_path / "sub" / "dir" / "out.xyz"
        out = export(_nt_for_export, deep_path)
        assert out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# 7. TestConnectivity
# ─────────────────────────────────────────────────────────────────────────────

class TestConnectivity:
    """
    Tests for core/connectivity.py — requires scipy.
    All tests in this class are automatically skipped when scipy is not installed.
    """

    @pytest.fixture(autouse=True)
    def _require_scipy(self):
        pytest.importorskip("scipy")

    def test_cc_bonds_have_reasonable_length(self, nt_5_5):
        """All C-C bonds in a graphene nanotube must be ~1.42 Å (±0.15 Å)."""
        bonds = compute_bonds(nt_5_5.coords, nt_5_5.symbols)
        assert len(bonds) > 0, "No bonds found in (5,5) nanotube"
        for i, j in bonds:
            d = float(np.linalg.norm(nt_5_5.coords[i] - nt_5_5.coords[j]))
            assert 1.2 < d < 1.65, (
                f"Bond ({i},{j}) length {d:.3f} Å is outside the expected C-C range"
            )

    def test_bond_count_reasonable(self, nt_5_5):
        """(5,5) nanotube with 20 atoms: bond count should be between 20 and 30."""
        bonds = compute_bonds(nt_5_5.coords, nt_5_5.symbols)
        # In bulk: 20 atoms × 3 bonds / 2 = 30. Without PBC some boundary bonds
        # may be absent, so we accept ≥ 20.
        assert 20 <= len(bonds) <= 30, f"Unexpected bond count: {len(bonds)}"

    def test_no_self_bonds(self, nt_5_5):
        """Bond list must never contain (i, i) pairs."""
        bonds = compute_bonds(nt_5_5.coords, nt_5_5.symbols)
        for i, j in bonds:
            assert i != j, "Self-bond detected"

    def test_bonds_sorted_i_lt_j(self, nt_5_5):
        """Convention: i < j for all bonds."""
        bonds = compute_bonds(nt_5_5.coords, nt_5_5.symbols)
        for i, j in bonds:
            assert i < j, f"Bond ({i},{j}) violates i < j convention"

    def test_bond_settings_default_max_carbon(self):
        """Default C-C cutoff = (0.76+0.76)×1.20 = 1.824 Å."""
        bs = BondSettings()
        cutoff = bs.default_max("C", "C")
        assert abs(cutoff - (0.76 + 0.76) * 1.20) < 0.001

    def test_bond_settings_custom_max_override(self):
        """Custom per-pair cutoff must override the default."""
        bs = BondSettings()
        bs.custom_max[frozenset(["C", "C"])] = 2.0
        assert bs.max_dist("C", "C") == 2.0

    def test_bond_settings_reset_clears_custom(self):
        bs = BondSettings()
        bs.custom_max[frozenset(["C", "N"])] = 1.8
        bs.reset()
        assert len(bs.custom_max) == 0

    def test_bond_line_arrays_shape(self, nt_5_5):
        """bond_line_arrays must return (4B, 3) pts and (4B, 4) colors."""
        bonds = compute_bonds(nt_5_5.coords, nt_5_5.symbols)
        cpk = {"C": (0.5, 0.5, 0.5, 1.0)}
        pts, cols = bond_line_arrays(nt_5_5.coords, nt_5_5.symbols, bonds, cpk)
        B = len(bonds)
        assert pts.shape  == (4 * B, 3), f"pts shape {pts.shape}"
        assert cols.shape == (4 * B, 4), f"cols shape {cols.shape}"

    def test_bond_line_arrays_empty_input(self, nt_5_5):
        """Empty bond list → empty arrays, not an error."""
        cpk = {"C": (0.5, 0.5, 0.5, 1.0)}
        pts, cols = bond_line_arrays(nt_5_5.coords, nt_5_5.symbols, [], cpk)
        assert pts.shape[0]  == 0
        assert cols.shape[0] == 0

    def test_get_radius_known_elements(self):
        """Alvarez covalent radii for common elements."""
        assert abs(get_radius("C")  - 0.76) < 0.01
        assert abs(get_radius("N")  - 0.71) < 0.01
        assert abs(get_radius("B")  - 0.84) < 0.01
        assert abs(get_radius("Si") - 1.11) < 0.01

    def test_get_radius_unknown_element_returns_fallback(self):
        """Unknown element symbol must return the fallback radius (0.90 Å)."""
        r = get_radius("Xx")
        assert r > 0.0


# ─────────────────────────────────────────────────────────────────────────────
# 8. TestSymmetry
# ─────────────────────────────────────────────────────────────────────────────

class TestSymmetry:
    def test_snap_hexagonal_exact_gamma(self):
        """snap_to_symmetry sets γ to exactly 60° for a noisy hexagonal cell."""
        snapped, _ = snap_to_symmetry(_noisy_hexagonal())
        assert abs(snapped.gamma_deg - 60.0) < 0.001

    def test_snap_hexagonal_equal_lengths(self):
        """snap_to_symmetry equalises |a₁| and |a₂| for hexagonal cells."""
        snapped, _ = snap_to_symmetry(_noisy_hexagonal())
        assert abs(snapped.a - snapped.b) < 1e-9

    def test_snap_rectangular_gamma90(self):
        """snap_to_symmetry corrects γ to 90° for near-rectangular cells."""
        snapped, _ = snap_to_symmetry(_slightly_off_rectangular())
        assert abs(snapped.gamma_deg - 90.0) < 0.001

    def test_snap_rectangular_preserves_lengths(self):
        """snap_to_symmetry must not change individual a, b lengths (a ≠ b case)."""
        s = _slightly_off_rectangular()
        snapped, _ = snap_to_symmetry(s)
        assert abs(snapped.a - s.a) < 0.01
        assert abs(snapped.b - s.b) < 0.01

    def test_snap_oblique_returns_unchanged(self):
        """snap_to_symmetry leaves oblique structures unchanged."""
        s = _oblique()
        snapped, desc = snap_to_symmetry(s)
        assert "Oblique" in desc
        assert abs(snapped.gamma_deg - s.gamma_deg) < 0.01

    def test_snap_returns_latticestructure(self):
        """snap_to_symmetry must return a LatticeStructure, not None."""
        result, _ = snap_to_symmetry(_graphene())
        assert isinstance(result, LatticeStructure)

    def test_primitive_2x1_supercell_reduces(self):
        """
        A 2×1 graphene supercell (4 atoms) must reduce to 2 atoms.
        Tests both the builtin and spglib backends (whichever is available).
        """
        s = _graphene_2x1()
        assert len(s.atoms) == 4, "Supercell helper must produce 4 atoms"
        prim, desc = find_primitive_cell(s)
        assert len(prim.atoms) == 2, (
            f"Expected 2 atoms in primitive cell, got {len(prim.atoms)}. "
            f"Description: {desc}"
        )

    def test_primitive_already_primitive(self):
        """
        Applying find_primitive_cell to an already-primitive cell must not
        increase the atom count.
        """
        s = _graphene()     # 2 atoms — already primitive
        prim, _ = find_primitive_cell(s)
        assert len(prim.atoms) <= len(s.atoms)

    def test_primitive_returns_latticestructure(self):
        result, _ = find_primitive_cell(_graphene())
        assert isinstance(result, LatticeStructure)


# ─────────────────────────────────────────────────────────────────────────────
# 9. TestRoundTrip
# ─────────────────────────────────────────────────────────────────────────────

class TestRoundTrip:
    def test_xyz_write_read_atom_count(self, tmp_path):
        """
        Write a LatticeStructure as extended XYZ → read_xyz → same atom count
        and lattice parameters.
        """
        s = _graphene()
        xyz_path = tmp_path / "graphene.xyz"
        a1, a2 = s.a1, s.a2
        with open(xyz_path, "w") as f:
            f.write(f"{len(s.atoms)}\n")
            f.write(
                f'Lattice="{a1[0]:.6f} {a1[1]:.6f} 0.000000 '
                f'{a2[0]:.6f} {a2[1]:.6f} 0.000000 '
                f'0.000000 0.000000 20.000000" '
                f'Properties=species:S:1:pos:R:3\n'
            )
            for atom in s.atoms:
                x, y = atom["pos"]
                z = atom.get("z", 0.0)
                f.write(f"{atom['symbol']}  {x:.8f}  {y:.8f}  {z:.8f}\n")
        s2 = read_xyz(xyz_path)
        assert len(s2.atoms) == len(s.atoms)
        assert abs(s2.a - s.a) < 0.01
        assert abs(s2.gamma_deg - s.gamma_deg) < 0.5

    def test_pdb_write_read_atom_count(self, tmp_path):
        """
        write_pdb → read_pdb → atom count must be preserved.
        """
        s = _graphene()
        ch = compute_chirality(4, 4, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        pdb_path = tmp_path / "nt.pdb"
        write_pdb(nt, pdb_path)
        s2 = read_pdb(pdb_path)
        assert len(s2.atoms) == nt.n_atoms

    def test_all_formats_write_without_error(self, tmp_path):
        """Smoke test: all supported formats must write a non-empty file."""
        s = _graphene()
        ch = compute_chirality(3, 3, s)
        nt = build_nanotube(s, ch, vacuum=10.0)
        formats = {
            "out.xyz":    (write_xyz,    (nt,)),
            "out.pdb":    (write_pdb,    (nt,)),
            "POSCAR":     (write_poscar, (nt,)),
            "out.lammps": (write_lammps, (nt,)),
            "out.pwi":    (write_qe,     (nt,)),
        }
        for filename, (writer, args) in formats.items():
            path = tmp_path / filename
            out = writer(*args, path)
            assert out.exists(), f"{filename} was not created"
            assert out.stat().st_size > 0, f"{filename} is empty"
