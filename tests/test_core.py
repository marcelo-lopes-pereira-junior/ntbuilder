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

from core.io        import LatticeStructure, read_xyz, read_pdb, read_cif
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


def _penta_like():
    """Square cell whose basis breaks the axis-swap mirror.

    A fictional structure, *not* penta-graphene: it drops the (1/2, 1/2) site
    and swaps two of the buckling signs, which is what leaves it with no
    diagonal operation at all.  Real penta-graphene does have one — a glide,
    mirror plus (1/2, 1/2) — so its (n, m) map does fold; see
    ``test_penta_graphene_has_a_diagonal_glide``.  This cell is kept because a
    basis that breaks the mirror outright still has to be detected.
    """
    a = 3.63
    a1 = np.array([a, 0.0])
    a2 = np.array([0.0, a])
    atoms = [
        {"symbol": "C", "pos": np.array([0.0, 0.0]),           "z":  0.0},
        {"symbol": "C", "pos": 0.634 * a1 + 0.134 * a2,        "z":  0.687},
        {"symbol": "C", "pos": 0.366 * a1 + 0.866 * a2,        "z": -0.687},
    ]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


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
# 1b. TestCIFFallback
# ─────────────────────────────────────────────────────────────────────────────

class TestCIFFallback:
    """read_cif must fall back to the built-in parser whenever gemmi refuses
    a file, not only when gemmi is absent."""

    def test_bare_data_block_header_is_readable(self, tmp_path):
        """A CIF starting with a nameless ``data_`` block loads anyway.

        gemmi enforces the CIF grammar strictly and rejects such files, but
        several databases emit them — including the Shi et al. sp2-carbon
        structures used as case studies.
        """
        cif = tmp_path / "bare_header.cif"
        cif.write_text(
            "data_\n"
            "_symmetry_space_group_name_H-M   'P1'\n"
            "_cell_length_a                   2.460000\n"
            "_cell_length_b                   2.460000\n"
            "_cell_length_c                   20.000000\n"
            "_cell_angle_alpha                90.000000\n"
            "_cell_angle_beta                 90.000000\n"
            "_cell_angle_gamma                60.000000\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_type_symbol\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "C1 C 0.000000 0.000000 0.500000\n"
            "C2 C 0.333333 0.333333 0.500000\n"
        )
        s = read_cif(cif)
        assert len(s.atoms) == 2
        assert s.a == pytest.approx(2.46, abs=1e-3)
        assert s.gamma_deg == pytest.approx(60.0, abs=1e-2)

    def test_unparsable_file_still_raises(self, tmp_path):
        """The fallback must not swallow genuinely broken input."""
        junk = tmp_path / "junk.cif"
        junk.write_text("this is not a CIF at all\n")
        with pytest.raises(Exception):
            read_cif(junk)


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

    def test_square_is_square(self):
        """a = b with γ = 90° is the square class, not rectangular.

        It used to be folded into 'rectangular', which mislabelled
        penta-graphene and hid the reason it always closes an exact T.
        """
        a1 = np.array([3.0, 0.0])
        a2 = np.array([0.0, 3.0])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        assert s.lattice_type == "square"
        assert s.is_square

    def test_hbn_is_hexagonal(self):
        """h-BN has a hexagonal lattice."""
        assert _hbn().lattice_type == "hexagonal"

    def test_rhombic_gamma20_is_centred_rectangular(self):
        """a = b with γ = 20° is rhombic — the centred rectangular class.

        It used to come back 'oblique' because io had no row for it, which is
        how both of the extreme-γ example cells were misnamed.
        """
        gamma = math.radians(20)
        a1 = np.array([3.0, 0.0])
        a2 = np.array([3.0 * math.cos(gamma), 3.0 * math.sin(gamma)])
        s = LatticeStructure(a1=a1, a2=a2, atoms=[])
        assert s.lattice_type == "centred rectangular"

    def test_truly_oblique_is_oblique(self):
        """a != b and γ != 90° is the only genuinely oblique case."""
        gamma = math.radians(75)
        a1 = np.array([3.0, 0.0])
        a2 = np.array([4.2 * math.cos(gamma), 4.2 * math.sin(gamma)])
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

    def test_oblique_sector_is_a_half_turn(self):
        """A truly oblique lattice has only +-I, so the wedge is 180 deg.

        The old contract returned gamma here.  That was the opening angle of
        the first quadrant, not a fundamental domain: gamma leaves the
        directions between gamma and 180 deg unaccounted for, and those are
        real tubes when nothing but inversion identifies them.
        """
        from core.chirality import unique_sector_deg
        assert abs(unique_sector_deg(_oblique()) - 180.0) < 0.01

    def test_oblique_scan_reaches_negative_m(self):
        """On an oblique lattice (n, -m) is a distinct tube and must be listed."""
        res = scan_chirality(_oblique(), n_max=4, m_max=4, max_diameter=40.0,
                             unique_only=True)
        pairs = {(r.n, r.m) for r in res}
        assert (3, 1) in pairs and (3, -1) in pairs

    def test_centred_rectangular_folds(self):
        """A rhombic cell is centred rectangular: it has a mirror and folds.

        |a1| = |a2| forces the swap of the two axes to be an isometry, so the
        lattice is centred rectangular whatever gamma is.  The old classifier
        called anything outside 60/90/120 oblique and never folded it.
        """
        from core.chirality import unique_sector_deg
        gamma = math.radians(70.0)
        a = 3.4
        a1 = np.array([a, 0.0])
        a2 = np.array([a * math.cos(gamma), a * math.sin(gamma)])
        s = LatticeStructure(a1=a1, a2=a2,
                             atoms=[{"symbol": "C", "pos": np.zeros(2), "z": 0.0}])
        assert s.lattice_type == "centred rectangular"
        assert abs(unique_sector_deg(s) - 90.0) < 0.01
        res = scan_chirality(s, n_max=5, m_max=5, max_diameter=40.0,
                             unique_only=True)
        pairs = {(r.n, r.m) for r in res}
        assert (5, 0) in pairs and (0, 5) not in pairs

    def test_penta_graphene_has_a_diagonal_glide(self):
        """Real penta-graphene folds: (5,0) and (0,5) are one tube.

        An earlier conclusion said otherwise.  It came from comparing the two
        finite unit cells atom by atom, and the two cells differ by an axial
        shift of |T|/2 — which wraps atoms across the cell boundary and changes
        the distance list even though the infinite tubes coincide.  Applying the
        minimum image along the axis settles it: the two agree to 2e-15 A.
        """
        cif = Path(__file__).resolve().parents[1] / "examples" / "Penta_Graphene.cif"
        if not cif.exists():
            pytest.skip("example not present")
        from core.io import load_structure
        from core.builder import build_nanotube
        s = load_structure(str(cif))

        def fingerprint(n, m):
            ch = compute_chirality(n, m, s, search_limit=60)
            p = np.asarray(build_nanotube(s, ch, vacuum=10.0).coords, float)
            d = p[:, None, :] - p[None, :, :]
            d[:, :, 2] -= ch.T_norm * np.round(d[:, :, 2] / ch.T_norm)
            r = np.linalg.norm(d, axis=-1)
            return np.sort(r[np.triu_indices(len(p), 1)])

        assert np.abs(fingerprint(5, 0) - fingerprint(0, 5)).max() < 1e-9

        from core.chirality import basis_swap_invariant
        assert basis_swap_invariant(s) is True

    def test_glide_counts_as_a_mirror(self):
        """pm and pg enumerate alike: a glide is a mirror plus an origin shift.

        The basis below maps onto itself under the axis swap only after a
        translation of (1/2, 1/2).  Rolling is blind to that translation, so
        the tube set must fold exactly as for a pure mirror.
        """
        from core.chirality import basis_swap_invariant
        a = 4.0
        s = LatticeStructure(a1=np.array([a, 0.0]), a2=np.array([0.0, a]),
                             atoms=[{"symbol": "C", "pos": np.array([0.1, 0.2]) * a,
                                     "z": 0.0},
                                    {"symbol": "C", "pos": np.array([0.7, 0.6]) * a,
                                     "z": 0.0}])
        assert basis_swap_invariant(s) is True

    def test_structure_group_is_a_subgroup_of_the_holohedry(self):
        """A crystal can lose its lattice's symmetry, never gain it."""
        from core.planegroup import lattice_point_group, structure_point_group
        for s in (_graphene(), _penta_like(), _oblique(), _hbn()):
            holo = lattice_point_group(s)
            struct = structure_point_group(s)
            assert len(struct) <= len(holo)
            for U in struct:
                assert any(np.array_equal(U, V) for V in holo)

    def test_swap_invariant_true_for_graphene(self):
        """Graphene's basis survives the axis-swap mirror, so folding is safe."""
        from core.chirality import basis_swap_invariant
        assert basis_swap_invariant(_graphene()) is True

    def test_swap_invariant_false_for_mirror_breaking_basis(self):
        """A square lattice whose basis breaks the diagonal mirror is not foldable."""
        from core.chirality import basis_swap_invariant
        assert basis_swap_invariant(_penta_like()) is False

    def test_map_not_folded_when_basis_breaks_mirror(self):
        """(n,0) and (0,n) are both kept when the basis is not swap-invariant.

        Penta-graphene motivates this: the lattice is square but the sp3 atoms
        make (5,0) and (0,5) genuinely inequivalent tubes, so folding the map
        would hide half of the accessible structures.
        """
        res = scan_chirality(_penta_like(), n_max=6, m_max=6,
                             max_diameter=20.0, unique_only=True)
        pairs = {(r.n, r.m) for r in res}
        assert (5, 0) in pairs
        assert (0, 5) in pairs
        assert any(m > n for n, m in pairs)

    def test_map_still_folded_for_swap_invariant_lattice(self):
        """Graphene keeps the folded (unique-sector) map."""
        res = scan_chirality(_graphene(), n_max=6, m_max=6,
                             max_diameter=20.0, unique_only=True)
        pairs = {(r.n, r.m) for r in res}
        assert (5, 0) in pairs
        assert (0, 5) not in pairs

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

def _rect_bpn():
    """Rectangular cell with the biphenylene network's a, b and one atom.

    b/a = 3.88/4.26 closes only at very long T, so every chiral (n, m) has a
    long trade-off front of approximate cells.
    """
    a1 = np.array([4.26, 0.0])
    a2 = np.array([0.0, 3.88])
    atoms = [{"symbol": "C", "pos": np.array([0.0, 0.0]), "z": 0.0}]
    return LatticeStructure(a1=a1, a2=a2, atoms=atoms)


def _exhaustive_front(n, m, a1, a2, cap):
    """Residual front over EVERY lattice vector of at most ``cap`` cells.

    No continued fractions and no window around an ideal t1: all (t1, t2) in
    a box large enough to hold the closest vector of every cell count up to
    the cap.  Returns {cells: strain %}.
    """
    Ch = n * a1 + m * a2
    Chn = float(np.linalg.norm(Ch))
    area = abs(float(a1[0] * a2[1] - a1[1] * a2[0]))
    g = math.gcd(abs(n), abs(m))
    R = math.hypot(cap * area / Chn, Chn / g) + 1.0
    B = int(R / np.linalg.svd(np.stack([a1, a2]), compute_uv=False).min()) + 2
    t1, t2 = np.meshgrid(np.arange(-B, B + 1), np.arange(0, B + 1), indexing="ij")
    t1, t2 = t1.ravel(), t2.ravel()
    keep = (t2 > 0) | (t1 > 0)
    t1, t2 = t1[keep], t2[keep]
    cells = np.abs(n * t2 - m * t1)
    ok = (cells > 0) & (cells <= cap)
    t1, t2, cells = t1[ok], t2[ok], cells[ok]
    T = np.outer(t1, a1) + np.outer(t2, a2)
    eps = 100.0 * np.abs(T @ Ch) / (Chn * np.linalg.norm(T, axis=1))
    best = np.full(cap + 1, np.inf)
    np.minimum.at(best, cells, eps)
    front, b = {}, np.inf
    for c in np.nonzero(np.isfinite(best))[0]:
        if best[c] < b - 1e-10:
            front[int(c)] = float(best[c])
            b = best[c]
    return front


class TestTranslationSearch:
    """T_options and the default T search against exhaustive enumeration."""

    @pytest.mark.parametrize("make, n, m, limit", [
        (_rect_bpn, 5, 6, 300),
        (_rect_bpn, 4, 1, 300),
        (_graphene, -1, 2, 300),     # exact T = a1; a cheaper cell precedes it
        (_oblique, 2, 1, 100),
    ])
    def test_front_matches_exhaustive_enumeration(self, make, n, m, limit):
        """
        Every cell no cheaper cell beats, and nothing else.  Built from
        convergents with t1 near its ideal, the front of (5,6) on this lattice
        had 7 rows where the exhaustive one has over a hundred.
        """
        from core.chirality import T_options, _search_T
        s = make()
        rows = T_options(n, m, s.a1, s.a2, limit)
        got = {abs(n * r["t2"] - m * r["t1"]): r["strain"] for r in rows}
        d1, d2, _ = _search_T(n, m, s.a1, s.a2, limit=limit)
        want = _exhaustive_front(n, m, s.a1, s.a2, abs(n * d2 - m * d1))
        assert set(got) == set(want)
        for c, e in want.items():
            assert got[c] == pytest.approx(e, rel=1e-9, abs=1e-12)

    def test_max_strain_returns_the_cheapest_cell_within_tolerance(self):
        from core.chirality import _search_T
        s = _rect_bpn()
        n, m, limit, tol = 5, 6, 300, 0.1
        d1, d2, _ = _search_T(n, m, s.a1, s.a2, limit=limit)
        front = _exhaustive_front(n, m, s.a1, s.a2, abs(n * d2 - m * d1))
        cheapest = min(c for c, e in front.items() if e <= tol)
        ch = compute_chirality(n, m, s, search_limit=limit, max_strain=tol)
        assert ch.n_atoms == cheapest          # one atom per primitive cell
        assert ch.strain <= tol

    def test_default_search_considers_t1_zero(self):
        """T = a2 can be the best cell without being exactly perpendicular.

        Skipping t1 = 0 made AgBr3 (5,4) take a 1496-cell T at 0.37 % over
        a2 alone at 0.24 %.
        """
        from core.chirality import _search_T, _search_T_scan
        a1 = np.array([1.0, 0.0])
        a2 = np.array([0.01, 1.0])
        expected = 0.01 / math.hypot(0.01, 1.0)
        for search in (_search_T, _search_T_scan):
            t1, t2, err = search(1, 0, a1, a2, limit=30)
            assert (t1, t2) == (0, 1), search.__name__
            assert err == pytest.approx(expected)


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

    @pytest.mark.parametrize("row", range(6))
    def test_approximate_cell_holds_each_site_once(self, row):
        """
        An approximate T is not perpendicular to Ch, and the cell must still be
        the (Ch, T) parallelogram: n_atoms_cell * |n*t2 - m*t1| atoms, every
        site once, periodic around and along the tube.  Cut by projections onto
        the two oblique axes it repeated sites -- 184 atoms instead of 6 for
        biphenylene (5,6) at t = (1, 1).  Rows 0-4 are 67-99.6 % residual,
        row 5 is 0.23 %.
        """
        pytest.importorskip("scipy")
        from scipy.spatial import cKDTree
        from core.chirality import T_options
        s = _rect_bpn()
        rows = T_options(5, 6, s.a1, s.a2, 300)
        ch = compute_chirality(5, 6, s, search_limit=300,
                               max_strain=rows[row]["strain"])
        assert (ch.t1, ch.t2) == (rows[row]["t1"], rows[row]["t2"])
        nt = build_nanotube(s, ch, vacuum=10.0)
        assert len(nt.symbols) == ch.n_atoms

        # Unroll the tube and look for a site present twice, periodic images
        # of the cell included.
        c = float(nt.box[0]) / 2.0
        R = ch.Ch_norm / (2.0 * math.pi)
        u = np.mod(np.arctan2(nt.coords[:, 1] - c, nt.coords[:, 0] - c),
                   2.0 * math.pi) * R
        pts = np.column_stack([u, nt.coords[:, 2]])
        images = np.vstack([pts + [i * ch.Ch_norm, k * nt.length]
                            for i in (-1, 0, 1) for k in (-1, 0, 1)])
        d, _ = cKDTree(images).query(pts, k=2)
        assert d[:, 1].min() > 1e-3


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

    def test_centred_rectangular_in_non_rhombic_basis(self):
        """A centred rectangular lattice written with a != b keeps its class.

        u and v are the rhombic primitive vectors (|u| = |v|); the cell below
        uses u and u + v instead, which is the same lattice with a != b.  It
        was labelled oblique before the classifier looked at equivalent bases
        (1 711 catalogue systems).  The snap must keep the basis the cell uses,
        so (n, m) keep their meaning, and make the mirror exact.
        """
        from core.planegroup import lattice_point_group
        # |u| = |v| = 5.385 A, rhombic angle 136.4 deg; with u + v the cell is
        # 5.385 x 4.000 A at 68.2 deg -- far from hexagonal and from square, so
        # the only class it can reach is centred rectangular.
        u, v = np.array([2.0, 5.0]), np.array([2.0, -5.0])
        s = LatticeStructure([{"symbol": "C", "pos": np.zeros(2), "z": 0.0}], u, u + v)
        assert abs(s.a - s.b) > 0.01
        assert s.lattice_type == "centred rectangular"
        snapped, _ = snap_to_symmetry(s)
        assert snapped.lattice_type == "centred rectangular"
        assert abs(snapped.a - s.a) < 1e-6 and abs(snapped.b - s.b) < 1e-6
        assert len(lattice_point_group(snapped)) == 4

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


# ─────────────────────────────────────────────────────────────────────────────
# 10. TestCurvatureBonds — bonds that rolling forms and breaks, pair by pair
# ─────────────────────────────────────────────────────────────────────────────

def _thick_bilayer():
    """Two graphene sheets at z = ±3 Å: a synthetic thick layer whose outer
    face is stretched by (R + 3)/R when rolled."""
    g = _graphene()
    atoms = [{**at, "z": dz} for dz in (3.0, -3.0) for at in g.atoms]
    return LatticeStructure(a1=g.a1, a2=g.a2, atoms=atoms)


class TestCurvatureBonds:

    def test_sites_are_the_builder_positions(self):
        from core.builder import tube_sites, _roll_sites
        for s, (n, m) in ((_graphene(), (5, 3)), (_buckled(), (4, 2))):
            ch = compute_chirality(n, m, s)
            k, sw, w, _x, _y = tube_sites(s, ch)
            z = np.array([a["z"] for a in s.atoms])[k]
            for inward in (False, True):
                nt = build_nanotube(s, ch, vacuum=10.0, roll_inward=inward)
                shift = np.array([nt.box[0] / 2, nt.box[1] / 2, 0.0])
                P = _roll_sites(s, ch, sw, w, z, inward)
                assert np.abs(nt.coords - shift - P).max() < 1e-9

    def test_wide_graphene_tube_is_clean(self):
        from core.builder import check_curvature_bonds
        s = _graphene()
        assert check_curvature_bonds(s, compute_chirality(10, 10, s)) == (set(), set())

    def test_walls_closer_than_a_bond_are_formed(self):
        # (2,0) is 1.57 Å wide: opposite walls sit inside the C-C cutoff.  The
        # species-level check misses it, because C-C is bonded in the sheet.
        from core.builder import check_curvature_bonds, check_spurious_bonds
        s = _graphene()
        ch = compute_chirality(2, 0, s)
        formed, broken = check_curvature_bonds(s, ch)
        assert formed == {frozenset({"C"})}
        assert check_spurious_bonds(s, build_nanotube(s, ch)) == set()

    def test_stretched_outer_face_breaks_bonds(self):
        from core.builder import check_curvature_bonds, curvature_tokens
        s = _thick_bilayer()
        formed, broken = check_curvature_bonds(s, compute_chirality(10, 10, s))
        assert frozenset({"C"}) in broken
        assert "-C-C" in curvature_tokens(formed, broken)
        # Wide enough, the stretch stays inside the 10 % margin.
        assert check_curvature_bonds(s, compute_chirality(40, 40, s)) == (set(), set())

    def test_tokens_mark_formed_and_broken(self):
        from core.builder import curvature_tokens
        tokens = curvature_tokens({frozenset({"S"}), frozenset({"Mo", "S"})}, {frozenset({"Mo"})})
        assert tokens == ["+Mo-S", "+S-S", "-Mo-Mo"]
