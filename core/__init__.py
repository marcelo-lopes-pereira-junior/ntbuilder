"""
nanotube_builder.core
---------------------
Core scientific modules for nanotube generation from arbitrary 2D lattices.

Quick usage
-----------
>>> from nanotube_builder.core import load_structure, compute_chirality, build_nanotube, export
>>> struct = load_structure("graphene.xyz", a1=[2.46, 0], a2=[1.23, 2.13])
>>> chiral = compute_chirality(6, 6, struct)
>>> tube   = build_nanotube(struct, chiral, vacuum=10.0)
>>> export(tube, "armchair_6_6.pdb")
"""

from .io           import load_structure, LatticeStructure
from .chirality    import compute_chirality, scan_chirality, ChiralityResult
from .builder      import build_nanotube, NanotubeStructure
from .exporters    import export
from .connectivity import compute_bonds, bond_line_arrays, COVALENT_RADII, BondSettings
from .symmetry     import snap_to_symmetry, find_primitive_cell

__all__ = [
    "load_structure", "LatticeStructure",
    "compute_chirality", "scan_chirality", "ChiralityResult",
    "build_nanotube", "NanotubeStructure",
    "export",
    "compute_bonds", "bond_line_arrays", "COVALENT_RADII", "BondSettings",
    "snap_to_symmetry", "find_primitive_cell",
]
