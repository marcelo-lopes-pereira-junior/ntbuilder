# NTBuilder — Nanotube Structure Generator

NTBuilder generates nanotube structures from arbitrary 2D crystal inputs. It supports planar, buckled, and Janus materials across all 2D Bravais lattice types (hexagonal, rectangular, oblique) and exports ready-to-use files for VASP, Quantum ESPRESSO, LAMMPS, XYZ, and PDB workflows.

**Web interface (no installation):** [nanoeng.unb.br/ntbuilder](https://nanoeng.unb.br/ntbuilder)

---

## Features

- Load 2D crystal structures from **CIF** files (PDB and XYZ also supported in the desktop app)
- Interactive **polar Hamada chirality map** coloured by strain or atom count, with per-point spurious-bond markers for buckled/Janus lattices and real-time **roll inward** toggle
- **Multi-walled nanotubes (MWNT)**: automatic shell stacking with target interlayer spacing
- **Bundle builder**: periodic supercells in linear, triangle, square (2×2), and hexagonal (1+6) arrangements
- **Deformations**: axial strain and uniform torsion (twist rate in °/Å), applied post-build
- **Axial supercell replication** (1–12×): replicated geometry written to file on export
- **Batch export**: full (n, m) family with diameter and atom-count filters — individual files (desktop) or ZIP archive (web)
- **CLI** (`cli.py`): headless build, batch, polar scan, and MWNT generation
- **Post-build analysis**: bond statistics, electronic character (metallic/semiconducting), line-group symmetry info, and auto-generated Methods paragraph for publications
- Smart large-tube handling: tubes > 500 000 atoms bypass the viewer and export directly via a memory-efficient analytic iterator
- Export to **VASP POSCAR**, **Quantum ESPRESSO** (`.pwi`), **LAMMPS** (`.lammps`), **XYZ**, and **PDB**

---

## Quick Start — Web Interface

No installation needed. Go to [nanoeng.unb.br/ntbuilder](https://nanoeng.unb.br/ntbuilder), upload a CIF file (or pick one of the bundled examples), select (n, m) on the chirality map, and download in your preferred format.

---

## Desktop Installation

```bash
git clone https://github.com/marcelo-lopes-pereira-junior/ntbuilder.git
cd ntbuilder
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python main.py
```

**Requirements:** Python ≥ 3.10, PyQt6 ≥ 6.5, NumPy ≥ 1.24, SciPy ≥ 1.10, matplotlib ≥ 3.7, pyqtgraph ≥ 0.13, PyOpenGL ≥ 3.1, gemmi ≥ 0.6 *(optional, recommended for CIF parsing)*.

---

## CLI

```bash
python cli.py build  graphene.cif --n 10 --m 10 --output nt.xyz
python cli.py batch  graphene.cif --type armchair --nfrom 5 --nto 20 --fmt xyz --outdir ./batch
python cli.py polar  graphene.cif --dmax 30 --output map.csv
python cli.py mwnt   graphene.cif --n 5 --m 5 --walls 3 --spacing 3.4 --output mwnt.xyz
```

---

## Example Structures

| File | Material | Lattice | Type |
|------|----------|---------|------|
| `Graphene.cif` | Graphene | Hexagonal | Planar |
| `Biphenylene_Network.cif` | Biphenylene network | Rectangular | Planar |
| `Penta_Graphene.cif` | Penta-graphene | Rectangular | Buckled |
| `MoS2.cif` | Molybdenum disulfide | Hexagonal | Buckled (TMD) |
| `MoSSe.cif` | Janus MoSSe | Hexagonal | Buckled (Janus) |

---

## Project Structure

```
ntbuilder/
├── core/                # Core library — GUI-free, importable independently
│   ├── io.py            # Structure reader (CIF, PDB, XYZ)
│   ├── chirality.py     # Chirality vectors, T-vector search, Hamada map scan
│   ├── builder.py       # Nanotube rolling algorithm (analytic-bounds iterator)
│   ├── mwnt.py          # Multi-walled nanotube builder
│   ├── bundles.py       # Periodic bundle supercell builder
│   ├── deformations.py  # Axial strain and torsion transformations
│   ├── analysis.py      # Bond statistics, electronic character, Methods text
│   ├── exporters.py     # File writers (VASP, QE, LAMMPS, XYZ, PDB)
│   ├── connectivity.py  # Bond detection and rendering arrays
│   └── symmetry.py      # Lattice snapping and primitive-cell reduction
├── gui/                 # PyQt6 desktop interface
│   ├── main_window.py
│   ├── style.py
│   ├── dialogs/         # Bond settings, advanced operations (MWNT, bundle, deform)
│   └── panels/          # Input, polar map, 3D viewer
├── web/                 # FastAPI web interface
│   ├── api/             # Backend (FastAPI + uvicorn)
│   ├── static/          # Frontend (HTML/CSS/JS — Plotly.js + 3Dmol.js)
│   └── deploy/          # nginx config, systemd service, cleanup cron
├── examples/            # Example CIF files (5 canonical + 4 in extra/)
├── tests/               # pytest suite
├── assets/              # Logos and icons
├── cli.py               # Command-line entry point
├── main.py              # Desktop GUI entry point
└── requirements.txt
```

---

## Citation

If you use NTBuilder in your research, please cite:

> **Marcelo Lopes Pereira Junior**
> *NTBuilder — Nanotube Structure Generator from Arbitrary 2D Lattices* (2026)
> University of Brasília (UnB)
> GitHub: https://github.com/marcelo-lopes-pereira-junior/ntbuilder
> Web: https://nanoeng.unb.br/ntbuilder

---

## Author

**Prof. Dr. Marcelo Lopes Pereira Junior**
University of Brasília (UnB) — NanoEng
marcelo.lopes@unb.br | https://nanoeng.unb.br

---

## License

MIT — see [LICENSE](LICENSE).
