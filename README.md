# NTBuilder — Nanotube Structure Generator

NTBuilder generates nanotube structures from arbitrary 2D crystal inputs. It supports planar, buckled, and Janus materials across all five 2D Bravais lattices (hexagonal, square, rectangular, centred rectangular, oblique) and exports ready-to-use files for VASP, Quantum ESPRESSO, LAMMPS, XYZ, and PDB workflows.

- **Web builder (no installation):** [nanoeng.unb.br/ntbuilder](https://nanoeng.unb.br/ntbuilder)
- **Nanotube catalogue:** [nanoeng.unb.br/ntbuilder/catalogo](https://nanoeng.unb.br/ntbuilder/catalogo) — 20.3 million nanotubes from 46 403 2D systems of eight public 2D databases, with filters and downloads

---

## Features

- Load 2D crystal structures from **CIF** files (PDB and XYZ also supported in the desktop app)
- Interactive **polar chirality map** coloured by strain or atom count. The map shows one (n, m) per family of indices that the lattice symmetry makes equivalent, over the wedge the plane group actually gives, with a **roll inward** toggle for buckled and Janus layers
- **Curvature bond check**: every pair of atoms of the tube is compared with the same pair in the flat sheet (cutoff 1.2 × the sum of covalent radii, 10 % margin). A **formed** bond is a pair squeezed below the cutoff, a **broken** bond a pair stretched past it (for example on the outer face of a thick layer). Both are marked with an X on the map and reported after a build
- **Exact or tolerant cells**: the exact translational vector T can need huge cells on low-symmetry lattices; a maximum strain (e.g. 0.5 %) returns the shortest cell within that tolerance
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

No installation needed. Go to [nanoeng.unb.br/ntbuilder](https://nanoeng.unb.br/ntbuilder), upload a CIF file (or pick one of the bundled examples), select (n, m) on the chirality map, and download in your preferred format. To start from a known 2D material, search it in the [catalogue](https://nanoeng.unb.br/ntbuilder/catalogo) and use **Open in NTBuilder**.

---

## Desktop Installation

Requires **Python ≥ 3.10** and `git`. Install inside a virtual environment: a system-wide matplotlib or PyQt mixed with packages in `~/.local` breaks the 3D projection and can crash the viewer.

**Linux and macOS**

```bash
git clone https://github.com/marcelo-lopes-pereira-junior/ntbuilder.git
cd ntbuilder
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python main.py
```

**Windows (PowerShell)**

```powershell
git clone https://github.com/marcelo-lopes-pereira-junior/ntbuilder.git
cd ntbuilder
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
python main.py
```

**Updating an existing copy**

The published history was rewritten on 13 September 2026. A copy cloned before that date does not update with `git pull`; bring it to the current version once with `git fetch origin && git reset --hard origin/main` (this discards local changes to tracked files), and use `git pull` afterwards.

```bash
cd ntbuilder
git pull
source .venv/bin/activate      # Windows: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

**Check the installation** (optional):

```bash
pytest tests/
python cli.py build examples/Graphene.cif --n 10 --m 10 --output nt.xyz
```

**Requirements:** PyQt6 ≥ 6.5, NumPy ≥ 1.24, SciPy ≥ 1.10, spglib ≥ 2.0, matplotlib ≥ 3.7, pyqtgraph ≥ 0.13, PyOpenGL ≥ 3.1, gemmi ≥ 0.6, pytest ≥ 7.4 (all installed by `requirements.txt`).

### Troubleshooting

- **The window closes at start on Linux with a Wayland message** (`xdg_surface buffer … does not match the configured maximized state`): update to the current version, where the window fits any screen. On an older copy, run `QT_QPA_PLATFORM=xcb python main.py`.
- **`Unable to import Axes3D`**: matplotlib is installed twice (system and `pip --user`). Use the virtual environment above.
- **Blank 3D viewer**: the viewer needs OpenGL. Update the graphics driver, or on a remote session run with `QT_QPA_PLATFORM=xcb` and a working X server.

---

## CLI

```bash
python cli.py build  examples/Graphene.cif --n 10 --m 10 --output nt.xyz
python cli.py batch  examples/Graphene.cif --type armchair --nfrom 5 --nto 20 --fmt xyz --outdir ./batch
python cli.py polar  examples/Graphene.cif --dmax 30 --output map.csv
python cli.py mwnt   examples/Graphene.cif --n 5 --m 5 --walls 3 --spacing 3.4 --output mwnt.xyz
```

---

## Example Structures

| File | Material | Lattice | Type |
|------|----------|---------|------|
| `Graphene.cif` | Graphene | Hexagonal | Planar |
| `Biphenylene_Network.cif` | Biphenylene network | Rectangular | Planar |
| `Penta_Graphene.cif` | Penta-graphene | Square | Buckled |
| `MoS2.cif` | Molybdenum disulfide | Hexagonal | Buckled (TMD) |
| `MoSSe.cif` | Janus MoSSe | Hexagonal | Buckled (Janus) |

---

## Running the Web Interface Locally

The builder page runs with the packages of `web/requirements.txt`:

```bash
source .venv/bin/activate
pip install -r web/requirements.txt
cd web
uvicorn api.main:app --host 127.0.0.1 --port 8765
```

Then open http://127.0.0.1:8765. The catalogue page needs its SQLite database (`web/data/catalogue.db`, several GB), which is built by the cluster campaign described in `docs/CATALOGUE_DESIGN.md` and is not distributed with the repository. Production deployment files (nginx, systemd, Slurm job) are in `web/deploy/`.

---

## Project Structure

```
ntbuilder/
├── core/                # Core library — GUI-free, importable independently
│   ├── io.py            # Structure reader (CIF, PDB, XYZ), Bravais classification
│   ├── chirality.py     # Chirality vectors, T-vector search, chirality map scan
│   ├── planegroup.py    # Plane-group operations on (n, m) and the irreducible wedge
│   ├── symmetry.py      # Lattice snapping, primitive cell, chirality group
│   ├── builder.py       # Rolling algorithm and the curvature bond check
│   ├── mwnt.py          # Multi-walled nanotube builder
│   ├── bundles.py       # Periodic bundle supercell builder
│   ├── deformations.py  # Axial strain and torsion transformations
│   ├── analysis.py      # Bond statistics, electronic character, Methods text
│   ├── exporters.py     # File writers (VASP, QE, LAMMPS, XYZ, PDB)
│   └── connectivity.py  # Bond detection and rendering arrays
├── gui/                 # PyQt6 desktop interface
│   ├── main_window.py
│   ├── style.py
│   ├── dialogs/         # Bond settings, advanced operations (MWNT, bundle, deform)
│   └── panels/          # Input, polar map, 3D viewer
├── web/                 # FastAPI web interface
│   ├── api/             # Builder API, catalogue API (catalogue.py, cat_query.py, cat_build.py), mail
│   ├── static/          # Builder and catalogue pages (HTML/CSS/JS — Plotly.js + 3Dmol.js)
│   └── deploy/          # nginx, systemd, cleanup cron, catalogue install and Slurm request job
├── docs/                # Catalogue design notes
├── examples/            # Example CIF files
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
