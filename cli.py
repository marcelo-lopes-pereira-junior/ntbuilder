#!/usr/bin/env python3
"""
NTBuilder CLI — headless nanotube generation from the command line.

Usage
-----
  python cli.py build  graphene.cif --n 10 --m 10 --output nt.xyz
  python cli.py batch  graphene.cif --type armchair --nfrom 5 --nto 20 --fmt xyz --outdir ./batch
  python cli.py polar  graphene.cif --dmax 30 --output map.csv
  python cli.py mwnt   graphene.cif --n 5 --m 5 --walls 3 --spacing 3.4 --output mwnt.xyz
  python cli.py deform nt.xyz --strain 0.05 --twist 1.0 --output nt_deformed.xyz
  python cli.py bundle nt.xyz --geometry hexagonal7 --spacing 3.4 --output bundle.xyz
  python cli.py prim   graphene.cif
  python cli.py query  --formula MoS2 --source cod

Run ``python cli.py <command> --help`` for per-command options.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# ── Bootstrap core importability ──────────────────────────────────────────────
_ROOT = Path(__file__).parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _load(path: str):
    from core.io import load_structure
    return load_structure(path)


def _snap(structure):
    try:
        from core.symmetry import snap_to_symmetry
        structure, desc = snap_to_symmetry(structure)
        if desc:
            print(f"[snap]  {desc}", file=sys.stderr)
    except Exception:
        pass
    return structure


def _chirality(n, m, structure):
    from core.chirality import compute_chirality
    ch = compute_chirality(n, m, structure)
    if ch is None:
        sys.exit(f"Error: (n,m) = ({n},{m}) is degenerate (n=m=0).")
    return ch


def _build(structure, chirality, vacuum, roll_inward):
    from core.builder import build_nanotube
    return build_nanotube(structure, chirality,
                          vacuum=vacuum, roll_inward=roll_inward)


def _export(nt, path: str):
    from core.exporters import export
    out = export(nt, path)
    print(f"Saved → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Commands
# ─────────────────────────────────────────────────────────────────────────────

def cmd_build(args):
    """Build a single nanotube and export."""
    structure  = _snap(_load(args.input))
    chirality  = _chirality(args.n, args.m, structure)

    print(f"Building ({args.n},{args.m})  D={chirality.diameter:.4f} Å  "
          f"L={chirality.T_norm:.4f} Å  atoms={chirality.n_atoms}")

    nt = _build(structure, chirality, args.vacuum, args.roll_inward)

    # Axial replication
    if args.repeat > 1:
        from core.builder import NanotubeStructure
        import numpy as np
        Lz  = float(nt.box[2])
        new_coords = np.vstack([nt.coords + np.array([0,0,i*Lz]) for i in range(args.repeat)])
        nt = NanotubeStructure(
            chirality=nt.chirality, symbols=list(nt.symbols)*args.repeat,
            coords=new_coords, box=np.array([nt.box[0], nt.box[1], Lz*args.repeat]),
            vacuum=nt.vacuum,
        )
        print(f"Replicated ×{args.repeat}  →  {nt.n_atoms} atoms")

    _export(nt, args.output)


def cmd_mwnt(args):
    """Build a multi-walled nanotube."""
    from core.mwnt import build_mwnt, mwnt_summary

    structure  = _snap(_load(args.input))
    chirality  = _chirality(args.n, args.m, structure)

    print(f"Building MWNT  inner=({args.n},{args.m})  "
          f"walls={args.walls}  spacing={args.spacing} Å ...")

    result = build_mwnt(
        structure, chirality,
        n_walls=args.walls,
        interlayer_spacing=args.spacing,
        vacuum=args.vacuum,
        roll_inward=args.roll_inward,
    )

    print(mwnt_summary(result))
    _export(result.nanotube, args.output)


def cmd_deform(args):
    """Apply deformations (strain, torsion) to an existing structure file."""
    from core.deformations import apply_axial_strain, apply_torsion
    from core.io import load_structure

    # Load the nanotube from file — must be xyz or similar
    # We re-wrap it as a NanotubeStructure with dummy chirality
    print("Loading structure for deformation...")
    structure  = _load(args.input)

    # We need a NanotubeStructure; for raw XYZ/POSCAR we create a stub
    try:
        from core.builder import NanotubeStructure
        from core.chirality import ChiralityResult
        import numpy as np

        atoms   = structure.atoms
        symbols = [a["symbol"] for a in atoms]
        coords  = np.array([[a["pos"][0], a["pos"][1], a.get("z", 0.0)]
                             for a in atoms])

        # Dummy chirality for export purposes
        dummy_ch = ChiralityResult(n=0, m=0, a1=structure.a1, a2=structure.a2)

        box_xy = float(np.ptp(coords[:, 0])) + 2.0 * 10.0
        box_z  = float(np.ptp(coords[:, 2]))

        nt = NanotubeStructure(
            chirality=dummy_ch,
            symbols=symbols,
            coords=coords,
            box=np.array([box_xy, box_xy, box_z]),
            vacuum=10.0,
        )
    except Exception as exc:
        sys.exit(f"Error wrapping structure: {exc}")

    if args.strain != 0.0:
        nt = apply_axial_strain(nt, args.strain)
        print(f"Applied axial strain {args.strain*100:+.2f}%")
    if args.twist != 0.0:
        nt = apply_torsion(nt, args.twist)
        print(f"Applied torsion {args.twist:+.4f} °/Å")

    _export(nt, args.output)


def cmd_bundle(args):
    """Build a nanotube bundle from an existing nanotube file."""
    from core.bundles import build_bundle, GEOMETRY_LABELS

    print(f"Building bundle  geometry={args.geometry}  spacing={args.spacing} Å ...")

    # Load source nanotube (same stub approach as deform)
    structure = _load(args.input)
    try:
        from core.builder import NanotubeStructure
        from core.chirality import ChiralityResult
        import numpy as np

        atoms   = structure.atoms
        symbols = [a["symbol"] for a in atoms]
        coords  = np.array([[a["pos"][0], a["pos"][1], a.get("z", 0.0)] for a in atoms])
        dummy_ch = ChiralityResult(n=0, m=0, a1=structure.a1, a2=structure.a2)
        box_xy   = float(np.ptp(coords[:, :2].max(0) - coords[:, :2].min(0))) + 2.0 * 10.0
        box_z    = float(np.ptp(coords[:, 2]))
        nt = NanotubeStructure(
            chirality=dummy_ch, symbols=symbols, coords=coords,
            box=np.array([box_xy, box_xy, box_z]), vacuum=10.0,
        )
    except Exception as exc:
        sys.exit(f"Error wrapping structure: {exc}")

    result = build_bundle(nt, geometry=args.geometry, spacing=args.spacing,
                          nx=args.nx, ny=args.ny)
    label = GEOMETRY_LABELS.get(args.geometry, args.geometry)
    print(f"{label}  pitch={result.pitch:.4f} Å  {result.n_tubes} tubes  "
          f"→ {result.nanotube.n_atoms} atoms total")
    _export(result.nanotube, args.output)


def cmd_batch(args):
    """Batch-build a series of nanotubes and export each to a directory."""
    from core.chirality import scan_chirality, compute_chirality
    from core.builder import build_nanotube

    structure = _snap(_load(args.input))
    outdir    = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Build list of (n, m) pairs
    if args.type == "all":
        results = scan_chirality(structure, n_max=args.nmax,
                                 max_diameter=args.dmax, unique_only=True)
        pairs = [(r.n, r.m) for r in results]
    else:
        pairs = []
        for n_idx in range(args.nfrom, args.nto + 1):
            if args.type == "armchair":
                pairs.append((n_idx, n_idx))
            elif args.type == "zigzag":
                pairs.append((n_idx, 0))

    print(f"Batch: {len(pairs)} nanotubes → {outdir}/")

    ok, fail = 0, 0
    for n, m in pairs:
        try:
            ch = compute_chirality(n, m, structure)
            if ch is None:
                continue
            nt = build_nanotube(structure, ch, vacuum=args.vacuum)
            name = outdir / f"nt_{n}_{m}.{args.fmt.lstrip('.')}"
            _export(nt, str(name))
            ok += 1
        except Exception as exc:
            print(f"  ({n},{m}) FAILED: {exc}", file=sys.stderr)
            fail += 1

    print(f"Done: {ok} OK, {fail} failed.")


def cmd_polar(args):
    """Print or export the chirality map as CSV."""
    from core.chirality import scan_chirality

    structure = _snap(_load(args.input))
    results   = scan_chirality(structure, n_max=args.nmax,
                               max_diameter=args.dmax, unique_only=True)

    rows = [["n", "m", "diameter_A", "theta_deg", "strain_pct", "n_atoms"]]
    for r in results:
        rows.append([r.n, r.m, f"{r.diameter:.4f}", f"{r.theta_deg:.4f}",
                     f"{r.strain:.6f}", r.n_atoms])

    if args.output:
        with open(args.output, "w") as f:
            for row in rows:
                f.write(",".join(str(x) for x in row) + "\n")
        print(f"Chirality map → {args.output}  ({len(results)} points)")
    else:
        for row in rows:
            print(",".join(str(x) for x in row))


def cmd_prim(args):
    """Find and report the primitive cell of a structure."""
    from core.symmetry import find_primitive_cell

    structure = _load(args.input)
    prim, desc = find_primitive_cell(structure)
    print(desc)
    print(f"  Original : {len(structure.atoms)} atoms  a={structure.a:.4f} Å  "
          f"b={structure.b:.4f} Å  γ={structure.gamma_deg:.2f}°")
    print(f"  Primitive: {len(prim.atoms)} atoms  a={prim.a:.4f} Å  "
          f"b={prim.b:.4f} Å  γ={prim.gamma_deg:.2f}°  ({prim.lattice_type})")

    if args.output:
        from core.exporters import export
        # Write primitive cell as a flat XYZ
        import numpy as np
        from core.builder import NanotubeStructure
        from core.chirality import ChiralityResult
        atoms   = prim.atoms
        symbols = [a["symbol"] for a in atoms]
        coords  = np.array([[a["pos"][0], a["pos"][1], a.get("z", 0.0)] for a in atoms])
        dummy   = ChiralityResult(n=0, m=0, a1=prim.a1, a2=prim.a2)
        nt_stub = NanotubeStructure(
            chirality=dummy, symbols=symbols, coords=coords,
            box=np.array([prim.a, prim.b, 20.0]), vacuum=10.0,
        )
        export(nt_stub, args.output)
        print(f"  Primitive cell → {args.output}")


def cmd_query(args):
    """Query an online structural database."""
    from core.analysis import query_cod, query_mp

    print(f"Querying {args.source} for '{args.formula}' ...")

    if args.source.lower() == "cod":
        results = query_cod(args.formula, max_results=args.max)
    elif args.source.lower() == "mp":
        if not args.api_key:
            sys.exit("Materials Project requires --api-key.")
        results = query_mp(args.formula, args.api_key, max_results=args.max)
    else:
        sys.exit(f"Unknown source '{args.source}'. Use 'cod' or 'mp'.")

    if not results:
        print("No results found.")
        return

    for r in results:
        url = r.get("file_url") or r.get("mp_url", "")
        print(f"  [{r['id']}] {r['formula']}  sg={r['sg']}  "
              f"a={r.get('a','')} b={r.get('b','')}  {url}")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults → {args.output}")


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────────────────

def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ntbuilder",
        description="NTBuilder — headless nanotube generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    # ── build ─────────────────────────────────────────────────────────────────
    pb = sub.add_parser("build", help="Build a single nanotube")
    pb.add_argument("input",           help="Input structure file (CIF, POSCAR, XSF, …)")
    pb.add_argument("--n",    type=int, required=True,  help="Chiral index n")
    pb.add_argument("--m",    type=int, required=True,  help="Chiral index m")
    pb.add_argument("--output", "-o",  default="nanotube.xyz", help="Output file")
    pb.add_argument("--vacuum",        type=float, default=10.0, help="Vacuum padding (Å)")
    pb.add_argument("--repeat",        type=int,   default=1,    help="Axial replications")
    pb.add_argument("--roll-inward",   action="store_true",      help="Roll buckled structures inward")

    # ── mwnt ──────────────────────────────────────────────────────────────────
    pm = sub.add_parser("mwnt", help="Build a multi-walled nanotube")
    pm.add_argument("input")
    pm.add_argument("--n",       type=int,   required=True)
    pm.add_argument("--m",       type=int,   required=True)
    pm.add_argument("--walls",   type=int,   default=2,    help="Number of walls")
    pm.add_argument("--spacing", type=float, default=3.4,  help="Interlayer spacing (Å)")
    pm.add_argument("--vacuum",  type=float, default=10.0)
    pm.add_argument("--roll-inward", action="store_true")
    pm.add_argument("--output",  "-o", default="mwnt.xyz")

    # ── deform ────────────────────────────────────────────────────────────────
    pd = sub.add_parser("deform", help="Apply strain/torsion to an existing file")
    pd.add_argument("input")
    pd.add_argument("--strain",  type=float, default=0.0, help="Axial strain fraction (e.g. 0.05)")
    pd.add_argument("--twist",   type=float, default=0.0, help="Torsion rate (°/Å)")
    pd.add_argument("--output",  "-o", default="deformed.xyz")

    # ── bundle ────────────────────────────────────────────────────────────────
    pbu = sub.add_parser("bundle", help="Build a nanotube bundle")
    pbu.add_argument("input")
    pbu.add_argument("--geometry", default="hexagonal7",
                     choices=["linear","triangle","square4","hexagonal7","grid"])
    pbu.add_argument("--spacing",  type=float, default=3.4, help="Surface-to-surface gap (Å)")
    pbu.add_argument("--nx",       type=int,   default=2,   help="Grid columns (geometry=grid)")
    pbu.add_argument("--ny",       type=int,   default=2,   help="Grid rows (geometry=grid)")
    pbu.add_argument("--output",   "-o", default="bundle.xyz")

    # ── batch ─────────────────────────────────────────────────────────────────
    pba = sub.add_parser("batch", help="Batch-build a nanotube series")
    pba.add_argument("input")
    pba.add_argument("--type",   choices=["armchair","zigzag","all"], default="armchair")
    pba.add_argument("--nfrom",  type=int,   default=5,    help="n start (armchair/zigzag)")
    pba.add_argument("--nto",    type=int,   default=20,   help="n end")
    pba.add_argument("--nmax",   type=int,   default=60,   help="Max index for 'all'")
    pba.add_argument("--dmax",   type=float, default=30.0, help="Max diameter for 'all' (Å)")
    pba.add_argument("--vacuum", type=float, default=10.0)
    pba.add_argument("--fmt",    default="xyz",            help="Output format extension")
    pba.add_argument("--outdir", default="./batch_output", help="Output directory")

    # ── polar ─────────────────────────────────────────────────────────────────
    pp = sub.add_parser("polar", help="Export chirality map as CSV")
    pp.add_argument("input")
    pp.add_argument("--dmax",   type=float, default=25.0)
    pp.add_argument("--nmax",   type=int,   default=60)
    pp.add_argument("--output", "-o", default=None, help="CSV file (default: stdout)")

    # ── prim ──────────────────────────────────────────────────────────────────
    prm = sub.add_parser("prim", help="Find primitive cell")
    prm.add_argument("input")
    prm.add_argument("--output", "-o", default=None, help="Optional output file")

    # ── query ─────────────────────────────────────────────────────────────────
    pq = sub.add_parser("query", help="Query an online structure database")
    pq.add_argument("--formula",  required=True, help="Chemical formula, e.g. MoS2")
    pq.add_argument("--source",   default="cod", choices=["cod","mp"],
                    help="Database: 'cod' (free) or 'mp' (needs API key)")
    pq.add_argument("--api-key",  default="",    help="Materials Project API key (mp only)")
    pq.add_argument("--max",      type=int, default=20, help="Max results")
    pq.add_argument("--output",   "-o", default=None,   help="Save results as JSON")

    return p


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

_COMMANDS = {
    "build":  cmd_build,
    "mwnt":   cmd_mwnt,
    "deform": cmd_deform,
    "bundle": cmd_bundle,
    "batch":  cmd_batch,
    "polar":  cmd_polar,
    "prim":   cmd_prim,
    "query":  cmd_query,
}

def main(argv=None):
    parser = make_parser()
    args   = parser.parse_args(argv)

    # Normalise attribute names (argparse converts - to _)
    if hasattr(args, "roll_inward"):
        pass  # already snake_case
    if hasattr(args, "api_key") and hasattr(args, "api-key"):
        args.api_key = getattr(args, "api-key", "")

    fn = _COMMANDS.get(args.command)
    if fn is None:
        parser.print_help()
        sys.exit(1)

    try:
        fn(args)
    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        if "--debug" in sys.argv:
            raise
        sys.exit(1)


if __name__ == "__main__":
    main()
