"""
NTBuilder Web API — FastAPI backend
------------------------------------
Endpoints
    GET  /                            → serve SPA
    GET  /api/examples                → list bundled examples
    POST /api/upload                  → upload a structure file, returns file_id
    POST /api/polar                   → compute polar chirality map (JSON for Plotly)
    POST /api/build                   → build nanotube, returns job_id + metadata + XYZ
    GET  /api/xyz/{job_id}            → XYZ string for 3D viewer
    POST /api/primitive               → find primitive cell of uploaded structure
    POST /api/mwnt                    → build a scaled MWNT (k·(n,m) per wall)
    POST /api/bundle                  → build a nanotube bundle supercell
    POST /api/deform                  → apply axial strain / torsion / radial strain
    POST /api/analysis                → bond histogram, electronic, symmetry
    GET  /api/download/{job_id}/{fmt} → download output file
    POST /api/batch                   → batch-build a family of (n,m) tubes
    GET  /api/stats                   → usage statistics (JSON, for admins)

The Query Database endpoint is intentionally absent: the upstream
(C2DB / COD / Materials Project) APIs are unstable, and the feature is
currently listed as a planned extension in the manuscript.
"""
from __future__ import annotations

import datetime
import io
import json
import logging
import math
import shutil
import sys
import time
import uuid
import zipfile
from pathlib import Path
from typing import Optional

# Module-level logger so the server admin can see per-request phase
# timings without altering production code.  Defaults to INFO so the
# uvicorn output is reasonably terse.
_log = logging.getLogger("ntbuilder")
if not _log.handlers:
    _h = logging.StreamHandler(sys.stderr)
    _h.setFormatter(logging.Formatter("[%(asctime)s] ntbuilder %(message)s",
                                       datefmt="%H:%M:%S"))
    _log.addHandler(_h)
    _log.setLevel(logging.INFO)


def _timer():
    """Context manager-style stopwatch used by /api/build, /api/deform, …

    Usage:
        t = _timer()
        ...do work...
        t.lap("load")
        ...
        t.lap("write")
        t.report("deform N=29400")
    """
    class _T:
        def __init__(self):
            self.t0     = time.perf_counter()
            self.last   = self.t0
            self.phases = []
        def lap(self, name):
            now = time.perf_counter()
            self.phases.append((name, (now - self.last) * 1000.0))
            self.last = now
        def report(self, prefix):
            total = (time.perf_counter() - self.t0) * 1000.0
            parts = " ".join(f"{n}={ms:.0f}ms" for n, ms in self.phases)
            _log.info("%s total=%.0fms %s", prefix, total, parts)
    return _T()

import numpy as np
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

# ── Bootstrap Python path so 'core' is importable ────────────────────────────
_HERE = Path(__file__).parent          # web/api/
_ROOT = _HERE.parent.parent            # nanotube_builder/
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from core.io import load_structure, LatticeStructure                          # noqa: E402
from core.chirality import scan_chirality, compute_chirality, unique_sector_deg  # noqa: E402
from core.builder import build_nanotube, check_spurious_bonds                 # noqa: E402
from core.exporters import (                                                   # noqa: E402
    export, write_xyz, write_pdb, write_lammps, write_poscar, write_qe,
    write_xsf, write_cp2k, write_siesta, write_cif,
)
from core.symmetry import snap_to_symmetry, find_primitive_cell               # noqa: E402

from .models import (                                                          # noqa: E402
    BatchRequest, BuildRequest, PolarRequest,
    MWNTRequest, BundleRequest, DeformRequest, AnalysisRequest,
    MethodsRequest, DFTInputRequest,
)

# ── Paths ────────────────────────────────────────────────────────────────────
_WEB      = _HERE.parent                          # web/
_STATIC   = _WEB / "static"
_EXAMPLES = _ROOT / "examples"

# Job scratch directory.
#
# When the project lives inside a sync'd folder (OneDrive, Dropbox, iCloud
# Drive, …) every write under web/tmp/ is intercepted by the cloud agent and
# locked/uploaded immediately.  For a single deform call we used to write
# 9 file formats + a pickle per job — that means dozens of sync events per
# user interaction, which can make even small structures (~20 k atoms) take
# minutes.  We therefore default _TMP to a *local* temp directory; the user
# can opt back into the original `web/tmp/` by setting NTBUILDER_TMP=web.
import os, tempfile as _tempfile
_FORCE_WEBTMP = os.environ.get("NTBUILDER_TMP", "").lower() in ("web", "1", "yes")
_TMP = (_WEB / "tmp") if _FORCE_WEBTMP \
       else (Path(_tempfile.gettempdir()) / "ntbuilder")
_TMP.mkdir(parents=True, exist_ok=True)
# Logged at module import time so the user can confirm in the uvicorn console
# that the scratch directory is outside any sync'd folder (OneDrive, Dropbox,
# iCloud Drive) — otherwise file-system events bottleneck every operation.
_log.info("NTBuilder scratch dir: %s", _TMP)
if any(k in str(_TMP) for k in ("OneDrive", "iCloud", "Dropbox")):
    _log.warning("scratch dir is inside a cloud-sync folder — performance "
                 "will degrade.  Unset NTBUILDER_TMP or set it to '' to use "
                 "the system temp directory instead.")


def _ensure_tmp() -> Path:
    """Guarantee the scratch directory exists and return it.

    ``_TMP`` is created once at import time, but on a long-running server
    that defaults to the system temp dir a periodic /tmp cleaner
    (e.g. systemd-tmpfiles) can delete it after a few days, silently
    pulling the rug out from under every write.  Each endpoint that
    persists a job calls this first so a vanished scratch dir is
    transparently recreated instead of surfacing to the user as an
    HTTP 500 ("Build failed").
    """
    _TMP.mkdir(parents=True, exist_ok=True)
    return _TMP


def _new_job_dir() -> Path:
    """Create and return a fresh ``_TMP/<uuid>`` job directory.

    Recreates ``_TMP`` first (see :func:`_ensure_tmp`) so the server is
    self-healing if the scratch dir was cleaned away underneath it.
    """
    job_dir = _ensure_tmp() / str(uuid.uuid4())
    job_dir.mkdir()
    return job_dir

# ── Accepted upload extensions ────────────────────────────────────────────────
_ACCEPTED_EXTS = {
    ".cif", ".pdb", ".xyz",
    ".poscar", ".contcar", ".vasp",
    ".xsf",
    ".lammps", ".data",
    ".in", ".pwi",
}

# ── Anonymous stats log ───────────────────────────────────────────────────────
# The stats file is intentionally kept persistent (``web/stats.jsonl``) so
# the /api/stats endpoint can aggregate usage across server restarts.
#
# ``_log_stat`` is fire-and-forget telemetry — it must never block a request.
# The synchronous ``open("a")`` used previously could stall the FastAPI
# worker whenever the underlying file system stalled: typical culprits are
# a cloud-sync agent (OneDrive / Dropbox / iCloud Drive) holding a transient
# lock while uploading, antivirus on-access scans on Windows, or a paused
# network mount on a remote deployment.  We dispatch the write to a daemon
# thread so the request returns immediately regardless of what happens to
# the underlying file handle; exceptions inside the thread are swallowed,
# matching the previous best-effort contract.  ``NTBUILDER_STATS`` lets a
# server admin redirect the log to a path of their choice without editing
# code.
import threading as _threading

_stats_override = os.environ.get("NTBUILDER_STATS", "").strip()
_STATS_FILE = Path(_stats_override) if _stats_override else (_WEB / "stats.jsonl")


def _log_stat(event: str, **kwargs) -> None:
    """Append one JSON line to the stats log (fire-and-forget, never raises)."""
    record = {
        "ts":    datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "event": event,
        **kwargs,
    }
    payload = json.dumps(record) + "\n"

    def _flush() -> None:
        try:
            with _STATS_FILE.open("a", encoding="utf-8") as f:
                f.write(payload)
        except Exception:
            pass

    _threading.Thread(target=_flush, daemon=True).start()

# ── App ──────────────────────────────────────────────────────────────────────
app = FastAPI(
    title="NTBuilder Web API",
    version="1.1.0",
    description="REST API for nanotube generation from arbitrary 2D lattices.",
)


# ── Static files + SPA root ──────────────────────────────────────────────────
app.mount("/static", StaticFiles(directory=str(_STATIC)), name="static")


@app.get("/", include_in_schema=False)
async def root():
    return FileResponse(str(_STATIC / "index.html"))


# ── Examples ─────────────────────────────────────────────────────────────────
#
# The bundled ``examples/`` directory ships the five canonical lattices
# at its top level (Graphene, Biphenylene Network, Penta-Graphene, MoS₂,
# MoSSe); ancillary / research lattices (irida-graphene,
# rectangular-graphene, the oblique carbon allotropes, …) live in
# ``examples/extra/`` so they remain accessible to power users without
# cluttering the default UI listing.  ``/api/examples`` therefore enums
# the *root* directory by default and recurses into ``extra/`` only when
# the caller opts in with ``?all=true``.
@app.get("/api/examples")
async def list_examples(all: bool = False):
    """Return bundled structure examples as {name, filename} list.

    With ``all=false`` (default) only the canonical lattices at the
    top level of ``examples/`` are returned; with ``all=true`` the
    contents of ``examples/extra/`` are appended as well.  The
    ``filename`` field is path-relative to ``examples/`` so callers
    can pass it back unchanged as the ``example`` parameter to
    ``/api/build`` etc.
    """
    files: list[Path] = []
    for ext in _ACCEPTED_EXTS:
        files.extend(_EXAMPLES.glob(f"*{ext}"))
        if all:
            files.extend(_EXAMPLES.glob(f"extra/*{ext}"))
    files = sorted(set(files), key=lambda p: p.stem)
    return [
        {
            "name":     f.stem.replace("_", " "),
            "filename": str(f.relative_to(_EXAMPLES)).replace("\\", "/"),
        }
        for f in files
    ]


# ── Upload ───────────────────────────────────────────────────────────────────
@app.post("/api/upload")
async def upload_structure(file: UploadFile = File(...)):
    """Accept a structure file upload and return a temporary file_id."""
    fname = (file.filename or "").lower()
    ext = Path(fname).suffix
    # Accept files named POSCAR/CONTCAR without extension too
    if ext not in _ACCEPTED_EXTS and Path(fname).name not in ("poscar", "contcar"):
        raise HTTPException(
            400,
            f"Unsupported file type '{ext}'. Accepted: CIF, PDB, XYZ, POSCAR/CONTCAR, "
            "XSF, LAMMPS data (.lammps/.data), Quantum ESPRESSO (.in/.pwi)."
        )
    data = await file.read()
    if len(data) > 10 * 1024 * 1024:   # 10 MB hard limit
        raise HTTPException(413, "File too large (max 10 MB).")
    file_id = str(uuid.uuid4())
    # Store with original extension so load_structure can auto-detect format
    save_ext = ext if ext else ""
    (_ensure_tmp() / f"{file_id}{save_ext}").write_bytes(data)
    return {"file_id": file_id, "filename": file.filename}


# ── Polar map ─────────────────────────────────────────────────────────────────
@app.post("/api/polar")
async def polar_map(req: PolarRequest):
    """Compute all valid (n,m) chiralities and return Plotly-ready scatter data."""
    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc

    # Auto snap to symmetry
    try:
        structure, snap_desc = snap_to_symmetry(structure)
    except Exception:
        snap_desc = ""

    theta_max = unique_sector_deg(structure)

    results = scan_chirality(
        structure,
        n_max=req.n_max,
        max_diameter=req.max_diameter,
        unique_only=True,
    )

    # Detect curvature-induced spurious bonds for each chirality.  The
    # decision is *informational only* — the user remains free to build
    # any (n, m) they want (see /api/build); the polar map merely renders
    # an X marker on those points so the operator can spot incompatible
    # chiralities at a glance without having to build them first.  For a
    # flat (non-buckled) structure there can be no curvature-induced
    # spurious bonds, so we skip the per-point construction entirely.
    needs_spurious_check = structure.has_buckling
    if needs_spurious_check:
        # Local imports keep the cold-start cost of /api/polar low for the
        # common flat-structure case (graphene, biphenylene, …).
        from core.builder import build_nanotube as _bn
        from core.builder import check_spurious_bonds as _csb

    points: list[dict] = []
    for r in results:
        # Server-side strain filter
        if req.strain_max is not None and r.strain > req.strain_max:
            continue
        theta_rad = math.radians(r.theta_deg)
        x = r.diameter * math.cos(theta_rad)
        y = r.diameter * math.sin(theta_rad)

        point = {
            "n":          r.n,
            "m":          r.m,
            "diameter":   round(r.diameter,   4),
            "theta_deg":  round(r.theta_deg,  4),
            "strain":     round(r.strain,     6),
            "n_atoms":    r.n_atoms,
            "x":          round(x, 4),
            "y":          round(y, 4),
        }
        # Only attach the ``spurious`` field for structures where it is
        # physically meaningful — flat (non-buckled) lattices cannot
        # develop curvature-induced bonds, so omitting the field there
        # tells the frontend to suppress the legend entirely.
        if needs_spurious_check:
            spurious_pairs: list[str] = []
            try:
                _nt = _bn(structure, r, vacuum=0.0,
                          roll_inward=req.roll_inward)
                sp  = _csb(structure, _nt)
                spurious_pairs = sorted("-".join(sorted(p)) for p in sp)
            except Exception:
                # Construction failures (e.g. (n,m) collapsed by snap) are
                # not curvature issues; leave the marker as a regular dot.
                spurious_pairs = []
            point["spurious"] = spurious_pairs
        points.append(point)

    return {
        "points":        points,
        "dmax":          req.max_diameter,
        "theta_max":     theta_max,
        "lattice_type":  structure.lattice_type,
        "a":             round(structure.a, 4),
        "b":             round(structure.b, 4),
        "gamma_deg":     round(structure.gamma_deg, 2),
        "n_species":     len({a["symbol"] for a in structure.atoms}),
        "species":       list({a["symbol"] for a in structure.atoms}),
        "d_min":         round(structure.d_min, 4),
        "snap_desc":     snap_desc,
    }


# ── Build ─────────────────────────────────────────────────────────────────────
@app.post("/api/build")
async def build(req: BuildRequest):
    """Build a nanotube, save all output formats, return metadata + XYZ preview."""
    timer = _timer()
    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc

    # Auto snap to symmetry
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    chirality = compute_chirality(req.n, req.m, structure)
    if chirality is None:
        raise HTTPException(400, "Invalid chirality: n=m=0 is degenerate.")

    # The legacy "diameter < structure.d_min" rejection has been removed:
    # ``d_min`` was simply ``2·|z_max|`` (the trivial geometric lower bound
    # so that the innermost shell has positive radius after rolling), which
    # is far too permissive for buckled / Janus structures — e.g. MoSSe
    # (11,11) easily clears the 3.26 Å threshold yet still develops
    # curvature-induced Se–Se bonds.  A scientifically motivated user may
    # also want to study sub-threshold tubes deliberately, so we no longer
    # block the build at this layer.  Curvature-induced spurious bonds are
    # detected post-build below and surfaced as a warning, and the polar
    # map renders an X marker on chiralities that would trigger them so
    # the choice is informed from the start.

    build_kwargs: dict = dict(vacuum=req.vacuum, roll_inward=req.roll_inward)
    if req.bond_cutoffs:
        build_kwargs["bond_cutoffs"] = req.bond_cutoffs

    try:
        nt = build_nanotube(structure, chirality, **build_kwargs)
    except Exception as exc:
        raise HTTPException(500, f"Build failed: {exc}") from exc

    # Axial replication
    if req.n_repeat > 1:
        nt = _replicate_z(nt, req.n_repeat)

    # Spurious-bond check (requires scipy; silently skip if not installed)
    warning: Optional[str] = None
    try:
        spurious = check_spurious_bonds(structure, nt)
        if spurious:
            pairs = ", ".join("-".join(sorted(p)) for p in spurious)
            warning = (
                f"Curvature-induced spurious bonds detected ({pairs}). "
                "Consider using a larger diameter (higher n,m)."
            )
    except ImportError:
        pass  # scipy not available — skip check silently

    # Save all output formats
    job_dir = _new_job_dir()
    job_id  = job_dir.name
    _write_all(job_dir, nt)

    # Metadata JSON
    meta = {
        "n":          req.n,
        "m":          req.m,
        "diameter":   round(chirality.diameter, 4),
        "length":     round(float(nt.box[2]), 4),
        "n_atoms":    nt.n_atoms,
        "strain":     round(chirality.strain, 6),
        "theta_deg":  round(chirality.theta_deg, 2),
        "t1":         chirality.t1,
        "t2":         chirality.t2,
        "n_repeat":   req.n_repeat,
        "vacuum":     req.vacuum,
        "roll_inward": req.roll_inward,
        "box":        [round(float(x), 4) for x in nt.box],  # [Lx, Ly, Lz]
    }
    (job_dir / "meta.json").write_text(json.dumps(meta, indent=2))

    # XYZ string for 3D viewer (inline in response)
    xyz_str = (job_dir / "nanotube.xyz").read_text()

    timer.lap("write")
    _log_stat(
        "build",
        n=req.n, m=req.m,
        diameter=meta["diameter"],
        n_atoms=meta["n_atoms"],
        n_repeat=req.n_repeat,
        vacuum=req.vacuum,
    )
    timer.report(f"build ({req.n},{req.m}) reps={req.n_repeat} N={nt.n_atoms}")

    return {
        "job_id":    job_id,
        "warning":   warning,
        "xyz":       xyz_str,
        **meta,
    }


# ── XYZ for 3D viewer ────────────────────────────────────────────────────────
@app.get("/api/xyz/{job_id}")
async def get_xyz(job_id: str):
    """Return the XYZ file for the 3D viewer."""
    path = _TMP / _safe_id(job_id) / "nanotube.xyz"
    if not path.exists():
        raise HTTPException(404, "Job not found or expired.")
    return PlainTextResponse(path.read_text())


# ── Primitive cell ────────────────────────────────────────────────────────────
def _write_lattice_xyz(struct: LatticeStructure, path: Path) -> None:
    """Serialise a 2D LatticeStructure as an extended-XYZ file.

    The lattice vectors are embedded in the comment line as
    ``Lattice="a1x a1y 0 a2x a2y 0 0 0 100"`` so that
    :func:`core.io.read_xyz` can round-trip the structure.
    """
    a1 = struct.a1
    a2 = struct.a2
    lattice_str = (
        f'{a1[0]:.10f} {a1[1]:.10f} 0.0 '
        f'{a2[0]:.10f} {a2[1]:.10f} 0.0 '
        f'0.0 0.0 100.0'
    )
    lines = [
        str(len(struct.atoms)),
        f'Lattice="{lattice_str}" Properties=species:S:1:pos:R:3',
    ]
    for atom in struct.atoms:
        sym = atom["symbol"]
        x, y = float(atom["pos"][0]), float(atom["pos"][1])
        z    = float(atom.get("z", 0.0))
        lines.append(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@app.post("/api/primitive")
async def primitive(file_id: Optional[str] = None, example: Optional[str] = None):
    """Find the primitive cell of an uploaded structure and persist it as
    a new file so subsequent ``/api/polar``, ``/api/build``, … calls use
    the primitive lattice.  Returns the new ``file_id`` along with the
    diagnostic metadata."""
    struct_path = _resolve_struct(file_id, example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc

    try:
        prim, desc = find_primitive_cell(structure)
    except Exception as exc:
        raise HTTPException(500, f"find_primitive_cell failed: {exc}") from exc

    # Persist the primitive lattice as a new uploaded file so the rest
    # of the workflow can address it.
    new_id   = str(uuid.uuid4())
    new_path = _ensure_tmp() / f"{new_id}.xyz"
    try:
        _write_lattice_xyz(prim, new_path)
    except Exception as exc:
        raise HTTPException(500,
            f"Could not persist primitive cell: {exc}") from exc

    return {
        "file_id":     new_id,        # ← frontend should switch to this id
        "description": desc,
        "n_atoms_orig": len(structure.atoms),
        "n_atoms_prim": len(prim.atoms),
        "a_orig":  round(structure.a, 4),
        "b_orig":  round(structure.b, 4),
        "a_prim":  round(prim.a, 4),
        "b_prim":  round(prim.b, 4),
        "lattice_type": prim.lattice_type,
        "gamma_deg": round(prim.gamma_deg, 2),
        "species":   list({a["symbol"] for a in prim.atoms}),
    }


# ── Download ──────────────────────────────────────────────────────────────────
_FMT_MAP = {
    "xyz":    ("nanotube.xyz",    "text/plain"),
    "pdb":    ("nanotube.pdb",    "text/plain"),
    "lammps": ("nanotube.lammps", "text/plain"),
    "poscar": ("POSCAR",          "text/plain"),
    "qe":     ("nanotube.pwi",    "text/plain"),
    "xsf":    ("nanotube.xsf",    "text/plain"),
    "cp2k":   ("nanotube.inp",    "text/plain"),
    "siesta": ("nanotube.fdf",    "text/plain"),
    "cif":    ("nanotube.cif",    "text/plain"),
}


@app.get("/api/download/{job_id}/{fmt}")
async def download(job_id: str, fmt: str):
    """Serve an output file.

    Files other than the (eagerly-written) XYZ are generated lazily from the
    persisted pickle on first access — this keeps the build/deform/bundle
    response time independent of how many output formats the user actually
    needs.
    """
    if fmt not in _FMT_MAP:
        raise HTTPException(
            400, f"Unknown format '{fmt}'. Choose from: {', '.join(_FMT_MAP)}."
        )
    job_dir = _TMP / _safe_id(job_id)
    if not job_dir.exists():
        raise HTTPException(404, "Job not found or expired (files are kept 24 h).")

    filename, media_type = _FMT_MAP[fmt]
    file_path = job_dir / filename

    if not file_path.exists():
        # Lazy materialisation: load pickle, write the requested format.
        pkl = job_dir / "nanotube.pkl"
        if not pkl.exists():
            raise HTTPException(500, "Job pickle missing — please rebuild.")
        try:
            import pickle
            with pkl.open("rb") as f:
                nt = pickle.load(f)
            _write_one(fmt, nt, file_path)
        except Exception as exc:
            raise HTTPException(500, f"Could not generate {fmt}: {exc}") from exc

    _log_stat("download", fmt=fmt, job_id=job_id)

    return FileResponse(
        str(file_path),
        filename=filename,
        media_type=media_type,
    )


# ── Health & version (diagnostics) ────────────────────────────────────────────
@app.get("/api/health")
async def health():
    """Return server diagnostics.  Use this to verify that the running
    instance has the latest scratch-dir configuration (i.e., that ``_TMP``
    points outside a sync'd folder) — critical for deform-on-bundle
    performance.  Also reports the existing job count.
    """
    try:
        n_jobs = sum(1 for _ in _ensure_tmp().iterdir() if _.is_dir())
    except Exception:
        n_jobs = -1
    in_onedrive = "OneDrive" in str(_TMP) or "iCloud" in str(_TMP) \
                  or "Dropbox" in str(_TMP)
    return {
        "version":      app.version,
        "tmp":          str(_TMP),
        "tmp_in_cloud": in_onedrive,
        "jobs":         n_jobs,
        "max_deform":   _MAX_DEFORM_ATOMS,
        "examples":     str(_EXAMPLES),
    }


# ── Stats viewer ──────────────────────────────────────────────────────────────
@app.get("/api/stats")
async def get_stats():
    """Return aggregated usage stats from the local log (admin use)."""
    if not _STATS_FILE.exists():
        return {"total_builds": 0, "total_downloads": 0, "by_format": {}, "recent": []}

    builds, downloads, by_fmt = 0, 0, {}
    recent = []
    for line in _STATS_FILE.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
        except Exception:
            continue
        if rec.get("event") == "build":
            builds += 1
        elif rec.get("event") == "download":
            downloads += 1
            fmt = rec.get("fmt", "?")
            by_fmt[fmt] = by_fmt.get(fmt, 0) + 1
        recent.append(rec)

    return {
        "total_builds":    builds,
        "total_downloads": downloads,
        "by_format":       by_fmt,
        "recent":          recent[-20:],   # last 20 events
    }


# ── Batch build ───────────────────────────────────────────────────────────────
@app.post("/api/batch")
async def batch_build(req: BatchRequest):
    """Build multiple nanotubes and return a ZIP archive."""
    if len(req.chiralities) > 100:
        raise HTTPException(400, "Maximum 100 chiralities per batch.")
    if not req.chiralities:
        raise HTTPException(400, "No chiralities specified.")

    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc

    # Auto snap to symmetry
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    results = []
    errors  = []

    for pair in req.chiralities:
        try:
            n, m = int(pair[0]), int(pair[1])
        except Exception:
            errors.append(f"Invalid pair {pair}")
            continue
        try:
            chirality = compute_chirality(n, m, structure)
            if chirality is None:
                errors.append(f"({n},{m}): degenerate (n=m=0)")
                continue
            # The legacy ``D < d_min`` filter was dropped here for the same
            # reason as in /api/build: ``d_min`` is the trivial geometric
            # bound (``2·|z_max|``), not a real bond-spuriousness criterion,
            # and the user may legitimately want to include below-threshold
            # tubes in a batch (e.g. to map curvature-induced bonding
            # systematically).  Spurious bonds, when present, are reported
            # per-tube and the user can decide whether to keep or drop them.
            nt = build_nanotube(structure, chirality,
                                vacuum=req.vacuum, roll_inward=req.roll_inward)
            if req.n_repeat > 1:
                nt = _replicate_z(nt, req.n_repeat)
            results.append((n, m, nt))
        except Exception as exc:
            errors.append(f"({n},{m}): {exc}")

    if not results:
        raise HTTPException(400, f"No valid nanotubes built. Errors: {'; '.join(errors)}")

    # Package into ZIP.  Write each format directly into the zip stream to
    # avoid hitting OneDrive-watched scratch paths and to keep memory usage
    # bounded.
    _BATCH_FMTS = [
        ("xyz",    "nanotube.xyz"),    ("pdb",    "nanotube.pdb"),
        ("lammps", "nanotube.lammps"), ("poscar", "POSCAR"),
        ("qe",     "nanotube.pwi"),    ("xsf",    "nanotube.xsf"),
        ("cp2k",   "nanotube.inp"),    ("siesta", "nanotube.fdf"),
        ("cif",    "nanotube.cif"),
    ]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for n, m, nt in results:
            folder = f"nt_{n}_{m}"
            tmp = _ensure_tmp() / str(uuid.uuid4())
            tmp.mkdir()
            try:
                for fmt, fname in _BATCH_FMTS:
                    fp = tmp / fname
                    _write_one(fmt, nt, fp)
                    if fp.exists():
                        zf.write(fp, f"{folder}/{fname}")
                meta = {
                    "n": n, "m": m,
                    "diameter": round(nt.chirality.diameter, 4),
                    "length":   round(float(nt.box[2]), 4),
                    "n_atoms":  nt.n_atoms,
                    "strain":   round(nt.chirality.strain, 6),
                }
                zf.writestr(f"{folder}/meta.json", json.dumps(meta, indent=2))
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

        if errors:
            zf.writestr("errors.txt", "\n".join(errors))

    buf.seek(0)
    return StreamingResponse(
        buf,
        media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="nanotubes.zip"'},
    )


# ── Scaled MWNT ───────────────────────────────────────────────────────────────
@app.post("/api/mwnt")
async def mwnt(req: MWNTRequest):
    """Build a multi-walled nanotube via integer scaling of (n,m).

    Returns the wall plan plus a job_id so the merged structure can be
    downloaded via the standard /api/download/{job_id}/{fmt} endpoint.
    """
    from core.mwnt import plan_scaled_walls, build_mwnt_scaled

    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    inner = compute_chirality(req.n, req.m, structure)
    if inner is None:
        raise HTTPException(400, "(n=0, m=0) is degenerate.")

    plans = plan_scaled_walls(inner, req.n_walls,
                              interlayer_spacing=req.interlayer_spacing)
    try:
        result = build_mwnt_scaled(
            structure, inner,
            n_walls            = req.n_walls,
            interlayer_spacing = req.interlayer_spacing,
            vacuum             = req.vacuum,
            roll_inward        = req.roll_inward,
        )
    except Exception as exc:
        raise HTTPException(500, f"MWNT build failed: {exc}") from exc

    job_dir = _new_job_dir()
    job_id  = job_dir.name
    _t = _timer()
    _write_all(job_dir, result.nanotube)
    _t.lap("write")
    _t.report(f"mwnt ({req.n},{req.m}) walls={req.n_walls} N={result.nanotube.n_atoms}")
    _log_stat("mwnt_build", n=req.n, m=req.m, walls=req.n_walls,
              atoms=result.nanotube.n_atoms)

    return {
        "job_id": job_id,
        "walls": [
            {"index": p.index, "k": p.k, "n": p.n, "m": p.m,
             "diameter": round(p.diameter, 4),
             "target_diameter": round(p.target_diameter, 4),
             "actual_spacing": None if p.index == 0
                               else round(p.actual_spacing, 4)}
            for p in plans
        ],
        "requested_spacing": req.interlayer_spacing,
        "mean_spacing":      round(result.mean_spacing, 4),
        "n_atoms":           result.nanotube.n_atoms,
        "box":               [float(x) for x in result.nanotube.box],
    }


# ── Nanotube bundle ───────────────────────────────────────────────────────────
@app.post("/api/bundle")
async def bundle(req: BundleRequest):
    """Replicate a nanotube into a periodic bundle supercell.

    When ``from_job_id`` is provided, the bundle is built on top of the
    previously generated structure (e.g. an MWNT or a deformed tube) —
    this is the path used by the web UI when the user clicks "Bundle"
    after another operation.  Otherwise a fresh SWNT is built from
    ``file_id`` + ``(n, m)``.
    """
    from core.bundles import build_bundle

    if req.from_job_id:
        nt = _load_nanotube_from_job(req.from_job_id)
        if req.n_repeat > 1:
            nt = _replicate_z(nt, req.n_repeat)
    else:
        struct_path = _resolve_struct(req.file_id, req.example)
        try:
            structure = load_structure(str(struct_path))
        except Exception as exc:
            raise HTTPException(400, f"Could not read structure file: {exc}") from exc
        try:
            structure, _ = snap_to_symmetry(structure)
        except Exception:
            pass

        chirality = compute_chirality(req.n, req.m, structure)
        if chirality is None:
            raise HTTPException(400, "(n=0, m=0) is degenerate.")
        nt = build_nanotube(structure, chirality, vacuum=req.vacuum)
        if req.n_repeat > 1:
            nt = _replicate_z(nt, req.n_repeat)

    try:
        result = build_bundle(
            nt,
            geometry = req.geometry,
            spacing  = req.spacing,
            vacuum   = req.vacuum,
            nx       = req.nx,
            ny       = req.ny,
        )
    except Exception as exc:
        raise HTTPException(400, f"Bundle build failed: {exc}") from exc

    job_dir = _new_job_dir()
    job_id  = job_dir.name
    _t = _timer()
    _write_all(job_dir, result.nanotube)
    _t.lap("write")
    _t.report(f"bundle {req.geometry} ({req.n},{req.m}) "
              f"N={result.nanotube.n_atoms}")
    _log_stat("bundle_build", n=req.n, m=req.m, geometry=req.geometry,
              tubes=result.n_tubes, atoms=result.nanotube.n_atoms)

    return {
        "job_id":   job_id,
        "geometry": result.geometry,
        "n_tubes":  result.n_tubes,
        "pitch":    round(result.pitch, 4),
        "spacing":  round(result.spacing, 4),
        "vacuum":   req.vacuum,
        "n_atoms":  result.nanotube.n_atoms,
        "box":      [float(x) for x in result.nanotube.box],
    }


# ── Deformations ──────────────────────────────────────────────────────────────
# Hard cap: writing 9 file formats for very large structures (>~200 k atoms)
# becomes the dominant cost — a bundle of 7 tubes replicated 10× axially can
# easily hit this for low (n,m).  We bail out before the heavy phase so the
# request returns a clean error instead of hanging the worker.
_MAX_DEFORM_ATOMS = 250_000


@app.post("/api/deform")
async def deform(req: DeformRequest):
    """Apply axial strain, torsion and / or radial strain to a nanotube."""
    from core.deformations import (
        apply_axial_strain, apply_torsion, apply_radial_strain,
        deformation_description, torsion_warning,
    )

    timer = _timer()

    if req.from_job_id:
        nt = _load_nanotube_from_job(req.from_job_id)
        timer.lap("load")
        # Pre-check the post-replication atom count so we can refuse cheaply
        # instead of crashing inside the exporters.
        projected = int(nt.n_atoms) * max(1, int(req.n_repeat))
        if projected > _MAX_DEFORM_ATOMS:
            raise HTTPException(
                413,
                f"Estrutura projetada com {projected:,} átomos "
                f"({nt.n_atoms:,} × {req.n_repeat} repetições) excede o "
                f"limite de {_MAX_DEFORM_ATOMS:,} para /api/deform.  "
                "Reduza as repetições axiais ou aplique a torção antes do bundle."
            )
        if req.n_repeat > 1:
            nt = _replicate_z(nt, req.n_repeat)
            timer.lap("replicate")
    else:
        struct_path = _resolve_struct(req.file_id, req.example)
        try:
            structure = load_structure(str(struct_path))
        except Exception as exc:
            raise HTTPException(400, f"Could not read structure file: {exc}") from exc
        try:
            structure, _ = snap_to_symmetry(structure)
        except Exception:
            pass

        chirality = compute_chirality(req.n, req.m, structure)
        if chirality is None:
            raise HTTPException(400, "(n=0, m=0) is degenerate.")
        nt = build_nanotube(structure, chirality, vacuum=req.vacuum)
        if req.n_repeat > 1:
            nt = _replicate_z(nt, req.n_repeat)

    try:
        if abs(req.axial_strain) > 1e-9:
            nt = apply_axial_strain(nt, req.axial_strain)
        if abs(req.twist_rate) > 1e-9:
            nt = apply_torsion(nt, req.twist_rate, z_vacuum=req.z_vacuum)
        if abs(req.radial_strain) > 1e-9:
            nt = apply_radial_strain(nt, req.radial_strain)
    except Exception as exc:
        raise HTTPException(400, f"Deformation failed: {exc}") from exc
    timer.lap("deform")

    job_dir = _new_job_dir()
    job_id  = job_dir.name
    _write_all(job_dir, nt)
    timer.lap("write")
    _log_stat("deform_build", n=req.n, m=req.m,
              strain=req.axial_strain, twist=req.twist_rate,
              radial=req.radial_strain, atoms=nt.n_atoms)
    timer.report(f"deform N={nt.n_atoms} twist={req.twist_rate:+.3f}°/Å reps={req.n_repeat}")

    return {
        "job_id":      job_id,
        "description": deformation_description(
            axial_strain=req.axial_strain,
            twist_rate=req.twist_rate,
            radial_strain=req.radial_strain,
        ),
        "warning":     torsion_warning(req.twist_rate, req.z_vacuum),
        "n_atoms":     nt.n_atoms,
        "box":         [float(x) for x in nt.box],
    }


# ── Structural analysis ───────────────────────────────────────────────────────
@app.post("/api/analysis")
async def analysis(req: AnalysisRequest):
    """Return bond-length statistics, electronic character (zone-folding)
    and tube-symmetry information for a nanotube."""
    from core.analysis import (
        bond_analysis, electronic_character_label, tube_symmetry_info,
    )

    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    chirality = compute_chirality(req.n, req.m, structure)
    if chirality is None:
        raise HTTPException(400, "(n=0, m=0) is degenerate.")
    nt = build_nanotube(structure, chirality, vacuum=req.vacuum)

    ba   = bond_analysis(nt, cutoff=req.bond_cutoff)
    species = list(getattr(nt, "symbols", []))
    elec = electronic_character_label(
        req.n, req.m, structure.lattice_type, species=species,
    )
    sym  = tube_symmetry_info(req.n, req.m)
    return {
        "bond_analysis": {
            "n_bonds":      int(ba["n_bonds"]),
            "mean":         float(ba["mean"]),
            "std":          float(ba["std"]),
            "min":          float(ba["min"]),
            "max":          float(ba["max"]),
            # distances and per-pair labels are returned as plain lists
            # (FastAPI/Pydantic handles JSON serialisation transparently)
            "distances":    [float(d) for d in ba["distances"]],
            "pairs":        [str(p) for p in ba["pairs"]],
            "species":      sorted(ba["species_set"]),
        },
        "electronic_label": elec,
        "symmetry":         sym,
    }


# ── Methods text ──────────────────────────────────────────────────────────────
@app.post("/api/methods")
async def methods(req: MethodsRequest):
    """Return a ready-to-paste Methods-section paragraph for the nanotube.

    When ``from_job_id`` is provided the routine re-uses the persisted
    NanotubeStructure (so any axial / radial / torsional deformations are
    preserved in the box-size sentences), otherwise the SWNT is rebuilt
    from the parent lattice + ``(n, m)``.
    """
    from core.analysis import generate_methods_text

    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    if req.from_job_id:
        nt = _load_nanotube_from_job(req.from_job_id)
    else:
        chirality = compute_chirality(req.n, req.m, structure)
        if chirality is None:
            raise HTTPException(400, "(n=0, m=0) is degenerate.")
        nt = build_nanotube(structure, chirality, vacuum=req.vacuum)

    try:
        text = generate_methods_text(
            nt,
            structure   = structure,
            deform_desc = req.deform_desc or "",
            software    = req.software,
            version     = req.version,
            cite_key    = req.cite_key,
            n_walls     = req.n_walls,
            wall_info   = req.wall_info or "",
        )
    except Exception as exc:
        raise HTTPException(500, f"Methods text generation failed: {exc}") from exc

    return {"text": text}


# ── DFT input files ───────────────────────────────────────────────────────────
@app.post("/api/dft_inputs")
async def dft_inputs(req: DFTInputRequest):
    """Generate DFT input files (VASP / QE / CP2K / SIESTA) for the nanotube.

    The structure is rebuilt from ``(n, m)`` unless ``from_job_id`` is
    supplied — in that case the current persisted NanotubeStructure is
    used, preserving any deformations or multi-wall geometry.
    """
    from core.analysis import (
        generate_vasp_inputs, generate_qe_input, generate_cp2k_input,
    )
    from core.exporters import write_siesta

    struct_path = _resolve_struct(req.file_id, req.example)
    try:
        structure = load_structure(str(struct_path))
    except Exception as exc:
        raise HTTPException(400, f"Could not read structure file: {exc}") from exc
    try:
        structure, _ = snap_to_symmetry(structure)
    except Exception:
        pass

    if req.from_job_id:
        nt = _load_nanotube_from_job(req.from_job_id)
    else:
        chirality = compute_chirality(req.n, req.m, structure)
        if chirality is None:
            raise HTTPException(400, "(n=0, m=0) is degenerate.")
        nt = build_nanotube(structure, chirality, vacuum=req.vacuum)

    code = (req.code or "vasp").strip().lower()
    try:
        if code == "vasp":
            from core.exporters import write_poscar
            incar, kpoints = generate_vasp_inputs(nt)
            tmp_dir = _ensure_tmp() / f"_inp_{uuid.uuid4()}"
            tmp_dir.mkdir()
            try:
                write_poscar(nt, tmp_dir / "POSCAR")
                poscar = (tmp_dir / "POSCAR").read_text()
            finally:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            return {"code": "vasp",
                    "files": {"INCAR": incar, "KPOINTS": kpoints, "POSCAR": poscar}}

        if code == "qe":
            text = generate_qe_input(nt, structure)
            return {"code": "qe", "files": {"nanotube.pwi": text}}

        if code == "cp2k":
            text = generate_cp2k_input(nt)
            return {"code": "cp2k", "files": {"nanotube.inp": text}}

        if code == "siesta":
            tmp_dir = _ensure_tmp() / f"_inp_{uuid.uuid4()}"
            tmp_dir.mkdir()
            try:
                write_siesta(nt, tmp_dir / "nanotube.fdf")
                text = (tmp_dir / "nanotube.fdf").read_text()
            finally:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            return {"code": "siesta", "files": {"nanotube.fdf": text}}

    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(500, f"DFT input generation failed: {exc}") from exc

    raise HTTPException(400, f"Unknown DFT code '{req.code}'. "
                              "Use vasp | qe | cp2k | siesta.")


# ── Internal helpers ──────────────────────────────────────────────────────────
def _resolve_struct(file_id: Optional[str], example: Optional[str]) -> Path:
    """Locate an uploaded structure file or a bundled example."""
    if file_id:
        safe = _safe_id(file_id)
        # Try every accepted extension (upload saves with original ext)
        candidates = list(_TMP.glob(f"{safe}.*")) + [_TMP / safe]
        for p in candidates:
            if p.exists():
                return p
        raise HTTPException(404, "Uploaded file not found. Please re-upload.")
    if example:
        # Resolve relative to ``_EXAMPLES`` and re-anchor to its real path
        # so that a forged ``example=../../etc/passwd`` cannot escape the
        # bundled examples sandbox.  ``examples/extra/<file>.cif`` style
        # paths are accepted because they resolve cleanly underneath
        # ``_EXAMPLES``.
        candidate = (_EXAMPLES / example).resolve()
        base      = _EXAMPLES.resolve()
        try:
            candidate.relative_to(base)
        except ValueError:
            raise HTTPException(400, f"Invalid example path: {example!r}")
        if not candidate.exists():
            raise HTTPException(404, f"Example '{example}' not found.")
        return candidate
    raise HTTPException(400, "Provide either 'file_id' (uploaded) or 'example' (name).")


def _safe_id(value: str) -> str:
    """Strip any path separators to prevent directory traversal."""
    return Path(value).name


def _replicate_z(nt, n: int):
    from core.builder import NanotubeStructure
    Lz      = float(nt.box[2])
    offsets = [np.array([0.0, 0.0, i * Lz]) for i in range(n)]
    new_coords = np.vstack([nt.coords + off for off in offsets])
    return NanotubeStructure(
        chirality = nt.chirality,
        symbols   = list(nt.symbols) * n,
        coords    = new_coords,
        box       = np.array([nt.box[0], nt.box[1], Lz * n]),
        vacuum    = nt.vacuum,
    )


def _write_all(job_dir: Path, nt) -> None:
    """Persist a job to disk.

    Only the XYZ (needed by the 3D viewer) and the pickle (needed for
    operation chaining) are written eagerly.  The remaining 8 output
    formats are generated on demand by /api/download to keep build /
    deform / bundle latency low.

    Both writes use the fastest available binary protocol and a single
    ``write_bytes`` / ``write_text`` flush, so even on a sync'd folder
    the cloud agent fires only two events per operation.
    """
    write_xyz(nt, job_dir / "nanotube.xyz")
    import pickle
    pkl_path = job_dir / "nanotube.pkl"
    pkl_path.write_bytes(pickle.dumps(nt, protocol=pickle.HIGHEST_PROTOCOL))


# Map of download format → (filename, writer).  Used by the download endpoint
# to materialise files lazily from the persisted pickle.
def _write_one(fmt: str, nt, dest: Path) -> None:
    if   fmt == "xyz":    write_xyz(nt,    dest)
    elif fmt == "pdb":    write_pdb(nt,    dest)
    elif fmt == "lammps": write_lammps(nt, dest)
    elif fmt == "poscar": write_poscar(nt, dest)
    elif fmt == "qe":     write_qe(nt,     dest)
    elif fmt == "xsf":    write_xsf(nt,    dest)
    elif fmt == "cp2k":   write_cp2k(nt,   dest)
    elif fmt == "siesta": write_siesta(nt, dest)
    elif fmt == "cif":    write_cif(nt,    dest)
    else:
        raise ValueError(fmt)


def _load_nanotube_from_job(job_id: str):
    """Reload a previously built NanotubeStructure from its pickle blob.

    Used by /api/mwnt, /api/bundle, /api/deform when the request carries a
    ``from_job_id`` — that way the operation is applied to the current
    on-screen structure (e.g. an MWNT) instead of restarting from the
    (n, m) of the parent 2D lattice.
    """
    import pickle
    pkl = _TMP / _safe_id(job_id) / "nanotube.pkl"
    if not pkl.exists():
        raise HTTPException(
            404,
            "Previous job not found or expired — rebuild the nanotube first.",
        )
    with pkl.open("rb") as f:
        return pickle.load(f)
