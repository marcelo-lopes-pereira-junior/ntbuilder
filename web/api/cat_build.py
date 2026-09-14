"""
Construção dos tubos do catálogo — usada pelo site e pelo job do cluster.

Este módulo não importa FastAPI de propósito: o mesmo código roda no processo
do site (um tubo, na hora) e num nó do cluster com o python de módulo (uma
seleção inteira, dentro de um ZIP), onde não há FastAPI instalado.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

# No repositório, core/ está dois níveis acima (web/api/ -> raiz).  No nó do
# cluster o mesmo arquivo roda em projects/ntbuilder_cat/site/, com core/ ao
# lado, porque /home não existe lá.  Os dois caminhos entram, e o que tiver
# core/ ganha.
for _cand in (Path(__file__).resolve().parents[2], Path(__file__).resolve().parent):
    if (_cand / "core" / "io.py").exists() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

from core.io import LatticeStructure                                  # noqa: E402
from core.symmetry import snap_to_symmetry                            # noqa: E402
from core.chirality import compute_chirality                          # noqa: E402
from core.builder import build_nanotube                               # noqa: E402
from core.exporters import (write_xyz, write_cif, write_poscar,        # noqa: E402
                            write_lammps)

WRITERS = {"xyz": (write_xyz, "xyz"), "cif": (write_cif, "cif"),
           "poscar": (write_poscar, "POSCAR"), "lammps": (write_lammps, "lmp")}


def _render(nt, fmt: str, scratch: Path) -> tuple[str, str]:
    """O arquivo do tubo como texto.

    Os escritores de ``core.exporters`` recebem um caminho, nao um arquivo
    aberto -- e escrevem o corpo de uma vez, que e o que os deixa rapidos --,
    entao o texto passa por um arquivo de rascunho.  Em disco local do no isso
    nao custa nada e o ZIP sai do texto, sem um segundo caminho de escrita.
    """
    writer, ext = WRITERS[fmt]
    tmp = scratch / f"tube.{ext}"
    writer(nt, tmp)
    text = tmp.read_text()
    tmp.unlink(missing_ok=True)
    return text, ext

# Os mesmos parâmetros da varredura que gerou o catálogo; gravados também na
# tabela `meta` do banco, para o arquivo entregue ser o mesmo que a linha diz.
EPS_MAX = 0.5
NMAX = 50_000
SEARCH = 5000

CITATIONS = {
    "c2db": "S. Haastrup et al., 2D Mater. 5, 042002 (2018);\n"
            "M. N. Gjerding et al., 2D Mater. 8, 044002 (2021).",
    "2dmatpedia": "J. Zhou et al., Sci. Data 6, 86 (2019).",
    "jz1c03193": "High-Throughput Screening of Two-Dimensional Planar sp2 Carbon Space "
                 "Associated with a Labeled Quotient Graph, J. Phys. Chem. Lett. 12, "
                 "11511 (2021), doi:10.1021/acs.jpclett.1c03193 (material suplementar aberto).",
    "mc2d": "N. Mounet et al., Nat. Nanotechnol. 13, 246 (2018), doi:10.1038/s41565-017-0035-5;\n"
            "D. Campi et al., ACS Nano 17, 11268 (2023), doi:10.1021/acsnano.2c11510.",
    "jarvis": "K. Choudhary et al., Sci. Rep. 7, 5179 (2017), doi:10.1038/s41598-017-05402-0;\n"
              "K. Choudhary et al., npj Comput. Mater. 6, 173 (2020), doi:10.1038/s41524-020-00440-1.",
    "mxene": "D. Ontiveros et al., ACS Catal. 15, 14403 (2025), doi:10.1021/acscatal.5c04191.",
    "mathub2d": "M. Yao et al., Sci. China Mater. 66, 2768 (2023), doi:10.1007/s40843-022-2401-3.",
    "alexandria": "H.-C. Wang et al., 2D Mater. (2023), doi:10.1088/2053-1583/accc43;\n"
                  "Th. Cavignac et al., J. Phys. Mater. 9, 025014 (2026), doi:10.1088/2515-7639/ae6620.",
}


def structure_of(rec: dict) -> LatticeStructure:
    """A camada 2D como o manifesto a guardou, já na simetria da ferramenta."""
    atoms = [{"symbol": s, "pos": np.array([x, y], float), "z": z}
             for s, x, y, z in rec["atoms"]]
    s = LatticeStructure(atoms, np.array(rec["a1"], float),
                         np.array(rec["a2"], float),
                         source=f'{rec["source"]}:{rec["uid"]}')
    out, _ = snap_to_symmetry(s)
    return out


def chirality_of(s: LatticeStructure, n: int, m: int):
    """A mesma regra do catálogo: a menor célula com ε ≤ 0,5 %."""
    return compute_chirality(n, m, s, search_limit=SEARCH, max_strain=EPS_MAX,
                             max_cells=max(1, NMAX // len(s.atoms)))


def tube_text(rec: dict, n: int, m: int, fmt: str = "cif",
              vacuum: float = 10.0, inward: bool = False) -> tuple[str, str]:
    """O arquivo do tubo como texto, e o nome que ele deve ter."""
    if fmt not in WRITERS:
        raise ValueError(f"formato desconhecido: {fmt}")
    import tempfile
    s = structure_of(rec)
    ch = chirality_of(s, n, m)
    nt = build_nanotube(s, ch, vacuum=vacuum, roll_inward=inward)
    with tempfile.TemporaryDirectory(prefix="ntbcat_") as d:
        text, ext = _render(nt, fmt, Path(d))
    sense = "cw" if inward else "ccw"
    name = f'{rec["source"]}_{rec["uid"]}_{n}_{m}_{sense}.{ext}'
    return text, name


def build_one(rec: dict, n: int, m: int, fmt: str = "cif", vacuum: float = 10.0,
              inward: bool = False, out_dir: Path | None = None) -> tuple[Path, str]:
    """Escreve um tubo num arquivo temporário e devolve (caminho, nome)."""
    import tempfile
    text, name = tube_text(rec, n, m, fmt=fmt, vacuum=vacuum, inward=inward)
    d = Path(out_dir) if out_dir else Path(tempfile.mkdtemp(prefix="ntbcat_"))
    d.mkdir(parents=True, exist_ok=True)
    path = d / name
    path.write_text(text)
    return path, name


MANIFEST_COLS = ["source", "uid", "formula", "lattice_type", "n", "m", "sense",
                 "diameter", "d_inner", "thickness", "t_norm", "atoms", "eps",
                 "exact", "theta", "chirality", "spurious_ccw", "spurious_cw", "file"]


def sense_plan(row: dict, clean_only: bool) -> list[tuple[str, bool]]:
    """Os tubos que uma linha do catálogo entrega, como (sentido, para dentro).

    Quando os dois sentidos dão o mesmo tubo, entra um só, rotulado
    "equivalent".  Quando diferem, entram os dois, cada um só se passar no
    filtro de ligações espúrias.
    """
    if int(row.get("senses", 2)) == 1:
        return [("equivalent", False)] if (not clean_only or row.get("clean_ccw", 1)) else []
    return [(s, s == "cw") for s in ("ccw", "cw")
            if not clean_only or row.get(f"clean_{s}", 1)]


def write_selection(rows, structures, zip_path: Path, formats=("cif",),
                    vacuum: float = 10.0, clean_only: bool = True, progress=None,
                    rule: str = "") -> dict:
    """Constrói a seleção e a escreve num ZIP, com manifesto e citações.

    ``rows`` são as linhas do catálogo, ``structures`` um mapeamento
    material_id -> camada.  Cada linha entrega os tubos de :func:`sense_plan`.
    Um tubo que falhar entra em ``errors.txt`` e o resto segue.
    """
    zip_path = Path(zip_path)
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    plans = [sense_plan(r, clean_only) for r in rows]
    total = sum(len(pl) for pl in plans)
    done = n_ok = 0
    used_sources, errors = set(), []
    man = io.StringIO()
    w = csv.writer(man)
    w.writerow(MANIFEST_COLS)

    import tempfile
    scratch = Path(tempfile.mkdtemp(prefix="ntbcat_job_"))
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        cache: dict[int, LatticeStructure] = {}
        for row, plan in zip(rows, plans):
            if not plan:
                continue
            mid = row["material_id"]
            rec = structures[mid]
            try:
                s = cache.get(mid)
                if s is None:
                    s = cache[mid] = structure_of(rec)
                    if len(cache) > 64:
                        cache.pop(next(iter(cache)))
                ch = chirality_of(s, row["n"], row["m"])
            except Exception as e:
                errors.append(f'{rec["source"]}:{rec["uid"]} ({row["n"]},{row["m"]}): '
                              f"{type(e).__name__}: {e}")
                done += len(plan)
                if progress:
                    progress(done, total)
                continue
            for sense, inward in plan:
                try:
                    nt = build_nanotube(s, ch, vacuum=vacuum, roll_inward=inward)
                    suffix = "" if sense == "equivalent" else f"_{sense}"
                    names = []
                    for fmt in formats:
                        text, ext = _render(nt, fmt, scratch)
                        name = (f'{rec["source"]}/{rec["uid"]}_'
                                f'{row["n"]}_{row["m"]}{suffix}.{ext}')
                        z.writestr(name, text)
                        names.append(name)
                    w.writerow([rec["source"], rec["uid"], rec["formula"],
                                row["lattice_type"], row["n"], row["m"], sense,
                                round(row["diameter"], 4), round(row["d_inner"], 4),
                                round(row.get("thickness", 0.0), 4),
                                round(row["t_norm"], 4), row["atoms"],
                                row["eps"], row["exact"], round(row["theta"], 3),
                                "achiral" if row.get("kind") == 0 else "chiral",
                                row["spurious_ccw"], row["spurious_cw"],
                                ";".join(names)])
                    used_sources.add(rec["source"])
                    n_ok += 1
                except Exception as e:
                    errors.append(f'{rec["source"]}:{rec["uid"]} ({row["n"]},'
                                  f'{row["m"]}) {sense}: {type(e).__name__}: {e}')
                done += 1
                if progress and done % 200 == 0:
                    progress(done, total)
        z.writestr("manifest.csv", man.getvalue())
        cite = ["Nanotubos construídos com o NTBuilder.",
                f"Regra do catálogo: {rule}" if rule else "", "",
                "Cite as origens das camadas usadas nesta seleção:"]
        cite += [f"\n[{k}] {CITATIONS[k]}" for k in sorted(used_sources) if k in CITATIONS]
        z.writestr("CITATION.txt", "\n".join(cite) + "\n")
        if errors:
            z.writestr("errors.txt", "\n".join(errors) + "\n")
    import shutil
    shutil.rmtree(scratch, ignore_errors=True)
    if progress:
        progress(total, total)
    return {"n_files": n_ok, "n_errors": len(errors),
            "bytes": zip_path.stat().st_size}


def load_filter(path: Path) -> dict:
    return json.loads(Path(path).read_text())
