"""
Atende um pedido do catálogo — roda num nó do cluster, nunca no headnode.

Lê o filtro que o site gravou, seleciona as linhas no catálogo, constrói os
tubos, escreve o ZIP e vai avisando o progresso num arquivo ao lado.  O site
não fica esperando: ele lê ``progress.json`` quando alguém pergunta pelo
protocolo, e manda o e-mail final quando ``done.json`` aparece.

    NTB_PROTO=NTB-20260912-0001 NTB_OUT=<dir no NFS> python3 request_job.py

O progresso vai num arquivo, e não numa tabela: ``/home`` é local de cada nó,
então o job não alcança o SQLite de pedidos do site, e escrever SQLite
concorrente por NFS é justamente o que não se deve fazer.  Um JSON pequeno,
reescrito de forma atômica, atravessa o NFS sem drama.
"""
from __future__ import annotations

import json
import os
import resource
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
# Teto de memória: sem contabilidade de memória no Slurm daqui, um pedido
# patológico levaria o nó e os jobs dos outros.
_cap = int(os.environ.get("NTB_MEM_CAP_GB", "8")) * 2**30
resource.setrlimit(resource.RLIMIT_AS, (_cap, _cap))

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))            # core/ e cat_*.py ficam ao lado

import cat_build                                                  # noqa: E402
from cat_query import (SELECT_ROWS, material_ids, normalise,      # noqa: E402
                       where_sql)

PROTO = os.environ["NTB_PROTO"]
OUT = Path(os.environ["NTB_OUT"])
DB = Path(os.environ.get("NTB_CATALOGUE_DB",
                         HERE.parent / "catalogue.db"))


def _write_atomic(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload))
    tmp.replace(path)


def main() -> int:
    t0 = time.time()
    f = normalise(json.loads((OUT / "filter.json").read_text()))
    prog = OUT / "progress.json"
    _write_atomic(prog, {"state": "running", "done": 0, "total": 0,
                         "started": time.time()})

    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    ids, drop = material_ids(conn, f)
    where, p = where_sql(f, ids=ids, exclude=drop)
    order = "t.material_id, t.atoms"
    limit = int(f["per_material_max"])
    if limit:
        # Teto por material: os menores de cada um, para uma seleção ampla não
        # virar mil tubos do mesmo sistema.
        sql = (SELECT_ROWS + where +
               f" AND t.id IN (SELECT id FROM tubes s WHERE s.material_id ="
               f" t.material_id AND {where} ORDER BY s.atoms LIMIT {limit})"
               f" ORDER BY {order}")
        rows = [dict(r) for r in conn.execute(sql, p + p)]
    else:
        rows = [dict(r) for r in conn.execute(SELECT_ROWS + where +
                                              f" ORDER BY {order}", p)]
    mids = sorted({r["material_id"] for r in rows})
    structures = {}
    for i in range(0, len(mids), 500):
        chunk = mids[i:i + 500]
        for r in conn.execute("SELECT material_id, json FROM structures WHERE "
                              f"material_id IN ({','.join('?' * len(chunk))})", chunk):
            structures[r["material_id"]] = json.loads(r["json"])
    rule = dict(conn.execute("SELECT key, value FROM meta")).get("rule", "")
    conn.close()
    print(f"{PROTO}: {len(rows)} tubos, {len(structures)} materiais, "
          f"selecao em {time.time() - t0:.1f} s", flush=True)

    clean_only = f["clean"] == "clean"
    total = sum(len(cat_build.sense_plan(r, clean_only)) for r in rows)
    _write_atomic(prog, {"state": "running", "done": 0, "total": total,
                         "started": time.time()})

    last = [time.time()]

    def on_progress(done, tot):
        now = time.time()
        if now - last[0] < 2.0 and done < tot:
            return                       # o NFS não precisa de um write por tubo
        last[0] = now
        _write_atomic(prog, {"state": "running", "done": done, "total": tot,
                             "seconds": round(now - t0, 1)})

    scratch = Path(os.environ.get("SLURM_TMPDIR", tempfile.gettempdir())) / f"ntb_{PROTO}"
    scratch.mkdir(parents=True, exist_ok=True)
    zip_local = scratch / f"{PROTO}.zip"
    try:
        info = cat_build.write_selection(
            rows, structures, zip_local, formats=list(f["formats"]),
            vacuum=float(f["vacuum"]), clean_only=clean_only, progress=on_progress,
            rule=rule)
        dest = OUT / f"{PROTO}.zip"
        shutil.copy2(zip_local, dest)
        payload = {"state": "done", "done": total, "total": total,
                   "bytes": dest.stat().st_size, "n_files": info["n_files"],
                   "n_errors": info["n_errors"], "seconds": round(time.time() - t0, 1),
                   "path": str(dest)}
        _write_atomic(prog, payload)
        _write_atomic(OUT / "done.json", payload)
        print(f"{PROTO}: {info['n_files']} arquivos, {info['n_errors']} erros, "
              f"{dest.stat().st_size / 2**20:.1f} MB em {time.time() - t0:.0f} s",
              flush=True)
        return 0
    except Exception as e:
        payload = {"state": "failed", "message": f"{type(e).__name__}: {e}"[:400],
                   "seconds": round(time.time() - t0, 1)}
        _write_atomic(prog, payload)
        _write_atomic(OUT / "done.json", payload)
        print(f"{PROTO}: FALHOU {payload['message']}", flush=True)
        return 1
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
