"""
Checagem de ligações espúrias do mapa de quiralidade, fora do processo do site.

O /api/polar devolve o mapa na hora (a busca dos (n,m) leva menos de 1 s) e
dispara este processo, que monta cada tubo e confere as ligações, do menor para
o maior, gravando o progresso num JSON reescrito de forma atômica.  A página
lê esse arquivo a cada segundo e vai marcando os X.

Antes a mesma checagem rodava dentro da requisição: numa camada ondulada de
rede oblíqua (C2DB 1BrCrMoSTeSe3-1) eram 337 tubos, 2,9 milhões de átomos e
170 s com D até 25 Å, e o nginx devolvia 504 aos 120 s.

    python3 polar_spurious.py <diretório do job>

O diretório tem ``spec.json`` (arquivo da estrutura, parâmetros da busca e a
lista de (n,m) em ordem) e recebe ``progress.json``.  Um arquivo ``cancel`` no
diretório interrompe a checagem: é o que um mapa novo do mesmo arquivo faz.
"""
from __future__ import annotations

import json
import os
import resource
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
# Teto de memória: roda no mesmo servidor do site.
_cap = int(os.environ.get("NTB_SPURIOUS_MEM_GB", "6")) * 2**30
resource.setrlimit(resource.RLIMIT_AS, (_cap, _cap))

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from core.builder import check_curvature_bonds, curvature_tokens   # noqa: E402
from core.chirality import scan_chirality                        # noqa: E402
from core.io import load_structure                               # noqa: E402
from core.symmetry import snap_to_symmetry                       # noqa: E402


def _write(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload))
    tmp.replace(path)


def main(job_dir: Path) -> int:
    spec = json.loads((job_dir / "spec.json").read_text())
    prog = job_dir / "progress.json"
    points = spec["points"]
    total = len(points)
    results: dict[str, list[str] | None] = {}
    _write(prog, {"state": "running", "done": 0, "total": total, "results": results})
    try:
        s = load_structure(spec["struct_path"])
        try:
            s, _ = snap_to_symmetry(s)          # o mesmo que o /api/polar faz
        except Exception:
            pass
        found = {(r.n, r.m): r for r in scan_chirality(
            s, n_max=spec["n_max"], max_diameter=spec["max_diameter"], unique_only=True)}
        last = 0.0
        for i, (n, m) in enumerate(points, 1):
            if (job_dir / "cancel").exists():
                _write(prog, {"state": "cancelled", "done": i - 1, "total": total, "results": results})
                return 0
            pairs: list[str] | None = []
            r = found.get((n, m))
            if r is not None:
                try:
                    # Ligações formadas (+A-B) e rompidas (-A-B), átomo a átomo,
                    # a mesma checagem do catálogo.
                    pairs = curvature_tokens(*check_curvature_bonds(
                        s, r, roll_inward=spec["roll_inward"]))
                except MemoryError:
                    pairs = None                 # não verificado, nunca "limpo"
                except Exception:
                    # Falha de construção não é ligação espúria: o ponto fica
                    # como um ponto comum, igual ao que o /api/polar fazia.
                    pairs = []
            results[f"{n},{m}"] = pairs
            now = time.time()
            if now - last > 0.8 or i == total:
                _write(prog, {"state": "running", "done": i, "total": total, "results": results})
                last = now
        _write(prog, {"state": "done", "done": total, "total": total, "results": results})
        return 0
    except Exception as e:
        _write(prog, {"state": "failed", "done": len(results), "total": total, "results": results,
                      "message": f"{type(e).__name__}: {e}"[:300]})
        return 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
