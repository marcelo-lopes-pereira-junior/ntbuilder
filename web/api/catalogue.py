"""
Catálogo de nanotubos — consulta, pedido e entrega.

O catálogo é um SQLite de 7,9 milhões de tubos (``web/data/catalogue.db``,
escrito pela varredura do cluster e só lido aqui) e uma tabela de pedidos
separada (``web/data/requests.db``), que é a única coisa que este módulo
escreve.  Um tubo sai na hora.  Uma seleção com construção estimada em até
1 minuto é montada no próprio servidor, com prioridade baixa e uma de cada vez
por worker, e o resto vai para a fila ``high`` do cluster.  As duas rodam o
mesmo ``request_job.py`` e avisam o progresso pelos mesmos arquivos.

Três regras de consulta vêm de medida, não de gosto (ver
docs/CATALOGUE_DESIGN.md):

* filtro de elemento passa pela lista de materiais, nunca por ``NOT EXISTS``
  tubo a tubo — 227 ms contra 42 s;
* ``clean_ccw``/``clean_cw`` são 0/1 porque ``spurious = ''`` não entra em
  índice;
* a estimativa que acompanha o cursor sai da tabela ``summary`` (2 ms), e o
  número exato do índice (centenas de ms).
"""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import subprocess
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field, field_validator

_HERE = Path(__file__).parent
_DATA = _HERE.parent / "data"
CATALOGUE_DB = Path(os.environ.get("NTB_CATALOGUE_DB", _DATA / "catalogue.db"))
REQUESTS_DB = Path(os.environ.get("NTB_REQUESTS_DB", _DATA / "requests.db"))
# Onde os ZIPs dos pedidos ficam.  Fora de web/, porque quem escreve é o job
# do cluster e o que o site faz aqui é só servir o arquivo.
REQUEST_OUT = Path(os.environ.get(
    "NTB_REQUEST_OUT", Path.home() / "projects" / "ntbuilder_cat" / "out"))
REQUEST_SBATCH = _HERE.parent / "deploy" / "request.sbatch"
# O mesmo diretório tem dois nomes: no headnode ele está montado em
# ~/projects, e nos nós em /mnt/home/$USER/projects.  O job roda num nó, então
# tudo que é caminho passado para ele precisa ser traduzido -- sem isso o job
# não achava o filtro e morria antes de começar.
LOCAL_ROOT = Path(os.environ.get("NTB_LOCAL_ROOT", Path.home() / "projects"))
CLUSTER_ROOT = Path(os.environ.get(
    "NTB_CLUSTER_ROOT", f"/mnt/home/{os.environ.get('USER', 'mpereira')}/projects"))
RETENTION_DAYS = 7
QUEUE = os.environ.get("NTB_REQUEST_QUEUE", "high")
# Montagem na hora, no servidor.  Reserva do headnode: prioridade baixa,
# INSTANT_SLOTS montagens por worker do uvicorn, no máximo INSTANT_QUEUE_MAX
# esperando (acima disso o pedido vai para o cluster) e teto de memória.
INSTANT_MAX_SECONDS = float(os.environ.get("NTB_INSTANT_MAX_SECONDS", "60"))
INSTANT_SLOTS = int(os.environ.get("NTB_INSTANT_SLOTS", "1"))
INSTANT_QUEUE_MAX = int(os.environ.get("NTB_INSTANT_QUEUE_MAX", "6"))
INSTANT_MEM_GB = os.environ.get("NTB_INSTANT_MEM_GB", "6")
SITE_JOB = LOCAL_ROOT / "ntbuilder_cat" / "site" / "request_job.py"

import logging                                                     # noqa: E402

from . import mail                                                # noqa: E402
from .cat_query import (FORMATS, LATTICES, SOURCES, FilterError,   # noqa: E402
                        estimate, material_ids, material_where, summary_sql,
                        where_sql, tube_columns, cube_columns, normalise, material_cube_sql,
                        cube_rows_column, SELECT_ROWS)

_log = logging.getLogger("ntbuilder.catalogue")
# Usado nos links dos e-mails; o servidor não sabe o nome público sozinho.
BASE_URL = os.environ.get("NTB_BASE_URL", "https://nanoeng.unb.br/ntbuilder")

router = APIRouter(prefix="/api/cat", tags=["catalogue"])


# ── conexões ────────────────────────────────────────────────────────────────

def _cat() -> sqlite3.Connection:
    """O catálogo, aberto somente para leitura.

    Modo read-only de propósito: se o arquivo não existir, o erro aparece aqui
    e não como uma tabela vazia, e nenhum caminho deste módulo pode escrever
    no banco que a varredura produziu.
    """
    if not CATALOGUE_DB.exists():
        raise HTTPException(503, "catálogo ainda não instalado neste servidor")
    conn = sqlite3.connect(f"file:{CATALOGUE_DB}?mode=ro", uri=True,
                           check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


def _req() -> sqlite3.Connection:
    REQUESTS_DB.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(REQUESTS_DB), check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("""CREATE TABLE IF NOT EXISTS requests (
        protocol TEXT PRIMARY KEY, email TEXT, name TEXT,
        filter_json TEXT, filter_hash TEXT,
        n_tubes INTEGER, n_atoms INTEGER, formats TEXT,
        state TEXT, slurm_id TEXT, progress INTEGER, bytes INTEGER,
        created TEXT, finished TEXT, expires TEXT, path TEXT, message TEXT)""")
    conn.execute("CREATE INDEX IF NOT EXISTS requests_hash "
                 "ON requests(filter_hash, state)")
    conn.execute("CREATE INDEX IF NOT EXISTS requests_email "
                 "ON requests(email, state)")
    return conn


# ── o filtro ────────────────────────────────────────────────────────────────

class CatFilter(BaseModel):
    """A seleção do usuário, em uma peça, para poder ser guardada e repetida."""

    sources:  Optional[list[str]] = None
    lattices: Optional[list[str]] = None
    # "contém todos estes" e "não contém nenhum destes"
    elements_all:  Optional[list[str]] = None
    elements_none: Optional[list[str]] = None
    # "contém uma classe": cada grupo pede ao menos um elemento dele
    elements_any: Optional[list[list[str]]] = None
    species_min: int = Field(1, ge=1, le=12)
    species_max: int = Field(12, ge=1, le=12)

    d_min: float = Field(0.0, ge=0.0, le=500.0)
    d_max: float = Field(80.0, ge=0.5, le=500.0)
    t_min: float = Field(0.0, ge=0.0, le=500.0)     # espessura da camada
    t_max: float = Field(60.0, ge=0.1, le=500.0)
    din_min: float = Field(0.0, ge=0.0, le=500.0)   # cavidade interna, D - t
    atoms_max: int = Field(50_000, ge=1, le=50_000)
    eps_max: float = Field(0.5, ge=0.0, le=0.5)
    exact_only: bool = False
    # achiral (Ch ao longo de um espelho ou deslizamento do sistema 2D), chiral
    chirality: Optional[list[str]] = None
    # 'clean' = só tubos sem ligação espúria (cada sentido distinto julgado
    # por conta própria), 'all' = todos
    clean: str = "clean"
    per_material_max: int = Field(0, ge=0, le=10_000)   # 0 = sem limite

    formats: list[str] = ["cif"]
    vacuum: float = Field(10.0, ge=1.0, le=50.0)

    @field_validator("sources", "lattices", "chirality", "formats",
                     "elements_all", "elements_none", mode="before")
    @classmethod
    def _clean_list(cls, v):
        if v is None:
            return None
        return [str(x).strip() for x in v if str(x).strip()]

    @field_validator("elements_any", mode="before")
    @classmethod
    def _clean_groups(cls, v):
        if v is None:
            return None
        groups = [[str(x).strip() for x in g if str(x).strip()][:120] for g in v][:12]
        return [g for g in groups if g] or None

    @field_validator("formats")
    @classmethod
    def _known_formats(cls, v):
        bad = [f for f in v if f not in FORMATS]
        if bad:
            raise ValueError(f"formato desconhecido: {', '.join(bad)}")
        return v or ["cif"]

    def key(self) -> str:
        """Impressão digital do filtro, para o cache de pedidos."""
        payload = json.dumps(self.model_dump(), sort_keys=True,
                             separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _sql(f: CatFilter, conn: sqlite3.Connection) -> tuple[str, list]:
    """O WHERE do filtro, com o erro de filtro virando 400.

    A lista de materiais é resolvida antes de montar a consulta grande: é o
    que escolhe o plano certo (ver :func:`cat_query.material_ids`).
    """
    try:
        d = f.model_dump()
        ids, drop = material_ids(conn, d)
        return where_sql(d, ids=ids, exclude=drop)
    except FilterError as e:
        raise HTTPException(400, str(e)) from e


# ── consulta ────────────────────────────────────────────────────────────────

@router.get("/meta")
async def meta():
    conn = _cat()
    try:
        m = {r["key"]: r["value"] for r in conn.execute("SELECT key, value FROM meta")}
        # Pelo cubo, não pela tabela de tubos: o GROUP BY sobre 10 milhões de
        # linhas levava 5 s a cada carga da página, e a soma do cubo é a mesma.
        by_source = [dict(r) for r in conn.execute(
            "SELECT source, SUM(n_tubes) AS n_tubes FROM summary GROUP BY source")]
        by_lattice = [dict(r) for r in conn.execute(
            "SELECT lattice_type, SUM(n_tubes) AS n_tubes FROM summary "
            "GROUP BY lattice_type ORDER BY 2 DESC")]
        # Quantos materiais têm cada elemento: é o que a tabela periódica da
        # página mostra ao passar o mouse, e o que apaga os 34 que não existem.
        element_counts = {r[0]: r[1] for r in conn.execute(
            "SELECT element, COUNT(*) FROM material_elements GROUP BY element")}
        # Por Base de Dados 2D: Sistemas 2D, tubos, exatos e limpos.  A página
        # soma as bases marcadas sem perguntar de novo ao servidor.
        n_mat = {r[0]: r[1] for r in conn.execute(
            "SELECT source, COUNT(*) FROM materials GROUP BY source")}
        by_source_stats = [
            {"source": r[0], "n_materials": n_mat.get(r[0], 0), "n_tubes": r[1],
             "n_exact": r[2], "n_clean": r[3]}
            for r in conn.execute(
                "SELECT source, SUM(n_tubes), "
                "SUM(CASE WHEN exact = 1 THEN n_tubes ELSE 0 END), "
                "SUM(n_clean) "
                "FROM summary GROUP BY source")]
        # Os limites reais do banco, para os controles da página não passarem
        # do que existe.
        limits = {"d_max": conn.execute("SELECT MAX(d_bin) + 1 FROM summary").fetchone()[0],
                  "t_max": conn.execute("SELECT MAX(thickness) FROM materials").fetchone()[0]}
        return {"meta": m, "by_source": by_source, "by_lattice": by_lattice,
                "by_source_stats": by_source_stats, "limits": limits,
                "elements": sorted(element_counts),
                "element_counts": element_counts,
                "formats": list(FORMATS), "retention_days": RETENTION_DAYS,
                # O que esta API sabe fazer.  A página esconde o que faltar, para
                # não mostrar contagem errada enquanto o serviço não reinicia.
                "features": ["elements_any", "pages", "instant", "layer", "systems", "polar"],
                "instant_max_seconds": INSTANT_MAX_SECONDS}
    finally:
        conn.close()


@router.post("/elements")
async def element_counts(f: CatFilter):
    """Quantos materiais de cada elemento sobram com a seleção atual.

    É o que apaga, na tabela periódica, os elementos que não combinam com o
    que já foi marcado: com Mo e S escolhidos, só fica aceso o que aparece
    junto com os dois em algum material.  Conta materiais, não tubos — a
    pergunta é de química, e responde em milissegundos sobre 24 586 linhas.
    """
    try:
        d = f.model_dump()
        where, p = material_where(d)
        cube_where, cp = material_cube_sql(d)
    except FilterError as e:
        raise HTTPException(400, str(e)) from e
    # Um Sistema 2D conta quando passa nos filtros de material E tem pelo menos
    # um tubo que passa nos filtros de tubo.  Antes só o primeiro contava, e a
    # tabela deixava acesos elementos que davam 0 tubos na seleção.
    #
    # EXISTS por Sistema 2D, e não IN sobre o cubo: o índice começa pelo
    # material e a busca para no primeiro tubo que passa.  Medido: "Mo, Fe e O,
    # sem ligação espúria" caiu de 2,2 s para 3 ms, o filtro padrão de 2,9 s
    # para 1 s.  O conjunto é materializado uma vez e serve às duas contagens.
    ok = ("SELECT m.id AS id FROM materials m WHERE " + where +
          " AND EXISTS (SELECT 1 FROM material_cube WHERE material_id = m.id AND "
          + cube_where + ")")
    params = p + cp
    conn = _cat()
    try:
        t0 = time.time()
        rows = conn.execute(
            "WITH ok AS MATERIALIZED (" + ok + ") "
            "SELECT e.element, COUNT(*), (SELECT COUNT(*) FROM ok) "
            "FROM material_elements e JOIN ok ON ok.id = e.material_id "
            "GROUP BY e.element", params).fetchall()
        counts = {r[0]: r[1] for r in rows}
        n = rows[0][2] if rows else 0
        return {"n_materials": n, "counts": counts,
                "ms": round((time.time() - t0) * 1e3, 1)}
    finally:
        conn.close()


@router.post("/count")
async def count(f: CatFilter, fast: bool = False):
    """Quantos tubos, quantos átomos, quantos bytes, quanto tempo.

    Soma ``n_clean`` ou ``n_all``, não conta linhas: um (n,m) de um Sistema 2D
    assimétrico vale dois tubos.  ``fast=true`` responde pelo cubo quando ele
    conhece a pergunta, e diz se o número é exato.
    """
    conn = _cat()
    try:
        t0 = time.time()
        try:
            d = f.model_dump()
            normalise(d)
        except FilterError as e:
            raise HTTPException(400, str(e)) from e
        if fast:
            cube = summary_sql(d)
            if cube is not None:
                where, p, cube_exact = cube
                c_t, c_a = cube_columns(d)
                row = conn.execute(
                    f"SELECT COALESCE(SUM({c_t}),0), COALESCE(SUM({c_a}),0) "
                    "FROM summary WHERE " + where, p).fetchone()
                out = estimate(int(row[0]), int(row[1]), f.formats)
                out.update(exact_count=cube_exact, source="cubo",
                           n_rows=_cube_rows(conn, d, where, p),
                           ms=round((time.time() - t0) * 1e3, 1))
                return out

        where, p = _sql(f, conn)
        c_t, c_a = tube_columns(d)
        row = conn.execute(f"SELECT COALESCE(SUM({c_t}),0), COALESCE(SUM({c_a}),0), "
                           "COUNT(*) FROM tubes WHERE " + where, p).fetchone()
        out = estimate(int(row[0]), int(row[1]), f.formats)
        out["n_rows"] = int(row[2])
        # COUNT(DISTINCT) sai do índice de cobertura para uma árvore temporária
        # (14 s no filtro padrão); só numa seleção pequena.
        if int(row[0]) <= 200_000:
            out["n_materials"] = conn.execute(
                "SELECT COUNT(DISTINCT material_id) FROM tubes WHERE " + where,
                p).fetchone()[0]
        out.update(exact_count=True, source="indice",
                   ms=round((time.time() - t0) * 1e3, 1))
        return out
    finally:
        conn.close()


def _cube_rows(conn: sqlite3.Connection, d: dict, where: str, p: list) -> int | None:
    """Linhas da seleção pelo cubo, para a amostra saber quantas páginas tem.

    None num banco de antes das colunas n_rows e n_rows_clean.
    """
    try:
        return int(conn.execute(f"SELECT COALESCE(SUM({cube_rows_column(d)}),0) "
                                "FROM summary WHERE " + where, p).fetchone()[0])
    except sqlite3.OperationalError:
        return None


@router.post("/sample")
async def sample(f: CatFilter, limit: int = 12, offset: int = 0,
                 material_id: Optional[int] = None):
    """Uma página da seleção, na ordem do banco (agrupada por Sistema 2D).

    Na ordem do banco e não por tamanho de célula: ORDER BY atoms ordena
    milhões de linhas numa árvore temporária (2,5 s na seleção padrão), e a
    ordem de gravação sai do índice na hora, em qualquer página.
    """
    limit = max(1, min(int(limit), 200))
    offset = max(0, int(offset))
    conn = _cat()
    where, p = _sql(f, conn)
    if material_id is not None:              # os tubos de um Sistema 2D da amostra
        where, p = where + " AND material_id = ?", p + [int(material_id)]
    try:
        # O filtro vale só para tubes e usa nomes sem prefixo (source e
        # lattice_type existem nas duas tabelas): filtra numa subconsulta e só
        # depois junta com materials.  Na junção direta, filtrar por base dava
        # "ambiguous column name: source".
        rows = conn.execute(
            "SELECT t.material_id, m.source, m.uid, m.formula, m.lattice_type, "
            "       ROUND(m.thickness,2) AS thickness, t.n, t.m, t.diameter, "
            "       t.d_inner, t.t_norm, t.atoms, t.eps, t.exact, t.theta, t.kind, "
            "       t.senses, t.clean_ccw, t.clean_cw, t.spurious_ccw, t.spurious_cw "
            "FROM (SELECT * FROM tubes WHERE " + where + " ORDER BY id LIMIT ? OFFSET ?) t "
            "JOIN materials m ON m.id = t.material_id ORDER BY t.id",
            p + [limit, offset]).fetchall()
        return {"tubes": [dict(r) for r in rows], "limit": limit, "offset": offset}
    finally:
        conn.close()


@router.post("/systems")
async def systems(f: CatFilter, limit: int = 12, offset: int = 0):
    """Uma página dos Sistemas 2D da seleção, cada um com os números dos seus tubos.

    A amostra por tubo mostrava dezenas de páginas do mesmo Sistema 2D (com Li
    marcado, 35 páginas de Li2VF6).  Aqui a unidade é o Sistema 2D: quais passam
    sai do mesmo teste do /elements (algum tubo do cubo por Sistema 2D passa no
    filtro), e o total de páginas é o n_materials dele.  Os números de cada um
    vêm do índice, só para os da página.
    """
    limit = max(1, min(int(limit), 100))
    offset = max(0, int(offset))
    try:
        d = f.model_dump()
        mw, mp = material_where(d)
        cw, cp = material_cube_sql(d)
    except FilterError as e:
        raise HTTPException(400, str(e)) from e
    conn = _cat()
    try:
        ids = [r[0] for r in conn.execute(
            "SELECT m.id FROM materials m WHERE " + mw +
            " AND EXISTS (SELECT 1 FROM material_cube WHERE material_id = m.id AND " + cw + ")"
            " ORDER BY m.id LIMIT ? OFFSET ?", mp + cp + [limit, offset])]
        if not ids:
            return {"systems": [], "limit": limit, "offset": offset}
        where, p = _sql(f, conn)
        c_t, _ = tube_columns(d)
        marks = ",".join("?" * len(ids))
        stats = {r[0]: r for r in conn.execute(
            f"SELECT material_id, SUM({c_t}), COUNT(*), MIN(diameter), MAX(diameter), MIN(atoms) "
            f"FROM tubes WHERE {where} AND material_id IN ({marks}) GROUP BY material_id", p + ids)}
        info = {r["id"]: r for r in conn.execute(
            "SELECT id, source, uid, formula, lattice_type, thickness FROM materials "
            f"WHERE id IN ({marks})", ids)}
        out = []
        for mid in ids:
            st, m = stats.get(mid), info[mid]
            if st is None:
                continue
            out.append({"material_id": mid, "source": m["source"], "uid": m["uid"],
                        "formula": m["formula"], "lattice_type": m["lattice_type"],
                        "thickness": round(m["thickness"], 2), "n_tubes": int(st[1]),
                        "n_rows": int(st[2]), "d_min": round(st[3], 2), "d_max": round(st[4], 2),
                        "atoms_min": int(st[5])})
        return {"systems": out, "limit": limit, "offset": offset}
    finally:
        conn.close()


@router.get("/tube/{material_id}/{n}/{m}")
async def tube(material_id: int, n: int, m: int, fmt: str = "cif",
               vacuum: float = 10.0, inward: bool = False):
    """Um tubo, construído na hora a partir da camada guardada no catálogo."""
    if fmt not in FORMATS:
        raise HTTPException(400, f"formato desconhecido: {fmt}")
    conn = _cat()
    try:
        row = conn.execute("SELECT json FROM structures WHERE material_id = ?",
                           (material_id,)).fetchone()
        if row is None:
            raise HTTPException(404, "material não está no catálogo")
        tube_row = conn.execute("SELECT eps, atoms FROM tubes WHERE material_id = ? "
                                "AND n = ? AND m = ?", (material_id, n, m)).fetchone()
    finally:
        conn.close()
    if tube_row is None:
        raise HTTPException(404, "este (n,m) não está no catálogo deste material")

    from .cat_build import build_one            # importa core sob demanda
    path, name = build_one(json.loads(row["json"]), n, m, fmt=fmt,
                           vacuum=vacuum, inward=inward)
    return FileResponse(str(path), filename=name, media_type="text/plain")


def layer_cif(rec: dict, vacuum: float = 20.0) -> str:
    """A camada 2D guardada no catálogo como CIF P1, para abrir no construtor.

    É a mesma camada que gerou os tubos (o registro do manifesto, já no plano
    xy), sem a simetria aplicada: o construtor faz a própria busca ao carregar.
    O eixo c leva a espessura mais o vácuo, com a camada no meio.
    """
    import math
    import numpy as np
    a1 = np.array(rec["a1"], float)
    a2 = np.array(rec["a2"], float)
    xyz = np.array([[at[1], at[2], at[3]] for at in rec["atoms"]], float)
    z = xyz[:, 2]
    c = float(z.max() - z.min()) + vacuum
    frac_xy = np.linalg.solve(np.array([a1, a2]).T, xyz[:, :2].T).T % 1.0
    frac_z = (z - 0.5 * (z.max() + z.min())) / c + 0.5
    a, b = float(np.linalg.norm(a1)), float(np.linalg.norm(a2))
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, float(a1 @ a2) / (a * b)))))
    name = f'{rec["source"]}_{rec["uid"]}'.replace(" ", "_")
    lines = [f"data_{name}",
             f"# {rec.get('formula', '')} — catálogo de nanotubos do NTBuilder, origem {rec['source']}:{rec['uid']}",
             f"_cell_length_a    {a:.6f}", f"_cell_length_b    {b:.6f}", f"_cell_length_c    {c:.6f}",
             "_cell_angle_alpha 90.000000", "_cell_angle_beta  90.000000", f"_cell_angle_gamma {gamma:.6f}",
             "_symmetry_space_group_name_H-M 'P 1'", "_symmetry_Int_Tables_number 1",
             "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
             "loop_", "_atom_site_label", "_atom_site_type_symbol",
             "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    for k, (at, (fx, fy), fz) in enumerate(zip(rec["atoms"], frac_xy, frac_z), 1):
        lines.append(f"{at[0]}{k} {at[0]} {fx:.8f} {fy:.8f} {fz:.8f}")
    return "\n".join(lines) + "\n"


@router.get("/polar/{material_id}")
async def polar_from_catalogue(material_id: int, roll_inward: bool = False):
    """O mapa de quiralidade de um Sistema 2D a partir do banco, sem calcular nada.

    Mesmo formato do /api/polar do construtor.  Os pontos são os tubos do
    catálogo (a menor célula com ε ≤ 0,5 %, até 50 000 átomos, na janela de
    diâmetro do catálogo), e os X vêm da checagem do sentido pedido.
    """
    import math
    conn = _cat()
    try:
        mat = conn.execute("SELECT source, uid, lattice_type, sector_deg, a, b, gamma, thickness "
                           "FROM materials WHERE id = ?", (material_id,)).fetchone()
        if mat is None:
            raise HTTPException(404, "Sistema 2D não está no catálogo")
        rows = conn.execute("SELECT n, m, diameter, theta, eps, exact, atoms, spurious_ccw, spurious_cw "
                            "FROM tubes WHERE material_id = ?", (material_id,)).fetchall()
        rec = json.loads(conn.execute("SELECT json FROM structures WHERE material_id = ?",
                                      (material_id,)).fetchone()["json"])
    finally:
        conn.close()
    if not rows:
        raise HTTPException(404, "Sistema 2D sem tubos no catálogo")
    # O mesmo início de setor do /api/polar: o complemento do maior vão angular.
    th = sorted(float(r["theta"]) % 180.0 for r in rows)
    gaps = [(th[(i + 1) % len(th)] - th[i]) % 180.0 for i in range(len(th))]
    k = max(range(len(gaps)), key=lambda i: gaps[i])
    start = th[(k + 1) % len(th)] if len(th) > 1 else 0.0
    col = "spurious_cw" if roll_inward else "spurious_ccw"
    buckled = float(mat["thickness"] or 0.0) > 1e-3
    points = []
    for r in rows:
        tp = math.radians((float(r["theta"]) - start) % 180.0)
        pt = {"n": r["n"], "m": r["m"], "diameter": round(r["diameter"], 4),
              "theta_deg": round(float(r["theta"]) % 180.0, 4),
              "strain": 0.0 if r["exact"] else round(float(r["eps"]), 6),
              "n_atoms": r["atoms"],
              "x": round(r["diameter"] * math.cos(tp), 4), "y": round(r["diameter"] * math.sin(tp), 4)}
        if buckled:
            pt["spurious"] = [x for x in (r[col] or "").split(",") if x]
        points.append(pt)
    # Menor distância entre átomos da camada, como no painel da estrutura.
    A = [rec["a1"], rec["a2"]]
    atoms = rec["atoms"]
    dmin = math.inf
    for i, ai in enumerate(atoms):
        for aj in atoms:
            for p in (-1, 0, 1):
                for q in (-1, 0, 1):
                    dx = aj[1] + p * A[0][0] + q * A[1][0] - ai[1]
                    dy = aj[2] + p * A[0][1] + q * A[1][1] - ai[2]
                    d = math.sqrt(dx * dx + dy * dy + (aj[3] - ai[3]) ** 2)
                    if d > 1e-6:
                        dmin = min(dmin, d)
    species = sorted({a[0] for a in atoms})
    return {"points": points, "dmax": round(max(r["diameter"] for r in rows), 3),
            "theta_max": float(mat["sector_deg"] or 180.0), "theta_start": float(start),
            "lattice_type": mat["lattice_type"], "a": round(mat["a"], 4), "b": round(mat["b"], 4),
            "gamma_deg": round(mat["gamma"], 2), "n_species": len(species), "species": species,
            "d_min": round(dmin, 4) if math.isfinite(dmin) else None, "snap_desc": "",
            "spurious_job": None, "source": "catalogue",
            "d_cover": [round(min(r["diameter"] for r in rows), 3), round(max(r["diameter"] for r in rows), 3)],
            "rule": "menor célula com ε ≤ 0,5 %, até 50 000 átomos"}


@router.get("/layer/{material_id}")
async def layer(material_id: int):
    """O CIF da camada 2D de um Sistema 2D do catálogo ("Abrir no NTBuilder")."""
    from fastapi.responses import Response
    conn = _cat()
    try:
        row = conn.execute("SELECT json FROM structures WHERE material_id = ?",
                           (material_id,)).fetchone()
    finally:
        conn.close()
    if row is None:
        raise HTTPException(404, "Sistema 2D não está no catálogo")
    rec = json.loads(row["json"])
    fname = f'{rec["source"]}_{rec["uid"]}.cif'.replace(" ", "_").replace("/", "_")
    return Response(layer_cif(rec), media_type="chemical/x-cif",
                    headers={"Content-Disposition": f'inline; filename="{fname}"',
                             "X-NTB-Filename": fname})


# ── pedido ──────────────────────────────────────────────────────────────────

class RequestIn(BaseModel):
    filter: CatFilter
    email: str
    name: str = ""

    @field_validator("email")
    @classmethod
    def _looks_like_email(cls, v):
        v = v.strip()
        if "@" not in v or "." not in v.split("@")[-1] or len(v) < 6:
            raise ValueError("e-mail inválido")
        return v


def _protocol(conn: sqlite3.Connection) -> str:
    """O próximo protocolo do dia cuja pasta de saída ainda não existe.

    A pasta decide, não só a contagem de linhas: um servidor de teste com outro
    banco de pedidos numerou um pedido como NTB-20260913-0001, a mesma pasta de
    um pedido real, e sobrescreveu o arquivo dele (13/09/2026).
    """
    day = datetime.now(timezone.utc).strftime("%Y%m%d")
    n = conn.execute("SELECT COUNT(*) FROM requests WHERE protocol LIKE ?",
                     (f"NTB-{day}-%",)).fetchone()[0]
    while True:
        n += 1
        proto = f"NTB-{day}-{n:04d}"
        taken = conn.execute("SELECT 1 FROM requests WHERE protocol = ?", (proto,)).fetchone()
        if not taken and not (REQUEST_OUT / proto).exists():
            return proto


@router.post("/request")
async def make_request(req: RequestIn):
    """Grava o pedido, manda para a fila do cluster e devolve o protocolo."""
    f = req.filter
    fhash = f.key()
    conn = _req()
    try:
        # Cache: pedido idêntico ainda no prazo devolve o mesmo arquivo.
        now = datetime.now(timezone.utc)
        hit = conn.execute(
            "SELECT protocol, expires, path FROM requests WHERE filter_hash = ? "
            "AND state = 'done' ORDER BY finished DESC LIMIT 1", (fhash,)).fetchone()
        if hit and hit["expires"] and hit["expires"] > now.isoformat():
            return {"protocol": hit["protocol"], "state": "done", "cached": True,
                    "mode": "cached",
                    "download": f"/api/cat/download/{hit['protocol']}"}

        est = await count(f)
        if est["n_tubes"] == 0:
            raise HTTPException(400, "a seleção está vazia")

        # Fila justa: sem teto de tamanho, mas um e-mail não enfileira o
        # cluster inteiro de uma vez.
        pend = conn.execute("SELECT COUNT(*) FROM requests WHERE email = ? AND "
                            "state IN ('queued','running')", (req.email,)).fetchone()[0]
        if pend >= 3:
            raise HTTPException(429, "você já tem 3 pedidos na fila; "
                                     "espere um terminar")

        # Até 1 minuto de construção estimada: monta aqui e o download sai na
        # tela.  O e-mail com o protocolo e o link de 7 dias vai do mesmo jeito.
        instant = (est["seconds"] <= INSTANT_MAX_SECONDS and SITE_JOB.exists()
                   and _instant_state["waiting"] < INSTANT_QUEUE_MAX)

        proto = _protocol(conn)
        expires = now + timedelta(days=RETENTION_DAYS)
        # resolve(): ~/projects é um link para /mnt/home/$USER/projects, e o nó
        # do cluster só vê o caminho do NFS -- /home é local de cada máquina.
        # O job recebia o caminho do headnode e não achava o filtro.
        out_dir = REQUEST_OUT / proto
        out_dir.mkdir(parents=True, exist_ok=False)     # nunca por cima de outro pedido
        (out_dir / "filter.json").write_text(
            json.dumps(f.model_dump(), indent=1, ensure_ascii=False))
        conn.execute(
            "INSERT INTO requests (protocol, email, name, filter_json, filter_hash,"
            " n_tubes, n_atoms, formats, state, progress, created, expires, path)"
            " VALUES (?,?,?,?,?,?,?,?,'queued',0,?,?,?)",
            (proto, req.email, req.name, json.dumps(f.model_dump()), fhash,
             est["n_tubes"], est["n_atoms"], ",".join(f.formats),
             now.isoformat(), expires.isoformat(),
             str(out_dir / f"{proto}.zip")))
        conn.commit()

        if instant:
            import asyncio
            conn.execute("UPDATE requests SET slurm_id = 'local:wait' WHERE protocol = ?",
                         (proto,))
            conn.commit()
            asyncio.get_running_loop().create_task(_run_local(proto, out_dir))
            return {"protocol": proto, "state": "queued", "mode": "instant",
                    "cached": False, "estimate": est, "expires": expires.isoformat(),
                    "status": f"/api/cat/status/{proto}",
                    "download": f"/api/cat/download/{proto}"}

        slurm_id, err = _submit(proto, out_dir)
        if slurm_id:
            conn.execute("UPDATE requests SET slurm_id = ? WHERE protocol = ?",
                         (slurm_id, proto))
        else:
            conn.execute("UPDATE requests SET state = 'failed', message = ? "
                         "WHERE protocol = ?", (err[:400], proto))
        conn.commit()
        if not slurm_id:
            raise HTTPException(503, f"não foi possível submeter o pedido: {err[:200]}")
        sent, why = mail.queued(req.email, req.name, proto, est,
                                _filter_text(f), BASE_URL)
        return {"protocol": proto, "state": "queued", "cached": False, "mode": "cluster",
                "estimate": est, "expires": expires.isoformat(),
                "mail_sent": sent, "mail_error": why,
                "status": f"/api/cat/status/{proto}"}
    finally:
        conn.close()


_instant_state = {"waiting": 0, "sem": None}


def _set_sid(proto: str, sid: str) -> None:
    conn = _req()
    try:
        conn.execute("UPDATE requests SET slurm_id = ? WHERE protocol = ?", (sid, proto))
        conn.commit()
    finally:
        conn.close()


async def _run_local(proto: str, out_dir: Path) -> None:
    """Monta um pedido pequeno no servidor, com o mesmo request_job.py do cluster.

    O job escreve progress.json e done.json, e o observador de pedidos manda o
    e-mail do link como faz com os do cluster.  Roda com nice 10 e teto de
    memória, uma montagem por vez em cada worker.
    """
    import asyncio
    import sys
    st = _instant_state
    if st["sem"] is None:
        st["sem"] = asyncio.Semaphore(INSTANT_SLOTS)
    st["waiting"] += 1
    acquired = False
    try:
        async with st["sem"]:
            st["waiting"] -= 1
            acquired = True
            env = dict(os.environ, NTB_PROTO=proto, NTB_OUT=str(out_dir),
                       NTB_CATALOGUE_DB=str(CATALOGUE_DB), NTB_MEM_CAP_GB=INSTANT_MEM_GB,
                       OMP_NUM_THREADS="1")
            with open(out_dir / "job.log", "ab") as log:
                proc = await asyncio.create_subprocess_exec(
                    "nice", "-n", "10", sys.executable, str(SITE_JOB),
                    env=env, stdout=log, stderr=subprocess.STDOUT, cwd=str(out_dir))
            _set_sid(proto, f"local:{proc.pid}")
            await proc.wait()
    except Exception as e:                      # o observador transforma em falha
        _log.exception("montagem local %s", proto)
        (out_dir / "done.json").write_text(json.dumps(
            {"state": "failed", "message": f"{type(e).__name__}: {e}"[:400]}))
    finally:
        if not acquired:
            st["waiting"] -= 1


def _fail(conn: sqlite3.Connection, row: sqlite3.Row, msg: str) -> None:
    """Marca a falha e manda o e-mail uma vez só, mesmo com dois workers olhando."""
    cur = conn.execute("UPDATE requests SET state = 'failed', finished = ?, message = ?"
                       " WHERE protocol = ? AND state IN ('queued','running')",
                       (datetime.now(timezone.utc).isoformat(), msg, row["protocol"]))
    conn.commit()
    if cur.rowcount:
        mail.failed(row["email"], row["name"] or "", row["protocol"], msg, BASE_URL)


def _check_local(conn: sqlite3.Connection, row: sqlite3.Row) -> None:
    """Montagem no servidor que morreu sem escrever o resultado (reinício, memória)."""
    sid = row["slurm_id"] or ""
    age = time.time() - _epoch(row["created"])
    if sid == "local:wait":
        if age < 3600:
            return
        _fail(conn, row, "o serviço reiniciou antes de montar o pedido")
        return
    try:
        pid = int(sid.split(":", 1)[1])
    except (IndexError, ValueError):
        return
    if age < 20 or Path(f"/proc/{pid}").exists():
        return
    if (Path(row["path"]).parent / "done.json").exists():
        return
    _fail(conn, row, "a montagem no servidor terminou sem gerar o arquivo")


def _cluster_path(p: Path) -> str:
    """O mesmo caminho, como o nó do cluster o vê."""
    p = Path(p)
    try:
        return str(CLUSTER_ROOT / p.relative_to(LOCAL_ROOT))
    except ValueError:
        return str(p)               # já é um caminho do NFS, ou está fora dele


def _submit(proto: str, out_dir: Path) -> tuple[str, str]:
    """``sbatch`` na fila rápida.  Construção de tubo nunca roda no headnode."""
    cmd = ["sbatch", "--parsable", "-p", QUEUE, "-J", f"ntb_req_{proto}",
           "-n", "1", "-c", "4", "-t", "01:00:00",
           f"--export=ALL,NTB_PROTO={proto},NTB_OUT={_cluster_path(out_dir)}",
           str(REQUEST_SBATCH)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except Exception as e:                       # sbatch ausente, fila fora
        return "", f"{type(e).__name__}: {e}"
    if r.returncode != 0:
        return "", (r.stderr or r.stdout or "sbatch falhou").strip()
    return r.stdout.strip().split(";")[0], ""


@router.get("/status/{protocol}")
async def status(protocol: str):
    conn = _req()
    try:
        row = conn.execute("SELECT * FROM requests WHERE protocol = ?",
                           (protocol,)).fetchone()
        if row is None:
            raise HTTPException(404, "protocolo não encontrado")
        _refresh(conn, row)
        row = conn.execute("SELECT * FROM requests WHERE protocol = ?",
                           (protocol,)).fetchone()
        d = dict(row)
        d.pop("email", None)                    # o protocolo não revela o e-mail
        total = d["n_tubes"] or 0
        d["percent"] = (round(100.0 * min(d["progress"] or 0, total) / total, 1)
                        if total else 0.0)
        d["mode"] = "instant" if (d.get("slurm_id") or "").startswith("local") else "cluster"
        if d["state"] == "done":
            d["download"] = f"/api/cat/download/{protocol}"
        return d
    finally:
        conn.close()


@router.get("/download/{protocol}")
async def download(protocol: str):
    conn = _req()
    try:
        row = conn.execute("SELECT state, path, expires FROM requests "
                           "WHERE protocol = ?", (protocol,)).fetchone()
    finally:
        conn.close()
    if row is None:
        raise HTTPException(404, "protocolo não encontrado")
    if row["state"] != "done":
        raise HTTPException(409, f"pedido está em '{row['state']}'")
    if row["expires"] and row["expires"] < datetime.now(timezone.utc).isoformat():
        raise HTTPException(410, f"o link expirou ({RETENTION_DAYS} dias); "
                                 "faça o pedido de novo")
    path = Path(row["path"])
    if not path.exists():
        raise HTTPException(410, "o arquivo já foi removido")
    return FileResponse(str(path), filename=path.name,
                        media_type="application/zip")

# ── o job do cluster falando com o site ─────────────────────────────────────
#
# O nó não alcança o SQLite de pedidos (``/home`` é local de cada máquina), e
# SQLite concorrente por NFS é justamente o que não se faz.  O job escreve
# ``progress.json`` e ``done.json`` no diretório do pedido, no NFS; o site lê
# esses arquivos, atualiza a linha e manda o e-mail.  Assim só um processo
# escreve o banco de pedidos: este.

def _filter_text(f: CatFilter) -> str:
    """O filtro em linhas legíveis, para o corpo do e-mail."""
    d = f.model_dump()
    out = []
    for k, v in d.items():
        if v in (None, [], "", 0) and k not in ("d_min",):
            continue
        out.append(f"  {k}: {v}")
    return "\n".join(out)


def _refresh(conn: sqlite3.Connection, row: sqlite3.Row) -> None:
    """Atualiza uma linha de pedido a partir do que o job deixou no NFS."""
    if row["state"] in ("done", "failed", "expired"):
        _expire_if_due(conn, row)
        return
    out_dir = Path(row["path"]).parent
    prog = out_dir / "progress.json"
    done = out_dir / "done.json"
    payload = None
    for path in (done, prog):
        if path.exists():
            try:
                payload = json.loads(path.read_text())
                break
            except (json.JSONDecodeError, OSError):
                continue
    if payload is None:
        # Nada escrito ainda: o job pode estar na fila, ou pode ter morrido
        # antes de começar.  Quem sabe disso é o Slurm.
        _check_slurm(conn, row)
        return

    state = payload.get("state", "running")
    if state == "running":
        conn.execute("UPDATE requests SET state = 'running', progress = ? "
                     "WHERE protocol = ?", (int(payload.get("done", 0)),
                                            row["protocol"]))
        conn.commit()
        if (row["slurm_id"] or "").startswith("local"):
            _check_local(conn, row)
        return

    now = datetime.now(timezone.utc).isoformat()
    if state == "done":
        # Só quem de fato muda o estado manda o e-mail: com dois workers, a
        # consulta do protocolo e a ronda podiam mandar o mesmo e-mail duas vezes.
        cur = conn.execute("UPDATE requests SET state = 'done', progress = ?, bytes = ?,"
                           " finished = ?, message = ? WHERE protocol = ?"
                           " AND state IN ('queued','running')",
                           (int(payload.get("done", 0)), int(payload.get("bytes", 0)),
                            now, f"{payload.get('n_files', 0)} arquivos, "
                                 f"{payload.get('n_errors', 0)} erros", row["protocol"]))
        conn.commit()
        if cur.rowcount:
            mail.ready(row["email"], row["name"] or "", row["protocol"],
                       int(payload.get("n_files", 0)), int(payload.get("bytes", 0)),
                       row["expires"] or "", BASE_URL)
    else:
        _fail(conn, row, payload.get("message", "falha no job"))


def _check_slurm(conn: sqlite3.Connection, row: sqlite3.Row) -> None:
    """Job que morreu sem escrever nada: o pedido não pode ficar 'na fila'."""
    if (row["slurm_id"] or "").startswith("local"):
        _check_local(conn, row)
        return
    if not row["slurm_id"]:
        return
    age = time.time() - _epoch(row["created"])
    if age < 120:                     # ainda é cedo para desconfiar
        return
    try:
        r = subprocess.run(["sacct", "-j", row["slurm_id"], "-n", "-X", "-o", "State"],
                           capture_output=True, text=True, timeout=15)
    except Exception:
        return
    states = {line.strip().split()[0] for line in r.stdout.splitlines() if line.strip()}
    if states & {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED"}:
        return
    if states and not states & {"COMPLETED"}:
        _fail(conn, row, f"o job terminou em {'/'.join(sorted(states))} sem gerar o arquivo")


def _epoch(iso: str | None) -> float:
    try:
        return datetime.fromisoformat(iso).timestamp()
    except (TypeError, ValueError):
        return 0.0


def _expire_if_due(conn: sqlite3.Connection, row: sqlite3.Row) -> None:
    """Sete dias: o arquivo sai do disco e a linha diz que expirou."""
    if row["state"] != "done" or not row["expires"]:
        return
    if row["expires"] > datetime.now(timezone.utc).isoformat():
        return
    path = Path(row["path"])
    try:
        if path.exists():
            path.unlink()
        for extra_file in ("progress.json", "done.json"):
            f = path.parent / extra_file
            if f.exists():
                f.unlink()
    except OSError:
        pass
    conn.execute("UPDATE requests SET state = 'expired' WHERE protocol = ?",
                 (row["protocol"],))
    conn.commit()


async def watch_requests(interval: float = 20.0) -> None:
    """Acompanha os pedidos abertos, sem ninguém precisar abrir a página.

    É o que faz o segundo e-mail sair na hora em que o job termina, em vez de
    quando alguém consulta o protocolo.  Uma volta custa a leitura de um JSON
    por pedido aberto.
    """
    import asyncio
    while True:
        try:
            conn = _req()
            try:
                rows = conn.execute(
                    "SELECT * FROM requests WHERE state IN "
                    "('queued','running') OR (state = 'done' AND expires < ?)",
                    (datetime.now(timezone.utc).isoformat(),)).fetchall()
                for row in rows:
                    try:
                        _refresh(conn, row)
                    except Exception:            # um pedido ruim não para os outros
                        _log.exception("pedido %s", row["protocol"])
            finally:
                conn.close()
        except Exception:
            _log.exception("ronda de pedidos")
        await asyncio.sleep(interval)
