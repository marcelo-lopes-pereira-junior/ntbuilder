"""
O filtro do catálogo como SQL — a mesma tradução no site e no cluster.

Vive num módulo só, sem FastAPI e sem Pydantic, porque quem consulta é tanto o
processo do site (para contar e mostrar a amostra) quanto o job que monta o
ZIP num nó, e os dois precisam selecionar exatamente os mesmos tubos.

Unidade de contagem.  Uma linha de ``tubes`` é um (n,m) de um Sistema 2D; ela
vale um tubo quando os dois sentidos de enrolamento dão o mesmo tubo e dois
quando dão tubos diferentes (Janus, penta-grafeno).  ``n_all`` guarda esse
número e ``n_clean`` quantos desses tubos não têm ligação espúria.  Toda
contagem soma uma das duas colunas, nunca conta linhas.
"""
from __future__ import annotations

LATTICES = ("hexagonal", "square", "rectangular",
            "centred rectangular", "oblique")
SOURCES = ("c2db", "2dmatpedia", "jz1c03193", "mc2d", "jarvis", "alexandria", "mxene", "mathub2d")
FORMATS = ("cif", "xyz", "poscar", "lammps")
# Faixas do cubo de agregados.  Um filtro cujos limites caem nas bordas delas é
# contado EXATAMENTE pelo cubo, em milissegundos; por isso os controles da
# página andam nesses valores.
ATOMS_BINS = (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
EPS_BINS = (0.0, 0.1, 0.2, 0.5)
DIN_BINS = (0.0, 5.0, 10.0, 20.0)
T_BINS = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 60.0)

# Medido nos arquivos que a ferramenta escreve.
BYTES_PER_ATOM = {"cif": 46, "xyz": 26, "poscar": 24, "lammps": 38}
ZIP_RATIO = 0.34
# Tempo de construção medido no headnode (13/09/2026, um núcleo, nice 10, CIF):
# 1,2 s para 257 tubos pequenos do C2DB (82 mil átomos), 3,8 s para 300 MXenes
# (351 mil) e 4,4 s para 209 tubos da seleção padrão (388 mil).  O ajuste
# t = 1,1e-5 s por átomo + 1,3e-3 s por tubo reproduz os três com erro < 10 %.
# Os 3 ms por tubo de antes subestimavam de 2 a 7 vezes.
SECONDS_PER_ATOM = 1.1e-5
SECONDS_PER_TUBE = 1.3e-3

DEFAULTS = {
    "sources": None, "lattices": None, "elements_all": None,
    "elements_any": None, "elements_none": None, "species_min": 1, "species_max": 12,
    "d_min": 0.0, "d_max": 80.0, "t_min": 0.0, "t_max": 60.0,
    "din_min": 0.0, "atoms_max": 50_000, "eps_max": 0.5, "exact_only": False,
    "chirality": None, "clean": "clean", "per_material_max": 0,
    "formats": ["cif"], "vacuum": 10.0,
}

# Aquiral: Ch paralelo ou perpendicular a um espelho ou deslizamento do Sistema
# 2D, calculado pelo grupo pontual na carga do banco.
KINDS = {"achiral": 0, "chiral": 1}


class FilterError(ValueError):
    """Filtro que não dá para atender — vira 400 na API."""


def normalise(f: dict) -> dict:
    out = dict(DEFAULTS)
    out.update({k: v for k, v in (f or {}).items() if k in DEFAULTS})
    # Valores de antes dos sentidos virarem tubos próprios.
    if out["clean"] in ("both", "any", "ccw", "cw"):
        out["clean"] = "clean"
    if out["clean"] not in ("clean", "all"):
        raise FilterError(f"clean deve ser clean ou all (veio {out['clean']!r})")
    return out


def clean_only(f: dict) -> bool:
    return normalise(f)["clean"] == "clean"


def tube_columns(f: dict) -> tuple[str, str]:
    """O que somar em ``tubes``: (tubos, átomos)."""
    return (("n_clean", "atoms * n_clean") if clean_only(f)
            else ("n_all", "atoms * n_all"))


def cube_columns(f: dict) -> tuple[str, str]:
    """O que somar em ``summary``: (tubos, átomos)."""
    return (("n_clean", "n_atoms_clean") if clean_only(f)
            else ("n_tubes", "n_atoms"))


def cube_rows_column(f: dict) -> str:
    """Quantas LINHAS (pares (n,m) de um Sistema 2D) o cubo tem por faixa.

    É a unidade da amostra paginada da página, uma linha por (n,m), com os
    sentidos que passam no filtro dentro dela.
    """
    return "n_rows_clean" if clean_only(f) else "n_rows"


def _kinds(f: dict) -> list[int] | None:
    want = f.get("chirality")
    if not want:
        return None
    bad = set(want) - set(KINDS)
    if bad:
        raise FilterError(f"quiralidade desconhecida: {', '.join(sorted(bad))}")
    ks = sorted({KINDS[k] for k in want})
    return None if len(ks) == len(KINDS) else ks


def where_sql(f: dict, ids: list[int] | None = None,
              exclude: list[int] | None = None) -> tuple[str, list]:
    """O filtro como WHERE sobre ``tubes``, com os valores em parâmetros.

    ``ids`` vem de :func:`material_ids`; sem ele, o filtro de material entra
    como subconsulta.
    """
    f = normalise(f)
    sql, p = ["1=1"], []
    for field, allowed, col in (("sources", SOURCES, "source"),
                                ("lattices", LATTICES, "lattice_type")):
        vals = f.get(field)
        if vals:
            bad = [v for v in vals if v not in allowed]
            if bad:
                raise FilterError(f"{field}: valor desconhecido {', '.join(bad)}")
            sql.append(f"{col} IN ({','.join('?' * len(vals))})")
            p += list(vals)
    # Intervalo meio-aberto: é o que faz o cubo (por angstrom inteiro) e o
    # índice darem o mesmo número.
    sql.append("diameter >= ? AND diameter < ?")
    p += [float(f["d_min"]), float(f["d_max"])]
    if float(f["din_min"]) > 0:
        sql.append("d_inner >= ?")
        p.append(float(f["din_min"]))
    sql.append("atoms <= ?")
    p.append(int(f["atoms_max"]))
    if f["exact_only"]:
        sql.append("exact = 1")
    else:
        sql.append("eps <= ?")
        p.append(float(f["eps_max"]))
    if f["clean"] == "clean":
        sql.append("n_clean > 0")
    kinds = _kinds(f)
    if kinds is not None:
        sql.append(f"kind IN ({','.join('?' * len(kinds))})")
        p += kinds
    if ids is not None:
        sql.append(f"material_id IN ({','.join('?' * len(ids))})")
        p += list(ids)
    else:
        clause, mp = _material_clause(f)
        if clause:
            sql.append("material_id IN (SELECT id FROM materials WHERE " + clause + ")")
            p += mp
    return " AND ".join(sql), p


def _material_clause(f: dict) -> tuple[str, list]:
    """O que é propriedade do Sistema 2D: elementos, espécies e espessura."""
    mat, mp = [], []
    for el in (f.get("elements_all") or []):
        mat.append("id IN (SELECT material_id FROM material_elements WHERE element = ?)")
        mp.append(el)
    # "Contém uma classe": cada grupo pede ao menos um elemento dele (metais de
    # transição, halogênios...), e grupos diferentes se somam com E.
    for group in (f.get("elements_any") or []):
        group = [e for e in group if e]
        if group:
            mat.append("id IN (SELECT material_id FROM material_elements WHERE element IN ("
                       + ",".join("?" * len(group)) + "))")
            mp += list(group)
    none = f.get("elements_none") or []
    if none:
        mat.append("id NOT IN (SELECT material_id FROM material_elements "
                   f"WHERE element IN ({','.join('?' * len(none))}))")
        mp += list(none)
    if int(f["species_min"]) > 1 or int(f["species_max"]) < 12:
        mat.append("n_species BETWEEN ? AND ?")
        mp += [int(f["species_min"]), int(f["species_max"])]
    if float(f["t_min"]) > 0 or float(f["t_max"]) < T_BINS[-1]:
        mat.append("thickness >= ? AND thickness < ?")
        mp += [float(f["t_min"]), float(f["t_max"])]
    return " AND ".join(mat), mp


def material_where(f: dict) -> tuple[str, list]:
    """O filtro que cabe na tabela ``materials`` — o que decide quais elementos
    ainda combinam com a seleção na tabela periódica."""
    f = normalise(f)
    sql, p = [], []
    for field, allowed, col in (("sources", SOURCES, "source"),
                                ("lattices", LATTICES, "lattice_type")):
        vals = f.get(field)
        if vals:
            bad = [v for v in vals if v not in allowed]
            if bad:
                raise FilterError(f"{field}: valor desconhecido {', '.join(bad)}")
            sql.append(f"{col} IN ({','.join('?' * len(vals))})")
            p += list(vals)
    clause, mp = _material_clause(f)
    if clause:
        sql.append(clause)
        p += mp
    return (" AND ".join(sql) or "1=1"), p


MATERIAL_IDS_MAX = 4000        # acima disto, a subconsulta é mais rápida


def material_ids(conn, f: dict) -> tuple[list[int] | None, None]:
    """Os Sistemas 2D que o filtro deixa passar, quando são poucos.

    Medido no banco de 10,4 M linhas: "contém Mo e S" (150 sistemas) como lista
    de ids leva 176 ms e como subconsulta 10 s; "sem metais de transição"
    (20 636 sistemas) como subconsulta leva 634 ms e como lista 85 s.
    """
    f = normalise(f)
    clause, params = _material_clause(f)
    if not clause:
        return None, None
    keep = [r[0] for r in conn.execute(
        "SELECT id FROM materials WHERE " + clause + " LIMIT ?",
        params + [MATERIAL_IDS_MAX + 1])]
    return (keep, None) if len(keep) <= MATERIAL_IDS_MAX else (None, None)


def _edge_down(x, edges):
    return max(e for e in edges if e <= x) if x >= edges[0] else edges[0]


def _edge_up(x, edges):
    ups = [e for e in edges if e >= x]
    return ups[0] if ups else edges[-1]


def summary_sql(f: dict) -> tuple[str, list, bool] | None:
    """A mesma pergunta no cubo, ou None se o cubo não a conhece.

    Devolve (sql, params, exato); ``exato`` diz se os limites caem nas bordas
    das faixas, e aí o número é o mesmo do índice.
    """
    f = normalise(f)
    if (f.get("elements_all") or f.get("elements_none") or f.get("elements_any")
            or int(f["species_min"]) > 1 or int(f["species_max"]) < 12
            or int(f["per_material_max"])):
        return None
    d_min, d_max = float(f["d_min"]), float(f["d_max"])
    exato = d_min == int(d_min) and d_max == int(d_max)
    sql = ["d_bin >= ?", "d_bin < ?"]
    p = [int(d_min), int(d_max)]

    t_min, t_max = float(f["t_min"]), float(f["t_max"])
    if t_min > 0 or t_max < T_BINS[-1]:
        exato = exato and t_min in T_BINS and t_max in T_BINS
        sql += ["t_bin >= ?", "t_bin < ?"]
        p += [_edge_down(t_min, T_BINS), _edge_up(t_max, T_BINS)]

    atoms_max = int(f["atoms_max"])
    bin_at = _edge_up(atoms_max, ATOMS_BINS)
    exato = exato and atoms_max == bin_at
    sql.append("atoms_bin <= ?")
    p.append(bin_at)

    if f["exact_only"]:
        sql.append("exact = 1")
    else:
        eps = float(f["eps_max"])
        exato = exato and eps in EPS_BINS
        sql.append("eps_bin <= ?")
        p.append(0 if eps <= 0 else 1 if eps <= 0.1 else 2 if eps <= 0.2 else 5)

    din = float(f["din_min"])
    if din > 0:
        exato = exato and din in DIN_BINS
        sql.append("din_bin >= ?")
        p.append(0 if din < 5 else 5 if din < 10 else 10 if din < 20 else 20)

    kinds = _kinds(f)
    if kinds is not None:
        sql.append(f"kind IN ({','.join('?' * len(kinds))})")
        p += kinds
    for field, col in (("sources", "source"), ("lattices", "lattice_type")):
        if f.get(field):
            sql.append(f"{col} IN ({','.join('?' * len(f[field]))})")
            p += list(f[field])
    return " AND ".join(sql), p, exato


def material_cube_sql(f: dict) -> tuple[str, list]:
    """Os filtros de TUBO sobre ``material_cube``, para a tabela periódica.

    Mesmas faixas do cubo geral: com os controles da página nas bordas delas,
    "este Sistema 2D tem algum tubo que passa?" é respondido exatamente.
    """
    f = normalise(f)
    sql = ["d_bin >= ?", "d_bin < ?", "atoms_bin <= ?"]
    p = [int(float(f["d_min"])), int(float(f["d_max"])),
         _edge_up(int(f["atoms_max"]), ATOMS_BINS)]
    if f["exact_only"]:
        sql.append("exact = 1")
    else:
        eps = float(f["eps_max"])
        sql.append("eps_bin <= ?")
        p.append(0 if eps <= 0 else 1 if eps <= 0.1 else 2 if eps <= 0.2 else 5)
    din = float(f["din_min"])
    if din > 0:
        sql.append("din_bin >= ?")
        p.append(0 if din < 5 else 5 if din < 10 else 10 if din < 20 else 20)
    kinds = _kinds(f)
    if kinds is not None:
        sql.append(f"kind IN ({','.join('?' * len(kinds))})")
        p += kinds
    if f["clean"] == "clean":
        sql.append("n_clean > 0")
    return " AND ".join(sql), p


def estimate(n_tubes: int, n_atoms: int, formats) -> dict:
    """Tamanho e tempo do pedido, dos números medidos."""
    formats = [x for x in (formats or []) if x in BYTES_PER_ATOM] or ["cif"]
    raw = sum(n_atoms * BYTES_PER_ATOM[x] for x in formats)
    return {"n_tubes": n_tubes, "n_atoms": n_atoms, "bytes_raw": raw,
            "bytes_zip": int(raw * ZIP_RATIO),
            "seconds": round((n_atoms * SECONDS_PER_ATOM + n_tubes * SECONDS_PER_TUBE)
                             * len(formats), 1)}


# A espessura vem por subconsulta, sem JOIN: o filtro de where_sql usa nomes sem
# prefixo, e source e lattice_type existem em tubes e em materials.  Com o JOIN,
# um pedido filtrado por Base de Dados 2D falhava com "ambiguous column name".
SELECT_ROWS = (
    "SELECT t.material_id, t.n, t.m, t.diameter, t.d_inner, t.t_norm, t.atoms, "
    "       t.eps, t.exact, t.theta, t.kind, t.senses, t.n_all, t.n_clean, "
    "       t.clean_ccw, t.clean_cw, t.spurious_ccw, t.spurious_cw, "
    "       t.lattice_type, "
    "       (SELECT thickness FROM materials WHERE id = t.material_id) AS thickness "
    "FROM tubes t WHERE ")
