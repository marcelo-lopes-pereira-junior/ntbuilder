"""
core/mwnt.py
------------
Multi-walled nanotube builder.

Given an innermost ChiralityResult and a number of walls N, automatically
finds the (n,m) for each additional shell that minimises the deviation from
the ideal interlayer spacing, then merges all walls into a single
NanotubeStructure ready for export.

Algorithm
---------
1.  Inner wall diameter d₀ is known from the user's (n,m) selection.
2.  For wall k, the target diameter is d_target = d₀ + 2k · spacing.
3.  We scan all (n,m) with diameter near d_target and pick the best match.
4.  Each wall is built independently (vacuum=0) so atoms are centred at
    the tube axis origin (0,0) before the box shift.
5.  Walls are re-centred, Z-tiled to match the inner wall's period, and
    merged into one NanotubeStructure with the final vacuum box.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .io import LatticeStructure
from .chirality import ChiralityResult, compute_chirality, scan_chirality
from .builder import NanotubeStructure, build_nanotube


# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class WallInfo:
    """Metadata for one wall of a MWNT."""
    index:              int     # 0 = innermost
    n:                  int
    m:                  int
    diameter:           float   # Å — actual diameter
    target_diameter:    float   # Å — ideal diameter from spacing rule
    diameter_error:     float   # Å — |actual − target|
    n_atoms_per_cell:   int     # atoms in the unit cell of this wall
    strain:             float   # periodicity strain %
    z_repetitions:      int     # how many times this wall was tiled in Z


@dataclass
class MWNTResult:
    """Complete result of build_mwnt()."""
    nanotube:           NanotubeStructure   # merged, export-ready
    walls:              list[WallInfo]
    n_walls:            int
    interlayer_spacing: float               # Å, as requested
    mean_spacing:       float               # Å, actual mean spacing


@dataclass
class WallCandidate:
    """One candidate chirality for an outer MWNT wall, with commensurability metrics."""
    chirality:       ChiralityResult
    target_diameter: float    # Å — ideal diameter we wanted
    delta_spacing:   float    # Å — signed: actual_spacing − requested_spacing
    T_ref:           float    # Å — T-period of the reference wall (inner)
    z_ratio:         float    # T_this / T_ref (float; near an integer = commensurate)
    z_n_rep:         int      # repetitions of the shorter wall to cover the longer
    z_strain_pct:    float    # % strain applied to whichever wall is shorter
    score:           float    # combined score (lower = better match on both axes)


# ─────────────────────────────────────────────────────────────────────────────
# Shell-finding helpers
# ─────────────────────────────────────────────────────────────────────────────

def find_shell_chirality(
    structure: LatticeStructure,
    target_diameter: float,
    search_window: float = 0.5,    # fraction of target_diameter
    n_max: int = 120,
) -> ChiralityResult | None:
    """
    Return the (n,m) ChiralityResult whose diameter is closest to
    *target_diameter*.

    Parameters
    ----------
    structure        : parent 2D lattice
    target_diameter  : ideal outer diameter in Å
    search_window    : scan range is [d*(1-w), d*(1+w)]
    n_max            : upper bound for index scanning
    """
    d_lo = target_diameter * (1.0 - search_window)
    d_hi = target_diameter * (1.0 + search_window)

    candidates = scan_chirality(
        structure,
        n_max=n_max,
        max_diameter=d_hi,
        unique_only=False,    # include all sectors to find the best physical match
        search_limit=300,
    )
    candidates = [c for c in candidates if c.diameter >= d_lo]
    if not candidates:
        return None
    return min(candidates, key=lambda c: abs(c.diameter - target_diameter))


def find_shell_candidates(
    structure:        LatticeStructure,
    target_diameter:  float,
    T_ref:            float,
    requested_spacing: float,
    n_candidates:     int   = 10,
    search_window:    float = 0.6,
    n_max:            int   = 150,
) -> list[WallCandidate]:
    """
    Return up to *n_candidates* WallCandidate objects for one MWNT outer wall.

    Each candidate scores the trade-off between:
      • spacing accuracy  — how close the actual interlayer gap is to the
                           requested value  (lower |Δspacing| = better)
      • Z commensurability — how much axial strain must be applied to make
                            the two wall T-periods fit the same supercell
                            (lower z_strain_pct = better)

    Candidates are returned sorted by combined score (equal weights).

    Parameters
    ----------
    structure         : 2D parent lattice
    target_diameter   : ideal outer-wall diameter (= inner_d + 2·k·spacing)
    T_ref             : T-vector length of the reference (inner/previous) wall
    requested_spacing : the spacing the user asked for (used in score normalisation)
    n_candidates      : how many top results to return
    search_window     : diameter scan window as fraction of target_diameter
    n_max             : upper index limit for scan_chirality
    """
    d_lo = target_diameter * (1.0 - search_window)
    d_hi = target_diameter * (1.0 + search_window)

    raw = scan_chirality(
        structure,
        n_max=n_max,
        max_diameter=d_hi,
        unique_only=False,
        search_limit=300,
    )
    raw = [c for c in raw if c.diameter >= d_lo]
    if not raw:
        return []

    results: list[WallCandidate] = []
    for c in raw:
        actual_spacing = (c.diameter - (target_diameter - 2.0 * requested_spacing)) / 2.0
        delta_spacing  = actual_spacing - requested_spacing

        # Z-commensurability: which wall is shorter?
        T_this = c.T_norm
        if T_this <= 0:
            continue

        T_long  = max(T_this, T_ref)
        T_short = min(T_this, T_ref)
        n_rep   = max(1, round(T_long / T_short))
        # strain on the shorter wall to fit the longer
        z_strain_pct = abs(n_rep * T_short - T_long) / T_long * 100.0

        # Combined score (equal weights, both normalised to ~[0,1])
        spacing_score = abs(delta_spacing) / max(requested_spacing, 1e-3)
        z_score       = z_strain_pct / 100.0
        score         = 0.5 * spacing_score + 0.5 * z_score

        results.append(WallCandidate(
            chirality       = c,
            target_diameter = target_diameter,
            delta_spacing   = delta_spacing,
            T_ref           = T_ref,
            z_ratio         = T_this / T_ref,
            z_n_rep         = n_rep,
            z_strain_pct    = z_strain_pct,
            score           = score,
        ))

    # Sort by combined score; deduplicate (n,m) — keep best per pair
    seen: set[tuple[int,int]] = set()
    unique: list[WallCandidate] = []
    for cand in sorted(results, key=lambda x: x.score):
        key = (cand.chirality.n, cand.chirality.m)
        if key not in seen:
            seen.add(key)
            unique.append(cand)
        if len(unique) >= n_candidates:
            break

    return unique


# ─────────────────────────────────────────────────────────────────────────────
# Main builder
# ─────────────────────────────────────────────────────────────────────────────

def build_mwnt_with_choices(
    structure:          LatticeStructure,
    chiralities:        list[ChiralityResult],
    interlayer_spacing: float = 3.4,
    vacuum:             float = 10.0,
    roll_inward:        bool  = False,
) -> MWNTResult:
    """
    Build a MWNT from a *user-selected* list of wall chiralities.

    Unlike :func:`build_mwnt` this function does not search for outer walls
    automatically — the caller supplies the full list (index 0 = innermost).
    All other parameters are the same.
    """
    if len(chiralities) < 1:
        raise ValueError("chiralities list must contain at least one entry.")

    raw_walls = [
        build_nanotube(structure, ch, vacuum=0.0, roll_inward=roll_inward)
        for ch in chiralities
    ]
    return _merge_walls(
        raw_walls, chiralities,
        inner_chirality    = chiralities[0],
        interlayer_spacing = interlayer_spacing,
        vacuum             = vacuum,
    )


def _merge_walls(
    raw_walls:          list[NanotubeStructure],
    chiralities:        list[ChiralityResult],
    inner_chirality:    ChiralityResult,
    interlayer_spacing: float,
    vacuum:             float,
) -> MWNTResult:
    """Internal: merge pre-built per-wall NanotubeStructures into one MWNT."""
    n_walls  = len(raw_walls)
    Lz_walls = [float(nt.box[2]) for nt in raw_walls]
    Lz_ref   = max(Lz_walls)

    all_symbols: list[str]         = []
    all_coords:  list[np.ndarray]  = []
    wall_infos:  list[WallInfo]    = []

    for k, (nt, ch) in enumerate(zip(raw_walls, chiralities)):
        cx_k     = float(nt.box[0]) / 2.0
        cy_k     = float(nt.box[1]) / 2.0
        coords_k = nt.coords.copy()
        coords_k[:, 0] -= cx_k
        coords_k[:, 1] -= cy_k

        Lz_k = Lz_walls[k]
        if Lz_k >= Lz_ref - 1e-6:
            syms_k = list(nt.symbols)
            n_rep  = 1
        else:
            n_rep    = max(1, round(Lz_ref / Lz_k))
            tiles    = [coords_k + np.array([0.0, 0.0, i * Lz_k])
                        for i in range(n_rep)]
            coords_k = np.vstack(tiles)
            syms_k   = list(nt.symbols) * n_rep
            actual_z = n_rep * Lz_k
            if abs(actual_z - Lz_ref) > 1e-6:
                coords_k[:, 2] *= Lz_ref / actual_z

        all_symbols.extend(syms_k)
        all_coords.append(coords_k)

        target_d = inner_chirality.diameter + 2.0 * k * interlayer_spacing
        wall_infos.append(WallInfo(
            index=k,
            n=ch.n, m=ch.m,
            diameter=ch.diameter,
            target_diameter=target_d,
            diameter_error=abs(ch.diameter - target_d),
            n_atoms_per_cell=ch.n_atoms,
            strain=ch.strain,
            z_repetitions=n_rep,
        ))

    merged_coords = np.vstack(all_coords)
    outer_r  = chiralities[-1].diameter / 2.0
    box_xy   = outer_r * 2.0 + vacuum
    merged_coords[:, 0] += box_xy / 2.0
    merged_coords[:, 1] += box_xy / 2.0

    merged = NanotubeStructure(
        chirality = inner_chirality,
        symbols   = all_symbols,
        coords    = merged_coords,
        box       = np.array([box_xy, box_xy, Lz_ref]),
        vacuum    = vacuum,
    )

    if n_walls > 1:
        spacings = [
            (chiralities[k].diameter - chiralities[k - 1].diameter) / 2.0
            for k in range(1, n_walls)
        ]
        mean_spacing = float(np.mean(spacings))
    else:
        mean_spacing = interlayer_spacing

    return MWNTResult(
        nanotube           = merged,
        walls              = wall_infos,
        n_walls            = n_walls,
        interlayer_spacing = interlayer_spacing,
        mean_spacing       = mean_spacing,
    )


def build_mwnt(
    structure:           LatticeStructure,
    inner_chirality:     ChiralityResult,
    n_walls:             int   = 2,
    interlayer_spacing:  float = 3.4,
    vacuum:              float = 10.0,
    roll_inward:         bool  = False,
) -> MWNTResult:
    """
    [LEGACY] Build a multi-walled nanotube via free (n,m) search.

    This routine is preserved for reference but **is no longer wired into
    the GUI/CLI** (see :func:`build_mwnt_scaled` for the current default).
    It performs an expensive scan over candidate chiralities to find the
    pair that best matches the requested interlayer spacing for each
    outer wall, which produces incommensurate walls along the axis and
    requires per-wall axial-strain remapping (handled in ``_merge_walls``).

    The mismatched/heterostructure MWNT this approach can in principle
    address is now listed as future work in the manuscript.

    Parameters
    ----------
    structure           : 2D parent lattice (any supported format).
    inner_chirality     : chirality of the innermost wall.
    n_walls             : total number of concentric walls (≥ 1).
    interlayer_spacing  : surface-to-surface gap between consecutive walls
                          in Å (default 3.4 Å, the graphite vdW value).
    vacuum              : XY-plane vacuum padding around the outermost wall.
    roll_inward         : passed through to :func:`build_nanotube`.

    Returns
    -------
    MWNTResult with ``.nanotube`` ready for export and ``.walls`` with
    per-wall metadata.
    """
    if n_walls < 1:
        raise ValueError("n_walls must be ≥ 1")

    # ── Step 1: determine chirality for each wall ─────────────────────────────
    chiralities: list[ChiralityResult] = [inner_chirality]
    for k in range(1, n_walls):
        target_d = inner_chirality.diameter + 2.0 * k * interlayer_spacing
        ch = find_shell_chirality(structure, target_d)
        if ch is None:
            raise ValueError(
                f"Could not find a valid (n,m) for wall {k + 1} "
                f"(target diameter {target_d:.2f} Å). "
                "Try fewer walls or a different interlayer spacing."
            )
        chiralities.append(ch)

    # ── Step 2: build each wall independently with zero vacuum ────────────────
    raw_walls: list[NanotubeStructure] = []
    for ch in chiralities:
        nt = build_nanotube(structure, ch, vacuum=0.0, roll_inward=roll_inward)
        raw_walls.append(nt)

    # ── Step 3 & 4: merge using shared helper ─────────────────────────────────
    return _merge_walls(
        raw_walls, chiralities,
        inner_chirality    = inner_chirality,
        interlayer_spacing = interlayer_spacing,
        vacuum             = vacuum,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Scaled MWNT (current default)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ScaledWallPlan:
    """Plan for one wall of a scaled MWNT, prior to construction."""
    index:           int      # 0 = innermost
    k:               int      # scaling factor; inner wall has k=1
    n:               int      # = k * n_inner
    m:               int      # = k * m_inner
    diameter:        float    # = k * d_inner (Å)
    target_diameter: float    # ideal diameter = d_inner + 2·index·spacing
    actual_spacing:  float    # gap to the previous wall, in Å


def plan_scaled_walls(
    inner_chirality:    ChiralityResult,
    n_walls:            int,
    interlayer_spacing: float = 3.4,
) -> list[ScaledWallPlan]:
    """
    Compute the integer-scaled (n,m) for each wall of a MWNT.

    The chiral angle is preserved by every wall, which guarantees that
    the T-vector of every wall is parallel to the inner one and the
    merged structure is exactly periodic along Z.  Two precursors share
    the same chiral angle iff their (n,m) lie on the same ray from the
    origin in the integer lattice — i.e. they are integer multiples of
    the *primitive direction* ``(n₀, m₀) = (n/g, m/g)`` with
    ``g = gcd(n, m)``.

    Examples
    --------
    * Inner (0, 14)  →  g = 14, primitive (0, 1).  Family:
      (0, 1), (0, 2), (0, 3), …, (0, 14), (0, 15), (0, 16), …
      Every option is a zigzag tube (same chiral angle 90°).
    * Inner (5, 5)   →  g = 5,  primitive (1, 1).  Family:
      (1, 1), (2, 2), …, (5, 5), (6, 6), (7, 7), … — all armchair.
    * Inner (5, 3)   →  g = 1,  primitive (5, 3).  Family:
      (5, 3), (10, 6), (15, 9), … — strict k·(n,m) scaling.

    Algorithm
    ---------
    1. Wall 0 (innermost) is the user-selected (n, m).
    2. For wall i ≥ 1, the target diameter is
       ``d_target = d_primitive · q_i`` where
       ``q_i = (d_inner + 2·i·spacing) / d_primitive``.
       ``q_i`` is rounded to the nearest integer ``k_i``; if monotonicity
       would be violated (``k_i ≤ k_{i-1}``), ``k_i`` is incremented.
    3. Wall i uses ``(k_i · n₀, k_i · m₀)``.

    The realised interlayer spacing is reported in
    :func:`scaled_mwnt_warning` so the user knows the deviation from
    the requested value.

    Parameters
    ----------
    inner_chirality    : chirality of the innermost wall.
    n_walls            : total number of walls (≥ 1).
    interlayer_spacing : requested surface-to-surface gap in Å.

    Returns
    -------
    List of ``ScaledWallPlan`` of length ``n_walls`` (innermost first).
    """
    if n_walls < 1:
        raise ValueError("n_walls must be ≥ 1")
    d_inner = float(inner_chirality.diameter)
    if d_inner <= 0:
        raise ValueError("inner diameter is non-positive")

    n_in, m_in = int(inner_chirality.n), int(inner_chirality.m)
    g  = math.gcd(abs(n_in), abs(m_in))
    if g == 0:
        raise ValueError("inner chirality is the degenerate (0,0) pair")
    n0, m0 = n_in // g, m_in // g                # primitive direction
    # Diameter scales linearly along the primitive ray:
    #   d(k·n₀, k·m₀) = (d_inner / g) · k          [for any 2D lattice]
    # ``d_primitive`` is the diameter the (n₀, m₀) tube would have.
    d_primitive = d_inner / g
    k_inner     = g                              # current wall is the g-th step

    plans: list[ScaledWallPlan] = []
    prev_k = 0
    for i in range(n_walls):
        target_d = d_inner + 2.0 * i * interlayer_spacing
        if i == 0:
            k = k_inner
        else:
            # Real-valued step on the primitive ray that would land on
            # ``target_d`` exactly.
            k_real = target_d / d_primitive
            k      = int(round(k_real))
            # Enforce strict monotonicity along the stack.
            if k <= prev_k:
                k = prev_k + 1
        diameter   = k * d_primitive
        prev_d     = plans[-1].diameter if plans else 0.0
        actual_gap = (diameter - prev_d) / 2.0 if i > 0 else float("nan")
        plans.append(ScaledWallPlan(
            index           = i,
            k               = k,
            n               = k * n0,
            m               = k * m0,
            diameter        = diameter,
            target_diameter = target_d,
            actual_spacing  = actual_gap,
        ))
        prev_k = k
    return plans


def scaled_mwnt_warning(plans: list[ScaledWallPlan],
                        requested_spacing: float) -> str | None:
    """
    Return a user-facing message describing the realised interlayer
    spacings, or ``None`` for a single-wall MWNT.
    """
    if len(plans) < 2:
        return None
    actuals = [p.actual_spacing for p in plans[1:]]
    mean_s  = sum(actuals) / len(actuals)
    max_dev = max(abs(s - requested_spacing) for s in actuals)
    lines = [
        f"Scaled MWNT — {len(plans)} walls",
        f"  Requested interlayer spacing : {requested_spacing:.3f} Å",
        f"  Mean realised spacing        : {mean_s:.3f} Å",
        f"  Max deviation from request   : {max_dev:.3f} Å",
        "",
        "Wall-by-wall realisation:",
    ]
    for p in plans:
        if p.index == 0:
            lines.append(
                f"  inner  k=1  ({p.n},{p.m})  d = {p.diameter:.3f} Å"
            )
        else:
            lines.append(
                f"  +{p.index:<5d} k={p.k}  ({p.n},{p.m})  "
                f"d = {p.diameter:.3f} Å   "
                f"gap = {p.actual_spacing:.3f} Å"
            )
    lines.append(
        "\nNote: walls share the (n,m) direction (k-scaling), so all "
        "T-vectors coincide — the merged MWNT is exactly periodic "
        "along Z and requires no axial-strain remapping."
    )
    return "\n".join(lines)


def build_mwnt_scaled(
    structure:           LatticeStructure,
    inner_chirality:     ChiralityResult,
    n_walls:             int   = 2,
    interlayer_spacing:  float = 3.4,
    vacuum:              float = 10.0,
    roll_inward:         bool  = False,
) -> MWNTResult:
    """
    Build a multi-walled nanotube via integer scaling of the inner (n,m).

    See :func:`plan_scaled_walls` for the algorithm.  In contrast to
    :func:`build_mwnt` this routine:

    * does not search the chirality space, so it is *O(n_walls)* instead
      of *O(n_walls · n_max²)* and runs in milliseconds even for many
      walls;
    * produces walls with strictly identical axial periodicity, so the
      merged supercell is exactly commensurate along Z;
    * may deviate from the requested interlayer spacing by up to half a
      bond length; the realised values are reported via
      :func:`scaled_mwnt_warning`.

    Heterostructure / mismatched MWNTs (different (n,m) directions per
    wall) are intentionally **out of scope** for this routine and are
    listed as future work in the manuscript.

    Parameters mirror :func:`build_mwnt`.
    """
    plans = plan_scaled_walls(
        inner_chirality, n_walls, interlayer_spacing=interlayer_spacing,
    )

    # ── Compute ChiralityResult for each scaled wall ────────────────────────
    chiralities: list[ChiralityResult] = [inner_chirality]
    for p in plans[1:]:
        ch = compute_chirality(p.n, p.m, structure)
        if ch is None:
            raise ValueError(
                f"compute_chirality returned None for ({p.n},{p.m}); "
                f"k={p.k} produced an invalid (0,0) index pair."
            )
        chiralities.append(ch)

    # ── Build each wall with zero vacuum (vacuum applied to merged result) ──
    raw_walls: list[NanotubeStructure] = []
    for ch in chiralities:
        nt = build_nanotube(structure, ch, vacuum=0.0, roll_inward=roll_inward)
        raw_walls.append(nt)

    # ── Merge: share helper, but spacings come from realised diameters ──────
    return _merge_walls(
        raw_walls, chiralities,
        inner_chirality    = inner_chirality,
        interlayer_spacing = interlayer_spacing,
        vacuum             = vacuum,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: wall summary string (for GUI / CLI display)
# ─────────────────────────────────────────────────────────────────────────────

def mwnt_summary(result: MWNTResult) -> str:
    """Return a human-readable summary of wall assignments."""
    lines = [
        f"MWNT — {result.n_walls} walls  "
        f"(requested spacing {result.interlayer_spacing:.2f} Å, "
        f"mean actual {result.mean_spacing:.2f} Å)",
        "",
        f"  {'Wall':<5} {'(n,m)':<10} {'Diameter (Å)':<15} "
        f"{'Error (Å)':<12} {'Atoms/cell':<12} {'Z reps':<8} {'Strain (%)':<10}",
        "  " + "-" * 74,
    ]
    for w in result.walls:
        label = "inner" if w.index == 0 else f"  +{w.index}"
        lines.append(
            f"  {label:<5} ({w.n},{w.m}){'':<6} {w.diameter:<15.4f} "
            f"{w.diameter_error:<12.4f} {w.n_atoms_per_cell:<12} "
            f"{w.z_repetitions:<8} {w.strain:<10.4f}"
        )
    lines.append(f"\n  Total atoms : {result.nanotube.n_atoms:,}")
    lines.append(f"  Box Z       : {result.nanotube.box[2]:.4f} Å  "
                 f"(= longest wall T-vector)")
    return "\n".join(lines)
