"""
Pydantic request/response models for NTBuilder Web API.
"""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, Field


class PolarRequest(BaseModel):
    file_id:      Optional[str]   = None
    example:      Optional[str]   = None
    n_max:        int             = Field(30, ge=1, le=60)
    max_diameter: float           = Field(25.0, ge=1.0, le=100.0)
    # Strain filter (applied server-side before returning points)
    strain_max:   Optional[float] = Field(None, ge=0.0, le=100.0)
    # Rolling direction — only meaningful for buckled / Janus structures.
    # Affects which species ends up on the concave (compressed) inner ring
    # and therefore which chiralities develop curvature-induced spurious
    # bonds.  Default False = conventional outward roll (z > 0 atoms move
    # to larger radius).
    roll_inward:  bool            = False


class BuildRequest(BaseModel):
    file_id:     Optional[str]  = None
    example:     Optional[str]  = None
    n:           int            = Field(..., ge=0, le=60)
    m:           int            = Field(..., ge=0, le=60)
    n_repeat:    int            = Field(1, ge=1, le=20)
    vacuum:      float          = Field(10.0, ge=1.0, le=50.0)
    roll_inward: bool           = False
    # Optional per-pair bond cutoffs: {"C-C": 1.8, "B-N": 1.65, ...}
    bond_cutoffs: Optional[dict] = None


class BatchRequest(BaseModel):
    file_id:     Optional[str]         = None
    example:     Optional[str]         = None
    chiralities: list                  # list of [n, m] pairs
    n_repeat:    int                   = Field(1,    ge=1,  le=20)
    vacuum:      float                 = Field(10.0, ge=1.0, le=50.0)
    roll_inward: bool                  = False


class MWNTRequest(BaseModel):
    """Scaled MWNT (k·(n,m) per wall).

    When ``from_job_id`` is supplied, the request is applied on top of an
    existing job (e.g. the user just built a single tube and wants to
    promote it to an MWNT without losing the current state).  Otherwise
    the nanotube is rebuilt from ``file_id`` + ``(n, m)``.
    """
    file_id:            Optional[str]  = None
    example:            Optional[str]  = None
    from_job_id:        Optional[str]  = None
    n:                  int            = Field(..., ge=0, le=60)
    m:                  int            = Field(..., ge=0, le=60)
    n_walls:            int            = Field(2, ge=1, le=10)
    interlayer_spacing: float          = Field(3.4, ge=1.0, le=10.0)
    vacuum:             float          = Field(10.0, ge=1.0, le=50.0)
    roll_inward:        bool           = False


class BundleRequest(BaseModel):
    """Build a periodic bundle from a freshly built nanotube.

    With ``from_job_id`` the bundle wraps the structure of that previous
    job (MWNT, deformed tube, …) instead of rebuilding a single SWNT.
    """
    file_id:      Optional[str] = None
    example:      Optional[str] = None
    from_job_id:  Optional[str] = None
    n:            int           = Field(..., ge=0, le=60)
    m:            int           = Field(..., ge=0, le=60)
    geometry:     str           = Field("hexagonal7",
        description="linear | triangle | square4 | hexagonal7 | grid")
    spacing:      float         = Field(3.4,  ge=0.5, le=30.0)
    vacuum:       float         = Field(10.0, ge=0.0, le=50.0)
    nx:           int           = Field(2,    ge=1,   le=10)
    ny:           int           = Field(2,    ge=1,   le=10)
    n_repeat:     int           = Field(1,    ge=1,   le=999)


class DeformRequest(BaseModel):
    """Apply axial strain and / or torsion to a nanotube.

    Pass ``from_job_id`` to deform an existing structure (MWNT, bundle,
    …) instead of rebuilding the SWNT from scratch.
    """
    file_id:      Optional[str] = None
    example:      Optional[str] = None
    from_job_id:  Optional[str] = None
    n:            int           = Field(..., ge=0, le=60)
    m:            int           = Field(..., ge=0, le=60)
    axial_strain: float         = Field(0.0,  ge=-0.5, le=2.0,
        description="Fractional axial strain (e.g. 0.05 = 5 %).")
    twist_rate:   float         = Field(0.0,  ge=-90.0, le=90.0,
        description="Torsion rate in degrees per Å.")
    radial_strain: float        = Field(0.0,  ge=-0.5, le=1.0)
    vacuum:       float         = Field(10.0, ge=1.0,  le=50.0)
    z_vacuum:     float         = Field(10.0, ge=0.0,  le=100.0,
        description="Z padding added when torsion ≠ 0.")
    n_repeat:     int           = Field(1,    ge=1,    le=999)


class AnalysisRequest(BaseModel):
    """Compute bond, electronic and symmetry analysis for a nanotube."""
    file_id:      Optional[str] = None
    example:      Optional[str] = None
    from_job_id:  Optional[str] = None
    n:            int           = Field(..., ge=0, le=60)
    m:            int           = Field(..., ge=0, le=60)
    vacuum:       float         = Field(10.0, ge=1.0, le=50.0)
    bond_cutoff:  float         = Field(2.0,  ge=0.5, le=5.0)


class MethodsRequest(BaseModel):
    """Generate a Methods-section paragraph for the current nanotube.

    The structure is rebuilt from ``file_id`` (or ``example``) + ``(n, m)``
    so that the lattice parameters used in the text are always consistent
    with the parent 2D lattice.  Optional ``deform_desc`` lets the caller
    inject the deformation description string produced server-side by the
    /api/deform endpoint.
    """
    file_id:      Optional[str] = None
    example:      Optional[str] = None
    from_job_id:  Optional[str] = None
    n:            int           = Field(..., ge=0, le=60)
    m:            int           = Field(..., ge=0, le=60)
    vacuum:       float         = Field(10.0, ge=1.0, le=50.0)
    deform_desc:  Optional[str] = None
    n_walls:      int           = Field(1,   ge=1, le=10)
    wall_info:    Optional[str] = None
    cite_key:     str           = "Pereira2026"
    software:     str           = "NTBuilder"
    version:      str           = "1.1"


class DFTInputRequest(BaseModel):
    """Generate DFT input files (VASP / QE / CP2K / SIESTA) for a nanotube.

    The ``code`` selector returns:
        "vasp"   → {"incar": ..., "kpoints": ..., "poscar": ...}
        "qe"     → {"input":  ...}
        "cp2k"   → {"input":  ...}
        "siesta" → {"input":  ...}    (re-uses core.exporters.write_siesta)
    """
    file_id:      Optional[str] = None
    example:      Optional[str] = None
    from_job_id:  Optional[str] = None
    n:            int           = Field(..., ge=0, le=60)
    m:            int           = Field(..., ge=0, le=60)
    vacuum:       float         = Field(10.0, ge=1.0, le=50.0)
    code:         str           = Field("vasp",
        description="vasp | qe | cp2k | siesta")
