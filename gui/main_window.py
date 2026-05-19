"""
gui/main_window.py
------------------
Main application window.

Layout (horizontal QSplitter):
  ┌──────────────┬────────────────────────┬────────────────────────┐
  │  InputPanel  │      PolarPanel        │     ViewerPanel        │
  │  (fixed)     │   (interactive map)    │  (3D viewer + export)  │
  └──────────────┴────────────────────────┴────────────────────────┘

Signal flow:
  InputPanel  ──structure_loaded──▶  PolarPanel  (redraw map)
  InputPanel  ──structure_loaded──▶  MainWindow  (store struct)
  PolarPanel  ──indices_selected──▶  InputPanel  (set n, m spinboxes)
  InputPanel  ──build_requested ──▶  MainWindow  (run core, emit result)
  MainWindow  ──nanotube_ready  ──▶  ViewerPanel (update 3D + info)
"""

from __future__ import annotations

import traceback
from pathlib import Path

import numpy as np
from PyQt6.QtCore    import Qt
from PyQt6.QtGui     import QPixmap, QAction, QIcon
from PyQt6.QtWidgets import (
    QMainWindow, QSplitter, QWidget, QVBoxLayout,
    QStatusBar, QMessageBox, QLabel, QDialog,
    QHBoxLayout, QPushButton, QFileDialog,
)

from .utils import load_pixmap as _load_pixmap, ScalablePixmapLabel

from core import load_structure, compute_chirality, build_nanotube, LatticeStructure
from core.builder import NanotubeStructure

from .panels.input_panel  import InputPanel
from .panels.polar_panel  import PolarPanel
from .panels.viewer_panel import ViewerPanel

# ── Version info ──────────────────────────────────────────────────────────────
_VERSION   = "1.1.0"
_DATE      = "2026"
_AUTHOR    = "Prof. Dr. Marcelo Lopes Pereira Junior"
_AFFIL     = "University of Brasília (UnB)"
_EMAIL     = "marcelo.lopes@unb.br"
_GITHUB    = "https://github.com/marcelo-lopes-pereira-junior/ntbuilder"
_ASSETS    = Path(__file__).parent.parent / "assets"



class MainWindow(QMainWindow):
    """Top-level application window."""

    # Fraction of window height used for the UnB logo in the status bar
    _UNB_FRAC = 0.048   # ≈ 40 px at 820px tall, ≈ 52 px when maximised at 1080px
    _UNB_MIN  = 20      # px — never smaller than this
    _UNB_MAX  = 70      # px — never larger than this

    def __init__(self):
        super().__init__()
        self.setWindowTitle("NTBuilder")

        # ── Declare ALL instance attributes before any Qt geometry calls ─────
        # showMaximized() / resize() can fire resizeEvent immediately, so every
        # attribute that resizeEvent reads must exist before those calls.
        self._structure:   LatticeStructure | None = None
        self._nanotube:    NanotubeStructure | None = None
        self._unb_lbl:     QLabel  | None = None
        self._unb_src_pix: QPixmap | None = None

        # Window icon
        win_icon = QIcon()
        icon_png = _ASSETS / "ntbuilder_icon.png"
        if icon_png.exists():
            src = QPixmap(str(icon_png))
            if not src.isNull():
                for sz in (16, 32, 48, 64, 128):
                    win_icon.addPixmap(
                        src.scaled(sz, sz,
                                   Qt.AspectRatioMode.KeepAspectRatio,
                                   Qt.TransformationMode.SmoothTransformation)
                    )
        if not win_icon.isNull():
            self.setWindowIcon(win_icon)

        self._build_ui()
        self._connect_signals()
        self.statusBar().showMessage("Ready — load a structure file to begin.")

        # Maximise after the UI is fully built so resizeEvent has valid state
        self.resize(1400, 820)
        self.showMaximized()

    # ─────────────────────────────────────────────────────────────────────────
    # UI construction
    # ─────────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(4)

        self.input_panel  = InputPanel()
        self.polar_panel  = PolarPanel()
        self.viewer_panel = ViewerPanel()

        splitter.addWidget(self.input_panel)
        splitter.addWidget(self.polar_panel)
        splitter.addWidget(self.viewer_panel)

        # Proportions: 22% | 38% | 40%
        splitter.setSizes([240, 420, 440])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)

        self.setCentralWidget(splitter)

        # ── Menu bar ─────────────────────────────────────────────────────────
        mb = self.menuBar()

        # Analysis / Methods / DFT Inputs as top-level menu items so they
        # share the menubar with Help (replacing the old QTabWidget at the
        # top of the right pane).
        self._act_analysis = QAction("Analysis", self)
        self._act_analysis.setEnabled(False)
        self._act_analysis.triggered.connect(self.viewer_panel.open_analysis_dialog)
        mb.addAction(self._act_analysis)

        self._act_methods = QAction("Methods", self)
        self._act_methods.setEnabled(False)
        self._act_methods.triggered.connect(self.viewer_panel.open_methods_dialog)
        mb.addAction(self._act_methods)

        self._act_dft = QAction("DFT Inputs", self)
        self._act_dft.setEnabled(False)
        self._act_dft.triggered.connect(self.viewer_panel.open_dft_inputs_dialog)
        mb.addAction(self._act_dft)

        # Enable the menubar items the first time a nanotube is loaded.
        def _enable_advanced(_loaded: bool) -> None:
            self._act_analysis.setEnabled(True)
            self._act_methods.setEnabled(True)
            self._act_dft.setEnabled(True)
        self.viewer_panel.nanotube_loaded.connect(_enable_advanced)

        help_menu = mb.addMenu("Help")
        act_about = QAction("About NTBuilder…", self)
        act_about.triggered.connect(self._on_about)
        help_menu.addAction(act_about)

        # ── Status bar with UnB logo ──────────────────────────────────────────
        sb = QStatusBar()
        self.setStatusBar(sb)

        # Load UnB logo from PNG once; resizeEvent will scale proportionally.
        _unb_h0 = int(self.height() * self._UNB_FRAC)
        _unb_h0 = max(self._UNB_MIN, min(self._UNB_MAX, _unb_h0))

        _unb_png = _ASSETS / "university_of_brasilia.png"
        _src = QPixmap(str(_unb_png)) if _unb_png.exists() else QPixmap()
        self._unb_src_pix = _src if not _src.isNull() else None

        unb_lbl = QLabel()
        unb_lbl.setToolTip(f"{_AFFIL}\n{_AUTHOR}\n{_EMAIL}")
        if self._unb_src_pix:
            display_pix = self._unb_src_pix.scaledToHeight(
                _unb_h0, Qt.TransformationMode.SmoothTransformation
            )
            unb_lbl.setPixmap(display_pix)
            unb_lbl.setFixedHeight(_unb_h0)
            unb_lbl.setFixedWidth(display_pix.width() + 6)
        else:
            # Text fallback: always readable
            unb_lbl.setText("University of Brasília")
            unb_lbl.setFixedWidth(190)
            unb_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            unb_lbl.setStyleSheet(
                "color: #003F8A; font-weight: bold; font-size: 14px;"
            )
        self._unb_lbl = unb_lbl
        sb.addPermanentWidget(unb_lbl)

    # ─────────────────────────────────────────────────────────────────────────
    # Signal wiring
    # ─────────────────────────────────────────────────────────────────────────

    def _connect_signals(self):
        # File loaded → store structure + update polar map
        self.input_panel.structure_loaded.connect(self._on_structure_loaded)

        # Polar map click → set n, m in input panel
        self.polar_panel.indices_selected.connect(self.input_panel.set_indices)

        # Roll-inward toggle → re-evaluate curvature-induced spurious bond
        # markers on the polar map.  Flipping the rolling direction swaps
        # the concave / convex face of a buckled monolayer (relevant for
        # Janus structures such as MoSSe where the two faces carry
        # different species), so chiralities that were "clean" may become
        # marked with × and vice-versa.
        self.input_panel.roll_direction_changed.connect(
            self.polar_panel.set_roll_inward
        )

        # Build button → run pipeline
        self.input_panel.build_requested.connect(self._on_build_requested)

        # MWNT from viewer Tools panel → run MWNT pipeline
        self.viewer_panel.mwnt_requested.connect(self._on_mwnt_requested)

        # Export requests from viewer panel (fmt, path, vacuum, n_rep)
        self.viewer_panel.export_requested.connect(self._on_export_requested)

        # Systematize request
        self.viewer_panel.systematize_requested.connect(self._on_systematize)

    # ─────────────────────────────────────────────────────────────────────────
    # Resize — keep UnB logo proportional to window height
    # ─────────────────────────────────────────────────────────────────────────

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not hasattr(self, "_unb_lbl") or self._unb_lbl is None or self._unb_src_pix is None:
            return
        new_h = int(self.height() * self._UNB_FRAC)
        new_h = max(self._UNB_MIN, min(self._UNB_MAX, new_h))
        pix = self._unb_src_pix.scaledToHeight(
            new_h, Qt.TransformationMode.SmoothTransformation
        )
        self._unb_lbl.setPixmap(pix)
        self._unb_lbl.setFixedHeight(new_h)
        self._unb_lbl.setFixedWidth(pix.width() + 6)

    # ─────────────────────────────────────────────────────────────────────────
    # Slots
    # ─────────────────────────────────────────────────────────────────────────

    def _on_about(self):
        _AboutDialog(self).exec()

    def _on_structure_loaded(self, structure: LatticeStructure):
        self._structure = structure
        self.polar_panel.set_structure(structure)
        # Hand the parent-lattice metadata to the viewer so AnalysisDialog
        # can apply the zone-folding rule only when appropriate (graphene).
        self.viewer_panel.set_parent_lattice(structure)
        self.input_panel.show_logo()   # reveal sidebar icon now that a file is loaded
        self.statusBar().showMessage(
            f"Loaded: {structure}  — click the polar map or set n, m and press Build."
        )

    def _on_build_requested(self, n: int, m: int, vacuum: float,
                             roll_inward: bool = False):
        if self._structure is None:
            QMessageBox.warning(self, "No structure",
                                "Please load a structure file first.")
            return

        try:
            self.statusBar().showMessage(f"Computing chirality for ({n},{m})…")
            chirality = compute_chirality(n, m, self._structure)
            if chirality is None:
                QMessageBox.warning(self, "Invalid indices",
                                    "n = m = 0 is not a valid nanotube.")
                return

            # ── Large-nanotube guard (> 500 000 atoms) ────────────────────────
            # The 3D viewer cannot handle tubes this large.  Offer to skip the
            # viewer and export directly to a file instead.
            _LARGE_ATOMS = 500_000
            if chirality.n_atoms > _LARGE_ATOMS:
                msg = (
                    f"The ({n},{m}) nanotube has <b>{chirality.n_atoms:,} atoms</b> "
                    f"per unit cell (tube length ≈ {chirality.T_norm:.0f} Å)."
                    f"<br><br>"
                    f"The 3D viewer cannot display a structure this large — it "
                    f"would likely freeze or crash the application."
                    f"<br><br>"
                    f"<b>Would you like to build and export it directly to a file?</b><br>"
                    f"<small>You can choose the format (VASP, QE, LAMMPS, PDB) "
                    f"in the next step.  The nanotube will <i>not</i> be shown in "
                    f"the 3D viewer.</small>"
                )
                reply = QMessageBox.question(
                    self, "Very large nanotube — export only?",
                    msg,
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if reply != QMessageBox.StandardButton.Yes:
                    self.statusBar().showMessage("Build cancelled.")
                    return
                self._build_and_export_large(chirality, vacuum, roll_inward)
                return
            # ── Normal path ───────────────────────────────────────────────────

            self.statusBar().showMessage(f"Building ({n},{m}) nanotube…")
            self._nanotube = build_nanotube(
                self._structure, chirality,
                vacuum=vacuum, roll_inward=roll_inward,
            )
            self.viewer_panel.set_nanotube(self._nanotube)
            self.polar_panel.highlight_point(n, m)

            direction = " (inward)" if roll_inward else ""
            self.statusBar().showMessage(
                f"Built ({n},{m}){direction}  D = {chirality.diameter:.3f} Å  "
                f"atoms = {chirality.n_atoms}  "
                f"strain = {chirality.strain:.4f}%"
            )

            # Bond validation for buckled/layered structures
            if self._structure.has_buckling:
                self._check_spurious_bonds()

        except Exception as exc:
            QMessageBox.critical(self, "Build error",
                                 f"{type(exc).__name__}: {exc}\n\n"
                                 + traceback.format_exc())

    def _on_mwnt_requested(
        self,
        n: int, m: int,
        vacuum: float, roll_inward: bool,
        n_walls: int, spacing: float,
    ):
        """
        Build a multi-walled nanotube via integer scaling of the inner (n,m).

        Each outer wall uses ``(k·n, k·m)`` with ``k`` chosen as the integer
        that best approximates the requested interlayer spacing.  All walls
        share the same T-vector so the resulting MWNT is exactly periodic
        along Z and requires no axial-strain remapping.  The realised
        spacings are reported to the user in a follow-up dialog.

        Heterostructure / mismatched-wall MWNTs (different (n,m) per shell)
        are out of scope for this routine and listed as future work.
        """
        if self._structure is None:
            QMessageBox.warning(self, "No structure",
                                "Please load a structure file first.")
            return
        try:
            from core import compute_chirality
            from core.mwnt import (
                plan_scaled_walls, build_mwnt_scaled, scaled_mwnt_warning,
            )

            inner_ch = compute_chirality(n, m, self._structure)
            if inner_ch is None:
                QMessageBox.warning(self, "Invalid indices",
                                    "n = m = 0 is not valid.")
                return

            # ── Plan first so we can show the user the scaling before building.
            plans = plan_scaled_walls(
                inner_ch, n_walls, interlayer_spacing=spacing,
            )

            self.statusBar().showMessage(
                f"Building scaled MWNT — {n_walls} walls "
                f"(k = {', '.join(str(p.k) for p in plans)})…"
            )

            result = build_mwnt_scaled(
                self._structure,
                inner_chirality    = inner_ch,
                n_walls            = n_walls,
                interlayer_spacing = spacing,
                vacuum             = vacuum,
                roll_inward        = roll_inward,
            )
            self._nanotube = result.nanotube
            self.viewer_panel.set_nanotube(self._nanotube, is_raw=False)

            d_outer = result.walls[-1].diameter if result.walls else 0.0
            self.statusBar().showMessage(
                f"MWNT ({n},{m}) — {result.n_walls} walls  "
                f"D_outer = {d_outer:.3f} Å  "
                f"atoms = {result.nanotube.n_atoms:,}  "
                f"mean spacing = {result.mean_spacing:.2f} Å"
            )

            # Show the realised wall plan so the user sees the spacing
            # deviation explicitly.
            warn = scaled_mwnt_warning(plans, requested_spacing=spacing)
            if warn is not None:
                QMessageBox.information(self, "MWNT — wall plan", warn)

        except Exception as exc:
            QMessageBox.critical(
                self, "MWNT build error",
                f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}"
            )

    def _build_and_export_large(self, chirality, vacuum: float,
                                 roll_inward: bool):
        """Build a very large nanotube and export it directly — skip the viewer."""
        from core import export as _export

        # Ask the user for the output file path + format
        fmt_filter = (
            "VASP POSCAR (*.vasp *.poscar);;"
            "Quantum ESPRESSO (*.in *.pwi);;"
            "LAMMPS data (*.lammps *.data);;"
            "PDB (*.pdb);;"
            "All files (*)"
        )
        path, selected = QFileDialog.getSaveFileName(
            self, "Export large nanotube", "", fmt_filter
        )
        if not path:
            self.statusBar().showMessage("Export cancelled.")
            return

        # Derive format + default extension from the selected filter
        ext = Path(path).suffix.lower()
        if "VASP" in selected or ext in (".vasp", ".poscar"):
            fmt, default_ext = ".poscar", ".poscar"
        elif "Quantum" in selected or ext in (".in", ".pwi"):
            fmt, default_ext = ".in", ".in"
        elif "LAMMPS" in selected or ext in (".lammps", ".data"):
            fmt, default_ext = ".lammps", ".lammps"
        elif ext == ".pdb" or "PDB" in selected:
            fmt, default_ext = ".pdb", ".pdb"
        elif ext == ".xyz":
            fmt, default_ext = ".xyz", ".xyz"
        else:
            fmt, default_ext = ext or ".poscar", ext or ".poscar"

        # Auto-append extension if the user omitted it
        if not Path(path).suffix:
            path = path + default_ext

        try:
            self.statusBar().showMessage(
                f"Building ({chirality.n},{chirality.m}) — "
                f"{chirality.n_atoms:,} atoms — this may take a while…"
            )
            nt = build_nanotube(
                self._structure, chirality,
                vacuum=vacuum, roll_inward=roll_inward,
            )
            out = _export(nt, path, fmt=fmt)
            self.statusBar().showMessage(
                f"Exported ({chirality.n},{chirality.m}) → {out}  "
                f"({chirality.n_atoms:,} atoms)"
            )
            QMessageBox.information(
                self, "Export complete",
                f"({chirality.n},{chirality.m}) nanotube written to:\n{out}\n\n"
                f"Atoms: {chirality.n_atoms:,}  |  "
                f"D = {chirality.diameter:.3f} Å  |  "
                f"L = {chirality.T_norm:.3f} Å"
            )
        except Exception as exc:
            QMessageBox.critical(self, "Export error",
                                 f"{type(exc).__name__}: {exc}\n\n"
                                 + traceback.format_exc())

    def _check_spurious_bonds(self):
        """Warn if the nanotube contains bonds absent from the flat 2D structure."""
        try:
            from core.builder import check_spurious_bonds
            spurious = check_spurious_bonds(
                self._structure, self._nanotube,
                settings=self.viewer_panel.bond_settings,
            )
            if not spurious:
                return

            pairs = sorted(
                "–".join(sorted(p)) for p in spurious
            )
            msg = (
                "<b>⚠ Spurious bonds detected</b><br><br>"
                "The following bond types appear in the 3D nanotube but do "
                "<i>not</i> exist in the flat 2D structure:<br><br>"
                f"&nbsp;&nbsp;<b>{',  '.join(pairs)}</b><br><br>"
                "This usually means the nanotube diameter is too small for "
                "this material — atoms from opposite sides of the tube, or "
                "from different layers, are getting close enough to appear "
                "bonded. Consider using larger (n,&nbsp;m) indices."
            )
            QMessageBox.warning(self, "Bond validation warning", msg)
        except Exception:
            pass   # don't crash the build if validation fails

    def _on_export_requested(self, fmt: str, path: str,
                             vacuum: float, n_rep: int, tube_axis: str = "Z"):
        if self._nanotube is None:
            QMessageBox.warning(self, "Nothing to export",
                                "Build a nanotube first.")
            return
        try:
            from core import export
            nt = _apply_vacuum_and_reps(self._nanotube, vacuum, n_rep,
                                        tube_axis=tube_axis)
            out = export(nt, path, fmt=fmt)
            self.statusBar().showMessage(f"Exported → {out}")
        except Exception as exc:
            QMessageBox.critical(self, "Export error", str(exc))

    def _on_systematize(self, params: dict):
        """Batch-generate a range of (n,m) tubes and export them."""
        if self._structure is None:
            QMessageBox.warning(self, "No structure", "Load a file first.")
            return
        try:
            from core import scan_chirality, build_nanotube, export
            from pathlib import Path

            out_dir  = Path(params["output_dir"])
            fmt      = params["format"]
            vacuum   = params["vacuum"]
            max_diam = params["max_diameter"]
            max_atm  = params.get("max_atoms")

            # Auto-compute per-axis index limits from D_max (same as polar panel).
            import math
            a = max(self._structure.a, 0.5)
            b = max(self._structure.b, 0.5)
            n_max = max(1, math.ceil(max_diam * math.pi / a))
            m_max = max(1, math.ceil(max_diam * math.pi / b))

            results = scan_chirality(
                self._structure,
                n_max=n_max,
                m_max=m_max,
                max_diameter=max_diam,
                max_atoms=max_atm,
            )

            # Apply minimum diameter filter
            min_diam = params.get("min_diameter", 0.0)
            if min_diam > 0:
                results = [r for r in results if r.diameter >= min_diam]

            self.statusBar().showMessage(
                f"Systematizing {len(results)} nanotubes…"
            )
            for ch in results:
                nt = build_nanotube(self._structure, ch, vacuum=vacuum)
                name = f"NT_n{ch.n}_m{ch.m}{fmt}"
                export(nt, out_dir / name, fmt=fmt)

            self.statusBar().showMessage(
                f"Done — {len(results)} files written to {out_dir}"
            )
            QMessageBox.information(
                self, "Systematization complete",
                f"{len(results)} nanotube files written to:\n{out_dir}"
            )
        except Exception as exc:
            QMessageBox.critical(self, "Systematization error",
                                 f"{exc}\n\n{traceback.format_exc()}")


# ─────────────────────────────────────────────────────────────────────────────
# About dialog
# ─────────────────────────────────────────────────────────────────────────────

class _AboutDialog(QDialog):
    """Credits / about dialog."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About NTBuilder")
        # Larger fixed size to avoid clipping the logo + author block
        # and to give the description label room to wrap onto two lines.
        self.setFixedSize(460, 470)
        self.setWindowFlags(
            self.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint
        )

        # ── Force solid, opaque background ───────────────────────────────────
        # Without this the dialog inherits a translucent paint from the global
        # stylesheet on some platforms, so the content panels show through.
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setAutoFillBackground(True)
        self.setStyleSheet(
            "QDialog { background-color: #F4F6FB; }"
            "QLabel  { background-color: transparent; color: #1A3A6B; }"
            "QPushButton {"
            "  background-color: #2851A3; color: white;"
            "  border: none; border-radius: 4px; padding: 6px 16px;"
            "}"
            "QPushButton:hover { background-color: #4A7FD4; }"
            "QPushButton:pressed { background-color: #1A3A6B; }"
        )

        root = QVBoxLayout(self)
        root.setSpacing(10)
        root.setContentsMargins(24, 20, 24, 20)

        # ── NTBuilder logo ───────────────────────────────────────────────────
        _logo_w, _logo_h = 310, 177   # display px
        logo_pix = _load_pixmap(
            _ASSETS / "ntbuilder_complete.png",
            _logo_w * 2, _logo_h * 2,   # load at 2× then smooth-scale
        )
        logo_lbl = QLabel()
        logo_lbl.setFixedSize(_logo_w, _logo_h)
        logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        if logo_pix:
            logo_lbl.setPixmap(
                logo_pix.scaled(
                    _logo_w, _logo_h,
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
        else:
            logo_lbl.setText(
                "<span style='font-size:26px; font-weight:bold; color:#003F8A;'>NT</span>"
                "<span style='font-size:22px; color:#003F8A;'>Builder</span>"
            )
            logo_lbl.setTextFormat(Qt.TextFormat.RichText)
        logo_row = QHBoxLayout()
        logo_row.addStretch()
        logo_row.addWidget(logo_lbl)
        logo_row.addStretch()
        root.addLayout(logo_row)

        # ── Version line ─────────────────────────────────────────────────────

        ver_lbl = QLabel(f"Version {_VERSION}  ·  {_DATE}")
        ver_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ver_lbl.setStyleSheet("color: #666666; font-size: 11px;")
        root.addWidget(ver_lbl)

        desc_lbl = QLabel(
            "Nanotube structure generator from arbitrary 2D crystal inputs. "
            "Supports CIF, PDB, XYZ, POSCAR, XSF, LAMMPS and QE input as readers; "
            "exports to VASP, QE, LAMMPS, CP2K, SIESTA, XSF, CIF, XYZ and PDB. "
            "Includes multi-walled tubes, bundles, axial strain and torsion."
        )
        desc_lbl.setWordWrap(True)
        desc_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        desc_lbl.setStyleSheet(
            "font-size: 10px; color: #444444; background-color: transparent;"
        )
        root.addWidget(desc_lbl)

        # ── Author block ─────────────────────────────────────────────────────
        info = QLabel(
            f"<b>{_AUTHOR}</b><br>"
            f"{_AFFIL}<br>"
            f"<a href='mailto:{_EMAIL}'>{_EMAIL}</a><br>"
            f"<a href='{_GITHUB}'>{_GITHUB}</a>"
        )
        info.setOpenExternalLinks(True)
        info.setTextFormat(Qt.TextFormat.RichText)
        info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        info.setStyleSheet("font-size: 10px; line-height: 1.6;")
        root.addWidget(info)

        root.addStretch()

        # ── Close button ─────────────────────────────────────────────────────
        close_btn = QPushButton("Close")
        close_btn.setFixedWidth(90)
        close_btn.clicked.connect(self.accept)
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        btn_row.addStretch()
        root.addLayout(btn_row)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _apply_vacuum_and_reps(
    nt: "NanotubeStructure",
    vacuum: float,
    n_rep: int,
    tube_axis: str = "Z",
) -> "NanotubeStructure":
    """
    Return a NanotubeStructure ready for export, applying:
      • vacuum padding around the cross-section
      • n_rep replications along the tube axis
      • optional axis permutation (tube_axis = "Z" | "X" | "Y")

    The nanotube is always built along Z internally.  When the user requests
    a different longitudinal axis the coordinate columns and box dimensions
    are permuted accordingly:

      Z (default): (x, y, z) → (x, y, z)   box = [D+vac, D+vac, L]
      X:           (x, y, z) → (z, x, y)   box = [L, D+vac, D+vac]
      Y:           (x, y, z) → (y, z, x)   box = [D+vac, L, D+vac]
    """
    from core.builder import NanotubeStructure

    box_xy     = float(nt.diameter) + vacuum
    old_centre = nt.box[0] / 2.0
    new_centre = box_xy    / 2.0
    delta_xy   = new_centre - old_centre

    # Shift atoms to new cross-section centre (still in Z-along frame)
    coords = nt.coords.copy()
    coords[:, 0] += delta_xy
    coords[:, 1] += delta_xy

    if n_rep > 1:
        L    = float(nt.length)
        reps = [coords + np.array([0.0, 0.0, i * L]) for i in range(n_rep)]
        coords  = np.vstack(reps)
        symbols = list(nt.symbols) * n_rep
        box_z   = L * n_rep
    else:
        symbols = list(nt.symbols)
        box_z   = float(nt.length)

    # Permute columns so the requested axis becomes longitudinal
    axis = tube_axis.upper()
    if axis == "X":
        # tube along X: new (x, y, z) = old (z, x, y)
        coords = coords[:, [2, 0, 1]]
        box    = np.array([box_z, box_xy, box_xy])
    elif axis == "Y":
        # tube along Y: new (x, y, z) = old (y, z, x)
        coords = coords[:, [1, 2, 0]]
        box    = np.array([box_xy, box_z, box_xy])
    else:
        # Z (default) — no permutation
        box    = np.array([box_xy, box_xy, box_z])

    return NanotubeStructure(
        chirality = nt.chirality,
        symbols   = symbols,
        coords    = coords,
        box       = box,
        vacuum    = vacuum,
    )
