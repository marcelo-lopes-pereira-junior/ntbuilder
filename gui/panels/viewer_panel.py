"""
gui/panels/viewer_panel.py
--------------------------
Right panel: 3D nanotube viewer + display options + info + export/systematize.

3D rendering uses pyqtgraph.opengl (hardware-accelerated).
Falls back to a matplotlib 3D canvas if OpenGL is unavailable.

Display strategy: bond-only stick model (bicolour VESTA convention).
No spheres are drawn; structure is inferred from connectivity.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PyQt6.QtCore    import pyqtSignal, Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QPushButton, QFileDialog, QComboBox,
    QFormLayout, QSizePolicy, QFrame, QDialog,
    QDialogButtonBox, QSpinBox, QDoubleSpinBox,
    QCheckBox, QLineEdit, QTabWidget, QGridLayout,
)
from PyQt6.QtGui import QFont

from core.builder       import NanotubeStructure
from core.connectivity  import BondSettings

# ─── Attempt to load pyqtgraph OpenGL ────────────────────────────────────────
try:
    import pyqtgraph.opengl as gl
    _HAS_GL = True
except Exception:
    _HAS_GL = False


class ViewerPanel(QWidget):
    """3D viewer + display options + export controls."""

    # fmt, path, vacuum, n_rep, tube_axis
    export_requested      = pyqtSignal(str, str, float, int, str)
    systematize_requested = pyqtSignal(dict)
    # n, m, vacuum, roll_inward, n_walls, spacing  — emitted by Multi-Wall tool
    mwnt_requested        = pyqtSignal(int, int, float, bool, int, float)
    # True the first time a nanotube becomes available (used by main_window
    # to enable the Analysis / Methods / DFT Inputs menubar items).
    nanotube_loaded       = pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._nanotube:     NanotubeStructure | None = None
        self._nanotube_raw: NanotubeStructure | None = None   # pre-deformation copy
        # Parent 2D lattice metadata — needed by AnalysisDialog to decide
        # whether the zone-folding electronic-character rule applies
        # (graphene only).  Set externally via set_parent_lattice().
        self._lattice_type:        str = "hexagonal"
        self._lattice_n_atoms_uc:  int | None   = None
        self._lattice_a:           float | None = None
        # Transformation labels associated with the current nanotube — used
        # to warn the user when chaining operations that don't compose
        # well (e.g. MWNT on top of a bundle).  Each undo/redo step keeps
        # its own label so the undo stack restores the correct context.
        self._transform_kind:  str = "raw"
        self._transform_stack: list[str] = []
        self._redo_kinds:      list[str] = []
        # Reps value that was active when each undo entry was pushed —
        # we restore it on undo so that a torsion applied at Reps=12
        # comes back to Reps=12 (not 1) when the user reverts.
        self._reps_stack:      list[int] = []
        self._redo_reps:       list[int] = []
        # Undo / redo stacks for transformation history (strain, torsion,
        # bundle, MWNT, etc.).  Each entry is a NanotubeStructure snapshot.
        self._undo_stack: list[NanotubeStructure] = []
        self._redo_stack: list[NanotubeStructure] = []
        self._UNDO_LIMIT: int = 32
        # Internal flag set during undo/redo so we don't push onto the wrong stack.
        self._suppress_history: bool = False
        self._bond_settings = BondSettings()
        self._build_ui()

    def set_parent_lattice(self, structure) -> None:
        """Record metadata of the parent 2D lattice (LatticeStructure).

        Used by AnalysisDialog to decide whether the zone-folding
        electronic-character rule applies (graphene-specific).
        """
        if structure is None:
            self._lattice_type       = "hexagonal"
            self._lattice_n_atoms_uc = None
            self._lattice_a          = None
            return
        self._lattice_type       = getattr(structure, "lattice_type", "hexagonal")
        try:
            self._lattice_n_atoms_uc = len(structure.atoms)
        except Exception:
            self._lattice_n_atoms_uc = None
        self._lattice_a = float(getattr(structure, "a", 0.0)) or None

    @property
    def bond_settings(self) -> BondSettings:
        """Shared BondSettings instance (used by viewer AND main_window)."""
        return self._bond_settings

    # ─────────────────────────────────────────────────────────────────────────
    # UI construction
    # ─────────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ── 3D Viewer ────────────────────────────────────────────────────────
        # Analysis, Methods and DFT Inputs are now available from the
        # application menubar (see main_window.py) rather than as tabs.
        if _HAS_GL:
            self._viewer = _GLViewer()
        else:
            self._viewer = _MPLViewer()
        self._viewer.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        root.addWidget(self._viewer, stretch=3)

        # ── Display Options ───────────────────────────────────────────────────
        disp_box = QGroupBox("Display")
        disp_lay = QVBoxLayout(disp_box)
        disp_lay.setSpacing(4)

        # Row 1 – simulation box toggle
        self._chk_box = QCheckBox("Show simulation box")
        self._chk_box.setChecked(True)
        self._chk_box.toggled.connect(self._on_view_changed)
        disp_lay.addWidget(self._chk_box)

        # Row 2 – vacuum + replications side by side
        row2 = QHBoxLayout()
        row2.setSpacing(8)
        row2.setAlignment(Qt.AlignmentFlag.AlignVCenter)

        row2.addWidget(QLabel("Vacuum:"))
        self._spin_vac = QDoubleSpinBox()
        self._spin_vac.setRange(0.0, 80.0)
        self._spin_vac.setValue(10.0)
        self._spin_vac.setSuffix(" Å")
        self._spin_vac.setSingleStep(2.0)
        self._spin_vac.setDecimals(1)
        self._spin_vac.setToolTip("Vacuum padding around nanotube (Å) — also applied on export")
        self._spin_vac.valueChanged.connect(self._on_view_changed)
        row2.addWidget(self._spin_vac)

        row2.addSpacing(16)
        self._lbl_reps = QLabel("Reps:")
        row2.addWidget(self._lbl_reps)
        self._spin_rep = QSpinBox()
        self._spin_rep.setRange(1, 999)
        self._spin_rep.setValue(1)
        self._spin_rep.setSuffix(" ×")
        self._spin_rep.setToolTip(
            "Replications along the tube axis (Z). "
            "Replicated structure is also exported when > 1.\n"
            "Large values may take noticeable time to render in the 3D viewer."
        )
        self._spin_rep.valueChanged.connect(self._on_view_changed)
        row2.addWidget(self._spin_rep)
        row2.addSpacing(16)
        btn_bonds = QPushButton("⚙  Bond Cutoffs…")
        btn_bonds.setFixedHeight(28)              # align with the spinboxes
        btn_bonds.setToolTip(
            "Edit per-species bond length cutoffs.\n"
            "Defaults use covalent radii from Alvarez (2008) × 1.20."
        )
        btn_bonds.clicked.connect(self._on_bond_cutoffs)
        row2.addWidget(btn_bonds)

        row2.addSpacing(8)

        # ── Undo / Redo with icon + visible "Undo" / "Redo" labels
        # (sharing row with Bond Cutoffs to avoid adding an extra row).
        # Qt's standard pixmaps guarantee glyph availability across
        # fonts where the unicode curly arrows ↶ ↷ are missing. ────────
        from PyQt6.QtWidgets import QStyle
        _std = self.style()
        self.btn_undo = QPushButton("Undo")
        self.btn_undo.setIcon(_std.standardIcon(QStyle.StandardPixmap.SP_ArrowBack))
        self.btn_undo.setFixedHeight(28)
        self.btn_undo.setMinimumWidth(72)
        self.btn_undo.setEnabled(False)
        self.btn_undo.setShortcut("Ctrl+Z")
        self.btn_undo.setToolTip(
            "Undo last structural change (strain, torsion, bundle, MWNT, …)."
            "\nShortcut: Ctrl+Z"
        )
        self.btn_undo.clicked.connect(self._on_undo)
        row2.addWidget(self.btn_undo)

        self.btn_redo = QPushButton("Redo")
        self.btn_redo.setIcon(_std.standardIcon(QStyle.StandardPixmap.SP_ArrowForward))
        self.btn_redo.setFixedHeight(28)
        self.btn_redo.setMinimumWidth(72)
        self.btn_redo.setEnabled(False)
        self.btn_redo.setShortcut("Ctrl+Y")
        self.btn_redo.setToolTip(
            "Redo last undone change.\nShortcut: Ctrl+Y"
        )
        self.btn_redo.clicked.connect(self._on_redo)
        row2.addWidget(self.btn_redo)

        row2.addStretch()

        disp_lay.addLayout(row2)
        root.addWidget(disp_box)

        # ── Info box (3-column grid to save vertical space) ──────────────────
        info_box = QGroupBox("Nanotube Properties")
        info_lay = QGridLayout(info_box)
        info_lay.setHorizontalSpacing(12)
        info_lay.setVerticalSpacing(2)
        info_lay.setContentsMargins(8, 14, 8, 6)

        self._info: dict[str, QLabel] = {}
        keys = ["Indices", "Diameter", "Chiral angle",
                "Atoms/cell", "Length |T|", "Strain",
                "Elements"]
        for i, key in enumerate(keys):
            row, col = divmod(i, 3)
            lbl_key = QLabel(f"{key}:")
            lbl_key.setAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            lbl_key.setStyleSheet("font-size: 10px; color: #1A3A6B;")
            lbl_val = QLabel("—")
            lbl_val.setStyleSheet("font-family: monospace; font-size: 10px;")
            if key == "Elements":
                lbl_val.setTextFormat(Qt.TextFormat.RichText)
            self._info[key] = lbl_val
            info_lay.addWidget(lbl_key, row, 2 * col,     Qt.AlignmentFlag.AlignRight)
            info_lay.addWidget(lbl_val, row, 2 * col + 1, Qt.AlignmentFlag.AlignLeft)
        # Three column-pairs: stretch the value columns so they take up
        # the spare horizontal space (keys stay tight).
        for col in (1, 3, 5):
            info_lay.setColumnStretch(col, 1)

        root.addWidget(info_box)

        # ── Export group ─────────────────────────────────────────────────────
        exp_box = QGroupBox("Export")
        exp_lay = QVBoxLayout(exp_box)

        fmt_row = QHBoxLayout()
        fmt_row.addWidget(QLabel("Format:"))
        self.combo_fmt = QComboBox()
        self.combo_fmt.addItems([
            ".pdb", ".xyz", "POSCAR", ".lammps", ".pwi (QE)",
            ".xsf", ".inp (CP2K)", ".fdf (SIESTA)", ".cif",
        ])
        fmt_row.addWidget(self.combo_fmt)
        fmt_row.addSpacing(12)
        fmt_row.addWidget(QLabel("Tube axis:"))
        self.combo_axis = QComboBox()
        self.combo_axis.addItems(["Z", "X", "Y"])
        self.combo_axis.setToolTip(
            "Longitudinal axis of the nanotube in the exported file.\n\n"
            "Z (default) — tube runs along Z; periodic box vector is Lz.\n"
            "X — coordinates are permuted so the tube runs along X.\n"
            "Y — coordinates are permuted so the tube runs along Y.\n\n"
            "Useful when the downstream code expects a specific orientation\n"
            "(e.g. LAMMPS fix deform along x, or QE ibrav conventions)."
        )
        fmt_row.addWidget(self.combo_axis)
        exp_lay.addLayout(fmt_row)

        # Unit-cell-only toggle — when checked, exports n_rep=1 regardless of the
        # Reps display spinbox so the user can view a supercell but export the
        # minimal unit cell.
        self._chk_unit_cell = QCheckBox("Export unit cell only  (ignore Reps)")
        self._chk_unit_cell.setToolTip(
            "When checked: export the single unit cell (1 repetition) even if\n"
            "Reps > 1 is selected for the 3D viewer.\n\n"
            "When unchecked: export exactly what you see (Reps repetitions)."
        )
        exp_lay.addWidget(self._chk_unit_cell)

        # Supercell warning — shown only when Reps > 1 AND unit-cell export is off
        self._lbl_rep_warn = QLabel()
        self._lbl_rep_warn.setStyleSheet(
            "color: #E69F00; font-size: 10px; font-style: italic;"
        )
        self._lbl_rep_warn.setVisible(False)
        exp_lay.addWidget(self._lbl_rep_warn)
        self._spin_rep.valueChanged.connect(self._update_rep_warning)
        self._chk_unit_cell.stateChanged.connect(self._update_rep_warning)

        btn_row = QHBoxLayout()

        self.btn_export = QPushButton("💾  Export …")
        self.btn_export.setObjectName("btn_export")
        self.btn_export.setFixedHeight(30)
        self.btn_export.setEnabled(False)
        self.btn_export.clicked.connect(self._on_export)
        btn_row.addWidget(self.btn_export)

        self.btn_sys = QPushButton("⚙  Batch …")
        self.btn_sys.setObjectName("btn_sys")
        self.btn_sys.setFixedHeight(30)
        self.btn_sys.setEnabled(False)
        self.btn_sys.clicked.connect(self._on_systematize)
        btn_row.addWidget(self.btn_sys)

        exp_lay.addLayout(btn_row)
        root.addWidget(exp_box)

        # ── Tools group ───────────────────────────────────────────────────────
        # All tools operate on whatever is currently visible in the viewer
        # ("what you see is what you get").  Each tool replaces _nanotube with
        # the result — the viewer updates immediately.
        tools_box = QGroupBox("Tools  (applied to current structure)")
        tools_lay = QVBoxLayout(tools_box)
        tools_lay.setSpacing(4)

        # All three transformation buttons in a single row.
        row_tools = QHBoxLayout()

        self.btn_mwnt = QPushButton("⧉  Multi-Wall…")
        self.btn_mwnt.setFixedHeight(28)
        self.btn_mwnt.setEnabled(False)
        self.btn_mwnt.setToolTip(
            "Wrap this nanotube as the innermost wall of a multi-walled structure.\n"
            "Wall k uses (k·n, k·m) of the primitive direction."
        )
        self.btn_mwnt.clicked.connect(self._on_mwnt)
        row_tools.addWidget(self.btn_mwnt)

        self.btn_deform = QPushButton("↔  Deform/Torsion…")
        self.btn_deform.setFixedHeight(28)
        self.btn_deform.setEnabled(False)
        self.btn_deform.setToolTip(
            "Apply axial strain and/or torsion to the current structure.\n"
            "Always starts from what is currently displayed."
        )
        self.btn_deform.clicked.connect(self._on_deform)
        row_tools.addWidget(self.btn_deform)

        self.btn_bundle = QPushButton("⬡  Bundle…")
        self.btn_bundle.setFixedHeight(28)
        self.btn_bundle.setEnabled(False)
        self.btn_bundle.setToolTip(
            "Replicate the current structure into a periodic bundle supercell.\n"
            "Always starts from what is currently displayed."
        )
        self.btn_bundle.clicked.connect(self._on_bundle)
        row_tools.addWidget(self.btn_bundle)

        tools_lay.addLayout(row_tools)

        # Analysis / Methods / DFT-Inputs used to be three separate buttons
        # here; they now live as tabs at the top of the right pane (see
        # _top_tabs above) so they don't crowd the toolbox.
        # We keep ``btn_analysis``, ``btn_methods`` and ``btn_dft`` as
        # disabled-no-op placeholders to satisfy any legacy ``setEnabled``
        # calls elsewhere in this class — the tabs themselves are the
        # source of truth for enablement.
        self.btn_analysis = QPushButton(); self.btn_analysis.hide()
        self.btn_methods  = QPushButton(); self.btn_methods.hide()
        self.btn_dft      = QPushButton(); self.btn_dft.hide()

        root.addWidget(tools_box)

    # ─────────────────────────────────────────────────────────────────────────
    # Public API
    # ─────────────────────────────────────────────────────────────────────────

    def _update_rep_warning(self, *_):
        """Show/hide the supercell warning based on Reps spinbox + unit-cell toggle."""
        n = self._spin_rep.value()
        unit_only = self._chk_unit_cell.isChecked()
        if n > 1 and not unit_only:
            nt    = self._nanotube
            atoms = (nt.n_atoms * n) if nt is not None else "?"
            self._lbl_rep_warn.setText(
                f"⚠  Supercell: {n} × unit cell along z  ({atoms} atoms)"
            )
            self._lbl_rep_warn.setVisible(True)
        else:
            self._lbl_rep_warn.setVisible(False)

    def set_nanotube(self, nt: NanotubeStructure, is_raw: bool = True,
                     kind: str | None = None):
        """Load a nanotube into the viewer.

        Parameters
        ----------
        nt      : the NanotubeStructure to display
        is_raw  : when True (default) also cache this as the pre-deformation
                  baseline so "Reset" inside DeformDialog can restore it.
                  Pass False when displaying a deformed or bundle supercell
                  to keep the previous raw reference intact.
        kind    : transformation label to associate with this state — one
                  of ``"raw"`` (freshly built SWNT), ``"mwnt"``,
                  ``"bundle"``, ``"strain"``, ``"torsion"``.  Used by the
                  operation-chaining warnings (e.g. MWNT-after-bundle).
                  Defaults to ``"raw"`` for is_raw=True calls and inherits
                  the previous kind otherwise.
        """
        # Resolve the transformation label for this state.
        if kind is None:
            kind = "raw" if is_raw else self._transform_kind

        # ── Push previous state onto the undo stack so the user can revert.
        # Skip when we are inside _on_undo / _on_redo (the stacks manage
        # themselves there) and when this is the initial build (is_raw=True
        # with no prior nanotube — there is nothing to revert to).
        if (not self._suppress_history) and self._nanotube is not None:
            self._undo_stack.append(self._nanotube)
            self._transform_stack.append(self._transform_kind)
            # Snapshot the Reps value too — torsion sets it to 1 and
            # hides the spinbox; on undo we need to bring back the
            # value the user had set before the torsion.
            self._reps_stack.append(int(self._spin_rep.value()))
            if len(self._undo_stack) > self._UNDO_LIMIT:
                self._undo_stack.pop(0)
                self._transform_stack.pop(0)
                self._reps_stack.pop(0)
            # A new branch invalidates the redo stack.
            self._redo_stack.clear()
            self._redo_kinds.clear()
            self._redo_reps.clear()

        self._nanotube = nt
        self._transform_kind = kind
        if is_raw:
            # A fresh build resets the entire history — the user explicitly
            # asked for a new structure, not a transformation of the old one.
            self._nanotube_raw = nt
            if not self._suppress_history:
                self._undo_stack.clear()
                self._redo_stack.clear()
                self._transform_stack.clear()
                self._redo_kinds.clear()
                self._reps_stack.clear()
                self._redo_reps.clear()
        # Reps spinbox is meaningful only when the structure is still a
        # periodic unit cell — i.e. anything other than a torqued state.
        # _set_reps_visible may not exist yet during the very first call
        # from __init__ → _build_ui (defensive ``hasattr`` check).
        if hasattr(self, "_lbl_reps"):
            self._set_reps_visible(kind != "torsion")
        self._update_undo_buttons()
        r = nt.chirality

        self._info["Indices"].setText(f"({r.n}, {r.m})")
        self._info["Diameter"].setText(f"{nt.diameter:.4f} Å")
        self._info["Chiral angle"].setText(f"{r.theta_deg:.2f}°")
        self._info["Atoms/cell"].setText(str(nt.n_atoms))
        self._info["Length |T|"].setText(f"{nt.length:.4f} Å")
        self._info["Strain"].setText(f"{r.strain:.6f} %")
        self._update_legend(nt)

        # Sync vacuum spinbox with built vacuum (block signal to avoid double render)
        self._spin_vac.blockSignals(True)
        self._spin_vac.setValue(nt.vacuum)
        self._spin_vac.blockSignals(False)

        self._refresh_viewer()
        self.btn_export.setEnabled(True)
        self.btn_sys.setEnabled(True)
        self.btn_mwnt.setEnabled(True)
        self.btn_deform.setEnabled(True)
        self.btn_bundle.setEnabled(True)
        # The Analysis / Methods / DFT actions live in the menubar now;
        # tell main_window to enable them.
        self.nanotube_loaded.emit(True)
        self._update_rep_warning()

    def _update_legend(self, nt: NanotubeStructure):
        """Show a compact colour legend for the elements present."""
        species = sorted(set(nt.symbols))
        parts = []
        for sym in species:
            col = _CPK_COLORS.get(sym, _CPK_DEFAULT)
            hex_col = "#{:02x}{:02x}{:02x}".format(
                int(col[0]*255), int(col[1]*255), int(col[2]*255))
            parts.append(
                f'<span style="color:{hex_col}; font-weight:bold;">■</span> {sym}'
            )
        self._info["Elements"].setText("  ".join(parts))

    # ─────────────────────────────────────────────────────────────────────────
    # Slots
    # ─────────────────────────────────────────────────────────────────────────

    def _on_view_changed(self):
        """Called whenever any display option changes."""
        if self._nanotube is None:
            return
        self._refresh_viewer()

    def _refresh_viewer(self):
        self._viewer.display(
            self._nanotube,
            vacuum   = self._spin_vac.value(),
            n_rep    = self._spin_rep.value(),
            show_box = self._chk_box.isChecked(),
            settings = self._bond_settings,
        )

    def _on_bond_cutoffs(self):
        """Open the Bond Cutoffs dialog; re-render if settings changed."""
        from gui.dialogs.bond_settings_dialog import BondSettingsDialog

        # Determine which species are present (use current nanotube or fall back)
        if self._nanotube is not None:
            species = list(set(self._nanotube.symbols))
        else:
            # No nanotube yet — open with empty species list (global params only)
            species = []

        dlg = BondSettingsDialog(self._bond_settings, species, parent=self)
        if dlg.exec() and self._nanotube is not None:
            # Settings were updated in-place by the dialog; re-render
            self._refresh_viewer()

    _FMT_MAP = {
        ".pdb":          ".pdb",
        ".xyz":          ".xyz",
        "POSCAR":        "poscar",
        ".lammps":       ".lammps",
        ".pwi (QE)":     ".pwi",
        ".xsf":          ".xsf",
        ".inp (CP2K)":   ".inp",
        ".fdf (SIESTA)": ".fdf",
        ".cif":          ".cif",
    }

    def _on_export(self):
        if self._nanotube is None:
            return
        fmt_key = self.combo_fmt.currentText()
        fmt     = self._FMT_MAP[fmt_key]
        r       = self._nanotube.chirality
        default = f"NT_n{r.n}_m{r.m}{fmt if fmt.startswith('.') else ''}"

        path, _ = QFileDialog.getSaveFileName(
            self, "Export nanotube", default, "All files (*)"
        )
        if path:
            # When "Export unit cell only" is ticked, always send n_rep=1
            # so the user can view a supercell but export the minimal cell.
            n_rep = 1 if self._chk_unit_cell.isChecked() else self._spin_rep.value()
            self.export_requested.emit(
                fmt, path,
                self._spin_vac.value(),
                n_rep,
                self.combo_axis.currentText(),   # "Z", "X", or "Y"
            )

    def _on_systematize(self):
        dlg = _SystematizeDialog(self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self.systematize_requested.emit(dlg.params())

    # ─────────────────────────────────────────────────────────────────────────
    # Operation-chaining confirmation
    # ─────────────────────────────────────────────────────────────────────────
    _CHAIN_WARNINGS = {
        # (next_op, current_kind) → user-facing message
        ("mwnt",   "bundle"):
            "You are about to build a multi-walled nanotube whose inner wall is "
            "a <b>bundle</b> ({n} tubes).  The MWNT pipeline treats the current "
            "structure as a single tube — the resulting geometry is not what "
            "you probably expect.<br><br>Continue anyway?",
        ("mwnt",   "mwnt"):
            "The current structure is already a multi-walled nanotube.  "
            "Re-running Multi-Wall will discard it and start from the inner "
            "wall's (n, m) again.<br><br>Continue?",
        ("mwnt",   "torsion"):
            "Building a MWNT on top of a torqued tube discards the torsion and "
            "starts from the unstrained inner (n, m).<br><br>Continue?",
        ("bundle", "torsion"):
            "Building a bundle on top of a torqued tube replicates the twisted "
            "geometry, which usually breaks the bundle periodicity.<br><br>"
            "Continue anyway?",
        ("strain", "torsion"):
            "Axial strain on top of a torqued tube produces a non-uniform "
            "twist density.  Use this only if you know what you are doing."
            "<br><br>Continue?",
    }

    def _confirm_chain(self, next_op: str, extra: dict | None = None) -> bool:
        """Show a QMessageBox.question when chaining a possibly-bogus operation.

        Returns True when the user wants to proceed (also when there is no
        warning for this combination).
        """
        msg = self._CHAIN_WARNINGS.get((next_op, self._transform_kind))
        if msg is None:
            return True
        if extra:
            msg = msg.format(**extra)
        from PyQt6.QtWidgets import QMessageBox
        reply = QMessageBox.question(
            self,
            f"Chained operation — {next_op}",
            msg,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return reply == QMessageBox.StandardButton.Yes

    def _on_mwnt(self):
        """Open MWNT dialog and emit mwnt_requested so main_window can run the pipeline."""
        if self._nanotube is None:
            return
        if not self._confirm_chain("mwnt", extra={"n": "?"}):
            return
        from gui.dialogs.advanced_dialogs import MWNTDialog
        r = self._nanotube.chirality
        dlg = MWNTDialog(r.n, r.m, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self.mwnt_requested.emit(
                r.n, r.m,
                dlg.vacuum,
                dlg.roll_inward,
                dlg.n_walls,
                dlg.spacing,
            )

    def _on_deform(self):
        """Apply strain/torsion to whatever is currently displayed."""
        if self._nanotube is None:
            return
        from gui.dialogs.advanced_dialogs import DeformDialog
        from PyQt6.QtWidgets import QMessageBox

        current_vacuum = float(getattr(self._nanotube, "vacuum", 10.0))
        dlg = DeformDialog(default_vacuum=current_vacuum, parent=self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        # Chain confirmation depends on which deformation is being applied.
        if abs(dlg.axial_strain) > 1e-9 and not self._confirm_chain("strain"):
            return

        # If the user has Reps > 1 set in the viewer and is applying torsion,
        # apply the twist to that supercell rather than the unit cell.
        n_rep_view = int(self._spin_rep.value())
        applies_torsion = abs(dlg.twist_rate) > 1e-9

        try:
            from core.deformations import (
                apply_axial_strain, apply_torsion, torsion_warning,
            )
            nt = self._nanotube   # "what you see is what you get"
            applied_kind = self._transform_kind
            if abs(dlg.axial_strain) > 1e-9:
                nt = apply_axial_strain(nt, dlg.axial_strain)
                applied_kind = "strain"
            if applies_torsion:
                nt = apply_torsion(
                    nt, dlg.twist_rate,
                    z_vacuum = dlg.z_vacuum,
                    n_rep    = n_rep_view,
                )
                applied_kind = "torsion"

                # Once a torsion is applied, the structure is no longer a
                # periodic unit cell — the Reps control has no meaning.
                # Reset it to 1 and hide it.  It will be restored by undo
                # when the user steps back to a periodic state.
                self._spin_rep.blockSignals(True)
                self._spin_rep.setValue(1)
                self._spin_rep.blockSignals(False)
                self._set_reps_visible(False)

            self.set_nanotube(nt, is_raw=False, kind=applied_kind)

            # Inform the user when torsion was applied, so the loss of
            # axial periodicity is explicit rather than implicit.
            warn = torsion_warning(
                dlg.twist_rate, dlg.z_vacuum, n_rep=n_rep_view,
            )
            if warn is not None:
                QMessageBox.information(self, "Torsion applied", warn)
        except Exception as exc:
            QMessageBox.critical(self, "Deform error", str(exc))

    def _set_reps_visible(self, visible: bool) -> None:
        """Show or hide the Reps spinbox and its label.

        Used to hide Reps after a torsion (the structure is no longer a
        periodic unit cell) and restore it when undo brings us back to
        a periodic state.
        """
        self._lbl_reps.setVisible(visible)
        self._spin_rep.setVisible(visible)

    def _on_bundle(self):
        if self._nanotube is None:
            return
        if not self._confirm_chain("bundle"):
            return
        from gui.dialogs.advanced_dialogs import BundleDialog
        # Inherit the current single-tube vacuum as the default so the user
        # gets the same lateral padding around the bundle unless they change it.
        current_vacuum = float(getattr(self._nanotube, "vacuum", 10.0))
        dlg = BundleDialog(self._nanotube.diameter,
                           vacuum=current_vacuum, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            try:
                from core.bundles import build_bundle
                result = build_bundle(
                    self._nanotube,
                    geometry=dlg.geometry,
                    spacing=dlg.spacing,
                    vacuum=dlg.vacuum,
                    nx=dlg.nx,
                    ny=dlg.ny,
                )
                self.set_nanotube(result.nanotube, is_raw=False, kind="bundle")
            except Exception as exc:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.critical(self, "Bundle error", str(exc))

    # ─────────────────────────────────────────────────────────────────────────
    # Public entry-points for menubar actions (Analysis / Methods / DFT)
    # ─────────────────────────────────────────────────────────────────────────

    def open_analysis_dialog(self) -> None:
        """Open the Analysis dialog with the appropriate lattice metadata."""
        if self._nanotube is None:
            return
        from gui.dialogs.advanced_dialogs import AnalysisDialog
        dlg = AnalysisDialog(
            self._nanotube,
            lattice_type       = self._lattice_type,
            n_atoms_unit_cell  = self._lattice_n_atoms_uc,
            lattice_constant_a = self._lattice_a,
            parent=self,
        )
        dlg.exec()

    def open_methods_dialog(self) -> None:
        if self._nanotube is None:
            return
        from gui.dialogs.advanced_dialogs import MethodsDialog
        MethodsDialog(self._nanotube, parent=self).exec()

    def open_dft_inputs_dialog(self) -> None:
        if self._nanotube is None:
            return
        from gui.dialogs.advanced_dialogs import DFTInputDialog
        DFTInputDialog(self._nanotube, parent=self).exec()

    # Legacy aliases for any code path that still uses the old names.
    _on_analysis   = open_analysis_dialog
    _on_methods    = open_methods_dialog
    _on_dft_inputs = open_dft_inputs_dialog

    # ─────────────────────────────────────────────────────────────────────────
    # Undo / Redo
    # ─────────────────────────────────────────────────────────────────────────

    def _update_undo_buttons(self) -> None:
        """Enable / disable the Undo and Redo buttons based on stack state."""
        if hasattr(self, "btn_undo"):
            self.btn_undo.setEnabled(bool(self._undo_stack))
        if hasattr(self, "btn_redo"):
            self.btn_redo.setEnabled(bool(self._redo_stack))

    def _on_undo(self) -> None:
        """Revert the last structural change."""
        if not self._undo_stack or self._nanotube is None:
            return
        # Push current state (incl. current Reps) to the redo side.
        self._redo_stack.append(self._nanotube)
        self._redo_kinds.append(self._transform_kind)
        self._redo_reps.append(int(self._spin_rep.value()))
        prev      = self._undo_stack.pop()
        prev_kind = self._transform_stack.pop() if self._transform_stack else "raw"
        prev_reps = self._reps_stack.pop()       if self._reps_stack       else 1
        self._suppress_history = True
        try:
            # is_raw=False so the pre-deformation reference is preserved.
            self.set_nanotube(prev, is_raw=False, kind=prev_kind)
            # Restore the Reps value that was active before the change.
            self._spin_rep.blockSignals(True)
            self._spin_rep.setValue(prev_reps)
            self._spin_rep.blockSignals(False)
            # _set_reps_visible was already invoked by set_nanotube via
            # the ``kind`` it received; an explicit refresh of the
            # supercell warning keeps the bottom-of-pane label honest.
            self._update_rep_warning()
            self._refresh_viewer()
        finally:
            self._suppress_history = False
        self._update_undo_buttons()

    def _on_redo(self) -> None:
        """Re-apply the last undone change."""
        if not self._redo_stack or self._nanotube is None:
            return
        self._undo_stack.append(self._nanotube)
        self._transform_stack.append(self._transform_kind)
        self._reps_stack.append(int(self._spin_rep.value()))
        nxt      = self._redo_stack.pop()
        nxt_kind = self._redo_kinds.pop() if self._redo_kinds else self._transform_kind
        nxt_reps = self._redo_reps.pop()  if self._redo_reps  else 1
        self._suppress_history = True
        try:
            self.set_nanotube(nxt, is_raw=False, kind=nxt_kind)
            self._spin_rep.blockSignals(True)
            self._spin_rep.setValue(nxt_reps)
            self._spin_rep.blockSignals(False)
            self._update_rep_warning()
            self._refresh_viewer()
        finally:
            self._suppress_history = False
        self._update_undo_buttons()


# ─────────────────────────────────────────────────────────────────────────────
# CPK colours + covalent radii
# ─────────────────────────────────────────────────────────────────────────────

_CPK_COLORS: dict[str, tuple] = {
    "H":  (0.95, 0.95, 0.95, 1.0),
    "C":  (0.40, 0.40, 0.40, 1.0),
    "N":  (0.20, 0.40, 0.90, 1.0),
    "O":  (0.90, 0.20, 0.20, 1.0),
    "F":  (0.70, 0.95, 0.70, 1.0),
    "P":  (1.00, 0.60, 0.00, 1.0),
    "S":  (0.95, 0.90, 0.10, 1.0),
    "Cl": (0.20, 0.85, 0.20, 1.0),
    "B":  (1.00, 0.72, 0.55, 1.0),
    "Si": (0.50, 0.60, 0.60, 1.0),
    "Fe": (0.88, 0.40, 0.20, 1.0),
    "Mo": (0.33, 0.33, 0.75, 1.0),
    "W":  (0.40, 0.60, 0.20, 1.0),
    "Se": (0.60, 0.30, 0.10, 1.0),
    "Te": (0.67, 0.40, 0.00, 1.0),
    "Bi": (0.62, 0.31, 0.71, 1.0),
    "Li": (0.80, 0.50, 1.00, 1.0),
    "Na": (0.67, 0.36, 0.95, 1.0),
    "K":  (0.56, 0.25, 0.83, 1.0),
    "Ca": (0.24, 1.00, 0.00, 1.0),
    "Ti": (0.75, 0.76, 0.78, 1.0),
    "Cr": (0.54, 0.60, 0.78, 1.0),
    "Mn": (0.61, 0.48, 0.78, 1.0),
    "Co": (0.94, 0.56, 0.63, 1.0),
    "Ni": (0.31, 0.82, 0.31, 1.0),
    "Cu": (0.78, 0.50, 0.20, 1.0),
    "Zn": (0.49, 0.50, 0.69, 1.0),
    "Ga": (0.76, 0.56, 0.56, 1.0),
    "Ge": (0.40, 0.56, 0.56, 1.0),
    "As": (0.74, 0.50, 0.89, 1.0),
    "Br": (0.65, 0.16, 0.16, 1.0),
    "Ru": (0.24, 0.61, 0.47, 1.0),
    "Rh": (0.04, 0.49, 0.55, 1.0),
    "Pd": (0.00, 0.41, 0.52, 1.0),
    "Ag": (0.75, 0.75, 0.75, 1.0),
    "In": (0.65, 0.46, 0.45, 1.0),
    "Sn": (0.40, 0.50, 0.50, 1.0),
    "Sb": (0.62, 0.39, 0.71, 1.0),
    "I":  (0.58, 0.00, 0.58, 1.0),
    "Pt": (0.82, 0.82, 0.88, 1.0),
    "Au": (1.00, 0.82, 0.14, 1.0),
    "Pb": (0.34, 0.35, 0.38, 1.0),
}
_CPK_DEFAULT = (0.80, 0.40, 0.80, 1.0)  # magenta for unknown elements


# ─────────────────────────────────────────────────────────────────────────────
# 3D viewer — pyqtgraph OpenGL
# ─────────────────────────────────────────────────────────────────────────────

class _GLViewer(QWidget if not _HAS_GL else gl.GLViewWidget):
    """
    Hardware-accelerated 3D molecular stick viewer.

    Rendering: bond-only model with bicolour VESTA convention.
    Each bond is split at its midpoint; each half inherits its endpoint's
    CPK colour.  No spheres are drawn.

    Optional overlays:
      • Unit cell wireframe box (semi-transparent light-blue lines)

    Replications stack copies along Z; inter-cell bonds are automatically
    detected because the replicated coords are used for the cKDTree search.
    """

    def __init__(self, parent=None):
        if _HAS_GL:
            super().__init__(parent)
            self.setMinimumHeight(280)
            self.opts["distance"] = 40
            self.opts["elevation"] = 20
            self.opts["azimuth"]   = 45
            # Dark background
            self.setBackgroundColor((15, 15, 30, 255))
        else:
            QWidget.__init__(self, parent)

    # ── public entry point ────────────────────────────────────────────────────

    def display(self, nt: NanotubeStructure,
                vacuum: float = 10.0,
                n_rep: int    = 1,
                show_box: bool = True,
                settings: "BondSettings | None" = None):
        if not _HAS_GL:
            return

        self.clear()

        self._draw_bonds(nt, n_rep, settings)

        if show_box:
            # Use the actual stored box (correct for both single tubes and bundles/MWNTs)
            box_xy = float(nt.box[0])
            box_z  = float(nt.length) * n_rep
            self._draw_box(box_xy, box_z)

        # ── Camera — base distance on the actual XY extent, not just diameter ──
        half_xy   = float(nt.box[0]) / 2.0
        total_z   = float(nt.length) * n_rep
        span      = max(total_z / 2.0, half_xy)
        self.opts["distance"]  = max(span * 2.8, half_xy * 5.5, 20.0)
        self.opts["elevation"] = 20
        self.opts["azimuth"]   = 45
        self.update()

    # ── internal helpers ──────────────────────────────────────────────────────

    def _draw_bonds(self, nt: NanotubeStructure, n_rep: int,
                    settings: "BondSettings | None" = None):
        """Build bond geometry for n_rep copies and add a GLLinePlotItem."""
        from core.connectivity import compute_bonds, bond_line_arrays

        base_coords  = nt.coords.astype(np.float32)
        base_coords -= base_coords.mean(axis=0)   # centre on origin
        base_syms    = list(nt.symbols)
        L            = float(nt.length)

        if n_rep > 1:
            all_coords_list: list[np.ndarray] = []
            all_syms: list[str] = []
            for i in range(n_rep):
                shift = np.array([0.0, 0.0, i * L], dtype=np.float32)
                all_coords_list.append(base_coords + shift)
                all_syms.extend(base_syms)
            disp_coords = np.vstack(all_coords_list)
            disp_coords -= disp_coords.mean(axis=0)   # re-centre
            disp_syms = all_syms
        else:
            disp_coords = base_coords
            disp_syms   = base_syms

        # cKDTree bond search on the full (replicated) coordinate set
        bonds = compute_bonds(disp_coords, disp_syms, settings=settings)
        if not bonds:
            return

        pts, bond_colors = bond_line_arrays(
            disp_coords, disp_syms, bonds, _CPK_COLORS, _CPK_DEFAULT
        )

        n_total = len(disp_syms)
        width   = 3.0 if n_total <= 600 else (2.0 if n_total <= 6_000 else 1.5)

        lines = gl.GLLinePlotItem(
            pos=pts, color=bond_colors,
            width=width, mode="lines", antialias=True,
        )
        self.addItem(lines)

    def _draw_box(self, box_xy: float, box_z: float):
        """Draw a semi-transparent wireframe parallelepiped (orthorhombic cell)."""
        h  = box_xy / 2.0
        z0 = -box_z  / 2.0
        z1 =  box_z  / 2.0

        corners = np.array([
            [-h, -h, z0], [ h, -h, z0], [ h,  h, z0], [-h,  h, z0],
            [-h, -h, z1], [ h, -h, z1], [ h,  h, z1], [-h,  h, z1],
        ], dtype=np.float32)

        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),   # bottom face
            (4, 5), (5, 6), (6, 7), (7, 4),   # top face
            (0, 4), (1, 5), (2, 6), (3, 7),   # vertical pillars
        ]

        pts = []
        for a, b in edges:
            pts.extend([corners[a], corners[b]])
        pts = np.array(pts, dtype=np.float32)

        # Light-blue, semi-transparent
        col = np.full((len(pts), 4), [0.55, 0.75, 1.0, 0.45], dtype=np.float32)

        box_item = gl.GLLinePlotItem(
            pos=pts, color=col, width=1.2, mode="lines", antialias=False
        )
        self.addItem(box_item)


# ─────────────────────────────────────────────────────────────────────────────
# 3D viewer — matplotlib fallback
# ─────────────────────────────────────────────────────────────────────────────

class _MPLViewer(QWidget):
    """Matplotlib 3D fallback viewer (no OpenGL required)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        import matplotlib
        matplotlib.use("QtAgg")
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

        self.fig = Figure(figsize=(4, 4), facecolor="#1A1A2E")
        self.ax  = self.fig.add_subplot(111, projection="3d")
        self.canvas = FigureCanvasQTAgg(self.fig)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas)

        lbl = QLabel("ℹ  Install PyOpenGL for interactive 3D rendering")
        lbl.setStyleSheet("color:#999; font-size:9px;")
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(lbl)

    def display(self, nt: NanotubeStructure,
                vacuum=10.0, n_rep=1, show_box=True,
                settings: "BondSettings | None" = None):
        from core.connectivity import compute_bonds

        self.ax.cla()
        self.ax.set_facecolor("#1A1A2E")
        self.fig.patch.set_facecolor("#1A1A2E")

        base_coords = nt.coords.astype(float)
        base_coords -= base_coords.mean(axis=0)
        base_syms = list(nt.symbols)
        L = float(nt.length)

        if n_rep > 1:
            parts = [base_coords + np.array([0, 0, i * L]) for i in range(n_rep)]
            all_coords = np.vstack(parts)
            all_syms   = base_syms * n_rep
            all_coords -= all_coords.mean(axis=0)
        else:
            all_coords = base_coords
            all_syms   = base_syms

        bonds = compute_bonds(all_coords, all_syms, settings=settings)
        _mpl_cpk = {
            k: "#{:02x}{:02x}{:02x}".format(int(v[0]*255), int(v[1]*255), int(v[2]*255))
            for k, v in _CPK_COLORS.items()
        }
        for (i, j) in bonds:
            mid = (all_coords[i] + all_coords[j]) / 2
            ci  = _mpl_cpk.get(all_syms[i], "#CC66CC")
            cj  = _mpl_cpk.get(all_syms[j], "#CC66CC")
            for start, end, col in [(all_coords[i], mid, ci), (mid, all_coords[j], cj)]:
                self.ax.plot(
                    [start[0], end[0]], [start[1], end[1]], [start[2], end[2]],
                    color=col, linewidth=0.8, alpha=0.9,
                )

        self.ax.set_axis_off()
        self.ax.set_title(
            f"({nt.chirality.n},{nt.chirality.m})  "
            f"D={nt.diameter:.2f} Å  atoms={nt.n_atoms}",
            color="white", fontsize=8,
        )
        self.fig.tight_layout()
        self.canvas.draw()


# ─────────────────────────────────────────────────────────────────────────────
# Systematize dialog
# ─────────────────────────────────────────────────────────────────────────────

class _SystematizeDialog(QDialog):
    """Dialog for batch nanotube generation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Batch Export")
        self.resize(360, 280)
        lay = QVBoxLayout(self)

        form = QFormLayout()

        self._dmin = QDoubleSpinBox()
        self._dmin.setRange(0.0, 99.0)
        self._dmin.setValue(0.0)
        self._dmin.setSuffix(" Å")
        self._dmin.setToolTip(
            "Minimum nanotube diameter to generate.\n"
            "For multi-layer structures this is typically 2 × max z-offset."
        )

        self._dmax = QDoubleSpinBox()
        self._dmax.setRange(1.0, 200.0)
        self._dmax.setValue(20.0)
        self._dmax.setSuffix(" Å")
        self._dmax.setToolTip(
            "Maximum nanotube diameter to generate.\n"
            "n_max and m_max are derived automatically from this value."
        )

        self._amax = QSpinBox()
        self._amax.setRange(100, 500000)
        self._amax.setValue(50000)
        self._amax.setToolTip("Skip nanotubes with more atoms than this limit.")

        self._vac  = QDoubleSpinBox()
        self._vac.setRange(0.0, 200.0)
        self._vac.setValue(50.0)
        self._vac.setSuffix(" Å")

        self._fmt = QComboBox()
        self._fmt.addItems([".pdb", ".xyz", "poscar", ".lammps", ".pwi"])

        self._dir_edit = QLineEdit()
        self._dir_edit.setPlaceholderText("Output folder…")
        btn_dir = QPushButton("Browse…")
        btn_dir.clicked.connect(self._browse)
        dir_row = QHBoxLayout()
        dir_row.addWidget(self._dir_edit)
        dir_row.addWidget(btn_dir)

        form.addRow("Min D (Å):",     self._dmin)
        form.addRow("Max D (Å):",     self._dmax)
        form.addRow("Max atoms:",     self._amax)
        form.addRow("Vacuum (Å):",    self._vac)
        form.addRow("Format:",        self._fmt)
        form.addRow("Output folder:", dir_row)

        lay.addLayout(form)

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Select output folder")
        if d:
            self._dir_edit.setText(d)

    def params(self) -> dict:
        return {
            "min_diameter": self._dmin.value(),
            "max_diameter": self._dmax.value(),
            "max_atoms":    self._amax.value(),
            "vacuum":       self._vac.value(),
            "format":       self._fmt.currentText(),
            "output_dir":   self._dir_edit.text() or ".",
        }
