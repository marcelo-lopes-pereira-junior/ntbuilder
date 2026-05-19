"""
gui/panels/input_panel.py
--------------------------
Left panel: file upload, structure info, n/m controls, build button.

Signals emitted
---------------
  structure_loaded(LatticeStructure)  — after successful file load
  build_requested(n, m, vacuum)       — when user clicks Build
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore    import pyqtSignal, Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QFormLayout,
    QPushButton, QSpinBox, QDoubleSpinBox, QLabel,
    QGroupBox, QFileDialog, QFrame, QCheckBox, QScrollArea,
)
from PyQt6.QtGui import QFont, QPixmap
from pathlib import Path

from core import load_structure, LatticeStructure
from core.symmetry import snap_to_symmetry, find_primitive_cell
from gui.utils import ScalablePixmapLabel

_ASSETS = Path(__file__).parent.parent.parent / "assets"


class InputPanel(QWidget):
    """Controls for loading a structure and specifying (n, m)."""

    structure_loaded = pyqtSignal(object)          # LatticeStructure
    build_requested  = pyqtSignal(int, int, float, bool)  # n, m, vacuum, roll_inward
    # n, m, vacuum, roll_inward, n_walls, interlayer_spacing
    mwnt_requested   = pyqtSignal(int, int, float, bool, int, float)
    # Emitted when the "Roll inward" checkbox toggles.  The polar map panel
    # listens to this so it can re-evaluate the curvature-induced spurious
    # bond check for every (n, m) — flipping the rolling direction swaps
    # which face of a buckled / Janus monolayer sits on the concave side,
    # and so changes which species pairs end up close enough to bond.
    roll_direction_changed = pyqtSignal(bool)      # new roll_inward state

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(230)
        self.setMaximumWidth(310)
        self._structure: LatticeStructure | None = None
        self._build_ui()

    # ─────────────────────────────────────────────────────────────────────────
    # UI
    # ─────────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        # ── Scroll wrapper so the panel works on small/laptop screens ─────────
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        inner = QWidget()
        scroll.setWidget(inner)

        outer_lay = QVBoxLayout(self)
        outer_lay.setContentsMargins(0, 0, 0, 0)
        outer_lay.addWidget(scroll)

        root = QVBoxLayout(inner)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(10)

        # ── Sidebar logo (hidden until first file is loaded) ─────────────────
        # Prefer the exported PNG (pixel-perfect, no render delay).
        # ScalablePixmapLabel auto-scales to fit the widget on resize.
        from PyQt6.QtWidgets import QHBoxLayout as _QHBox
        from PyQt6.QtWidgets import QSizePolicy as _QSP
        logo_png = _ASSETS / "ntbuilder.png"
        _tmp = QPixmap(str(logo_png)) if logo_png.exists() else QPixmap()
        logo_pix: QPixmap | None = _tmp if not _tmp.isNull() else None
        self._logo_widget = ScalablePixmapLabel(logo_pix)
        if logo_pix is None or logo_pix.isNull():
            self._logo_widget.setText("NTBuilder")
            self._logo_widget.setFont(QFont("Arial", 13, QFont.Weight.Bold))
        # Let the logo grow/shrink with the sidebar width, capped at a sensible height.
        self._logo_widget.setMinimumHeight(50)
        self._logo_widget.setMaximumHeight(160)
        self._logo_widget.setSizePolicy(_QSP.Policy.Expanding, _QSP.Policy.Preferred)
        self._logo_widget.hide()   # shown only after a file is loaded
        icon_row = _QHBox()
        icon_row.setContentsMargins(0, 4, 0, 4)
        icon_row.addStretch()
        icon_row.addWidget(self._logo_widget)
        icon_row.addStretch()
        root.addLayout(icon_row)

        _add_separator(root)

        # ── File group ───────────────────────────────────────────────────────
        file_box = QGroupBox("1. Structure File")
        file_lay = QVBoxLayout(file_box)

        self.btn_open = QPushButton("Open .cif / .pdb / .xyz …")
        self.btn_open.setToolTip(
            "Load a 2D crystal structure.\n"
            ".cif requires gemmi (pip install gemmi)."
        )
        self.btn_open.clicked.connect(self._on_open_file)
        file_lay.addWidget(self.btn_open)

        # ── Database query button (kept in code, hidden from UI) ─────────────
        # The COD/Materials Project/C2DB integration is functional but the
        # remote APIs are unstable enough that the feature is currently
        # disabled in the interface.  Listed as a planned extension in the
        # manuscript.  To re-enable, set _ENABLE_DB_QUERY=True below.
        _ENABLE_DB_QUERY = False
        self.btn_db = QPushButton("🔍  Query Database…")
        self.btn_db.setToolTip(
            "Search the Crystallography Open Database (COD) or Materials Project\n"
            "for structures and download them directly."
        )
        self.btn_db.clicked.connect(self._on_database_query)
        if _ENABLE_DB_QUERY:
            file_lay.addWidget(self.btn_db)
        else:
            self.btn_db.hide()

        self.lbl_file = QLabel("No file loaded")
        self.lbl_file.setWordWrap(True)
        self.lbl_file.setStyleSheet("color: #888888; font-size: 10px;")
        file_lay.addWidget(self.lbl_file)

        root.addWidget(file_box)

        # ── Structure info ───────────────────────────────────────────────────
        info_box = QGroupBox("Structure Info")
        info_lay = QFormLayout(info_box)
        info_lay.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        info_lay.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self.lbl_a      = _info_label()
        self.lbl_b      = _info_label()
        self.lbl_gamma  = _info_label()
        self.lbl_type   = _info_label()
        self.lbl_natoms = _info_label()

        info_lay.addRow("a (Å):",      self.lbl_a)
        info_lay.addRow("b (Å):",      self.lbl_b)
        info_lay.addRow("γ (°):",       self.lbl_gamma)
        info_lay.addRow("Lattice:",    self.lbl_type)
        info_lay.addRow("Atoms/cell:", self.lbl_natoms)
        root.addWidget(info_box)

        # ── Symmetry tools ───────────────────────────────────────────────────
        self._sym_box = QGroupBox("Symmetry Tools")
        sym_lay = QVBoxLayout(self._sym_box)
        sym_lay.setSpacing(4)

        self.btn_prim = QPushButton("Find primitive cell")
        self.btn_prim.setObjectName("btn_prim")   # targeted by global stylesheet
        self.btn_prim.setEnabled(False)
        self.btn_prim.setToolTip(
            "Find the smallest unit cell that reproduces this crystal.\n"
            "Useful when your input is a conventional or supercell\n"
            "(e.g. rectangular graphene → hexagonal primitive cell).\n"
            "Uses spglib if available, otherwise a built-in search."
        )
        self.btn_prim.clicked.connect(self._on_find_primitive)
        sym_lay.addWidget(self.btn_prim)

        self.lbl_sym_status = QLabel("")
        self.lbl_sym_status.setWordWrap(True)
        self.lbl_sym_status.setStyleSheet(
            "font-size: 9px; color: #555555; padding: 1px 2px;"
        )
        sym_lay.addWidget(self.lbl_sym_status)

        root.addWidget(self._sym_box)

        # ── Chiral indices ───────────────────────────────────────────────────
        nm_box = QGroupBox("2. Chiral Indices")
        nm_lay = QFormLayout(nm_box)
        nm_lay.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self.spin_n = QSpinBox()
        self.spin_n.setRange(0, 200)
        self.spin_n.setValue(5)
        self.spin_n.setToolTip("Chiral index n")
        self.spin_n.valueChanged.connect(self._on_indices_changed)

        self.spin_m = QSpinBox()
        self.spin_m.setRange(0, 200)
        self.spin_m.setValue(0)
        self.spin_m.setToolTip("Chiral index m")
        self.spin_m.valueChanged.connect(self._on_indices_changed)

        nm_lay.addRow("n:", self.spin_n)
        nm_lay.addRow("m:", self.spin_m)

        # Live preview
        self.lbl_preview = QLabel("—")
        self.lbl_preview.setStyleSheet(
            "font-size: 10px; color: #555555; padding: 2px;"
        )
        self.lbl_preview.setWordWrap(True)
        nm_lay.addRow("Preview:", self.lbl_preview)

        root.addWidget(nm_box)

        # ── Options ──────────────────────────────────────────────────────────
        opt_box = QGroupBox("3. Options")
        opt_lay = QFormLayout(opt_box)
        opt_lay.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self.spin_vacuum = QDoubleSpinBox()
        self.spin_vacuum.setRange(0.0, 200.0)
        self.spin_vacuum.setValue(50.0)
        self.spin_vacuum.setSingleStep(5.0)
        self.spin_vacuum.setSuffix(" Å")
        self.spin_vacuum.setToolTip(
            "Vacuum padding added around the nanotube in the xy plane.\n"
            "50 Å is a safe default for DFT calculations to avoid\n"
            "spurious inter-nanotube interactions."
        )
        opt_lay.addRow("Vacuum:", self.spin_vacuum)

        self.spin_search = QSpinBox()
        self.spin_search.setRange(10, 2000)
        self.spin_search.setValue(300)
        self.spin_search.setToolTip(
            "T-vector search limit.\n\n"
            "Every nanotube has a chiral vector Ch (along the\n"
            "circumference) and a translational vector T (along\n"
            "the tube axis). T must be perpendicular to Ch.\n\n"
            "For ideal lattices (exact γ) this solution is exact.\n"
            "For strained cases the algorithm searches integer\n"
            "pairs (t₁, t₂) up to this limit to minimise |Ch·T|.\n\n"
            "Higher limit → lower strain, but the nanotube unit\n"
            "cell may contain many more atoms and be slower to build.\n"
            "300 is sufficient for most real-world structures."
        )
        opt_lay.addRow("T search:", self.spin_search)

        # Roll direction — shown only for buckled structures
        self._roll_row_lbl = QLabel("Roll dir.:")
        self.chk_roll_in   = QCheckBox("Roll inward")
        self.chk_roll_in.setToolTip(
            "Only relevant for buckled/layered structures (e.g. pentagraphene, MoS₂).\n"
            "Checked  → atoms with z > 0 go to smaller radius (inner wall).\n"
            "Unchecked → atoms with z > 0 go to larger radius (outer wall).\n"
            "For Janus materials (MoSSe) this swaps which species faces inward."
        )
        self._roll_row_lbl.setVisible(False)
        self.chk_roll_in.setVisible(False)
        opt_lay.addRow(self._roll_row_lbl, self.chk_roll_in)
        # Forward toggle state to the rest of the application so the polar
        # map can re-evaluate which chiralities develop spurious bonds.
        self.chk_roll_in.toggled.connect(self.roll_direction_changed.emit)

        # Buckling info label — shown when structure has z-offsets.
        # MinimumExpanding vertical policy is required so that QFormLayout
        # allocates enough height for the wrapped multi-line text.
        self._lbl_buckling = QLabel()
        self._lbl_buckling.setStyleSheet(
            "color: #E69F00; font-size: 9px; padding: 2px;"
        )
        self._lbl_buckling.setWordWrap(True)
        from PyQt6.QtWidgets import QSizePolicy as _QSP2
        self._lbl_buckling.setSizePolicy(
            _QSP2.Policy.Expanding, _QSP2.Policy.MinimumExpanding
        )
        self._lbl_buckling.setVisible(False)
        opt_lay.addRow("", self._lbl_buckling)

        root.addWidget(opt_box)

        # ── Build button ─────────────────────────────────────────────────────
        self.btn_build = QPushButton("▶  Build Nanotube")
        self.btn_build.setObjectName("btn_build")   # targeted by global stylesheet
        self.btn_build.setFixedHeight(38)
        self.btn_build.setEnabled(False)
        self.btn_build.clicked.connect(self._on_build)
        root.addWidget(self.btn_build)

        root.addStretch()

        # Footer
        foot = QLabel("Tip: click any point on the polar map\nto select (n, m) automatically.")
        foot.setStyleSheet("color:#999999; font-size:9px;")
        foot.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(foot)

    # ─────────────────────────────────────────────────────────────────────────
    # Slots
    # ─────────────────────────────────────────────────────────────────────────

    def _on_open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Structure File", "",
            "Structure files (*.cif *.pdb *.xyz *.poscar *.contcar *.vasp "
            "*.xsf *.lammps *.data *.in *.pwi);;"
            "CIF (*.cif);;"
            "VASP POSCAR (*.poscar *.contcar *.vasp);;"
            "XSF / XCrysDen (*.xsf);;"
            "LAMMPS data (*.lammps *.data);;"
            "Quantum ESPRESSO (*.in *.pwi);;"
            "PDB (*.pdb);;"
            "XYZ (*.xyz);;"
            "All files (*)"
        )
        if not path:
            return
        try:
            from pathlib import Path
            ext = Path(path).suffix.lower()

            # For plain XYZ we need lattice vectors — ask core to try, then prompt
            if ext == ".xyz":
                struct = self._load_xyz_with_dialog(path)
            else:
                struct = load_structure(path)

            if struct is None:
                return

            # Auto-snap to ideal symmetry on every import.
            # This silently fixes floating-point deviations (e.g. γ = 60.0001°)
            # so downstream chirality computations always use exact parameters.
            try:
                snapped, snap_desc = snap_to_symmetry(struct)
                if "nothing to snap" not in snap_desc:
                    struct = snapped
                    self.lbl_sym_status.setText(f"↺ {snap_desc}")
                else:
                    self.lbl_sym_status.setText("")
            except Exception:
                pass   # if snap fails for any reason, use the original

            self._structure = struct
            self._update_info(struct)
            self.lbl_file.setText(Path(path).name)
            self.btn_build.setEnabled(True)
            self.structure_loaded.emit(struct)

        except Exception as exc:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Load error", str(exc))

    def _load_xyz_with_dialog(self, path: str):
        """Try extended XYZ first; if that fails, ask the user for a, b."""
        try:
            return load_structure(path)
        except ValueError:
            pass

        # Ask for lattice parameters
        dlg = _LatticeParamsDialog(self)
        if dlg.exec():
            a1, a2 = dlg.get_vectors()
            return load_structure(path, a1=a1, a2=a2)
        return None

    def _update_info(self, struct: LatticeStructure):
        self.lbl_a.setText(f"{struct.a:.4f}")
        self.lbl_b.setText(f"{struct.b:.4f}")
        self.lbl_gamma.setText(f"{struct.gamma_deg:.2f}")
        self.lbl_type.setText(struct.lattice_type.capitalize())
        self.lbl_natoms.setText(str(len(struct.atoms)))
        self.lbl_sym_status.setText("")

        # Enable the primitive-cell button for every loaded structure.
        # The actual search only runs when the user clicks the button.
        self.btn_prim.setEnabled(True)
        self.btn_prim.setToolTip(
            "Find the smallest unit cell that reproduces this crystal.\n"
            "Useful when your input is a conventional or supercell\n"
            "(e.g. rectangular graphene → hexagonal primitive cell).\n"
            "Uses spglib if available, otherwise a built-in search."
        )

        # Show/hide rolling direction controls
        buckled = struct.has_buckling
        self._roll_row_lbl.setVisible(buckled)
        self.chk_roll_in.setVisible(buckled)
        if buckled:
            dz = struct.max_z_offset
            species = sorted(set(a["symbol"] for a in struct.atoms))
            self._lbl_buckling.setText(
                f"Buckled structure detected (Δz ≈ {dz:.2f} Å). "
                "Choose rolling direction above."
            )
        self._lbl_buckling.setVisible(buckled)

        self._on_indices_changed()

    def _on_indices_changed(self):
        if self._structure is None:
            self.lbl_preview.setText("—")
            return
        n, m = self.spin_n.value(), self.spin_m.value()
        if n == 0 and m == 0:
            self.lbl_preview.setText("Invalid: n = m = 0")
            return
        try:
            from core import compute_chirality
            ch = compute_chirality(
                n, m, self._structure,
                search_limit=self.spin_search.value()
            )
            self.lbl_preview.setText(
                f"D = {ch.diameter:.3f} Å\n"
                f"θ = {ch.theta_deg:.2f}°\n"
                f"atoms = {ch.n_atoms}\n"
                f"strain = {ch.strain:.4f}%"
            )
        except Exception:
            self.lbl_preview.setText("—")

    def _on_find_primitive(self):
        if self._structure is None:
            return
        try:
            orig_type = self._structure.lattice_type
            new_struct, desc = find_primitive_cell(self._structure)

            # Sanity check: primitive-cell reduction should never change the
            # lattice type from oblique to rectangular (or vice versa) unless
            # the atom count also decreased, which implies a genuine reduction.
            # When atom count is unchanged and the type changed, spglib has
            # likely mis-classified the cell; we keep the original in that case.
            new_type = new_struct.lattice_type
            atoms_changed = len(new_struct.atoms) != len(self._structure.atoms)
            unexpected_type_change = (
                not atoms_changed
                and orig_type == "oblique"
                and new_type != "oblique"
            )
            if unexpected_type_change:
                self.lbl_sym_status.setText(
                    f"⚠ Primitive cell search returned an unexpected lattice "
                    f"type ({orig_type} → {new_type}) with unchanged atom count. "
                    f"Original structure kept."
                )
                return

            self._structure = new_struct
            self._update_info(new_struct)
            self.lbl_sym_status.setText(f"✓ {desc}")
            self.structure_loaded.emit(new_struct)
        except Exception as exc:
            self.lbl_sym_status.setText(f"⚠ {exc}")

    def _on_build(self):
        n           = self.spin_n.value()
        m           = self.spin_m.value()
        vacuum      = self.spin_vacuum.value()
        roll_inward = self.chk_roll_in.isChecked()
        self.build_requested.emit(n, m, vacuum, roll_inward)

    def _on_database_query(self):
        """Open the database query dialog; load downloaded structure if accepted."""
        from PyQt6.QtWidgets import QMessageBox, QDialog, QVBoxLayout, QLabel, QDialogButtonBox, QFrame
        try:
            from gui.dialogs.advanced_dialogs import DatabaseQueryDialog
        except ImportError as e:
            QMessageBox.warning(self, "Not available", str(e))
            return

        dlg = DatabaseQueryDialog(parent=self)
        if dlg.exec() != dlg.DialogCode.Accepted:
            return

        result_path = dlg.downloaded_path()
        if not result_path:
            return

        try:
            struct = load_structure(str(result_path))
            try:
                from core.symmetry import snap_to_symmetry
                snapped, snap_desc = snap_to_symmetry(struct)
                if "nothing to snap" not in snap_desc:
                    struct = snapped
                    self.lbl_sym_status.setText(f"↺ {snap_desc}")
            except Exception:
                pass

            # ── Structure validation ──────────────────────────────────────────
            # Show a quick inspection dialog so the user can decide whether the
            # imported structure makes sense for nanotube rolling.
            proceed = self._validate_imported_structure(struct, Path(result_path).name)
            if not proceed:
                return   # user cancelled

            self._structure = struct
            self._update_info(struct)
            self.lbl_file.setText(Path(result_path).name)
            self.btn_build.setEnabled(True)
            self.structure_loaded.emit(struct)
        except Exception as exc:
            QMessageBox.critical(self, "Load error", str(exc))

    def _validate_imported_structure(self, struct, filename: str) -> bool:
        """
        Show a structure-metrics dialog after importing from the database.

        Returns True if the user wants to proceed, False to cancel.
        Always shows metrics; prominently warns when the structure looks 3D-like.
        """
        from PyQt6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QLabel,
            QDialogButtonBox, QFrame, QGridLayout,
        )
        from PyQt6.QtGui import QFont, QColor
        from PyQt6.QtCore import Qt

        n_atoms   = len(struct.atoms)
        max_z     = struct.max_z_offset
        a_val     = struct.a
        b_val     = struct.b
        gamma_val = struct.gamma_deg

        # Detect suspicious properties
        warnings: list[str] = []
        if n_atoms > 20:
            warnings.append(
                f"⚠  {n_atoms} atoms per unit cell — this may be a bulk crystal or supercell.\n"
                "   For rolling into a nanotube you usually want ≤ 8 atoms/cell.\n"
                "   Consider finding the primitive cell (Symmetry Tools → Find primitive cell)."
            )
        if max_z > 2.5:
            warnings.append(
                f"⚠  Large out-of-plane offset ({max_z:.2f} Å) — this looks like a 3D bulk\n"
                "   structure, not a 2D layer.  NTBuilder extracted one z-slice but the\n"
                "   result may not be physically meaningful."
            )
        elif max_z > 0.5:
            warnings.append(
                f"ℹ  Buckled/layered structure detected (max Δz = {max_z:.2f} Å).\n"
                "   This is normal for MoS₂, silicene, etc.  Use 'Roll inward' option\n"
                "   to control which species faces toward the tube axis."
            )

        dlg = QDialog(self)
        dlg.setWindowTitle("Imported Structure — Inspection")
        dlg.setMinimumWidth(480)
        dlg.setStyleSheet("""
            QDialog { background-color: #F4F6FB; }
            QLabel  { color: #1A3A6B; background: transparent; }
            QFrame[frameShape="4"] { color: #C4CDE0; }
        """)
        lay = QVBoxLayout(dlg)
        lay.setSpacing(10)
        lay.setContentsMargins(18, 14, 18, 14)

        # Title
        title = QLabel(f"Structure imported: <b>{filename}</b>")
        title.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(title)

        # Separator
        sep = QFrame(); sep.setFrameShape(QFrame.Shape.HLine)
        lay.addWidget(sep)

        # Metrics grid
        grid = QGridLayout()
        grid.setHorizontalSpacing(20)
        grid.setVerticalSpacing(4)
        metrics = [
            ("Atoms / cell",    str(n_atoms)),
            ("a (Å)",           f"{a_val:.4f}"),
            ("b (Å)",           f"{b_val:.4f}"),
            ("γ (°)",           f"{gamma_val:.2f}"),
            ("Lattice type",    struct.lattice_type.capitalize()),
            ("Max Δz (Å)",      f"{max_z:.3f}"),
        ]
        for i, (lbl_txt, val_txt) in enumerate(metrics):
            lbl = QLabel(lbl_txt + ":")
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight)
            val = QLabel(f"<b>{val_txt}</b>")
            val.setTextFormat(Qt.TextFormat.RichText)
            grid.addWidget(lbl, i, 0)
            grid.addWidget(val, i, 1)
        lay.addLayout(grid)

        # Warning box (only if there are warnings)
        if warnings:
            sep2 = QFrame(); sep2.setFrameShape(QFrame.Shape.HLine)
            lay.addWidget(sep2)
            bg_color = "#F8D7DA" if any("⚠" in w for w in warnings) else "#FFF3CD"
            for w in warnings:
                warn_lbl = QLabel(w)
                warn_lbl.setWordWrap(True)
                warn_lbl.setStyleSheet(
                    f"background:{bg_color}; color:#7B2D00; "
                    "border:1px solid #EEC0C0; border-radius:5px; padding:6px 8px; font-size:10px;"
                )
                lay.addWidget(warn_lbl)

        # Buttons
        if warnings:
            btn_box = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            )
            btn_box.button(QDialogButtonBox.StandardButton.Ok).setText("Proceed anyway")
            btn_box.button(QDialogButtonBox.StandardButton.Cancel).setText("Cancel import")
        else:
            btn_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
            btn_box.button(QDialogButtonBox.StandardButton.Ok).setText("Use this structure")
        btn_box.accepted.connect(dlg.accept)
        btn_box.rejected.connect(dlg.reject)
        lay.addWidget(btn_box)

        return dlg.exec() == QDialog.DialogCode.Accepted

    # ─────────────────────────────────────────────────────────────────────────
    # Public API
    # ─────────────────────────────────────────────────────────────────────────

    def set_indices(self, n: int, m: int):
        """Called when the user clicks a point on the polar map."""
        self.spin_n.setValue(n)
        self.spin_m.setValue(m)

    def show_logo(self):
        """Reveal the sidebar icon (called after the first file is loaded)."""
        self._logo_widget.show()


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _info_label() -> QLabel:
    lbl = QLabel("—")
    lbl.setStyleSheet("font-family: monospace; font-size: 10px;")
    return lbl


def _add_separator(layout: QVBoxLayout):
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    layout.addWidget(line)


# ─────────────────────────────────────────────────────────────────────────────
# Dialog: ask for lattice vectors when XYZ has no cell info
# ─────────────────────────────────────────────────────────────────────────────

class _LatticeParamsDialog(QWidget):
    """Simple dialog for entering a, b and γ lattice parameters.

    Uses plain QLineEdit fields so that:
    - No graphene default values are pre-filled
    - The decimal separator is always a dot (no locale-dependent comma)
    - The dialog background is explicitly styled, independent of the parent theme
    """

    def __init__(self, parent=None):
        from PyQt6.QtWidgets import (
            QDialog, QDialogButtonBox, QFormLayout, QVBoxLayout, QLineEdit,
            QMessageBox,
        )
        from PyQt6.QtGui import QDoubleValidator
        from PyQt6.QtCore import QLocale

        self._dlg = QDialog(parent)
        self._dlg.setWindowTitle("Lattice Parameters")
        self._dlg.setFixedSize(300, 180)

        # ── Force a clean light look regardless of the parent palette ──────────
        self._dlg.setStyleSheet("""
            QDialog {
                background-color: #F4F6FB;
            }
            QLabel {
                color: #1A3A6B;
                font-weight: 500;
            }
            QLineEdit {
                background-color: #FFFFFF;
                color: #111827;
                border: 1px solid #C4CDE0;
                border-radius: 5px;
                padding: 4px 8px;
                font-size: 13px;
            }
            QLineEdit:focus {
                border: 1.5px solid #2851A3;
            }
            QLineEdit[invalid="true"] {
                border: 1.5px solid #C0392B;
                background-color: #FFF0EE;
            }
        """)

        # ── Locale-neutral validator: always uses dot as decimal separator ──────
        _c_locale = QLocale(QLocale.Language.C)

        def _make_field(placeholder: str, bottom: float, top: float) -> QLineEdit:
            f = QLineEdit()
            f.setPlaceholderText(placeholder)
            v = QDoubleValidator(bottom, top, 6, f)
            v.setLocale(_c_locale)
            v.setNotation(QDoubleValidator.Notation.StandardNotation)
            f.setValidator(v)
            return f

        self._a     = _make_field("e.g.  2.4600",  0.1, 100.0)
        self._b     = _make_field("e.g.  2.4600",  0.1, 100.0)
        self._gamma = _make_field("e.g.  60.00",   1.0, 179.0)

        lay = QVBoxLayout(self._dlg)
        lay.setContentsMargins(16, 16, 16, 12)
        lay.setSpacing(8)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        form.setHorizontalSpacing(12)
        form.addRow("a (Å):", self._a)
        form.addRow("b (Å):", self._b)
        form.addRow("γ (°):", self._gamma)
        lay.addLayout(form)

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self._try_accept)
        btns.rejected.connect(self._dlg.reject)
        lay.addWidget(btns)

        self._QMessageBox = QMessageBox

    def _parse(self, field) -> float | None:
        """Return float from field text (accepts dot or comma), or None."""
        txt = field.text().strip().replace(",", ".")
        try:
            v = float(txt)
            return v if v > 0 else None
        except ValueError:
            return None

    def _try_accept(self):
        ok = True
        for field, lo, hi in [
            (self._a,     0.1, 100.0),
            (self._b,     0.1, 100.0),
            (self._gamma, 1.0, 179.0),
        ]:
            v = self._parse(field)
            bad = v is None or not (lo <= v <= hi)
            field.setProperty("invalid", "true" if bad else "false")
            field.style().unpolish(field)
            field.style().polish(field)
            if bad:
                ok = False
        if not ok:
            self._QMessageBox.warning(
                self._dlg, "Invalid input",
                "Please enter valid numbers (use a dot as decimal separator).\n"
                "a, b: 0.1 – 100 Å   |   γ: 1 – 179°"
            )
            return
        self._dlg.accept()

    def exec(self) -> bool:
        from PyQt6.QtWidgets import QDialog
        return self._dlg.exec() == QDialog.DialogCode.Accepted

    def get_vectors(self):
        import math
        a = self._parse(self._a)
        b = self._parse(self._b)
        g = math.radians(self._parse(self._gamma))
        a1 = np.array([a, 0.0])
        a2 = np.array([b * math.cos(g), b * math.sin(g)])
        return a1, a2
