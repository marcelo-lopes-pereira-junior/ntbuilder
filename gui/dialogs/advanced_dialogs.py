"""
gui/dialogs/advanced_dialogs.py
--------------------------------
All advanced feature dialogs: MWNT, Bundle, Deformations, Analysis,
Methods text, DFT input generator, and Database query.
"""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import numpy as np
from PyQt6.QtCore    import Qt, QThread, pyqtSignal as Signal
from PyQt6.QtGui     import QFont, QClipboard, QGuiApplication, QColor
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QPushButton, QSpinBox, QDoubleSpinBox,
    QComboBox, QTextEdit, QGroupBox, QTabWidget,
    QDialogButtonBox, QCheckBox, QLineEdit, QFileDialog,
    QProgressBar, QListWidget, QListWidgetItem, QMessageBox,
    QSizePolicy, QFrame, QScrollArea, QWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
)

from core.builder import NanotubeStructure


# ─────────────────────────────────────────────────────────────────────────────
# Shared dialog stylesheet
# ─────────────────────────────────────────────────────────────────────────────

_DIALOG_STYLE = """
    QDialog {
        background-color: #F4F6FB;
    }
    QGroupBox {
        background-color: #F4F6FB;
        border: 1px solid #C4CDE0;
        border-radius: 6px;
        margin-top: 8px;
        font-weight: 600;
        color: #1A3A6B;
    }
    QGroupBox::title {
        subcontrol-origin: margin;
        left: 8px;
        padding: 0 4px;
    }
    QLabel {
        color: #1A3A6B;
        background-color: transparent;
    }
    QLineEdit, QTextEdit {
        background-color: #FFFFFF;
        color: #111827;
        border: 1px solid #C4CDE0;
        border-radius: 5px;
        padding: 3px 7px;
    }
    QLineEdit:focus, QTextEdit:focus {
        border: 1.5px solid #2851A3;
    }
    QSpinBox, QDoubleSpinBox {
        background-color: #FFFFFF;
        color: #111827;
        border: 1px solid #C4CDE0;
        border-radius: 5px;
        padding: 2px 6px;
    }
    QComboBox {
        background-color: #FFFFFF;
        color: #111827;
        border: 1px solid #C4CDE0;
        border-radius: 5px;
        padding: 2px 6px;
    }
    QComboBox QAbstractItemView {
        background-color: #FFFFFF;
        color: #111827;
        selection-background-color: #D0E0FF;
    }
    QListWidget {
        background-color: #FFFFFF;
        color: #111827;
        border: 1px solid #C4CDE0;
        border-radius: 5px;
        alternate-background-color: #F0F4FA;
    }
    QListWidget::item:selected {
        background-color: #2851A3;
        color: #FFFFFF;
    }
    QTabWidget::pane {
        background-color: #F4F6FB;
        border: 1px solid #C4CDE0;
    }
    QTabBar::tab {
        background-color: #E2E9F5;
        color: #1A3A6B;
        padding: 5px 14px;
        border-top-left-radius: 4px;
        border-top-right-radius: 4px;
        border: 1px solid #C4CDE0;
        border-bottom: none;
    }
    QTabBar::tab:selected {
        background-color: #FFFFFF;
        font-weight: 600;
    }
    QProgressBar {
        background-color: #E2E9F5;
        border: 1px solid #C4CDE0;
        border-radius: 4px;
        text-align: center;
        color: #1A3A6B;
    }
    QProgressBar::chunk {
        background-color: #2851A3;
        border-radius: 3px;
    }
    QCheckBox {
        color: #1A3A6B;
        background-color: transparent;
    }
    QScrollArea {
        background-color: #F4F6FB;
        border: none;
    }
"""


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _divider():
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


def _copy_btn(text_edit: QTextEdit) -> QPushButton:
    btn = QPushButton("📋 Copy to clipboard")
    btn.clicked.connect(lambda: QGuiApplication.clipboard().setText(text_edit.toPlainText()))
    return btn


def _save_btn(text_edit: QTextEdit, default_name: str, parent=None) -> QPushButton:
    btn = QPushButton("💾 Save…")
    def _save():
        path, _ = QFileDialog.getSaveFileName(parent, "Save file", default_name)
        if path:
            Path(path).write_text(text_edit.toPlainText(), encoding="utf-8")
    btn.clicked.connect(_save)
    return btn


# ─────────────────────────────────────────────────────────────────────────────
# Wall picker dialog  (shown once per outer wall when building a MWNT)
# ─────────────────────────────────────────────────────────────────────────────

class WallPickerDialog(QDialog):
    """
    Let the user choose among ~10 candidate chiralities for one MWNT wall.

    The table shows, for each candidate:
      (n, m) | Diameter | Δ spacing | T period | Z-ratio | Z-strain | Score

    Columns are sortable.  The strain cell is colour-coded:
      green  < 5 %    excellent commensurability
      yellow 5–20 %   acceptable
      red    > 20 %   heavy strain — consider a different candidate
    """

    # ── column indices ────────────────────────────────────────────────────────
    _COL_NM      = 0
    _COL_DIAM    = 1
    _COL_DELTA   = 2
    _COL_T       = 3
    _COL_ZRATIO  = 4
    _COL_STRAIN  = 5
    _COL_SCORE   = 6

    _HEADERS = ["(n, m)", "Diameter (Å)", "Δ spacing (Å)",
                "T period (Å)", "Z-ratio", "Z-strain (%)", "Score"]

    def __init__(self, wall_index: int,
                 inner_n: int, inner_m: int,
                 inner_T: float,
                 candidates,          # list[WallCandidate]
                 requested_spacing: float,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Select Wall {wall_index + 1} chirality")
        self.resize(740, 420)
        self.setStyleSheet(_DIALOG_STYLE)
        self._candidates = candidates
        self._selected   = candidates[0] if candidates else None
        self._build_ui(wall_index, inner_n, inner_m, inner_T, requested_spacing)

    def _build_ui(self, wall_index, inner_n, inner_m, inner_T, req_spacing):
        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        lay.setContentsMargins(16, 12, 16, 12)

        # ── Header ────────────────────────────────────────────────────────────
        hdr = QLabel(
            f"<b>Wall {wall_index + 1}</b> — inner wall is "
            f"<b>({inner_n}, {inner_m})</b>  T = {inner_T:.3f} Å<br>"
            f"<small style='color:#555;'>Requested spacing: {req_spacing:.2f} Å  ·  "
            f"Click a row to select  ·  Click a column header to sort</small>"
        )
        hdr.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(hdr)

        # ── Table ─────────────────────────────────────────────────────────────
        self._table = QTableWidget(len(self._candidates), len(self._HEADERS))
        self._table.setHorizontalHeaderLabels(self._HEADERS)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        self._table.verticalHeader().setVisible(False)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(self._COL_NM, QHeaderView.ResizeMode.ResizeToContents)

        self._populate_table()
        self._table.selectRow(0)
        self._table.itemSelectionChanged.connect(self._on_selection)
        lay.addWidget(self._table)

        # ── Legend ────────────────────────────────────────────────────────────
        legend = QLabel(
            "<small>"
            "<span style='background:#C8F0C8; padding:2px 6px;'>■</span> Z-strain &lt; 5 % (excellent) &nbsp;"
            "<span style='background:#FFF3CD; padding:2px 6px;'>■</span> 5–20 % (acceptable) &nbsp;"
            "<span style='background:#F8D7DA; padding:2px 6px;'>■</span> &gt; 20 % (high — prefer another candidate)"
            "</small>"
        )
        legend.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(legend)

        # ── Buttons ───────────────────────────────────────────────────────────
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.button(QDialogButtonBox.StandardButton.Ok).setText("Select this wall")
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    def _populate_table(self):
        from PyQt6.QtGui import QColor, QBrush

        def _num(val, decimals=3):
            item = QTableWidgetItem(f"{val:.{decimals}f}")
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            return item

        for row, cand in enumerate(self._candidates):
            ch = cand.chirality

            nm_item = QTableWidgetItem(f"({ch.n}, {ch.m})")
            nm_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            nm_item.setFont(_bold_font())
            self._table.setItem(row, self._COL_NM,   nm_item)
            self._table.setItem(row, self._COL_DIAM,  _num(ch.diameter, 4))

            delta_item = _num(cand.delta_spacing, 3)
            if abs(cand.delta_spacing) < 0.2:
                delta_item.setForeground(QBrush(QColor("#1A6B1A")))  # green text
            elif abs(cand.delta_spacing) > 1.0:
                delta_item.setForeground(QBrush(QColor("#8B0000")))  # red text
            self._table.setItem(row, self._COL_DELTA, delta_item)

            self._table.setItem(row, self._COL_T,      _num(ch.T_norm, 3))

            ratio_item = _num(cand.z_ratio, 3)
            # How close to an integer?
            frac = abs(cand.z_ratio - round(cand.z_ratio))
            if frac < 0.05:
                ratio_item.setForeground(QBrush(QColor("#1A6B1A")))
            self._table.setItem(row, self._COL_ZRATIO, ratio_item)

            strain_item = _num(cand.z_strain_pct, 2)
            strain_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            if cand.z_strain_pct < 5.0:
                strain_item.setBackground(QBrush(QColor("#C8F0C8")))
            elif cand.z_strain_pct < 20.0:
                strain_item.setBackground(QBrush(QColor("#FFF3CD")))
            else:
                strain_item.setBackground(QBrush(QColor("#F8D7DA")))
            self._table.setItem(row, self._COL_STRAIN, strain_item)

            self._table.setItem(row, self._COL_SCORE, _num(cand.score, 4))

    def _on_selection(self):
        rows = self._table.selectionModel().selectedRows()
        if rows:
            # Map visual row → original candidate (table may be sorted)
            visual_row = rows[0].row()
            nm_text = self._table.item(visual_row, self._COL_NM).text()
            nm_text = nm_text.strip("()")
            n_str, m_str = [x.strip() for x in nm_text.split(",")]
            n, m = int(n_str), int(m_str)
            for cand in self._candidates:
                if cand.chirality.n == n and cand.chirality.m == m:
                    self._selected = cand
                    break

    @property
    def selected(self):
        """The WallCandidate the user selected (or None if cancelled)."""
        return self._selected


def _bold_font():
    from PyQt6.QtGui import QFont as _QFont
    f = _QFont()
    f.setBold(True)
    return f


# ─────────────────────────────────────────────────────────────────────────────
# MWNT dialog
# ─────────────────────────────────────────────────────────────────────────────

class MWNTDialog(QDialog):
    """Configure multi-walled nanotube parameters."""

    def __init__(self, inner_n: int, inner_m: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Build Multi-Walled Nanotube")
        self.setFixedSize(420, 320)
        self.setStyleSheet(_DIALOG_STYLE)
        self._inner_n = inner_n
        self._inner_m = inner_m
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setSpacing(10)
        lay.setContentsMargins(20, 16, 20, 16)

        # Header
        hdr = QLabel(
            f"<b>Innermost wall:</b> ({self._inner_n}, {self._inner_m})<br>"
            "<small style='color:#666;'>Additional walls are found automatically "
            "to match the requested interlayer spacing.</small>"
        )
        hdr.setTextFormat(Qt.TextFormat.RichText)
        hdr.setWordWrap(True)
        lay.addWidget(hdr)
        lay.addWidget(_divider())

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self._spin_walls = QSpinBox()
        self._spin_walls.setRange(2, 10)
        self._spin_walls.setValue(2)
        self._spin_walls.setToolTip(
            "Total number of concentric walls (2 = DWNT, 3 = TWNT, …)."
        )
        form.addRow("Number of walls:", self._spin_walls)

        self._spin_spacing = QDoubleSpinBox()
        self._spin_spacing.setRange(1.5, 20.0)
        self._spin_spacing.setValue(3.40)
        self._spin_spacing.setSuffix(" Å")
        self._spin_spacing.setSingleStep(0.05)
        self._spin_spacing.setDecimals(2)
        self._spin_spacing.setToolTip(
            "Surface-to-surface interlayer gap (Å).\n"
            "Default 3.40 Å is the graphene/C vdW spacing.\n"
            "Use 3.15 Å for h-BN, ~6.2 Å for MoS₂."
        )
        form.addRow("Interlayer spacing:", self._spin_spacing)

        self._spin_vacuum = QDoubleSpinBox()
        self._spin_vacuum.setRange(0.0, 100.0)
        self._spin_vacuum.setValue(10.0)
        self._spin_vacuum.setSuffix(" Å")
        self._spin_vacuum.setSingleStep(5.0)
        form.addRow("Vacuum padding:", self._spin_vacuum)

        self._chk_roll = QCheckBox("Roll inward (for buckled structures)")
        form.addRow("", self._chk_roll)

        lay.addLayout(form)
        lay.addStretch()

        info = QLabel(
            "<small style='color:#888;'>"
            "ℹ  The Z-period of the MWNT equals the inner wall's translational "
            "vector length. Outer walls are tiled to match and trimmed.</small>"
        )
        info.setTextFormat(Qt.TextFormat.RichText)
        info.setWordWrap(True)
        lay.addWidget(info)

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    @property
    def n_walls(self)    -> int:   return self._spin_walls.value()
    @property
    def spacing(self)    -> float: return self._spin_spacing.value()
    @property
    def vacuum(self)     -> float: return self._spin_vacuum.value()
    @property
    def roll_inward(self) -> bool: return self._chk_roll.isChecked()


# ─────────────────────────────────────────────────────────────────────────────
# Bundle dialog
# ─────────────────────────────────────────────────────────────────────────────

class BundleDialog(QDialog):
    """Configure nanotube bundle geometry."""

    _LABELS = {
        "linear":     "Linear (2 tubes along X)",
        "triangle":   "Equilateral triangle (3 tubes)",
        "square4":    "Square 2×2 (4 tubes)",
        "hexagonal7": "Hexagonal — 1+6 (7 tubes)",
        "grid":       "Custom grid (N×M)",
    }

    def __init__(self, diameter: float, vacuum: float = 10.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Build Nanotube Bundle")
        self.setFixedSize(420, 340)
        self.setStyleSheet(_DIALOG_STYLE)
        self._diameter = diameter
        self._default_vacuum = float(vacuum)
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setSpacing(10)
        lay.setContentsMargins(20, 16, 20, 16)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self._cmb_geom = QComboBox()
        for key, label in self._LABELS.items():
            self._cmb_geom.addItem(label, key)
        self._cmb_geom.setCurrentIndex(3)   # hexagonal7 default
        self._cmb_geom.currentIndexChanged.connect(self._on_geom_changed)
        form.addRow("Geometry:", self._cmb_geom)

        self._spin_spacing = QDoubleSpinBox()
        self._spin_spacing.setRange(0.5, 30.0)
        self._spin_spacing.setValue(3.40)
        self._spin_spacing.setSuffix(" Å")
        self._spin_spacing.setSingleStep(0.1)
        self._spin_spacing.setDecimals(2)
        self._spin_spacing.setToolTip("Surface-to-surface intertube gap (Å).")
        form.addRow("Intertube spacing:", self._spin_spacing)

        self._spin_vacuum = QDoubleSpinBox()
        self._spin_vacuum.setRange(0.0, 50.0)
        self._spin_vacuum.setValue(self._default_vacuum)
        self._spin_vacuum.setSuffix(" Å")
        self._spin_vacuum.setSingleStep(0.5)
        self._spin_vacuum.setDecimals(2)
        self._spin_vacuum.setToolTip(
            "Lateral vacuum measured from the surface of the outermost\n"
            "tube to the simulation-box boundary.\n"
            "Use 0 to obtain a tightly packed periodic bundle."
        )
        form.addRow("Outer vacuum:", self._spin_vacuum)

        self._spin_nx = QSpinBox()
        self._spin_nx.setRange(1, 10)
        self._spin_nx.setValue(3)
        self._spin_nx.setEnabled(False)
        form.addRow("Grid columns (nx):", self._spin_nx)

        self._spin_ny = QSpinBox()
        self._spin_ny.setRange(1, 10)
        self._spin_ny.setValue(3)
        self._spin_ny.setEnabled(False)
        form.addRow("Grid rows (ny):", self._spin_ny)

        lay.addLayout(form)

        self._lbl_info = QLabel()
        self._lbl_info.setStyleSheet("color:#555; font-size:10px;")
        self._lbl_info.setWordWrap(True)
        lay.addWidget(self._lbl_info)
        self._update_info()

        lay.addStretch()
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

        self._spin_spacing.valueChanged.connect(self._update_info)

    def _on_geom_changed(self):
        grid = self._cmb_geom.currentData() == "grid"
        self._spin_nx.setEnabled(grid)
        self._spin_ny.setEnabled(grid)
        self._update_info()

    def _update_info(self):
        pitch = self._diameter + self._spin_spacing.value()
        geom  = self._cmb_geom.currentData()
        n_map = {"linear":2,"triangle":3,"square4":4,"hexagonal7":7,"grid":
                 self._spin_nx.value() * self._spin_ny.value()}
        n = n_map.get(geom, "?")
        self._lbl_info.setText(
            f"Tube diameter: {self._diameter:.3f} Å  |  "
            f"Centre-to-centre pitch: {pitch:.3f} Å  |  {n} tubes total"
        )

    @property
    def geometry(self) -> str:   return self._cmb_geom.currentData()
    @property
    def spacing(self)  -> float: return self._spin_spacing.value()
    @property
    def vacuum(self)   -> float: return self._spin_vacuum.value()
    @property
    def nx(self)       -> int:   return self._spin_nx.value()
    @property
    def ny(self)       -> int:   return self._spin_ny.value()


# ─────────────────────────────────────────────────────────────────────────────
# Deformations dialog
# ─────────────────────────────────────────────────────────────────────────────

class DeformDialog(QDialog):
    """Configure and preview axial strain + torsion."""

    def __init__(self, default_vacuum: float = 10.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Apply Deformations")
        self.setFixedSize(420, 350)
        self.setStyleSheet(_DIALOG_STYLE)
        self._default_vacuum = float(default_vacuum)
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setSpacing(10)
        lay.setContentsMargins(20, 16, 20, 16)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self._spin_strain = QDoubleSpinBox()
        self._spin_strain.setRange(-50.0, 200.0)
        self._spin_strain.setValue(0.0)
        self._spin_strain.setSuffix(" %")
        self._spin_strain.setSingleStep(1.0)
        self._spin_strain.setDecimals(2)
        self._spin_strain.setToolTip(
            "Axial strain applied along the tube axis.\n"
            "Positive = tensile (stretching); negative = compressive."
        )
        form.addRow("Axial strain:", self._spin_strain)

        self._spin_twist = QDoubleSpinBox()
        self._spin_twist.setRange(-90.0, 90.0)
        self._spin_twist.setValue(0.0)
        self._spin_twist.setSuffix(" °/Å")
        self._spin_twist.setSingleStep(0.1)
        self._spin_twist.setDecimals(3)
        self._spin_twist.setToolTip(
            "Torsion applied as a rotation rate around the tube axis.\n"
            "Positive = right-hand (conventional) twist.\n"
            "Typical values: 0.1 – 2.0 °/Å for moderate twist."
        )
        form.addRow("Torsion rate:", self._spin_twist)

        self._spin_zvac = QDoubleSpinBox()
        self._spin_zvac.setRange(0.0, 100.0)
        self._spin_zvac.setValue(self._default_vacuum)
        self._spin_zvac.setSuffix(" Å")
        self._spin_zvac.setSingleStep(1.0)
        self._spin_zvac.setDecimals(2)
        self._spin_zvac.setEnabled(False)   # enabled when twist != 0
        self._spin_zvac.setToolTip(
            "Vacuum padding added to each end of the simulation box along\n"
            "the tube axis (Z).  Required when torsion is applied because\n"
            "the twist breaks axial periodicity."
        )
        form.addRow("Z vacuum (torsion):", self._spin_zvac)

        lay.addLayout(form)

        # Keep the Z-vacuum field disabled unless the user actually
        # applies a torsion — purely cosmetic, but avoids confusion.
        self._spin_twist.valueChanged.connect(
            lambda v: self._spin_zvac.setEnabled(abs(v) > 1e-9)
        )

        self._note = QLabel()
        self._note.setTextFormat(Qt.TextFormat.RichText)
        self._note.setWordWrap(True)
        self._note.setStyleSheet("color:#888; font-size:10px;")
        self._update_note()
        lay.addWidget(self._note)
        self._spin_twist.valueChanged.connect(self._update_note)

        lay.addStretch()
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    def _update_note(self):
        if abs(self._spin_twist.value()) > 1e-9:
            self._note.setText(
                "<b>ℹ Torsion breaks axial periodicity.</b>  A vacuum slab "
                "will be added along Z so the twisted segment is treated "
                "as a finite cluster in periodic codes."
            )
        else:
            self._note.setText(
                "ℹ  Deformations are applied to the current nanotube. "
                "Use Reset to return to the undeformed structure before "
                "applying a different deformation."
            )

    @property
    def axial_strain(self) -> float: return self._spin_strain.value() / 100.0
    @property
    def twist_rate(self)   -> float: return self._spin_twist.value()
    @property
    def z_vacuum(self)     -> float: return self._spin_zvac.value()


# ─────────────────────────────────────────────────────────────────────────────
# Analysis dialog (bond histogram + symmetry info)
# ─────────────────────────────────────────────────────────────────────────────

class AnalysisDialog(QDialog):
    """Bond length histogram + electronic character + symmetry."""

    def __init__(self, nt: NanotubeStructure,
                 lattice_type:        str = "hexagonal",
                 n_atoms_unit_cell:   int | None = None,
                 lattice_constant_a:  float | None = None,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Structure Analysis")
        self.resize(600, 480)
        self.setStyleSheet(_DIALOG_STYLE)
        self._nt = nt
        self._lattice_type       = lattice_type
        self._n_atoms_unit_cell  = n_atoms_unit_cell
        self._lattice_constant_a = lattice_constant_a
        self._build_ui()

    def _build_ui(self):
        from core.analysis import bond_analysis, electronic_character_label, tube_symmetry_info
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

        lay = QVBoxLayout(self)
        tabs = QTabWidget()

        # ── Tab 1: Bond histogram ─────────────────────────────────────────────
        tab_bonds = QWidget()
        tb_lay    = QVBoxLayout(tab_bonds)

        analysis = bond_analysis(self._nt)

        if analysis["n_bonds"] == 0:
            tb_lay.addWidget(QLabel("No bonds found within cutoff."))
        else:
            fig = Figure(figsize=(5, 3.5), facecolor="white")
            ax  = fig.add_subplot(111)
            canvas = FigureCanvasQTAgg(fig)

            dists = analysis["distances"]
            pairs = analysis["pairs"]
            species_set = sorted(analysis["species_set"])

            # Plot per-pair histogram
            import numpy as np
            colors = ["#1A3A6B","#E07000","#1C7A3E","#9B2335","#7B3FB0","#00838F"]
            for k, sp in enumerate(species_set):
                d_k = [d for d, p in zip(dists, pairs) if p == sp]
                if d_k:
                    ax.hist(d_k, bins=40, alpha=0.75,
                            label=sp, color=colors[k % len(colors)])

            ax.set_xlabel("Bond length (Å)")
            ax.set_ylabel("Count")
            ax.set_title(f"Bond length distribution  ({analysis['n_bonds']} bonds)")
            ax.legend(fontsize=8)
            fig.tight_layout()
            tb_lay.addWidget(canvas)

            # Stats table
            stats = QLabel(
                f"<b>Statistics:</b>  "
                f"mean = {analysis['mean']:.4f} Å  |  "
                f"σ = {analysis['std']:.4f} Å  |  "
                f"min = {analysis['min']:.4f} Å  |  "
                f"max = {analysis['max']:.4f} Å"
            )
            stats.setTextFormat(Qt.TextFormat.RichText)
            tb_lay.addWidget(stats)

        tabs.addTab(tab_bonds, "Bond Lengths")

        # ── Tab 2: Symmetry + electronic character ────────────────────────────
        tab_sym = QWidget()
        ts_lay  = QVBoxLayout(tab_sym)

        ch = self._nt.chirality
        n, m = ch.n, ch.m

        sym_info   = tube_symmetry_info(n, m)
        # Pass the parent 2D lattice markers so the zone-folding rule is
        # only applied when the precursor is *pristine graphene*
        # (hexagonal C with a 2-atom basis and a ≈ 2.46 Å).  The dialog's
        # caller supplies these via the constructor; if absent, the
        # function falls back to 'requires DFT' to prevent silent
        # misclassification of e.g. irida-graphene or h-BN.
        species_uc = None
        nt_syms = getattr(self._nt, "symbols", None)
        if nt_syms:
            species_uc = sorted(set(nt_syms))
        elec_label = electronic_character_label(
            n, m, self._lattice_type,
            species             = species_uc,
            n_atoms_unit_cell   = self._n_atoms_unit_cell,
            lattice_constant_a  = self._lattice_constant_a,
        )

        sym_text = QTextEdit()
        sym_text.setReadOnly(True)
        sym_text.setFontFamily("monospace")
        sym_text.setPlainText(
            f"Chiral indices : ({n}, {m})\n"
            f"Type           : {sym_info['type'].capitalize()}\n"
            f"Description    : {sym_info['description']}\n"
            f"gcd(n,m)       : {sym_info['d_nm']}\n"
            f"\n"
            f"Electronic character (zone-folding)\n"
            f"  {elec_label}\n"
            f"\n"
            f"Nanotube geometry\n"
            f"  Diameter     : {ch.diameter:.4f} Å\n"
            f"  |T|          : {ch.T_norm:.4f} Å\n"
            f"  θ (chiral)   : {ch.theta_deg:.4f}°\n"
            f"  Atoms/cell   : {self._nt.n_atoms}\n"
            f"  Strain       : {ch.strain:.6f} %\n"
            f"  Box (Å)      : {self._nt.box[0]:.3f} × "
            f"{self._nt.box[1]:.3f} × {self._nt.box[2]:.3f}\n"
        )
        ts_lay.addWidget(sym_text)
        ts_lay.addWidget(_copy_btn(sym_text))

        tabs.addTab(tab_sym, "Symmetry & Electronics")
        lay.addWidget(tabs)

        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        lay.addWidget(close)


# ─────────────────────────────────────────────────────────────────────────────
# Methods text dialog
# ─────────────────────────────────────────────────────────────────────────────

class MethodsDialog(QDialog):
    """Auto-generated methods-section paragraph."""

    def __init__(self, nt: NanotubeStructure,
                 structure=None, deform_desc: str = "",
                 n_walls: int = 1, wall_info: str = "",
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Methods Section Generator")
        self.setStyleSheet(_DIALOG_STYLE)
        self.resize(620, 380)
        self._nt          = nt
        self._structure   = structure
        self._deform_desc = deform_desc
        self._n_walls     = n_walls
        self._wall_info   = wall_info
        self._build_ui()

    def _build_ui(self):
        from core.analysis import generate_methods_text

        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        lay.setContentsMargins(16, 14, 16, 14)

        # Options row
        opt_row = QHBoxLayout()
        opt_row.addWidget(QLabel("Software name:"))
        self._edit_sw = QLineEdit("NTBuilder")
        self._edit_sw.setMaximumWidth(120)
        opt_row.addWidget(self._edit_sw)
        opt_row.addSpacing(12)
        opt_row.addWidget(QLabel("Version:"))
        self._edit_ver = QLineEdit("1.1")
        self._edit_ver.setMaximumWidth(60)
        opt_row.addWidget(self._edit_ver)
        opt_row.addSpacing(12)
        opt_row.addWidget(QLabel("Cite key:"))
        self._edit_cite = QLineEdit("[CITE]")
        self._edit_cite.setMaximumWidth(100)
        opt_row.addWidget(self._edit_cite)
        btn_regen = QPushButton("↻ Regenerate")
        btn_regen.clicked.connect(self._regen)
        opt_row.addWidget(btn_regen)
        opt_row.addStretch()
        lay.addLayout(opt_row)

        self._text = QTextEdit()
        self._text.setReadOnly(False)  # allow user edits
        self._text.setFont(QFont("Segoe UI", 10))
        lay.addWidget(self._text)

        btn_row = QHBoxLayout()
        btn_row.addWidget(_copy_btn(self._text))
        btn_row.addWidget(_save_btn(self._text, "methods.txt", self))
        btn_row.addStretch()
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        btn_row.addWidget(close)
        lay.addLayout(btn_row)

        self._regen()

    def _regen(self):
        from core.analysis import generate_methods_text
        txt = generate_methods_text(
            self._nt,
            structure=self._structure,
            deform_desc=self._deform_desc,
            software=self._edit_sw.text() or "NTBuilder",
            version=self._edit_ver.text() or "1.1",
            cite_key=self._edit_cite.text() or "[CITE]",
            n_walls=self._n_walls,
            wall_info=self._wall_info,
        )
        self._text.setPlainText(txt)


# ─────────────────────────────────────────────────────────────────────────────
# DFT input dialog
# ─────────────────────────────────────────────────────────────────────────────

class DFTInputDialog(QDialog):
    """Generate complete DFT input files (VASP, QE, CP2K)."""

    def __init__(self, nt: NanotubeStructure, structure=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("DFT Input Generator")
        self.resize(700, 500)
        self.setStyleSheet(_DIALOG_STYLE)
        self._nt        = nt
        self._structure = structure
        self._build_ui()

    def _build_ui(self):
        from core.analysis import generate_vasp_inputs, generate_qe_input, generate_cp2k_input

        lay  = QVBoxLayout(self)
        tabs = QTabWidget()

        def _make_tab(title: str, content: str, default_name: str):
            w   = QWidget()
            wl  = QVBoxLayout(w)
            txt = QTextEdit()
            txt.setFont(QFont("Courier New", 9))
            txt.setPlainText(content)
            wl.addWidget(txt)
            br = QHBoxLayout()
            br.addWidget(_copy_btn(txt))
            br.addWidget(_save_btn(txt, default_name, self))
            br.addStretch()
            wl.addLayout(br)
            tabs.addTab(w, title)

        # VASP
        incar, kpoints = generate_vasp_inputs(self._nt)
        _make_tab("VASP — INCAR", incar, "INCAR")
        _make_tab("VASP — KPOINTS", kpoints, "KPOINTS")

        # QE
        qe = generate_qe_input(self._nt, self._structure)
        _make_tab("Quantum ESPRESSO", qe, f"nt_{self._nt.chirality.n}_{self._nt.chirality.m}.in")

        # CP2K
        cp2k = generate_cp2k_input(self._nt)
        _make_tab("CP2K", cp2k, f"nt_{self._nt.chirality.n}_{self._nt.chirality.m}.inp")

        lay.addWidget(tabs)

        note = QLabel(
            "<small style='color:#888;'>ℹ  These inputs provide a sensible starting point. "
            "Always review pseudopotential choices, k-point density, and energy cutoffs "
            "before production calculations.</small>"
        )
        note.setTextFormat(Qt.TextFormat.RichText)
        note.setWordWrap(True)
        lay.addWidget(note)

        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        lay.addWidget(close)


# ─────────────────────────────────────────────────────────────────────────────
# C2DB download thread
# ─────────────────────────────────────────────────────────────────────────────

class _C2DBDownloadThread(QThread):
    """Background thread that downloads the C2DB SQLite database."""

    progress = Signal(int, int)   # bytes_done, bytes_total
    finished = Signal(str)        # "" on success, error message on failure

    def __init__(self):
        super().__init__()
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        from core.analysis import C2DBCache
        import urllib.request, shutil, tempfile
        from pathlib import Path

        target = C2DBCache._db_path()

        try:
            url = C2DBCache._find_download_url(timeout=10)
            hdr = {
                "User-Agent": "NTBuilder/1.1 (research use)",
                "Accept":     "*/*",
            }
            req = urllib.request.Request(url, headers=hdr)
            with urllib.request.urlopen(req, timeout=60) as resp:
                total = int(resp.headers.get("Content-Length", 0) or 0)
                done  = 0
                chunk = 1 << 17   # 128 KB

                tmp = tempfile.NamedTemporaryFile(
                    dir=target.parent, prefix=".c2db_dl_",
                    suffix=".db", delete=False,
                )
                try:
                    while not self._cancelled:
                        buf = resp.read(chunk)
                        if not buf:
                            break
                        tmp.write(buf)
                        done += len(buf)
                        self.progress.emit(done, total)
                    tmp.close()
                    if self._cancelled:
                        Path(tmp.name).unlink(missing_ok=True)
                        self.finished.emit("Download cancelled.")
                        return
                    shutil.move(tmp.name, str(target))
                except Exception:
                    tmp.close()
                    Path(tmp.name).unlink(missing_ok=True)
                    raise

            self.finished.emit("")
        except Exception as exc:
            self.finished.emit(
                f"{exc}\n\nYou can download the file manually from:\n"
                f"{C2DBCache._DOWNLOAD_URLS[0]}\n"
                f"and save it to:\n{C2DBCache._db_path()}"
            )


# ─────────────────────────────────────────────────────────────────────────────
# Database query dialog
# ─────────────────────────────────────────────────────────────────────────────

class DatabaseQueryDialog(QDialog):
    """Query C2DB (2DHHub) for 2D structures."""

    # ── C2DB column layout ────────────────────────────────────────────────────
    _C2DB_HEADERS = ["UID", "Formula", "Layer group", "E-hull (eV/at)", "Gap PBE (eV)", "Magnetic", "Source"]
    _C2DB_COL_UID        = 0
    _C2DB_COL_FORMULA    = 1
    _C2DB_COL_LAYERGROUP = 2
    _C2DB_COL_EHULL      = 3
    _C2DB_COL_GAP        = 4
    _C2DB_COL_MAG        = 5
    _C2DB_COL_SOURCE     = 6

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Database Structure Search")
        self.resize(780, 520)
        self.setStyleSheet(_DIALOG_STYLE)
        self._results: list[dict] = []
        self._downloaded_path: str | None = None
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        lay.setContentsMargins(16, 12, 16, 12)

        # ── Search row ────────────────────────────────────────────────────────
        src_row = QHBoxLayout()
        src_row.addWidget(QLabel("Formula:"))
        self._edit_formula = QLineEdit()
        self._edit_formula.setPlaceholderText("e.g.  C,  BN,  MoS2,  C2")
        self._edit_formula.returnPressed.connect(self._do_search)
        src_row.addWidget(self._edit_formula)

        src_row.addSpacing(8)
        _src_lbl = QLabel("Source: <b>C2DB — 2DHHub</b>")
        _src_lbl.setTextFormat(Qt.TextFormat.RichText)
        _src_lbl.setStyleSheet("color:#2851A3;")
        src_row.addWidget(_src_lbl)

        self._btn_search = QPushButton("Search")
        self._btn_search.clicked.connect(self._do_search)
        src_row.addWidget(self._btn_search)
        lay.addLayout(src_row)

        # ── C2DB stability filter row ─────────────────────────────────────────
        c2db_row = QHBoxLayout()
        c2db_row.addWidget(QLabel("Max E-hull (eV/atom):"))
        self._spin_ehull = QDoubleSpinBox()
        self._spin_ehull.setRange(0.0, 2.0)
        self._spin_ehull.setValue(0.5)
        self._spin_ehull.setSingleStep(0.1)
        self._spin_ehull.setDecimals(2)
        self._spin_ehull.setMaximumWidth(80)
        self._spin_ehull.setToolTip(
            "Only return materials with energy above the convex hull\n"
            "below this value.  0 = exactly on hull (most stable);\n"
            "0.2–0.5 eV/atom includes metastable but synthesisable materials."
        )
        c2db_row.addWidget(self._spin_ehull)
        _c2db_hint = QLabel("  0 = on hull · 0.2 = moderately stable · 0.5 = metastable")
        _c2db_hint.setStyleSheet("color:#666; font-size:10px;")
        c2db_row.addWidget(_c2db_hint)
        c2db_row.addStretch()
        self._c2db_widget = QWidget()
        self._c2db_widget.setLayout(c2db_row)
        lay.addWidget(self._c2db_widget)

        # ── C2DB local database status / download row ─────────────────────────
        self._db_row_widget = self._build_db_row()
        lay.addWidget(self._db_row_widget)

        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setVisible(False)
        lay.addWidget(self._progress)

        # ── Results table ─────────────────────────────────────────────────────
        _TABLE_STYLE = """
            QTableWidget {
                background-color: #FFFFFF;
                color: #111827;
                gridline-color: #D0D8EC;
                border: 1px solid #C4CDE0;
                border-radius: 5px;
                font-size: 11px;
            }
            QTableWidget::item:selected {
                background-color: #2851A3;
                color: #FFFFFF;
            }
            QTableWidget QHeaderView::section {
                background-color: #E2E9F5;
                color: #1A3A6B;
                font-weight: 600;
                border: none;
                border-right: 1px solid #C4CDE0;
                padding: 4px 6px;
            }
            QTableWidget::item:alternate {
                background-color: #F0F4FA;
            }
        """
        self._table = QTableWidget(0, len(self._C2DB_HEADERS))
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setAlternatingRowColors(True)
        self._table.setSortingEnabled(True)
        self._table.verticalHeader().setVisible(False)
        self._table.setStyleSheet(_TABLE_STYLE)
        lay.addWidget(self._table)

        # ── Legend (dynamically updated by _setup_table_columns) ─────────────
        self._legend_widget = QWidget()
        self._legend_lay    = QHBoxLayout(self._legend_widget)
        self._legend_lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self._legend_widget)

        # ── Status label ──────────────────────────────────────────────────────
        self._lbl_status = QLabel("Enter a chemical formula and click Search.")
        self._lbl_status.setStyleSheet("color:#666; font-size:10px;")
        lay.addWidget(self._lbl_status)

        # Initialise columns for the default source (C2DB)
        self._setup_table_columns()

        # ── Bottom buttons ────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._btn_open = QPushButton("Open selected in NTBuilder")
        self._btn_open.setEnabled(False)
        self._btn_open.clicked.connect(self._do_open)
        btn_row.addWidget(self._btn_open)
        btn_row.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        btn_row.addWidget(close_btn)
        lay.addLayout(btn_row)

        self._table.itemSelectionChanged.connect(
            lambda: self._btn_open.setEnabled(bool(self._table.selectedItems()))
        )

    # ─────────────────────────────────────────────────────────────────────────
    # C2DB local database: status bar + download
    # ─────────────────────────────────────────────────────────────────────────

    def _build_db_row(self) -> QWidget:
        from core.analysis import C2DBCache

        w   = QWidget()
        lay = QHBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)

        _icon = QLabel("🗄")
        _icon.setStyleSheet("font-size:13px;")
        lay.addWidget(_icon)

        self._lbl_db_status = QLabel()
        self._lbl_db_status.setStyleSheet("font-size:10px; color:#444;")
        lay.addWidget(self._lbl_db_status)

        self._btn_dl = QPushButton("⬇  Download full C2DB")
        self._btn_dl.setStyleSheet(
            "font-size:10px; padding:2px 8px; "
            "background:#2851A3; color:#FFF; border-radius:4px;"
        )
        self._btn_dl.clicked.connect(self._do_download)
        lay.addWidget(self._btn_dl)

        self._btn_dl_cancel = QPushButton("Cancel")
        self._btn_dl_cancel.setStyleSheet("font-size:10px; padding:2px 8px;")
        self._btn_dl_cancel.setVisible(False)
        self._btn_dl_cancel.clicked.connect(self._cancel_download)
        lay.addWidget(self._btn_dl_cancel)

        lay.addStretch()
        self._refresh_db_status()
        return w

    def _refresh_db_status(self):
        from core.analysis import C2DBCache
        if C2DBCache.is_available():
            mb = C2DBCache.db_size_mb()
            self._lbl_db_status.setText(
                f"Local C2DB: <b style='color:#1A6B1A;'>ready</b> "
                f"({mb:.0f} MB — full database, ~16 000 materials)"
            )
            self._lbl_db_status.setTextFormat(Qt.TextFormat.RichText)
            self._btn_dl.setVisible(False)
        else:
            self._lbl_db_status.setText(
                "Local C2DB: <b style='color:#8B4513;'>not downloaded</b> "
                "— searching in 30 curated materials"
            )
            self._lbl_db_status.setTextFormat(Qt.TextFormat.RichText)
            self._btn_dl.setVisible(True)

    def _do_download(self):
        reply = QMessageBox.question(
            self, "Download C2DB database",
            "The full C2DB database will be downloaded once and cached at\n"
            "~/.ntbuilder/c2db.db\n\n"
            "Estimated size: 100 MB – 2 GB (depends on server version).\n"
            "This may take several minutes on a slow connection.\n\n"
            "Proceed?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        self._btn_dl.setEnabled(False)
        self._btn_dl_cancel.setVisible(True)
        self._progress.setRange(0, 100)
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._lbl_db_status.setText("Downloading C2DB database…")

        self._dl_thread = _C2DBDownloadThread()
        self._dl_thread.progress.connect(self._on_dl_progress)
        self._dl_thread.finished.connect(self._on_dl_finished)
        self._dl_thread.start()

    def _cancel_download(self):
        if hasattr(self, "_dl_thread") and self._dl_thread.isRunning():
            self._dl_thread.cancel()
        self._progress.setVisible(False)
        self._btn_dl_cancel.setVisible(False)
        self._btn_dl.setEnabled(True)
        self._refresh_db_status()

    def _on_dl_progress(self, done: int, total: int):
        if total > 0:
            pct = int(done * 100 / total)
            self._progress.setRange(0, 100)
            self._progress.setValue(pct)
            self._lbl_db_status.setText(
                f"Downloading… {done/1_048_576:.1f} / {total/1_048_576:.1f} MB"
            )
        else:
            # Unknown total — spin
            self._progress.setRange(0, 0)
            self._lbl_db_status.setText(
                f"Downloading… {done/1_048_576:.1f} MB received"
            )

    def _on_dl_finished(self, error: str):
        self._progress.setVisible(False)
        self._btn_dl_cancel.setVisible(False)
        self._btn_dl.setEnabled(True)
        if error:
            QMessageBox.warning(self, "Download failed", error)
        else:
            QMessageBox.information(
                self, "Download complete",
                "C2DB database downloaded successfully!\n"
                "You can now search all ~16 000 2D materials."
            )
        self._refresh_db_status()

    # ─────────────────────────────────────────────────────────────────────────
    # Table column management
    # ─────────────────────────────────────────────────────────────────────────

    def _setup_table_columns(self):
        """Set table headers and legend for C2DB."""
        # Clear legend
        while self._legend_lay.count():
            item = self._legend_lay.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        def _swatch(color: str) -> QLabel:
            s = QLabel("  ")
            s.setFixedWidth(18); s.setFixedHeight(14)
            s.setStyleSheet(
                f"background:{color}; border:1px solid #aaa; border-radius:2px;"
            )
            return s

        # Legend: E-hull colour bands
        self._legend_lay.addWidget(_swatch("#C8F0C8"))
        self._legend_lay.addWidget(QLabel("E-hull ≤ 0.05 eV/at (highly stable)"))
        self._legend_lay.addSpacing(10)
        self._legend_lay.addWidget(_swatch("#FFF3CD"))
        self._legend_lay.addWidget(QLabel("0.05–0.2 (metastable)"))
        self._legend_lay.addSpacing(10)
        self._legend_lay.addWidget(_swatch("#F8D7DA"))
        self._legend_lay.addWidget(QLabel("> 0.2 (unstable — use with care)"))
        note_lbl = QLabel("  ✦ All entries are confirmed 2D — no 3D-bulk risk")
        note_lbl.setStyleSheet("color:#2851A3; font-size:10px; font-style:italic;")
        self._legend_lay.addSpacing(12)
        self._legend_lay.addWidget(note_lbl)
        self._legend_lay.addStretch()

        # Set table columns
        self._table.setSortingEnabled(False)
        self._table.setColumnCount(len(self._C2DB_HEADERS))
        self._table.setHorizontalHeaderLabels(self._C2DB_HEADERS)
        hdr = self._table.horizontalHeader()
        for col in range(len(self._C2DB_HEADERS)):
            hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.Interactive)
        hdr.setSectionResizeMode(self._C2DB_COL_FORMULA, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(self._C2DB_COL_SOURCE,  QHeaderView.ResizeMode.Stretch)
        self._table.setSortingEnabled(True)

    # ─────────────────────────────────────────────────────────────────────────
    # Search
    # ─────────────────────────────────────────────────────────────────────────

    def _do_search(self):
        from core.analysis import query_c2db

        formula = self._edit_formula.text().strip()
        if not formula:
            return

        self._btn_search.setEnabled(False)
        self._progress.setVisible(True)
        self._table.setRowCount(0)
        self._results = []
        self._lbl_status.setText("Searching…")

        try:
            results = query_c2db(formula, stability_max=self._spin_ehull.value())

            self._results = results
            self._table.setSortingEnabled(False)
            self._table.setRowCount(len(results))

            def _item(text: str) -> QTableWidgetItem:
                it = QTableWidgetItem(str(text))
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                return it

            def _num_item(val) -> QTableWidgetItem:
                text = f"{val:.3f}" if val is not None else "—"
                it = QTableWidgetItem(text)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                it.setData(Qt.ItemDataRole.UserRole + 1,
                           float(val) if val is not None else 1e18)
                return it

            for row, r in enumerate(results):
                # ── C2DB row (or MP if implemented later) ─────────────────────
                uid_item = _item(r.get("uid", r.get("id", "?")))
                uid_item.setData(Qt.ItemDataRole.UserRole, r)
                self._table.setItem(row, self._C2DB_COL_UID,        uid_item)
                self._table.setItem(row, self._C2DB_COL_FORMULA,    _item(r.get("formula", "?")))
                self._table.setItem(row, self._C2DB_COL_LAYERGROUP, _item(r.get("layer_group", "—")))
                self._table.setItem(row, self._C2DB_COL_EHULL,      _num_item(r.get("ehull")))
                self._table.setItem(row, self._C2DB_COL_GAP,        _num_item(r.get("gap_pbe")))
                self._table.setItem(row, self._C2DB_COL_MAG,        _item(r.get("magnetic", "—")))
                self._table.setItem(row, self._C2DB_COL_SOURCE,     _item(r.get("source", "C2DB")))

                # Colour-code by E-hull (stability)
                ehull = r.get("ehull")
                if ehull is not None:
                    if ehull <= 0.05:
                        row_color = "#C8F0C8"   # green — highly stable
                    elif ehull <= 0.20:
                        row_color = "#FFF3CD"   # yellow — metastable
                    else:
                        row_color = "#F8D7DA"   # red — unstable
                    for col in range(self._table.columnCount()):
                        it = self._table.item(row, col)
                        if it:
                            it.setBackground(QColor(row_color))

            self._table.setSortingEnabled(True)
            self._table.resizeColumnsToContents()

            self._lbl_status.setText(
                f"{len(results)} result(s) found for '{formula}'."
                "  ✦ All entries are confirmed 2D monolayers (DFT-relaxed)."
            )

        except Exception as exc:
            QMessageBox.warning(self, "Search failed", str(exc))
            self._lbl_status.setText(f"Error: {exc}")
        finally:
            self._btn_search.setEnabled(True)
            self._progress.setVisible(False)

    # ─────────────────────────────────────────────────────────────────────────
    # Selection & download
    # ─────────────────────────────────────────────────────────────────────────

    def _selected_result(self) -> dict | None:
        """Return the full result dict for the currently selected row."""
        rows = self._table.selectionModel().selectedRows()
        if not rows:
            return None
        row = rows[0].row()
        # The result dict is stored on column 0 of every row
        id_item = self._table.item(row, 0)
        if id_item is None:
            return None
        return id_item.data(Qt.ItemDataRole.UserRole)

    def _do_open(self):
        """Download the selected CIF and close the dialog."""
        import urllib.request, tempfile
        r = self._selected_result()
        if r is None:
            return

        url = r.get("file_url") or r.get("mp_url", "")
        src = r.get("source", "")

        # COD: ends in .cif   |   C2DB: ends in /download/cif
        can_download = (
            url.endswith(".cif") or
            url.endswith("/download/cif") or
            url.endswith("/cif")
        )
        if not url or not can_download:
            QMessageBox.information(
                self, "Cannot download",
                f"Direct CIF download is not available for {src} entries.\n"
                f"Visit: {url}"
            )
            return

        try:
            uid = r.get("id") or r.get("uid") or "x"
            self._lbl_status.setText(f"Downloading from {src}…")
            tmp = tempfile.NamedTemporaryFile(
                suffix=".cif", delete=False,
                prefix=f"ntb_{uid}_"
            )
            req = urllib.request.Request(url, headers={
                "User-Agent": "NTBuilder/1.1 (research use)",
            })
            with urllib.request.urlopen(req, timeout=20) as resp:
                tmp.write(resp.read())
            tmp.close()
            self._downloaded_path = tmp.name
            self.accept()
        except Exception as exc:
            QMessageBox.warning(self, "Download failed", str(exc))
            self._lbl_status.setText(f"Download error: {exc}")

    def downloaded_path(self) -> str | None:
        return self._downloaded_path
