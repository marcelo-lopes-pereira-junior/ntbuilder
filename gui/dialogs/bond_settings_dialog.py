"""
gui/dialogs/bond_settings_dialog.py
------------------------------------
Dialog for editing per-species-pair bond cutoffs.

Layout
------
  Header: tolerance spinbox  |  min-dist spinbox  |  Reset button
  Table : one row per unique species pair
    columns: Pair | Default max (Å) | Custom max (Å) | Custom min (Å)

The "Custom max" and "Custom min" cells are double-spinboxes.
If left at 0.00 they mean "not set" (use the default formula).

Reference: Alvarez (2008) Dalton Trans. 2832–2838. DOI 10.1039/b801115j
Default tolerance: 1.20  (same as VESTA, Mercury, ASE)
"""

from __future__ import annotations

from PyQt6.QtCore    import Qt
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QDoubleSpinBox, QGroupBox,
    QHBoxLayout, QLabel, QPushButton,
    QSizePolicy, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from core.connectivity import BondSettings, COVALENT_RADII, get_radius


# ─────────────────────────────────────────────────────────────────────────────
# Helper: spinbox factory
# ─────────────────────────────────────────────────────────────────────────────

def _make_dspin(lo: float, hi: float, val: float, step: float = 0.01) -> QDoubleSpinBox:
    sp = QDoubleSpinBox()
    sp.setRange(lo, hi)
    sp.setDecimals(3)
    sp.setSingleStep(step)
    sp.setValue(val)
    sp.setSuffix(" Å")
    sp.setFixedWidth(90)
    return sp


# ─────────────────────────────────────────────────────────────────────────────
# Dialog
# ─────────────────────────────────────────────────────────────────────────────

_COL_PAIR    = 0
_COL_DEF_MAX = 1
_COL_CUS_MAX = 2
_COL_CUS_MIN = 3
_HEADERS = ["Pair", "Default max (Å)", "Custom max (Å)", "Custom min (Å)"]


class BondSettingsDialog(QDialog):
    """
    Modal dialog for inspecting and editing bond cutoffs.

    Parameters
    ----------
    settings : BondSettings  — edited in-place when the user clicks OK
    species  : list[str]     — unique element symbols present in the structure
    parent   : QWidget or None
    """

    def __init__(
        self,
        settings: BondSettings,
        species:  list[str],
        parent=None,
    ):
        super().__init__(parent)
        self._settings = settings
        self._species  = sorted(set(species))
        self._pairs    = settings.pairs_for(self._species)

        self.setWindowTitle("Bond Cutoffs")
        self.setMinimumWidth(520)
        self._build_ui()
        self._populate()

    # ─────────────────────────────────────────────────────────────────────────
    # UI
    # ─────────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── Global parameters ─────────────────────────────────────────────
        grp_global = QGroupBox("Global parameters")
        gl = QHBoxLayout(grp_global)

        gl.addWidget(QLabel("Tolerance factor:"))
        self.spin_tol = _make_dspin(0.5, 3.0, self._settings.tolerance, step=0.05)
        self.spin_tol.setToolTip(
            "Scale factor applied to (r_A + r_B) to get the default max bond length.\n"
            "Default 1.20 matches VESTA, Mercury, and ASE."
        )
        gl.addWidget(self.spin_tol)

        gl.addSpacing(16)
        gl.addWidget(QLabel("Global min (Å):"))
        self.spin_gmin = _make_dspin(0.0, 2.0, self._settings.min_dist, step=0.05)
        self.spin_gmin.setToolTip("Minimum bond length accepted for any pair (avoids self-bonds).")
        gl.addWidget(self.spin_gmin)

        gl.addStretch()
        btn_reset = QPushButton("Reset all to defaults")
        btn_reset.clicked.connect(self._reset)
        gl.addWidget(btn_reset)

        root.addWidget(grp_global)

        # ── Reference note ────────────────────────────────────────────────
        lbl_ref = QLabel(
            "Covalent radii: Alvarez (2008) "
            "<i>Dalton Trans.</i> 2832–2838. "
            "DOI: <tt>10.1039/b801115j</tt>"
        )
        lbl_ref.setTextFormat(Qt.TextFormat.RichText)
        lbl_ref.setStyleSheet("font-size: 10px; color: #555555;")
        root.addWidget(lbl_ref)

        # ── Per-pair table ────────────────────────────────────────────────
        self.table = QTableWidget(len(self._pairs), 4)
        self.table.setHorizontalHeaderLabels(_HEADERS)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        root.addWidget(self.table)

        # Tooltip under table
        lbl_hint = QLabel(
            "Custom max / min = 0.000 Å means 'use default' (not set)."
        )
        lbl_hint.setStyleSheet("font-size: 10px; color: #777777;")
        root.addWidget(lbl_hint)

        # ── Dialog buttons ────────────────────────────────────────────────
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self._apply_and_accept)
        btn_box.rejected.connect(self.reject)
        root.addWidget(btn_box)

        # Refresh defaults when tolerance changes
        self.spin_tol.valueChanged.connect(self._refresh_defaults)

    # ─────────────────────────────────────────────────────────────────────────
    # Populate
    # ─────────────────────────────────────────────────────────────────────────

    def _populate(self):
        """Fill the table rows from current settings."""
        self._row_widgets: list[tuple[QDoubleSpinBox, QDoubleSpinBox]] = []

        for row, (a, b) in enumerate(self._pairs):
            pair_key = frozenset([a, b])
            def_max  = self._settings.default_max(a, b)

            # Col 0 — pair label
            lbl = QTableWidgetItem(f"{a}–{b}")
            lbl.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.table.setItem(row, _COL_PAIR, lbl)

            # Col 1 — default max (read-only, refreshed on tolerance change)
            item_def = QTableWidgetItem(f"{def_max:.3f}")
            item_def.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            item_def.setForeground(Qt.GlobalColor.gray)
            self.table.setItem(row, _COL_DEF_MAX, item_def)

            # Col 2 — custom max spinbox
            cur_max = self._settings.custom_max.get(pair_key, 0.0)
            sp_max  = _make_dspin(0.0, 10.0, cur_max)
            sp_max.setToolTip(
                f"Custom maximum bond length for {a}–{b}.\n"
                f"Set to 0.000 to use the default ({def_max:.3f} Å)."
            )
            self.table.setCellWidget(row, _COL_CUS_MAX, sp_max)

            # Col 3 — custom min spinbox
            cur_min = self._settings.custom_min.get(pair_key, 0.0)
            sp_min  = _make_dspin(0.0, 5.0, cur_min)
            sp_min.setToolTip(
                f"Custom minimum bond length for {a}–{b}.\n"
                "Set to 0.000 to use the global minimum."
            )
            self.table.setCellWidget(row, _COL_CUS_MIN, sp_min)

            self._row_widgets.append((sp_max, sp_min))

        self.table.resizeColumnsToContents()

    # ─────────────────────────────────────────────────────────────────────────
    # Slots
    # ─────────────────────────────────────────────────────────────────────────

    def _refresh_defaults(self):
        """Update the 'Default max' column when the tolerance spinbox changes."""
        tol = self.spin_tol.value()
        for row, (a, b) in enumerate(self._pairs):
            def_max = (get_radius(a) + get_radius(b)) * tol
            item    = self.table.item(row, _COL_DEF_MAX)
            if item:
                item.setText(f"{def_max:.3f}")
            # Update tooltip on custom-max spinbox
            sp_max = self._row_widgets[row][0]
            sp_max.setToolTip(
                f"Custom maximum bond length for {a}–{b}.\n"
                f"Set to 0.000 to use the default ({def_max:.3f} Å)."
            )

    def _reset(self):
        """Reset global params and clear all custom overrides in the table."""
        self.spin_tol.setValue(1.20)
        self.spin_gmin.setValue(0.40)
        for sp_max, sp_min in self._row_widgets:
            sp_max.setValue(0.0)
            sp_min.setValue(0.0)

    def _apply_and_accept(self):
        """Write table values back into the BondSettings and close."""
        self._settings.tolerance = self.spin_tol.value()
        self._settings.min_dist  = self.spin_gmin.value()
        self._settings.custom_max.clear()
        self._settings.custom_min.clear()

        for (a, b), (sp_max, sp_min) in zip(self._pairs, self._row_widgets):
            key = frozenset([a, b])
            v_max = sp_max.value()
            v_min = sp_min.value()
            if v_max > 1e-6:
                self._settings.custom_max[key] = v_max
            if v_min > 1e-6:
                self._settings.custom_min[key] = v_min

        self.accept()
