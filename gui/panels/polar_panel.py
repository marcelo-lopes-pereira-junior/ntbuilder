"""
gui/panels/polar_panel.py
--------------------------
Interactive polar map panel.

  r     = nanotube diameter (Å)
  theta = chiral angle (degrees)

Clicking on a point emits indices_selected(n, m).
Points are colour-coded by log10(atoms/cell).
The currently selected (n, m) is highlighted with a ring.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PyQt6.QtCore    import pyqtSignal, Qt
from PyQt6.QtGui     import QPixmap
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QDoubleSpinBox, QCheckBox, QComboBox, QGroupBox,
    QSizePolicy, QStackedWidget, QPushButton, QFileDialog,
    QButtonGroup,
)

_ASSETS = Path(__file__).parent.parent.parent / "assets"

# Load the placeholder logo once at import time.
# All resize events scale this cached pixmap — no per-resize disk reads.
_PLACEHOLDER_SRC: QPixmap | None = None

def _get_placeholder_pixmap() -> "QPixmap | None":
    global _PLACEHOLDER_SRC
    if _PLACEHOLDER_SRC is None:
        _png = _ASSETS / "ntbuilder_complete.png"
        _tmp = QPixmap(str(_png)) if _png.exists() else QPixmap()
        _PLACEHOLDER_SRC = _tmp if not _tmp.isNull() else None
    return _PLACEHOLDER_SRC


class _LogoPlaceholder(QWidget):
    """
    Shows ntbuilder_complete.png centred and letterboxed (aspect ratio preserved).
    The image is loaded once and the cached QPixmap is smooth-scaled on every
    resize.
    """
    _RATIO = 350.0 / 200.0   # width / height (matches original artwork dimensions)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Inner label — no background, transparent
        self._inner = QLabel(self)
        self._inner.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._inner.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)

        src = _get_placeholder_pixmap()
        if src is None or src.isNull():
            # Text fallback if PNG is missing
            self._inner.setText("NTBuilder")
            self._inner.setStyleSheet("font-size: 24px; color: #003F8A;")
            self._src = None
        else:
            self._src = src
        self._inner.show()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        aw, ah = self.width(), self.height()
        # Fit to available width first
        tw = aw
        th = int(tw / self._RATIO)
        # If that's too tall, fit to height instead
        if th > ah:
            th = ah
            tw = int(th * self._RATIO)
        # Centre in the available space
        x = (aw - tw) // 2
        y = (ah - th) // 2
        self._inner.setGeometry(x, y, tw, th)
        # Scale the cached high-res pixmap to the new display size
        if self._src is not None:
            self._inner.setPixmap(
                self._src.scaled(
                    tw, th,
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )


import matplotlib
import matplotlib.ticker as mticker
matplotlib.use("QtAgg")
# STIX font closely resembles LaTeX's Computer Modern — applies to all
# text elements (labels, ticks, titles, colorbar) in this figure.
matplotlib.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size":        10,      # STIX strokes are thinner; bump base size up
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
})
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from core import LatticeStructure, scan_chirality
from core.chirality import unique_sector_deg
from core.builder   import build_nanotube, check_spurious_bonds



class PolarPanel(QWidget):
    """Matplotlib polar map embedded in a Qt widget."""

    indices_selected = pyqtSignal(int, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._structure:  LatticeStructure | None = None
        self._all_results = []
        self._highlight:  tuple[int, int] | None  = None
        # Current rolling direction (mirrors the "Roll inward" checkbox on
        # the input panel).  False by default so the initial polar map is
        # drawn with the conventional outward roll; ``set_roll_inward``
        # flips this and forces a re-draw + invalidates the spurious cache.
        self._roll_inward: bool = False

        self._build_ui()

    # ─────────────────────────────────────────────────────────────────────────
    # UI
    # ─────────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ── Map settings group ───────────────────────────────────────────────
        map_box = QGroupBox("Map settings")
        map_lay = QVBoxLayout(map_box)
        map_lay.setSpacing(4)
        map_lay.setContentsMargins(8, 4, 8, 6)

        # Row 1 — diameter range
        ctrl1 = QHBoxLayout()
        ctrl1.setSpacing(4)

        _tip_dmin = (
            "<b>Minimum nanotube diameter (Å)</b><br>"
            "Only nanotubes with D ≥ this value are shown on the map.<br>"
            "For buckled or multi-layer structures this is set automatically<br>"
            "to 2 × max z-offset so the inner wall does not collapse."
        )
        lbl_dmin = QLabel("Minimum diameter (Å):")
        lbl_dmin.setToolTip(_tip_dmin)
        ctrl1.addWidget(lbl_dmin)
        self.spin_dmin = QDoubleSpinBox()
        self.spin_dmin.setRange(0.0, 500.0)
        self.spin_dmin.setValue(5.0)
        self.spin_dmin.setSingleStep(0.5)
        self.spin_dmin.setToolTip(_tip_dmin)
        self.spin_dmin.valueChanged.connect(self._redraw)
        ctrl1.addWidget(self.spin_dmin)

        ctrl1.addSpacing(12)

        _tip_dmax = (
            "<b>Maximum nanotube diameter (Å)</b><br>"
            "Sets the visible range of the chiral map. All (n, m) pairs<br>"
            "with diameter in [Min, Max] are computed and plotted.<br>"
            "Increasing this value recomputes the full map."
        )
        lbl_dmax = QLabel("Maximum diameter (Å):")
        lbl_dmax.setToolTip(_tip_dmax)
        ctrl1.addWidget(lbl_dmax)
        self.spin_dmax = QDoubleSpinBox()
        self.spin_dmax.setRange(1.0, 500.0)
        self.spin_dmax.setValue(30.0)
        self.spin_dmax.setSingleStep(1.0)
        self.spin_dmax.setToolTip(_tip_dmax)
        self.spin_dmax.valueChanged.connect(self._recompute_and_redraw)
        ctrl1.addWidget(self.spin_dmax)

        ctrl1.addStretch()
        map_lay.addLayout(ctrl1)

        # Row 2 — strain filter + colour selector
        ctrl2 = QHBoxLayout()
        ctrl2.setSpacing(4)

        _tip_strain_filter = (
            "<b>Maximum strain allowed (%)</b><br>"
            "When checked, hides nanotubes whose periodicity strain<br>"
            "exceeds the given threshold.<br><br>"
            "<i>Strain</i> = how far the translational vector T deviates<br>"
            "from exact perpendicularity with the chiral vector Ch.<br>"
            "Strain = 0 % → perfect periodicity (always true for zigzag<br>"
            "and armchair). Chiral tubes on non-hexagonal lattices may<br>"
            "have strain > 0 — lower values give cleaner unit cells."
        )
        self.chk_strain_filter = QCheckBox("Max strain allowed (%):")
        self.chk_strain_filter.setToolTip(_tip_strain_filter)
        self.chk_strain_filter.stateChanged.connect(self._redraw)
        ctrl2.addWidget(self.chk_strain_filter)

        self.spin_strain_max = QDoubleSpinBox()
        self.spin_strain_max.setRange(0.0, 100.0)
        self.spin_strain_max.setValue(0.01)
        self.spin_strain_max.setSingleStep(0.001)
        self.spin_strain_max.setDecimals(4)
        self.spin_strain_max.setSuffix(" %")
        self.spin_strain_max.setEnabled(False)
        self.spin_strain_max.setToolTip(_tip_strain_filter)
        self.spin_strain_max.valueChanged.connect(self._redraw)
        self.chk_strain_filter.stateChanged.connect(
            lambda s: self.spin_strain_max.setEnabled(bool(s))
        )
        ctrl2.addWidget(self.spin_strain_max)

        ctrl2.addSpacing(16)

        _tip_colour = (
            "<b>Colour dots by:</b><br>"
            "<b>Atom count</b> — colour encodes log10(atoms per unit cell).<br>"
            "Useful to identify small-cell tubes (fast to compute or simulate).<br><br>"
            "<b>Intrinsic strain</b> — colour encodes periodicity strain (%).<br>"
            "Useful on non-hexagonal lattices to find low-strain chiralities."
        )
        lbl_colour = QLabel("Colour by:")
        lbl_colour.setToolTip(_tip_colour)
        ctrl2.addWidget(lbl_colour)

        # ── Segmented toggle switch (replaces the old QComboBox) ────────────
        # Two exclusive checkable buttons styled as a segmented control.
        _toggle_style = (
            "QPushButton {"
            "  background-color: #EEF1F8; color: #1A3A6B;"
            "  border: 1px solid #C4CDE0;"
            "  padding: 3px 10px; font-size: 10px;"
            "}"
            "QPushButton:checked {"
            "  background-color: #2851A3; color: white;"
            "  border: 1px solid #1A3A6B;"
            "}"
            "QPushButton:hover:!checked {"
            "  background-color: #DCE3F2;"
            "}"
        )
        self.btn_colour_atoms  = QPushButton("Atom count")
        self.btn_colour_strain = QPushButton("Intrinsic strain")
        for b in (self.btn_colour_atoms, self.btn_colour_strain):
            b.setCheckable(True)
            b.setFixedHeight(24)
            b.setStyleSheet(_toggle_style)
            b.setToolTip(_tip_colour)
        # Round only the outer corners so they look like a single segmented pill.
        self.btn_colour_atoms.setStyleSheet(
            _toggle_style.replace(
                "padding: 3px 10px;",
                "padding: 3px 10px; border-top-left-radius: 12px;"
                " border-bottom-left-radius: 12px; border-right: none;"
            )
        )
        self.btn_colour_strain.setStyleSheet(
            _toggle_style.replace(
                "padding: 3px 10px;",
                "padding: 3px 10px; border-top-right-radius: 12px;"
                " border-bottom-right-radius: 12px;"
            )
        )
        self.btn_colour_atoms.setChecked(True)

        self._colour_group = QButtonGroup(self)
        self._colour_group.setExclusive(True)
        self._colour_group.addButton(self.btn_colour_atoms,  0)
        self._colour_group.addButton(self.btn_colour_strain, 1)
        self._colour_group.idToggled.connect(
            lambda _id, checked: self._redraw() if checked else None
        )

        ctrl2.addWidget(self.btn_colour_atoms)
        ctrl2.addWidget(self.btn_colour_strain)

        # Backwards-compatibility shim: anything else in the codebase that
        # still queries ``self.cmb_colour.currentIndex()`` keeps working.
        class _CmbShim:
            def __init__(self, group):  self._g = group
            def currentIndex(self):     return self._g.checkedId()
        self.cmb_colour = _CmbShim(self._colour_group)

        ctrl2.addStretch()
        map_lay.addLayout(ctrl2)

        root.addWidget(map_box)

        # ── QStackedWidget: page 0 = logo placeholder, page 1 = polar map ──
        self._stack = QStackedWidget()
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )

        # Page 0 — logo placeholder (shown before any structure is loaded)
        self._stack.addWidget(_LogoPlaceholder())   # index 0

        # Page 1 — Matplotlib polar map
        self.fig = Figure(figsize=(5, 5), facecolor="white", constrained_layout=True)
        self.ax  = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self._stack.addWidget(self.canvas)  # index 1

        root.addWidget(self._stack)

        # Bottom row: info label + export button
        bottom_row = QHBoxLayout()

        self.lbl_info = QLabel("Load a structure to populate the map.")
        self.lbl_info.setStyleSheet("font-size: 10px; color: #666666;")
        bottom_row.addWidget(self.lbl_info, stretch=1)

        self.btn_export_map = QPushButton("💾  Export map…")
        self.btn_export_map.setFixedHeight(26)
        self.btn_export_map.setEnabled(False)
        self.btn_export_map.setToolTip(
            "Save the current chirality map as an SVG or PDF vector image."
        )
        self.btn_export_map.clicked.connect(self._on_export_map)
        bottom_row.addWidget(self.btn_export_map)

        root.addLayout(bottom_row)

    # ─────────────────────────────────────────────────────────────────────────
    # Public API
    # ─────────────────────────────────────────────────────────────────────────

    def set_structure(self, structure: LatticeStructure):
        self._structure = structure
        self._highlight = None

        # Switch from logo placeholder to polar map on first load
        self._stack.setCurrentIndex(1)
        self.btn_export_map.setEnabled(True)

        # The legacy auto-fill of D_min from ``structure.d_min`` (= 2·|z_max|)
        # was misleading: it is the trivial geometric lower bound that lets
        # the innermost shell sit at positive radius, but says nothing about
        # curvature-induced spurious bonds.  For Janus monolayers like MoSSe
        # the real threshold is far larger (e.g. ~29 Å for Se-Se).  Rather
        # than pre-filtering with a wrong number we surface the constraint
        # visually: chiralities that would develop spurious bonds are drawn
        # with an X marker on the polar map (see ``_redraw``).  The user
        # can still build any (n, m) they want — including X-marked ones —
        # for studies that deliberately exploit curvature-induced bonding.
        # ``spin_dmin`` keeps whatever value the user last chose; the widget
        # was constructed with the default of 5 Å, which gives a clean view
        # of the technologically interesting diameter window without
        # auto-filtering anything based on the structure.
        #
        # Spurious-bond status is memoised in ``self._spurious_cache``,
        # keyed by ``(n, m, roll_inward)`` because the rolling direction
        # selects which face of the buckled monolayer ends up on the
        # concave (compressed) side — flipping the toggle changes the
        # answer for asymmetric (Janus) structures like MoSSe.
        self._spurious_cache: dict[tuple[int, int, bool], list[str]] = {}

        self._recompute_and_redraw()

    def highlight_point(self, n: int, m: int):
        self._highlight = (n, m)
        self._redraw()

    def set_roll_inward(self, roll_inward: bool) -> None:
        """Update the rolling direction and refresh the spurious markers.

        Called by ``main_window`` when the user toggles the "Roll inward"
        checkbox on the input panel.  The cache is keyed by
        ``(n, m, roll_inward)`` so previously computed entries remain
        valid; only the *current* draw needs to use the new flag.
        """
        new_value = bool(roll_inward)
        if new_value == self._roll_inward:
            return
        self._roll_inward = new_value
        if self._all_results:
            self._redraw()

    # ─────────────────────────────────────────────────────────────────────────
    # Internal
    # ─────────────────────────────────────────────────────────────────────────

    def _recompute_and_redraw(self):
        if self._structure is None:
            return
        dmax = self.spin_dmax.value()
        # Derive n_max and m_max independently for each lattice axis.
        # Zigzag formula along a1: D = n·a/π  →  n_max = ⌈D·π/a⌉
        # Zigzag formula along a2: D = m·b/π  →  m_max = ⌈D·π/b⌉
        # This ensures all (n,m) pairs within the diameter window are found
        # regardless of whether a ≠ b (e.g. rectangular lattices).
        a = max(self._structure.a, 0.5)
        b = max(self._structure.b, 0.5)
        n_max = max(1, math.ceil(dmax * math.pi / a))
        m_max = max(1, math.ceil(dmax * math.pi / b))
        self._all_results = scan_chirality(
            self._structure,
            n_max=n_max,
            m_max=m_max,
            max_diameter=dmax + 2,  # small margin so points near edge are included
            search_limit=50,        # lighter search for map display (build uses 300)
        )
        self._redraw()

    def _redraw(self):
        if not self._all_results:
            return

        dmax = self.spin_dmax.value()
        dmin = self.spin_dmin.value()
        pts  = [r for r in self._all_results
                if dmin <= r.diameter <= dmax]

        # Apply optional strain filter
        if self.chk_strain_filter.isChecked():
            s_max = self.spin_strain_max.value()
            pts   = [r for r in pts if r.strain <= s_max]

        # ── Rebuild axes as a plain Cartesian axis ───────────────────────────
        # We do the polar→Cartesian projection manually:
        #   x = D·cos(θ),  y = D·sin(θ)
        # This avoids matplotlib's polar-axis boundary clipping, which
        # asymmetrically clips points at θ=thetamax (armchair) but not at
        # θ=thetamin (zigzag).  With a regular axis, both boundaries are just
        # interior regions of the plot and are treated identically.
        self.fig.clf()
        # Re-enable constrained_layout after clf() resets it
        self.fig.set_constrained_layout(True)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect("equal")
        self.ax.set_axis_off()
        self.ax.set_facecolor("white")
        self.fig.patch.set_facecolor("white")

        theta_max_deg = (
            unique_sector_deg(self._structure)
            if self._structure is not None else 90.0
        )
        theta_max_rad = math.radians(theta_max_deg)

        lt = self._structure.lattice_type if self._structure else "oblique"
        _is_sym = lt == "hexagonal" or (
            lt == "rectangular"
            and self._structure is not None
            and abs(self._structure.a - self._structure.b) < 1e-3
        )

        if not pts:
            self._draw_sector_frame(dmax, theta_max_deg, theta_max_rad, _is_sym)
            _r = dmax * 1.12
            self.ax.set_xlim(-_r * 0.08, _r)
            self.ax.set_ylim(-_r * 0.12, _r * math.sin(math.radians(
                unique_sector_deg(self._structure) if self._structure else 90.0
            )) + _r * 0.18)
            self.canvas.draw()
            return

        # Convert polar → Cartesian
        xs     = np.array([r.diameter * math.cos(math.radians(r.theta_deg)) for r in pts])
        ys     = np.array([r.diameter * math.sin(math.radians(r.theta_deg)) for r in pts])
        atoms  = np.array([r.n_atoms  for r in pts])
        strain = np.array([r.strain   for r in pts])

        # Curvature-induced spurious-bond check for buckled / Janus
        # structures.  Flat lattices (graphene, biphenylene, …) cannot
        # develop curvature-induced bonds, so we skip the per-chirality
        # construction entirely.  Results are memoised in
        # ``self._spurious_cache`` so that re-draws triggered by filter
        # / colour-mode changes pay the cost only once per (n, m).
        has_buckling = (
            self._structure is not None
            and getattr(self._structure, "has_buckling", False)
        )
        spurious_mask = np.zeros(len(pts), dtype=bool)
        if has_buckling:
            roll_inward = bool(self._roll_inward)
            for i, r in enumerate(pts):
                key = (r.n, r.m, roll_inward)
                if key not in self._spurious_cache:
                    try:
                        _nt = build_nanotube(
                            self._structure, r,
                            vacuum=0.0, roll_inward=roll_inward,
                        )
                        sp  = check_spurious_bonds(self._structure, _nt)
                        self._spurious_cache[key] = sorted(
                            "-".join(sorted(p)) for p in sp
                        )
                    except Exception:
                        self._spurious_cache[key] = []
                if self._spurious_cache[key]:
                    spurious_mask[i] = True

        # Colour metric
        colour_by_strain = self.cmb_colour.currentIndex() == 1
        if colour_by_strain:
            c_vals = np.log1p(strain)
            cmap   = "plasma"
            cb_lbl = "strain (%)"
        else:
            c_vals = np.log10(np.clip(atoms, 1, None))
            cmap   = "viridis"
            cb_lbl = "atoms / cell"

        # Two scatters with the same colour-mapping: dots for clean
        # chiralities, X markers for the curvature-spurious ones.  We
        # share the colour normalisation so the X markers are read on
        # the same colourbar.  The X marker is intentionally larger so
        # it is unambiguously distinguishable from a dot.
        clean = ~spurious_mask
        sc = self.ax.scatter(
            xs[clean], ys[clean],
            c=c_vals[clean], cmap=cmap, s=24,
            alpha=0.90, edgecolors="white", linewidths=0.4,
            zorder=3, picker=6,
        )
        if spurious_mask.any():
            # ``vmin``/``vmax`` from the dotted scatter so colours match.
            vmin, vmax = sc.get_clim()
            self.ax.scatter(
                xs[spurious_mask], ys[spurious_mask],
                c=c_vals[spurious_mask], cmap=cmap,
                vmin=vmin, vmax=vmax,
                marker="x", s=48, linewidths=1.4,
                zorder=4, picker=6,
            )

        # Marker legend — only shown when the spurious-bond check is
        # actually meaningful (i.e. the structure is buckled / Janus).
        # For flat lattices both symbols would carry the same meaning,
        # so the legend would be visual noise; we suppress it.
        if has_buckling:
            from matplotlib.lines import Line2D as _L2D
            legend_handles = [
                _L2D([0], [0], marker="o", linestyle="None",
                     markerfacecolor="0.5", markeredgecolor="white",
                     markersize=6, label="clean"),
                _L2D([0], [0], marker="x", linestyle="None",
                     color="0.25", markersize=8, markeredgewidth=1.6,
                     label="curvature-induced\nspurious bonds"),
            ]
            self.ax.legend(
                handles=legend_handles,
                loc="upper right", fontsize=8,
                framealpha=0.88, facecolor="white",
                edgecolor="0.7", borderpad=0.5,
                handletextpad=0.4, labelspacing=0.5,
            )

        # Highlight selected point
        if self._highlight:
            hn, hm = self._highlight
            for r in pts:
                if r.n == hn and r.m == hm:
                    hx = r.diameter * math.cos(math.radians(r.theta_deg))
                    hy = r.diameter * math.sin(math.radians(r.theta_deg))
                    self.ax.scatter(
                        hx, hy,
                        s=140, facecolors="none",
                        edgecolors="#FF0000", linewidths=2.2,
                        zorder=10,
                    )
                    break

        # Draw the sector frame (gridlines + boundary lines + labels)
        self._draw_sector_frame(dmax, theta_max_deg, theta_max_rad, _is_sym)

        # Index labels (n, m) along the two boundary lines
        self._draw_index_labels(pts, theta_max_deg, theta_max_rad, dmax)

        # Axis limits: leave a margin around the sector.
        # For oblique lattices with γ > 90° the a₂ boundary line extends into
        # negative-x territory; we must widen the left margin to show it.
        r_out  = dmax * 1.12
        x_left = -r_out * 0.08
        if theta_max_deg > 90.0:
            # Leftmost point of the sector is the a₂ boundary at r = dmax
            x_left = min(x_left, dmax * math.cos(theta_max_rad) - r_out * 0.08)
        self.ax.set_xlim(x_left, r_out)

        # IMPORTANT: for sectors wider than 90° the arc reaches its maximum y
        # at θ = 90° (where y = r), NOT at θ = θ_max where sin(θ_max) < 1.
        # Using sin(θ_max) for these cases clips the top of the arc.
        y_arc_max = r_out if theta_max_deg >= 90.0 else r_out * math.sin(theta_max_rad)
        self.ax.set_ylim(
            -r_out * 0.12,          # room for Zigzag label + diameter ticks
            y_arc_max + r_out * 0.18,
        )

        # Symmetry label (title)
        sym_labels = {
            "hexagonal": f"Hexagonal  (0° – {theta_max_deg:.0f}° unique)",
            "rectangular": (
                f"Square  (0° – {theta_max_deg:.0f}° unique)"
                if _is_sym else "Rectangular  (full quadrant)"
            ),
            "oblique": f"Oblique  (0° – {theta_max_deg:.1f}° sector, γ = a1∧a2)",
        }
        sym_str = sym_labels.get(lt, "")
        self.ax.set_title(
            f"Chiral Map — {sym_str}\n(click to select)",
            fontsize=11, pad=4,
        )

        # Colourbar — anchored to the axes so constrained_layout reserves space
        # for it automatically (no fixed figure-coordinate axes needed).
        # Ticks are shown in real units (atom counts or strain %) even though
        # the internal colour mapping uses log scale.
        try:
            cbar = self.fig.colorbar(
                sc, ax=self.ax, label=cb_lbl,
                orientation="horizontal", location="bottom",
                shrink=0.80, aspect=25, pad=0.04,
            )
            cbar.ax.tick_params(labelsize=9)
            cbar.ax.xaxis.label.set_size(9)
            if colour_by_strain:
                # strain: values are log1p(strain%) → show back-transformed
                cbar.ax.xaxis.set_major_formatter(
                    mticker.FuncFormatter(
                        lambda x, _: f"{np.expm1(x):.3g}"
                    )
                )
            else:
                # atoms: values are log10(n) → show as 10, 100, 1k, 10k …
                def _atoms_fmt(x, _):
                    v = 10 ** x
                    if v >= 1_000_000:
                        return f"{v/1_000_000:.0f}M"
                    if v >= 1_000:
                        return f"{v/1_000:.0f}k"
                    return f"{v:.0f}"
                cbar.ax.xaxis.set_major_formatter(
                    mticker.FuncFormatter(_atoms_fmt)
                )
        except Exception:
            pass

        self.canvas.draw()
        dmin_str = f"D ≥ {dmin:.1f} Å  " if dmin > 0 else ""
        n_sp = int(spurious_mask.sum())
        sp_str = (
            f"  ·  {n_sp} marked × (curvature-induced spurious bonds)"
            if n_sp > 0 else ""
        )
        self.lbl_info.setText(
            f"{len(pts)} nanotubes | {dmin_str}D ≤ {dmax:.1f} Å  ·  "
            f"click a point to select (n, m){sp_str}"
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Adaptive tick helper
    # ─────────────────────────────────────────────────────────────────────────

    @staticmethod
    def _tick_radii(dmax: float) -> list:
        """
        Return 4–6 'nice' diameter tick values (in Å) ≤ dmax.

        The step size is chosen from a fixed ladder so the number of ticks
        stays between 4 and 6 regardless of dmax:

          dmax ≤ 30  → step  5  (e.g. 5 10 15 20 25 30)
          dmax ~ 80  → step 15  (e.g. 15 30 45 60 75)
          dmax ~ 110 → step 20  (e.g. 20 40 60 80 100)
          dmax ~ 200 → step 40  (e.g. 40 80 120 160 200)
        """
        for step in (5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200):
            ticks = [step * i for i in range(1, 500) if step * i <= dmax]
            if len(ticks) <= 6:
                return ticks
        # Fallback: 5 evenly-spaced values
        return [round(dmax * i / 5) for i in range(1, 6)]

    # ─────────────────────────────────────────────────────────────────────────
    # Index-label helper
    # ─────────────────────────────────────────────────────────────────────────

    def _draw_index_labels(self, pts, theta_max_deg, theta_max_rad, dmax):
        """
        Draw (n, m) index labels at major diameter gridlines along the two
        sector boundary lines.

        • Zigzag boundary (θ ≈ 0°)  → label shows n value
        • Far boundary   (θ ≈ θmax) → label shows (n, m) pair
        """
        ax      = self.ax
        tol_ang = max(1.5, theta_max_deg * 0.04)    # angular tolerance (degrees)
        r_ticks = self._tick_radii(dmax)
        # Radial tolerance: accept a point within 40% of the tick spacing
        step    = ((r_ticks[-1] - r_ticks[0]) / max(len(r_ticks) - 1, 1)
                   if len(r_ticks) > 1 else 5.0)
        tol_r   = step * 0.45

        # ── Zigzag boundary ───────────────────────────────────────────────────
        # Skip the last (outermost) tick — it sits right under the "0° Zigzag"
        # annotation placed at the tip of the boundary line and would overlap it.
        zig_pts = [r for r in pts if r.theta_deg <= tol_ang]
        zig_ticks = r_ticks[:-1] if len(r_ticks) > 1 else r_ticks
        for rd in zig_ticks:
            if not zig_pts:
                break
            closest = min(zig_pts, key=lambda r: abs(r.diameter - rd))
            if abs(closest.diameter - rd) > tol_r:
                continue
            ax.text(
                closest.diameter, -dmax * 0.080,
                f"n = {closest.n}",
                fontsize=10, ha="center", va="top", color="#1A3A6B",
                fontweight="bold", clip_on=False,
            )

        # ── Far boundary (θ ≈ θ_max) ─────────────────────────────────────────
        far_pts = [r for r in pts if abs(r.theta_deg - theta_max_deg) <= tol_ang]
        # Perpendicular direction (rotated +90° from boundary direction) for offset
        perp_x = -math.sin(theta_max_rad)
        perp_y =  math.cos(theta_max_rad)
        # Use a slightly larger offset when the boundary is the vertical axis
        # (θ_max ≈ 90°) to keep the "m = …" labels clear of the axis line.
        is_vertical_boundary = abs(theta_max_rad - math.pi / 2) < 0.05
        offset = dmax * (0.085 if is_vertical_boundary else 0.055)
        # Skip the last (outermost) tick — it would overlap the "Armchair" label
        # placed at the arc end just beyond dmax.
        far_ticks = r_ticks[:-1] if len(r_ticks) > 1 else r_ticks

        for rd in far_ticks:
            if not far_pts:
                break
            closest = min(far_pts, key=lambda r: abs(r.diameter - rd))
            if abs(closest.diameter - rd) > tol_r:
                continue
            bx = closest.diameter * math.cos(theta_max_rad)
            by = closest.diameter * math.sin(theta_max_rad)
            lbl = f"m = {closest.m}"
            # ha: for θ ≥ 90° the boundary is vertical or tilts left, so we
            # right-align the label so its text body stays on the outside of
            # the sector (negative-x side) and never crosses the axis line.
            ha = "right" if theta_max_rad >= math.pi / 2 - 1e-3 else "left"
            ax.text(
                bx + perp_x * offset,
                by + perp_y * offset,
                lbl,
                fontsize=10, ha=ha, va="center", color="#1A3A6B",
                fontweight="bold", clip_on=False,
            )

    # ─────────────────────────────────────────────────────────────────────────
    # Sector frame helper
    # ─────────────────────────────────────────────────────────────────────────

    def _draw_sector_frame(
        self,
        dmax: float,
        theta_max_deg: float,
        theta_max_rad: float,
        is_sym: bool,
    ):
        """Draw gridlines, arcs and boundary lines for the Cartesian polar view."""
        ax = self.ax
        arc_t = np.linspace(0, theta_max_rad, 200)

        # ── Concentric arc gridlines ─────────────────────────────────────────
        r_ticks = self._tick_radii(dmax)
        for r in r_ticks:
            ax.plot(r * np.cos(arc_t), r * np.sin(arc_t),
                    color="#CCCCCC", linewidth=0.6, zorder=1, linestyle=":")
            # Diameter labels sit just below the zigzag axis.
            # Skip the label at exactly dmax to avoid collision with the
            # "0° Zigzag" annotation placed at the end of that boundary line.
            if r < dmax:
                ax.text(r, -dmax * 0.045, f"{r} Å",
                        fontsize=10, ha="center", va="top", color="#333333",
                        fontweight="bold")

        # ── Radial (angular) gridlines ───────────────────────────────────────
        n_ticks = 4
        step    = theta_max_deg / n_ticks
        for i in range(n_ticks + 1):
            td  = i * step
            tr  = math.radians(td)
            cx  = dmax * math.cos(tr)
            cy  = dmax * math.sin(tr)

            if i == 0 or i == n_ticks:
                # Both boundary lines identical — solid black
                ax.plot([0, cx * 1.05], [0, cy * 1.05],
                        color="black", linewidth=1.2, zorder=2)
            else:
                # Interior angular gridlines — light dotted
                ax.plot([0, cx], [0, cy],
                        color="#CCCCCC", linewidth=0.6, zorder=1, linestyle=":")

            # Angle labels — placed just beyond the outer arc
            lx = (dmax * 1.08) * math.cos(tr)
            ly = (dmax * 1.08) * math.sin(tr)
            if is_sym and i == 0:
                # "Zigzag" label: right-aligned at x=dmax, just below the axis.
                # ha="right" means text extends LEFTWARD from that point — never cut.
                ax.text(dmax, -dmax * 0.07,
                        "0°  Zigzag", fontsize=11, ha="right", va="top",
                        fontweight="bold", clip_on=False)
            elif is_sym and i == n_ticks:
                # "Armchair" label: centred on the arc end point so it stays inside
                # the plot even at small figure widths (ha="left" was being cut).
                ax.text(lx, ly,
                        f"{td:.0f}°  Armchair", fontsize=11, ha="center", va="bottom",
                        fontweight="bold")
            else:
                lbl = f"{td:.1f}°"
                ha  = "left" if tr <= math.pi / 2 else "right"
                va  = "center"
                ax.text(lx, ly, lbl, fontsize=10, ha=ha, va=va)

        # ── Outer arc (boundary) ─────────────────────────────────────────────
        ax.plot(dmax * np.cos(arc_t), dmax * np.sin(arc_t),
                color="black", linewidth=1.2, zorder=2)

    # ─────────────────────────────────────────────────────────────────────────
    # Export map
    # ─────────────────────────────────────────────────────────────────────────

    def _on_export_map(self):
        """Save the current chirality map as SVG or PDF."""
        path, sel = QFileDialog.getSaveFileName(
            self, "Export chirality map",
            "chirality_map",
            "SVG vector image (*.svg);;PDF document (*.pdf);;PNG image (*.png)",
        )
        if not path:
            return
        # Auto-append extension if missing
        ext = Path(path).suffix.lower()
        if ext not in (".svg", ".pdf", ".png"):
            if "SVG" in sel:
                path += ".svg"
            elif "PDF" in sel:
                path += ".pdf"
            else:
                path += ".png"

        try:
            self.fig.savefig(path, dpi=200, bbox_inches="tight")
            self.lbl_info.setText(f"Map saved → {Path(path).name}")
        except Exception as exc:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Export error", str(exc))

    # ─────────────────────────────────────────────────────────────────────────
    # Click handler
    # ─────────────────────────────────────────────────────────────────────────

    def _on_click(self, event):
        if event.inaxes != self.ax or not self._all_results:
            return
        if event.xdata is None or event.ydata is None:
            return

        # Convert Cartesian click → polar (theta, r)
        click_r     = math.sqrt(event.xdata ** 2 + event.ydata ** 2)
        click_theta = math.atan2(event.ydata, event.xdata)   # radians

        dmax = self.spin_dmax.value()
        pts  = [r for r in self._all_results if r.diameter <= dmax]

        if not pts:
            return

        # Nearest point in (normalised_theta, normalised_r) space
        best, best_dist = None, float("inf")
        for r in pts:
            dth  = math.radians(r.theta_deg) - click_theta
            dr   = (r.diameter - click_r) / max(dmax, 1)
            dist = math.sqrt(dth ** 2 + dr ** 2)
            if dist < best_dist:
                best_dist, best = dist, r

        if best is not None and best_dist < 0.15:
            self.indices_selected.emit(best.n, best.m)
            self.highlight_point(best.n, best.m)
            self.lbl_info.setText(
                f"Selected ({best.n},{best.m})  "
                f"D={best.diameter:.3f} Å  "
                f"θ={best.theta_deg:.2f}°  "
                f"atoms={best.n_atoms}  "
                f"strain={best.strain:.4f}%"
            )
