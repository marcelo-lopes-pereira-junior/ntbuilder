"""
gui/utils.py
------------
Shared GUI helpers: PNG pixmap loading and a ScalablePixmapLabel widget.
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui  import QPixmap
from PyQt6.QtWidgets import QLabel, QSizePolicy


def load_pixmap(path: Path, width: int, height: int) -> "QPixmap | None":
    """
    Load a PNG (or any Qt-supported raster format) from *path* and scale it
    to fit within *width* × *height* pixels, preserving aspect ratio.

    Returns ``None`` when the file is missing or cannot be loaded.
    """
    if not path.exists():
        return None
    pix = QPixmap(str(path))
    if pix.isNull():
        return None
    return pix.scaled(
        width, height,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )


# ─────────────────────────────────────────────────────────────────────────────
# ScalablePixmapLabel
# ─────────────────────────────────────────────────────────────────────────────

class ScalablePixmapLabel(QLabel):
    """
    A QLabel that holds a high-resolution source pixmap and smooth-scales it
    to fill the widget's current size on every resize, preserving aspect ratio.

    Usage::

        lbl = ScalablePixmapLabel()
        lbl.set_source_pixmap(some_high_res_pixmap)
        # The label will auto-scale the pixmap whenever it is resized.
    """

    def __init__(self, src_pixmap: "QPixmap | None" = None, parent=None):
        super().__init__(parent)
        self._src: QPixmap | None = src_pixmap
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        if src_pixmap and not src_pixmap.isNull():
            self._refresh()

    def set_source_pixmap(self, pixmap: "QPixmap | None"):
        """Replace the source pixmap and immediately refresh the display."""
        self._src = pixmap
        self._refresh()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._refresh()

    def _refresh(self):
        """Scale the source pixmap to the widget's current size."""
        if self._src is None or self._src.isNull():
            return
        w, h = self.width(), self.height()
        if w <= 0 or h <= 0:
            return
        scaled = self._src.scaled(
            w, h,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        # Call QLabel.setPixmap directly to avoid infinite recursion
        QLabel.setPixmap(self, scaled)
