"""
main.py — entry point for the Nanotube Builder desktop application.

Usage:
    python main.py
"""

import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QFontDatabase
from gui.main_window import MainWindow
from gui.style import build_stylesheet


def main():
    # High-DPI support
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )
    app = QApplication(sys.argv)
    app.setApplicationName("NTBuilder")
    app.setOrganizationName("UnB")

    # Choose a modern UI font that is visibly different from the system default.
    # Priority: Calibri (bundled with Windows/Office, humanist and clean),
    # Candara, then Segoe UI variants, then macOS/Linux fallbacks.
    _families = set(QFontDatabase.families())
    _prefer   = [
        "Calibri",
        "Candara",
        "Segoe UI Variable",
        "Segoe UI",
        "Helvetica Neue",
        "Helvetica",
        "Arial",
    ]
    _chosen = next((f for f in _prefer if f in _families), None)
    _pt     = 10

    if _chosen:
        app.setFont(QFont(_chosen, _pt))
    else:
        _fb = app.font()
        _fb.setPointSize(_pt)
        app.setFont(_fb)
        _chosen = _fb.family()

    # Apply the full visual stylesheet (colours, borders, buttons, etc.)
    app.setStyleSheet(build_stylesheet(_chosen, _pt))

    win = MainWindow()
    # showNormal() first forces Qt to resolve the window geometry against the
    # correct screen before maximising — avoids the QWindowsWindow::setGeometry
    # warning that occurs when the initial placement lands off-screen on
    # multi-monitor setups with mixed DPI.
    win.showNormal()
    win.showMaximized()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
