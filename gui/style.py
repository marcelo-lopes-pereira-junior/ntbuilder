"""
gui/style.py
------------
Application-wide Qt stylesheet for NTBuilder.

Palette
-------
  Deep blue   #1A3A6B   (UnB / brand primary)
  Mid blue    #2851A3   (interactive / hover base)
  Sky blue    #4A7FD4   (hover highlight)
  Background  #F4F6FB   (window / panel background)
  Surface     #FFFFFF   (widget surface)
  Border      #C4CDE0   (subtle borders)
  Text        #1E1E2E   (primary text)
  Muted       #6B7A99   (secondary / placeholder text)
  Green       #1C7A3E   (Build / confirm actions)
  Green hover #25A356
  Red         #B03030   (destructive / warning)
"""

FONT_FAMILY = "Calibri"     # updated at runtime if unavailable
FONT_SIZE_PT = 10


def build_stylesheet(font_family: str = FONT_FAMILY,
                     font_size_pt: int = FONT_SIZE_PT) -> str:
    f  = font_family
    fs = font_size_pt
    return f"""
/* ── Global reset ────────────────────────────────────────────────────── */
* {{
    font-family: '{f}';
    font-size: {fs}pt;
    color: #1E1E2E;
}}

/* ── Window / container backgrounds ─────────────────────────────────── */
QMainWindow, QDialog {{
    background-color: #F4F6FB;
}}
/* Make ONLY decorative widgets inside the main window transparent — a
 * global ``QWidget {{ background-color: transparent; }}`` would propagate
 * to every dialog (QFileDialog, QMessageBox, …) and let the desktop
 * background show through, which is what we used to see on the file
 * dialog.  Targeting QLabel / QCheckBox / QRadioButton restores their
 * default see-through behaviour inside our panels without affecting
 * dialogs. */
QLabel, QCheckBox, QRadioButton {{
    background-color: transparent;
}}
QSplitter {{
    background-color: #F4F6FB;
}}
QSplitter::handle {{
    background-color: #C4CDE0;
    width: 3px;
}}

/* ── Menu bar ────────────────────────────────────────────────────────── */
QMenuBar {{
    background-color: #1A3A6B;
    color: #FFFFFF;
    padding: 2px 4px;
    spacing: 4px;
}}
QMenuBar::item {{
    background: transparent;
    padding: 4px 10px;
    border-radius: 4px;
}}
QMenuBar::item:selected {{
    background-color: #2851A3;
}}
QMenu {{
    background-color: #FFFFFF;
    border: 1px solid #C4CDE0;
    border-radius: 6px;
    padding: 4px 0;
}}
QMenu::item {{
    padding: 6px 22px 6px 14px;
    border-radius: 4px;
}}
QMenu::item:selected {{
    background-color: #EBF0FA;
    color: #1A3A6B;
}}
QMenu::separator {{
    height: 1px;
    background: #E0E5F0;
    margin: 3px 8px;
}}

/* ── Status bar ──────────────────────────────────────────────────────── */
QStatusBar {{
    background-color: #EEF1F8;
    border-top: 1px solid #C4CDE0;
    color: #4A5470;
    font-size: {fs - 1}pt;
}}

/* ── GroupBox ────────────────────────────────────────────────────────── */
QGroupBox {{
    background-color: #FFFFFF;
    border: 1.5px solid #C4CDE0;
    border-radius: 8px;
    margin-top: 10px;
    padding-top: 6px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 6px;
    color: #1A3A6B;
    font-weight: bold;
    font-size: {fs}pt;
}}

/* ── Push buttons ────────────────────────────────────────────────────── */
QPushButton {{
    background-color: #2851A3;
    color: #FFFFFF;
    border: none;
    border-radius: 6px;
    padding: 5px 16px;
    font-weight: bold;
    min-height: 24px;
}}
QPushButton:hover {{
    background-color: #4A7FD4;
}}
QPushButton:pressed {{
    background-color: #1A3A6B;
}}
QPushButton:disabled {{
    background-color: #B8C2D8;
    color: #E8EBF5;
}}

/* Build button — prominent green */
QPushButton#btn_build {{
    background-color: #1C7A3E;
    font-size: {fs + 1}pt;
    padding: 7px 20px;
    border-radius: 7px;
}}
QPushButton#btn_build:hover {{
    background-color: #25A356;
}}
QPushButton#btn_build:pressed {{
    background-color: #145C2E;
}}

/* Find primitive cell — secondary action */
QPushButton#btn_prim {{
    background-color: #FFFFFF;
    color: #2851A3;
    border: 1.5px solid #2851A3;
}}
QPushButton#btn_prim:hover {{
    background-color: #EBF0FA;
}}

/* Export button — teal green */
QPushButton#btn_export {{
    background-color: #1C7A5A;
}}
QPushButton#btn_export:hover {{
    background-color: #25A37A;
}}
QPushButton#btn_export:pressed {{
    background-color: #145C43;
}}

/* Batch button — amber */
QPushButton#btn_sys {{
    background-color: #8A6000;
}}
QPushButton#btn_sys:hover {{
    background-color: #B07D00;
}}
QPushButton#btn_sys:pressed {{
    background-color: #6A4800;
}}

/* ── Spin boxes ──────────────────────────────────────────────────────── */
QDoubleSpinBox, QSpinBox {{
    background-color: #FFFFFF;
    border: 1.5px solid #C4CDE0;
    border-radius: 5px;
    padding: 3px 6px;
    min-height: 22px;
    selection-background-color: #2851A3;
}}
QDoubleSpinBox:focus, QSpinBox:focus {{
    border-color: #2851A3;
}}
QDoubleSpinBox:disabled, QSpinBox:disabled {{
    background-color: #F0F2F8;
    color: #9BA8C0;
    border-color: #DDE2EE;
}}
QDoubleSpinBox::up-button, QSpinBox::up-button,
QDoubleSpinBox::down-button, QSpinBox::down-button {{
    width: 16px;
    border-radius: 3px;
}}

/* ── Combo box ───────────────────────────────────────────────────────── */
QComboBox {{
    background-color: #FFFFFF;
    border: 1.5px solid #C4CDE0;
    border-radius: 5px;
    padding: 3px 8px;
    min-height: 22px;
    selection-background-color: #2851A3;
}}
QComboBox:focus {{
    border-color: #2851A3;
}}
QComboBox::drop-down {{
    border: none;
    width: 20px;
}}
QComboBox QAbstractItemView {{
    background-color: #FFFFFF;
    border: 1px solid #C4CDE0;
    border-radius: 4px;
    selection-background-color: #EBF0FA;
    selection-color: #1A3A6B;
    outline: none;
}}

/* ── Check box ───────────────────────────────────────────────────────── */
QCheckBox {{
    spacing: 7px;
}}
QCheckBox::indicator {{
    width: 16px;
    height: 16px;
    border: 1.5px solid #C4CDE0;
    border-radius: 4px;
    background-color: #FFFFFF;
}}
QCheckBox::indicator:checked {{
    background-color: #2851A3;
    border-color: #2851A3;
    image: none;
}}
QCheckBox::indicator:hover {{
    border-color: #4A7FD4;
}}

/* ── Labels ──────────────────────────────────────────────────────────── */
QLabel {{
    background: transparent;
    color: #1E1E2E;
}}

/* ── Scroll area / bar ───────────────────────────────────────────────── */
QScrollArea {{
    border: none;
    background-color: transparent;
}}
QScrollBar:vertical {{
    background-color: #F0F2F8;
    width: 8px;
    border-radius: 4px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background-color: #B0BAD0;
    border-radius: 4px;
    min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{
    background-color: #2851A3;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
QScrollBar:horizontal {{
    background-color: #F0F2F8;
    height: 8px;
    border-radius: 4px;
}}
QScrollBar::handle:horizontal {{
    background-color: #B0BAD0;
    border-radius: 4px;
    min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{
    background-color: #2851A3;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}

/* ── Tab widget ──────────────────────────────────────────────────────── */
QTabWidget::pane {{
    border: 1.5px solid #C4CDE0;
    border-radius: 0 6px 6px 6px;
    background-color: #FFFFFF;
}}
QTabBar::tab {{
    background-color: #E8EDF8;
    border: 1.5px solid #C4CDE0;
    border-bottom: none;
    border-radius: 6px 6px 0 0;
    padding: 5px 14px;
    margin-right: 2px;
    color: #4A5470;
}}
QTabBar::tab:selected {{
    background-color: #FFFFFF;
    color: #1A3A6B;
    font-weight: bold;
}}
QTabBar::tab:hover:!selected {{
    background-color: #D5DDF5;
}}

/* ── Line edit ───────────────────────────────────────────────────────── */
QLineEdit {{
    background-color: #FFFFFF;
    border: 1.5px solid #C4CDE0;
    border-radius: 5px;
    padding: 3px 7px;
    selection-background-color: #2851A3;
}}
QLineEdit:focus {{
    border-color: #2851A3;
}}

/* ── Tooltip ─────────────────────────────────────────────────────────── */
QToolTip {{
    background-color: #1A3A6B;
    color: #FFFFFF;
    border: none;
    padding: 5px 8px;
    border-radius: 5px;
    font-size: {fs - 1}pt;
}}

/* ── Message / dialog boxes ──────────────────────────────────────────── */
QMessageBox {{
    background-color: #F4F6FB;
}}
QMessageBox QPushButton {{
    min-width: 80px;
}}
"""
