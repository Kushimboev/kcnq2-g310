#!/usr/bin/env python3
"""Colour dictionary of the manuscript figures: one colour has one meaning throughout the paper.

Usage:
    import sys; sys.path.insert(0, "<repository>/scripts/style")
    from palette import C, use_paper_style, panel_label, MM
    use_paper_style()

A reader learns a colour once and can then read every panel without a legend.
"""
import os
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
STYLE = os.path.join(_HERE, "paper.mplstyle")


def use_paper_style():
    plt.style.use(STYLE)


# --- journal widths (mm -> inch) ---
class MM:
    ONE_COL = 89 / 25.4       # 3.50"  one column
    ONE_HALF = 120 / 25.4     # 4.72"
    TWO_COL = 183 / 25.4      # 7.20"  two columns


class C:
    """Colour dictionary. Do not change these values: the figures rely on them being consistent."""

    # === Gate residues ===
    SER314 = "#D95F02"        # serine of the GSG motif (KCNQ2 Ser314)
    GLY310 = "#1B7EBD"        # glycine of the PAG motif (KCNQ2 Gly310)
    ALA306 = "#7570B3"        # Ala306

    # === Particles ===
    K = "#E7298A"             # K+ (pink, as in Mkrtchyan 2024, for easy comparison)
    KSF = "#66A61E"           # K+ placed in the filter at build time; kept distinct from the free ions
    CL = "#999999"            # Cl-
    WATER = "#92C5DE"         # water
    LIPID = "#7FCDBB"         # POPC membrane (35% transparent)

    # === Semantic colours ===
    ZERO = "#2B2B2B"          # our measurement (zero line)
    EXPECT = "#C7C7C7"        # expected from experiment (background, not data)
    EXPERIMENT = "#555555"    # experimental value / band
    CONDUCT = "#1A9850"       # conducting
    NOCONDUCT = "#D73027"     # not conducting
    HILITE = "#FDBF6F"        # highlight shading (filter zone etc.)

    # === Gate-width series: sequential quantity, sequential colours ===
    WIDTH = {
        "0.84": "#440154", "1.00": "#414487", "1.10": "#2A788E",
        "1.15": "#22A884", "1.25": "#7AD151", "1.30": "#BBDF27", "1.40": "#FDE725",
    }

    @staticmethod
    def width(nm):
        """'1.15' or 1.15 -> colour."""
        return C.WIDTH.get(f"{float(nm):.2f}", "#666666")


# --- data constants (kept in one place) ---
K_HYDRATED_R_NM = 0.31        # radius of hydrated K+
G_EXPERIMENT_PS = 5.8         # Selyanko 2001
G_EXPERIMENT_SD = 0.3
E_CHARGE = 1.602176634e-19


def panel_label(ax, letter, dx=-0.16, dy=1.04, size=11):
    """Panel letter: bold, at the top left of the panel and outside the axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=size, fontweight="bold", va="bottom", ha="left")


def poisson_ci(n, cl=0.95):
    """Poisson 95% confidence interval for an event count ((0, 3.0) for n = 0)."""
    from scipy.stats import chi2
    lo = 0.0 if n == 0 else chi2.ppf((1 - cl) / 2, 2 * n) / 2
    hi = chi2.ppf(1 - (1 - cl) / 2, 2 * n + 2) / 2
    return lo, hi
