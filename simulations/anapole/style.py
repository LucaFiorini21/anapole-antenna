"""Shared plot style: one colour per physical entity, everywhere."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = dict(dipole="#2a78d6", coil="#eb6834", ret="#1baf7a", total="#0b0b0b",
           rev="#e87ba4", theory="#898781", leads="#eda100", extra="#6250d6")
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.size": 10, "font.family": "DejaVu Sans",
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.titlelocation": "left",
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "legend.frameon": False, "legend.fontsize": 9,
    "figure.dpi": 110, "savefig.dpi": 170, "savefig.bbox": "tight",
})


def save(fig, name):
    import os
    d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig")
    fig.savefig(os.path.join(d, name + ".png"))
    plt.close(fig)
