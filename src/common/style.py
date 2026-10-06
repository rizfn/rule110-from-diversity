"""Shared figure style: fonts, sizes and the colour palette (rules: coral 36, yellow 37, teal 5; beats: BEATS)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
# Figures use Liberation Sans (metrically identical to Arial). If it is not installed, Arial or the default sans is used.
import glob as _glob
for _p in _glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
    try: font_manager.fontManager.addfont(_p)
    except Exception: pass
_FAMILY = ["Liberation Sans", "Arial", "DejaVu Sans"]
NAVY, CORAL, YELLOW, TEAL = "#2C4251", "#D16666", "#FFDD4A", "#58a4b0"
AGENT = {36: CORAL, 37: YELLOW, 5: TEAL}
def setup(base=14):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": _FAMILY, "font.size": base,
        "axes.labelsize": base + 2, "xtick.labelsize": base, "ytick.labelsize": base,
        "axes.linewidth": 1.6, "xtick.major.width": 1.6, "ytick.major.width": 1.6,
        "xtick.major.size": 6, "ytick.major.size": 6,
        "axes.edgecolor": NAVY, "axes.labelcolor": NAVY, "xtick.color": NAVY, "ytick.color": NAVY,
        "text.color": NAVY, "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False})

BEATS = ["#7088C4", "#7A4F8F", "#C384B5"]   # periwinkle, plum, orchid: the three beats (phases) of the code
