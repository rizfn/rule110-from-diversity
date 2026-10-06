"""Figure 1 (python src/construction/fig1.py): rule tables of rules 36, 37, 5 and 110 (two rows of four neighbourhoods each),
with zoomed space-time diagrams from a random start. -> plots/"""
import sys, pathlib; HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE.parent / "common"))
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from style import plt, setup, NAVY, CORAL, YELLOW, TEAL
setup(15)
RULES = [(36, CORAL), (37, YELLOW), (5, TEAL), (110, NAVY)]
N, T, SEED = 40, 24, 4

def eca_run(r, N, T, seed):
    t = np.array([(r >> k) & 1 for k in range(8)], np.uint8)
    x = np.random.default_rng(seed).integers(0, 2, N).astype(np.uint8); M = []
    for _ in range(T): M.append(x); x = t[(np.roll(x, 1) << 2) | (x << 1) | np.roll(x, -1)]
    return np.array(M)

def table(ax, r, col):
    s, gap, rowgap, lw = 1.0, 1.0, 2.9, 1.3
    for j, n in enumerate(range(7, -1, -1)):           # Wolfram order 111 ... 000
        x0 = (j % 4) * (3 * s + gap); y0 = rowgap if j < 4 else 0
        bits = [(n >> 2) & 1, (n >> 1) & 1, n & 1]
        for k, b in enumerate(bits):                   # fills without edges
            ax.add_patch(Rectangle((x0 + k * s, y0 + 1.2), s, s, fc=NAVY if b else "white", ec="none"))
        ax.add_patch(Rectangle((x0, y0 + 1.2), 3 * s, s, fc="none", ec=NAVY, lw=lw))   # one outline
        for k in (1, 2):                                # single dividing lines
            ax.plot([x0 + k * s] * 2, [y0 + 1.2 + 0.06, y0 + 1.2 + s - 0.06] if bits[k-1] and bits[k] else [y0 + 1.2, y0 + 1.2 + s], color="white" if bits[k-1] and bits[k] else NAVY, lw=lw, solid_capstyle="butt")
        ax.add_patch(Rectangle((x0 + s, y0), s, s, fc=col if (r >> n) & 1 else "white", ec=NAVY, lw=lw))
    ax.set_xlim(-0.2, 4 * (3 * s + gap) - gap + 0.2); ax.set_ylim(-0.2, rowgap + 2.4)
    ax.set_aspect("equal"); ax.axis("off")

fig = plt.figure(figsize=(7.0, 7.4))
gs = fig.add_gridspec(4, 3, width_ratios=[0.85, 2.0, 2.3], wspace=0.1, hspace=0.3,
                      left=0.01, right=0.99, top=0.985, bottom=0.015)
for i, (r, col) in enumerate(RULES):
    axl = fig.add_subplot(gs[i, 0]); axl.axis("off")
    axl.text(0.55, 0.5, f"Rule {r}", ha="center", va="center", fontsize=18, color=NAVY)
    axl.plot([0.12, 0.98], [0.33, 0.33], color=col, lw=5, solid_capstyle="round", transform=axl.transAxes)
    axl.set_xlim(0, 1); axl.set_ylim(0, 1)
    axl.text(0.0, 1.0, "ABCD"[i], transform=axl.transAxes, fontsize=20, fontweight="bold", va="top", ha="left")
    table(fig.add_subplot(gs[i, 1]), r, col)
    ax = fig.add_subplot(gs[i, 2])
    ax.imshow(eca_run(r, N, T, SEED), cmap=ListedColormap(["white", col]), interpolation="nearest", aspect="equal")
    ax.set_xticks(np.arange(-.5, N, 1), minor=True); ax.set_yticks(np.arange(-.5, T, 1), minor=True)
    ax.grid(which="minor", color="#d9dde0", lw=0.4); ax.tick_params(which="both", length=0)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(NAVY); sp.set_linewidth(1.6)
(HERE / "plots").mkdir(exist_ok=True)
fig.savefig(HERE / "plots" / f"fig1_rule_tables_{'_'.join(str(r) for r, _ in RULES)}_N{N}_T{T}_seed{SEED}.pdf")
