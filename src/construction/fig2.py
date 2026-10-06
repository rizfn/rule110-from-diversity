"""Figure 2 (python src/construction/fig2.py): layout of the rules, the rule-110 table rebuilt from the mixture,
and a decoded simulation. -> plots/"""
import sys, pathlib; HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE.parent / "common"))
import numpy as np
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import to_rgb, ListedColormap
from style import plt, setup, NAVY, CORAL, YELLOW, TEAL
from mixture import mixture_run, encode, eca_step
setup(15)
AG = [36, 37, 5]; COL = [CORAL, YELLOW, TEAL]
CODE = {0: (1, 0, 0), 1: (0, 0, 0)}          # logical bit -> physical block
def tint(c, a=0.22): return tuple(1 - a * (1 - v) for v in to_rgb(c))
def shade(c):  # black site: halfway between a saturated dark version and a navy mix of the rule colour
    import colorsys
    h, l, s_ = colorsys.rgb_to_hls(*to_rgb(c))
    sat = colorsys.hls_to_rgb(h, l * 0.62, min(1.0, s_ * 1.15))
    mix = tuple(0.5 * v + 0.5 * w for v, w in zip(to_rgb(c), to_rgb(NAVY)))
    return tuple(0.5 * a + 0.5 * b for a, b in zip(sat, mix))
def cell(ax, x, y, on, p=None, s=1.0, lw=1.4):
    fc = (shade(COL[p]) if on else tint(COL[p], 0.35)) if p is not None else (NAVY if on else "white")
    ax.add_patch(Rectangle((x, y), s, s, fc=fc, ec=NAVY, lw=lw))
def panel_label(ax, s, x=-0.02, y=1.0):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")

fig = plt.figure(figsize=(7.0, 8.2))
gs = fig.add_gridspec(3, 1, height_ratios=[0.95, 2.35, 2.6], hspace=0.22, left=0.12, right=0.99, top=0.98, bottom=0.07)

# (a) layout + code -----------------------------------------------------------
ax = fig.add_subplot(gs[0]); ax.axis("off"); ax.set_aspect("equal")
for i in range(12):
    ax.add_patch(Rectangle((i, 0), 1, 1, fc=COL[i % 3], ec=NAVY, lw=1.4))
for b in range(4):                                  # brackets grouping blocks
    ax.plot([3*b + .1, 3*b + .1, 3*b + 2.9, 3*b + 2.9], [-0.25, -0.45, -0.45, -0.25], color=NAVY, lw=1.8)
for k, (n, c) in enumerate(zip(AG, COL)):
    ax.text(k + 0.5, 1.25, str(n), ha="center", va="bottom", fontsize=15, fontweight="bold", color=NAVY)
x0 = 14.2
for row, bit in enumerate([0, 1]):
    y = 0.75 - 1.5 * row
    cell(ax, x0, y, bit)
    ax.add_patch(FancyArrowPatch((x0 + 1.25, y + .5), (x0 + 2.75, y + .5), arrowstyle="-|>", mutation_scale=18, color=NAVY, lw=1.8))
    for k in range(3): cell(ax, x0 + 3 + k, y, CODE[bit][k], k)
ax.set_xlim(-0.3, x0 + 6.3); ax.set_ylim(-1.0, 2.1)
panel_label(ax, "A", x=0.0)

# (b) light cones: one per logical neighbourhood ------------------------------
gsb = gs[1].subgridspec(2, 4, wspace=0.12, hspace=0.15)
for idx, n in enumerate(range(7, -1, -1)):
    ax = fig.add_subplot(gsb[idx // 4, idx % 4]); ax.axis("off"); ax.set_aspect("equal")
    lcr = [(n >> 2) & 1, (n >> 1) & 1, n & 1]
    for k, b in enumerate(lcr): cell(ax, 3 * k + 1, 6.4, b, s=1.0)         # logical neighbourhood
    w = np.array([v for b in lcr for v in CODE[b]], np.uint8)
    tabs = {r: [(r >> q) & 1 for q in range(8)] for r in AG}
    rows = [w.copy()]
    lo, hi = 0, 9
    for t in range(3):
        new = w.copy()
        for q in range(lo + 1, hi - 1):
            new[q] = tabs[AG[q % 3]][(w[q-1] << 2) | (w[q] << 1) | w[q+1]]
        w = new; lo, hi = lo + 1, hi - 1; rows.append(w.copy())
    for t, row in enumerate(rows):
        for q in range(t, 9 - t): cell(ax, q, 4.6 - t * 1.05, row[q], q % 3, lw=1.1)
    out_block = tuple(rows[3][3:6]); out = 1 if out_block == CODE[1] else 0
    assert out_block in CODE.values() and out == (110 >> n) & 1
    ax.add_patch(FancyArrowPatch((4.5, 1.25), (4.5, 0.25), arrowstyle="-|>", mutation_scale=14, color=NAVY, lw=1.6))
    cell(ax, 4.0, -0.95, out)
    ax.set_xlim(-0.2, 9.2); ax.set_ylim(-1.1, 7.5)
    if idx == 0: panel_label(ax, "B", x=-0.05)

# (c) space-time: mixture (coloured by rule) and decoded rule 110 ------------
gsc = gs[2].subgridspec(1, 2, wspace=0.42)
NL, TL, SEED = 60, 140, 13
T0, J0 = 104, 29   # window shown in panel C (logical units); axes relabelled from 0
x0 = np.random.default_rng(SEED).integers(0, 2, NL).astype(np.uint8)
P = encode(x0); phys = [P]
for _ in range(3 * TL - 1): P = mixture_run(P, 1); phys.append(P)
phys = np.array(phys)
img = np.zeros(phys.shape + (3,))
for p in range(3):
    on, off = np.array(shade(COL[p])), np.array(tint(COL[p], 0.35))
    img[:, p::3] = np.where(phys[:, p::3, None] == 1, on, off)
dec = (phys[::3].reshape(TL, NL, 3) == CODE[1]).all(-1).astype(np.uint8)
ref = [x0]
for _ in range(TL - 1): ref.append(eca_step(ref[-1], 110))
assert np.array_equal(dec, np.array(ref))
ax1 = fig.add_subplot(gsc[0]); ax1.imshow(np.roll(img, -3*J0, axis=1)[3*T0:3*T0+90, :90], interpolation="nearest", aspect="auto")
ax1.set_xlabel("Site"); ax1.set_ylabel("Time"); ax1.set_xticks([0, 45, 89]); ax1.set_yticks([0, 45, 89])
panel_label(ax1, "C", x=-0.2)
ax2 = fig.add_subplot(gsc[1]); ax2.imshow(np.roll(dec, -J0, axis=1)[T0:T0+30, :30], cmap=ListedColormap(["white", NAVY]), interpolation="nearest", aspect="auto")
ax2.set_xlabel("Logical site"); ax2.set_xticks([0, 15, 29]); ax2.set_yticks([0, 15, 29]); ax2.set_ylabel("Logical time")
for a in (ax1, ax2):
    for sp in a.spines.values(): sp.set_visible(False)
(HERE / "plots").mkdir(exist_ok=True)
fig.savefig(HERE / "plots" / f"fig2_layout_and_derivation_code100-000_tau3_NL{NL}_seed{SEED}_t{T0:04d}_j{J0:03d}.pdf")
