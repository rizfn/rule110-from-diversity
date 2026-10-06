"""Figure 3 (python src/random_starts/fig3.py): locking from random initial conditions. Reads outputs/, writes plots/.
Beat colours (periwinkle, plum, orchid): the time phase (mod 3) at which a block is a valid code word; dark: errors.
A1-A3: schematics of the three fates of the walls left after the fast stage (arrows: direction each wall moves).
B1-B3: hand-made starts on 256 blocks (perfect rule-110 domains whose beats differ by one step).
       All in the lab frame. B1 balanced; B2 excess left-moving walls; B3 excess right-moving walls.
C1: median lock time vs size, all runs from random starts (bars: quartiles).
C2: the same, split by the walls left at t = 3N (512 runs per size, lockcharge.c); lines: power-law fits
    (balanced 128-1,024; excess left-moving 512-16,384); hollow: points outside the fit range (balanced at 2,048 has 14 runs)."""
import sys, pathlib; HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE.parent / "common"))
OUT = HERE / "outputs"
import numpy as np, glob, re
from mixture import mixture_step, is_codeword
from beats import beat_rows, make_ic
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch
from style import plt, setup, NAVY, CORAL, YELLOW, TEAL
setup(15)
PERI, PLUM, ORCHID = "#7088C4", "#7A4F8F", "#C384B5"
CM = ListedColormap([NAVY, PERI, PLUM, ORCHID]); P_, L_, O_ = 0, 1, 2
fig = plt.figure(figsize=(7.0, 7.4))
gA = fig.add_gridspec(1, 3, left=0.145, right=0.97, top=0.965, bottom=0.81, wspace=0.5)
gB = fig.add_gridspec(1, 3, left=0.145, right=0.97, top=0.77, bottom=0.505, wspace=0.5)
gC = fig.add_gridspec(1, 2, left=0.145, right=0.97, top=0.39, bottom=0.095, wspace=0.42)
lab = lambda ax, s: ax.text(-0.40, 1.0, s, transform=ax.transAxes, fontsize=18, fontweight="bold", va="top", ha="left")
# ---------------- B1-B3: hand-made starts on 256 blocks, evolved and coloured by beat ----------------------------
def lock_time(x, tmax):
    for t in range(1, tmax + 1):
        x = mixture_step(x)
        if is_codeword(x).all(): return t
    return -1
x = make_ic(256, (56, 144, 56), (0, 2, 0), seed=1)                       # balanced
d = {"zero": {"R": beat_rows(x, lock_time(x.copy(), 3000) + 150)[1]}}
x = make_ic(256, (64, 20, 24, 148), (0, 1, 2, 0), seed=2)                # excess left-moving walls (locks at 3,326)
d["neg2"] = {"R": beat_rows(x, lock_time(x.copy(), 20000) + 500)[1]}
x = make_ic(256, (41, 66, 66, 83), (0, 2, 1, 0), seed=1)                 # excess right-moving walls (never locks)
d["pos2"] = {"R": beat_rows(x, 1600)[1]}

# ---------------- B1-B3 data and their colour assignment ----------------------------------------------------
lut = {"zero": {0: P_, 1: O_, 2: L_}, "neg": {0: L_, 2: O_, 1: P_}, "pos": {0: P_, 1: O_, 2: L_}}   # neg2/pos2 share these
t0 = {"zero": 9, "neg": 21, "pos": 11}
def colour(M, key):
    C = np.full(M.shape, -1)
    for b, c in lut[key].items(): C[M == b] = c
    return C

# ---------------- A1-A3: simple schematics (colours as in B1-B3), lab frame ---------------------------------
NX, NT = 400, 300
X = np.linspace(0, 1, NX)
def paint(edges_cols, t):
    """edges_cols: list of (x_left_edge, colour) pieces covering [0, 1] from left to right"""
    row = np.zeros(NX, int)
    for (x0, c), nxt in zip(edges_cols, edges_cols[1:] + [(1.0, None)]): row[(X >= x0) & (X < nxt[0])] = c
    return row
# A1 balanced: two walls close in and annihilate
def A1row(t):
    l, r = 0.2 + 0.4 * t, 0.8 - 0.4 * t
    return paint([(0, P_), (l, O_), (r, P_)], t) if t < 0.75 else paint([(0, P_)], t)
A1lines = [([0.2 + 0.4 * t for t in np.linspace(0, 0.75, 50)], np.linspace(0, 0.75, 50)),
           ([0.8 - 0.4 * t for t in np.linspace(0, 0.75, 50)], np.linspace(0, 0.75, 50))]
# A2 excess left-moving: a cluster drifts left, with a thin orchid band pinned to the leading wall and the third
#    wall closing in; the cluster collapses and orchid sweeps rightwards round the (periodic) lattice
BAND, TC, VF = 0.09, 0.6, 3.5
a = lambda t: 0.62 - 0.35 * t
c = lambda t: a(t) + BAND + 0.25 * max(0.0, 1 - t / TC)
front = lambda t: a(TC) + BAND + VF * (t - TC)                      # right edge of orchid after the collapse
T1 = TC + (1 - a(TC) - BAND) / VF                                   # front reaches the right edge and wraps
TW = (1 - a(TC) - BAND + VF * TC + 0.62) / (VF + 0.35)              # wrapped front meets the leading wall
def A2row(t):
    if t < TC: return paint([(0, L_), (a(t), O_), (a(t) + BAND, P_), (c(t), L_)], t)
    if t < T1: return paint([(0, L_), (a(t), O_), (front(t), L_)], t)
    if t < TW: return paint([(0, O_), (front(t) - 1, L_), (a(t), O_)], t)
    return paint([(0, O_)], t)
tt = np.linspace(0, 1, 600)
seg = lambda lo, hi: tt[(tt >= lo) & (tt <= hi)]
A2lines = [([a(t) for t in seg(0, TW)], seg(0, TW)), ([a(t) + BAND for t in seg(0, TC)], seg(0, TC)),
           ([c(t) for t in seg(0, TC)], seg(0, TC)), ([front(t) for t in seg(TC, T1)], seg(TC, T1)),
           ([front(t) - 1 for t in seg(T1, TW)], seg(T1, TW))]
# A3 excess right-moving: three walls drift right together and never meet
xs3 = (0.08, 0.34, 0.62)
def A3row(t): return paint([(0, P_), (xs3[0] + 0.3 * t, O_), (xs3[1] + 0.3 * t, L_), (xs3[2] + 0.3 * t, P_)], t)
A3lines = [([x0 + 0.3 * t for t in tt], tt) for x0 in xs3]
T = np.linspace(0, 1, NT)
schem = [(np.array([A1row(t) for t in T]), A1lines, [(0.2, +1), (0.8, -1)]),
         (np.array([A2row(t) for t in T]), A2lines, [(0.62, -1), (0.75, -1), (0.92, -1)]),
         (np.array([A3row(t) for t in T]), A3lines, [(0.08, +1), (0.34, +1), (0.62, +1)])]
for k, (F, lines, arr) in enumerate(schem):
    ax = fig.add_subplot(gA[0, k])
    ax.imshow(F, cmap=CM, vmin=-1, vmax=2, aspect="auto", interpolation="nearest", extent=[0, 1, 1, 0])
    for xl, tl_ in lines: ax.plot(xl, tl_, color=NAVY, lw=1.3, solid_capstyle="round")
    for x0, dd in arr:
        ax.add_patch(FancyArrowPatch((x0 - 0.07 * dd, 0.12), (x0 + 0.07 * dd, 0.12), arrowstyle="-|>", mutation_scale=10, color="white", lw=1.5, zorder=5))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(1, 0)
    for s_ in ax.spines.values(): s_.set_visible(True)
    lab(ax, f"A{k + 1}")

ext = {"zero": (9, 9 + len(d["zero"]["R"])), "neg": (9, 9 + len(d["neg2"]["R"])), "pos": (9, 9 + len(d["pos2"]["R"]))}
ticks = {"zero": [(200, "200"), (400, "400")], "neg": [(1000, "1k"), (2000, "2k"), (3000, "3k")], "pos": [(500, "500"), (1000, "1k"), (1500, "1.5k")]}
for k, key in enumerate(("zero", "neg", "pos")):
    M = d[{"zero": "zero", "neg": "neg2", "pos": "pos2"}[key]]["R"]; a, b = ext[key]
    ax = fig.add_subplot(gB[0, k])
    ax.imshow(colour(M, key), cmap=CM, vmin=-1, vmax=2, aspect="auto", interpolation="nearest", extent=[0, 256, b, a])
    ax.set_xticks([0, 256]); ax.set_xlabel("Block")
    ax.set_yticks([t for t, _ in ticks[key]]); ax.set_yticklabels([s for _, s in ticks[key]])
    if k == 0: ax.set_ylabel("Time")
    for s in ax.spines.values(): s.set_visible(True)
    lab(ax, f"B{k + 1}")

# ---------------- C1, C2: lock time vs size -----------------------------------------------------------------
Ns = [32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
med, q1, q3 = [], [], []
for n in Ns:
    v = np.concatenate([np.loadtxt(f, ndmin=1) for f in glob.glob(str(OUT / "lock_times" / f"lock_N{n}.txt"))]); v = np.where(v > 0, v, np.inf)
    med.append(np.median(v)); q1.append(np.percentile(v, 25)); q3.append(np.percentile(v, 75))
for n in (8192, 16384):                                                   # largest sizes: lock times from lockcharge.c
    v = np.concatenate([np.loadtxt(f, ndmin=2)[:, 1] for f in glob.glob(str(OUT / "lock_charge" / f"charge_N{n}_tmax*_s*.txt"))]); v = np.where(v > 0, v, np.inf)
    Ns.append(n); med.append(np.median(v)); q1.append(np.percentile(v, 25)); q3.append(np.percentile(v, 75))
Ns, med, q1, q3 = map(np.array, (Ns, med, q1, q3))
eb = dict(fmt="o", ms=7, mec=NAVY, mew=1.2, ecolor=NAVY, elinewidth=1.1, capsize=2.5, zorder=3)
def axes_style(ax):
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(25, 25000); ax.set_ylim(30, 4e9); ax.minorticks_off()
    ax.set_xticks([32, 512, 8192]); ax.set_xticklabels(["32", "512", "8192"]); ax.set_xlabel("Lattice size")
ax = fig.add_subplot(gC[0, 0])
ax.errorbar(Ns, med, yerr=[med - q1, q3 - med], color=TEAL, **eb)
ok = Ns >= 512; z, c0 = np.polyfit(np.log(Ns[ok]), np.log(med[ok]), 1)
xs = np.array([350, 25000]); ax.plot(xs, np.exp(c0) * xs ** z, color=CORAL, lw=2.2, zorder=2, label=f"Slope {z:.1f}")
axes_style(ax); ax.set_ylabel("Time to lock"); ax.legend(frameon=False, fontsize=12, loc="upper left", handlelength=1.2)
ax.text(-0.33, 1.08, "C1", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
ax = fig.add_subplot(gC[0, 1])
cls = {}
for f in glob.glob(str(OUT / "lock_charge" / "charge_N*_tmax*_s*.txt")):
    n = int(re.search(r"charge_N(\d+)", f).group(1)); dd_ = np.loadtxt(f, ndmin=2)
    for q, lt in dd_: cls.setdefault(n, []).append((int(q), lt))
fits = {}
for name, sel, col, fit_range in (("Balanced", lambda q: q == 0, YELLOW, (128, 1024)), ("Excess left-moving", lambda q: q < 0, CORAL, (512, 16384))):
    nn, m, a_, b_, cnt = [], [], [], [], []
    for n in sorted(cls):
        v = np.array([lt for q, lt in cls[n] if sel(q)])
        if len(v) < 10: continue
        v = np.where(v > 0, v, np.inf); nn.append(n); m.append(np.median(v)); a_.append(np.percentile(v, 25)); b_.append(np.percentile(v, 75)); cnt.append(len(v))
    nn, m, a_, b_, cnt = map(np.array, (nn, m, a_, b_, cnt))
    inr = (nn >= fit_range[0]) & (nn <= fit_range[1]); zz, cc = np.polyfit(np.log(nn[inr]), np.log(m[inr]), 1); fits[name] = zz
    ax.errorbar(nn[inr], m[inr], yerr=[m[inr] - a_[inr], b_[inr] - m[inr]], color=col, label=f"{name}\nslope {zz:.1f}", **eb)
    if (~inr).any():                                   # points outside the fit range (few runs): hollow
        ax.errorbar(nn[~inr], m[~inr], yerr=[m[~inr] - a_[~inr], b_[~inr] - m[~inr]], fmt="o", ms=7, mfc="white", mec=col, mew=1.4, ecolor=col, elinewidth=1.1, capsize=2.5, zorder=3)
    xs = np.array([fit_range[0] / 1.5, fit_range[1] * 1.5]); ax.plot(xs, np.exp(cc) * xs ** zz, color=col, lw=2.0, zorder=2)
    print(name, "slope %.3f over %s; runs per size:" % (zz, fit_range), dict(zip(nn.tolist(), cnt.tolist())))
axes_style(ax); ax.legend(frameon=False, fontsize=12, loc="upper left", handlelength=1.0, labelspacing=0.4, borderaxespad=0.2)
ax.text(-0.25, 1.08, "C2", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
(HERE / "plots").mkdir(exist_ok=True)
fig.savefig(HERE / "plots" / "fig3_random_starts.pdf")
print("z =", round(z, 3))
