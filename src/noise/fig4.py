"""Figure 4 (python src/noise/fig4.py), reads outputs/, writes plots/: random bit flips in the dynamics. Beat colours (periwinkle, plum, orchid): time phase of the code in each block; dark: errors.
A 160 blocks from a perfect start with noise 3e-4: each surviving flip opens a bubble of another beat whose two walls
move apart and vanish when they meet other walls; B probability that one flipped cell escapes healing (still unhealed
after 60 steps), by the rule of the flipped cell and the step of the 3-step cycle at which it is flipped (400 flips each);
dashed: the average, p = 0.23
C errors per block vs noise; D rescaled approach to the steady state.
Old docstring:
A median lock time vs size from random starts (bars: quartiles; 512 runs per size up to 1024, 256 for 1536-2048, 64 for 3072, 128 for 4096); line: fit for N >= 512
B space-time map of errors at noise 3e-4 from a perfect start (160 blocks)
C errors per block vs noise rate; line: rate-equation prediction n* = 2 sqrt(3 p eps / v_rel), no free parameters
D approach to the steady state, rescaled by the predicted n* and tau = (3 p eps v_rel)^(-1/2); dashed: n*tanh(t/tau)"""
import sys, pathlib; HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE.parent / "common"))
import numpy as np, glob
from mixture import mixture_step, is_codeword, code_state
from beats import beat_rows
from params import P_ESC, VR, VL
from style import plt, setup, NAVY, CORAL, YELLOW, TEAL, BEATS
setup(15)
VREL = VR + VL                                    # measured by measure_parameters.py
fig, axs = plt.subplots(2, 2, figsize=(7.0, 6.6), gridspec_kw={"hspace": 0.5, "wspace": 0.55})
lab = lambda ax, s: ax.text(-0.52, 1.07, s, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top")

from matplotlib.colors import ListedColormap
PCM = ListedColormap([NAVY] + BEATS)
rng = np.random.default_rng(5)                                       # A: 160 blocks, perfect start, noise 3e-4
NM = {"A": beat_rows(code_state(160, rng), 900, noise=3e-4, rng=rng)[1]}
def escape_grid(n_blocks=128, flips=400, seed=7):
    """probability that one flipped cell is still unhealed after 60 steps, by step after a code time (rows) and by the
    cell of the block that is flipped (columns: rules 36, 37, 5)"""
    rng = np.random.default_rng(seed)
    X = np.zeros((flips, 3 * n_blocks), np.uint8); X[:, 0::3] = 1 - rng.integers(0, 2, (flips, n_blocks))
    for _ in range(300): X = mixture_step(X)
    G = np.zeros((3, 3))
    for ph in range(3):
        Y0 = X.copy()
        for _ in range(ph): Y0 = mixture_step(Y0)
        for off in range(3):
            Y = Y0.copy(); j = rng.integers(0, n_blocks, flips); Y[np.arange(flips), 3 * j + off] ^= 1
            healed = np.zeros(flips, bool)
            for _ in range(60): Y = mixture_step(Y); healed |= is_codeword(Y).all(1)
            G[ph, off] = 1 - healed.mean()
    return G
for ax, M, s, yt in ((axs[0, 0], NM["A"], "A", [0, 400, 800]),):
    ax.imshow(M, cmap=PCM, vmin=-1, vmax=2, aspect="auto", interpolation="nearest", extent=[0, M.shape[1], len(M), 0])
    ax.set_xlabel("Block"); ax.set_ylabel("Time"); ax.set_yticks(yt); ax.set_xticks([0, M.shape[1] // 2, M.shape[1]])
    for sp in ax.spines.values(): sp.set_visible(True)
    lab(ax, s)
# B: probability that a single flip escapes healing, by the rule of the flipped cell and the step of the cycle
ax = axs[0, 1]; G = escape_grid()          # rows: steps after a code time; columns: cells with rules 36, 37, 5
w = 0.26
for j, (rule, col) in enumerate(((36, CORAL), (37, YELLOW), (5, TEAL))):
    ax.bar(np.arange(3) + (j - 1) * w, G[:, j], width=w * 0.92, color=col, edgecolor=NAVY, lw=1.0, label=f"Rule {rule}", zorder=3)
ax.axhline(G.mean(), color=NAVY, lw=1.6, ls=(0, (3, 2)), label=r"$p$", zorder=4)
ax.set_xticks(range(3)); ax.set_xticklabels(["0", "1", "2"]); ax.set_xlabel("Step after a code time")
ax.set_ylabel("Escape probability"); ax.set_ylim(0, 0.75)
h, l = ax.get_legend_handles_labels(); o = [l.index(k) for k in ("Rule 36", "Rule 37", "Rule 5", "$p$")]
ax.legend([h[i] for i in o], [l[i] for i in o], frameon=False, fontsize=13, loc="upper left", ncol=2, handlelength=1.0, columnspacing=0.8, labelspacing=0.2)
lab(ax, "B")

E, NS, curves = [], [], []
for fn in glob.glob(str(HERE / "outputs" / "noise_NB2048_T200000_seed41" / "noise_eps*.txt")):
    t, r = np.loadtxt(fn).T; e = float(fn.split("eps")[-1][:-4]); E.append(e); NS.append(r[t > t.max() / 4].mean()); curves.append((e, t, r))
o = np.argsort(E); E, NS = np.array(E)[o], np.array(NS)[o]
ax = axs[1, 0]
xs = np.logspace(-6.3, -1.8, 50); ax.plot(xs, 2 * np.sqrt(3 * P_ESC * xs / VREL), color=CORAL, lw=2.6, label=r"$n^\ast=\sqrt{6p\varepsilon/v}$")
ax.plot(E, NS, "o", ms=9, color=YELLOW, mec=NAVY, mew=1.4, zorder=3, label="Simulation")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"Noise $\varepsilon$"); ax.set_ylabel("Domain walls per block")
ax.set_xticks([1e-6, 1e-4, 1e-2]); ax.minorticks_off()
ax.legend(frameon=False, fontsize=13, loc="upper left", handlelength=1.5); lab(ax, "C")

ax = axs[1, 1]; cols = {1e-6: TEAL, 1e-5: YELLOW, 1e-4: CORAL}
xx = np.logspace(-1, 2, 300); ax.plot(xx, np.tanh(xx), color=NAVY, lw=2.2, ls=(0, (3, 2)), label="Rate equation")
for e, t, r in sorted(curves):
    if e not in cols: continue
    J = 3 * P_ESC * e; nst = 2 * np.sqrt(J / VREL); tau = 1 / np.sqrt(J * VREL)
    ax.plot(t / tau, r / nst, lw=2.6, color=cols[e], label=rf"$\varepsilon = 10^{{{int(round(np.log10(e)))}}}$")
ax.set_xscale("log"); ax.set_xlim(0.1, 100); ax.set_ylim(0, 1.3); ax.minorticks_off()
ax.set_xlabel(r"Time $\times\sqrt{6p\varepsilon v}$"); ax.set_ylabel(r"Domain walls / $n^\ast$")
ax.legend(frameon=False, fontsize=13, loc="lower right", handlelength=1.2, borderaxespad=0.1, labelspacing=0.1); ax.set_ylim(-0.12, 1.25); ax.set_yticks([0, 0.5, 1]); lab(ax, "D")
fig.subplots_adjust(left=0.18, right=0.965, top=0.95, bottom=0.125)
(HERE / "plots").mkdir(exist_ok=True)
fig.savefig(HERE / "plots" / "fig4_random_flips.pdf")
print("escape probabilities:\n", np.round(G, 3), "\naverage %.3f" % G.mean())
