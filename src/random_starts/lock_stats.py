"""Numbers quoted in Sec. 3: fraction of random starts that lock, median lock times, and the scaling exponents.
Reads outputs/lock_times and outputs/lock_charge."""
import numpy as np, glob, re, pathlib
OUT = pathlib.Path(__file__).resolve().parent / "outputs"
def load(n):
    lt = [np.loadtxt(f, ndmin=1) for f in glob.glob(str(OUT / "lock_times" / f"lock_N{n}.txt"))]
    if lt: return np.concatenate(lt)
    return np.concatenate([np.loadtxt(f, ndmin=2)[:, 1] for f in glob.glob(str(OUT / "lock_charge" / f"charge_N{n}_tmax*_s*.txt"))])
Ns = [32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 8192, 16384]
print(f"{'N':>6} {'runs':>5} {'locked':>7} {'median lock time':>17}")
med = []
for n in Ns:
    v = load(n); med.append(np.median(np.where(v > 0, v, np.inf)))
    print(f"{n:6d} {len(v):5d} {np.mean(v > 0):7.1%} {med[-1]:17.0f}")
Ns, med = np.array(Ns), np.array(med); ok = Ns >= 512
print("all runs, N >= 512: median ~ N^%.2f" % np.polyfit(np.log(Ns[ok]), np.log(med[ok]), 1)[0])
cls = {}
for f in glob.glob(str(OUT / "lock_charge" / "charge_N*_tmax*_s*.txt")):
    n = int(re.search(r"charge_N(\d+)", f).group(1)); cls.setdefault(n, []).append(np.loadtxt(f, ndmin=2))
print(f"\n{'N':>6} {'balanced':>9} {'excess left':>12} {'excess right':>13}   (runs out of 512; charge measured at t = 3N)")
for n in sorted(cls):
    q = np.vstack(cls[n])[:, 0]
    print(f"{n:6d} {np.sum(q == 0):9d} {np.sum(q < 0):12d} {np.sum((q > 0) & (q != 99)):13d}")
for name, sel, lo, hi in (("balanced", lambda q: q == 0, 128, 1024), ("excess left-moving", lambda q: q < 0, 512, 16384)):
    nn, mm = [], []
    for n in sorted(cls):
        d = np.vstack(cls[n]); v = d[sel(d[:, 0]), 1]
        if lo <= n <= hi and len(v) >= 10: nn.append(n); mm.append(np.median(np.where(v > 0, v, np.inf)))
    print("%s, N = %d-%d: median ~ N^%.2f" % (name, lo, hi, np.polyfit(np.log(nn), np.log(mm), 1)[0]))
