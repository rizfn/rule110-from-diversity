"""Measure the parameters of the rate equation (Sec. 4): the escape probability p and the wall speeds v_R, v_L.
Compiles and runs walls.cpp (in this folder); takes about a minute. Results go into params.py by hand.
  p     : fraction of single flips (random cell, flipped at each of the 3 steps of the cycle) still unhealed after 60 steps
  v_R,L : mean speed of the two walls created by one flipped cell, tracked for 60,000 steps on a 32,768-block ring"""
import subprocess, numpy as np, os
cdir = os.path.dirname(os.path.abspath(__file__))
subprocess.run(["g++", "-O3", "-march=native", "-o", "walls", "walls.cpp"], cwd=cdir, check=True)
out = subprocess.run(["./walls", "pairf", "2048", "16", "5"], cwd=cdir, capture_output=True, text=True, check=True).stdout
esc = [float(l.split()[3]) for l in out.splitlines()]
print(out.strip()); print("p = %.4f" % np.mean(esc))
out = subprocess.run(["./walls", "track", "32768", "4", "60000", "9"], cwd=cdir, capture_output=True, text=True, check=True).stdout
d = np.array([list(map(float, l.split())) for l in out.splitlines()])
rep, t, R, L = d[:, 0].astype(int), d[:, 1], d[:, 2], d[:, 3]
T = np.unique(t); ti = np.searchsorted(T, t); n = rep.max() + 1
Rm = np.full((n, len(T)), np.nan); Lm = Rm.copy()
Rm[rep, ti] = np.where(R > 0, R, np.nan); Lm[rep, ti] = np.where(L < 0, L, np.nan)
alive = (np.isfinite(Rm[:, -50:]).mean(1) > 0.5) & (np.isfinite(Lm[:, -50:]).mean(1) > 0.5)   # pairs that escaped
late = T > 5000
vR = np.polyfit(T[late], np.nanmean(Rm[alive], 0)[late], 1)[0]
vL = -np.polyfit(T[late], np.nanmean(Lm[alive], 0)[late], 1)[0]
print("tracked pairs: %d   v_R = %.4f   v_L = %.4f blocks per step" % (alive.sum(), vR, vL))
