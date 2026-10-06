"""Beat (time phase of the code) of every block, used for the colour maps of Figs. 3 and 4.

A block's candidate beats are the phases (time mod 3) at which it was a valid code word at every visit in the last
9 steps. One candidate -> that beat. Several -> the candidate its neighbours have. None -> an error (-1, drawn dark)."""
import numpy as np
from mixture import mixture_step, is_codeword

def resolve(full):
    nb = full.shape[1]; ncand = full.sum(0)
    ph = np.where(ncand == 1, np.argmax(full, 0), -1); amb = ncand >= 2
    for _ in range(nb):
        if not (amb & (ph < 0)).any(): break
        changed = False
        for nbr in (np.roll(ph, 1), np.roll(ph, -1)):
            ok = amb & (ph < 0) & (nbr >= 0)
            if ok.any():
                good = ok & full[np.clip(nbr, 0, 2), np.arange(nb)]
                if good.any(): ph = np.where(good, nbr, ph); changed = True
        if not changed: break
    return np.where(amb & (ph < 0), np.argmax(full, 0), ph)

def beat_rows(x, steps, every=1, window=9, noise=0.0, rng=None):
    """evolve x for `steps` steps (optionally flipping each bit with probability `noise` per step) and return the beat
    of every block, one row every `every` steps once the window is full"""
    rows, buf = [], []
    for s in range(steps):
        x = mixture_step(x)
        if noise: x = x ^ (rng.random(x.shape) < noise).astype(np.uint8)
        buf.append((is_codeword(x), (s + 1) % 3))
        if len(buf) > window: buf.pop(0)
        if len(buf) == window and s % every == 0:
            cnt = np.zeros((3, len(buf[0][0])), int)
            for cw, r in buf: cnt[r] += cw
            rows.append(resolve(cnt == window // 3))
    return x, np.array(rows)

def make_ic(n_blocks, widths, offsets, seed=0, burn=600):
    """hand-made start (Fig. 3B): a perfect rule-110 code lattice cut into domains advanced by 0, 1 or 2 extra steps,
    so that neighbouring domains run the code at different beats"""
    rng = np.random.default_rng(seed)
    x = np.zeros(3 * n_blocks, np.uint8); x[0::3] = 1 - rng.integers(0, 2, n_blocks)
    for _ in range(burn): x = mixture_step(x)
    states = [x]
    for _ in range(2): states.append(mixture_step(states[-1]))
    out = np.empty_like(x); s = 0
    for w, d in zip(widths, offsets):
        out[3 * s:3 * (s + w)] = states[d][3 * s:3 * (s + w)]; s += w
    assert s == n_blocks
    return out
