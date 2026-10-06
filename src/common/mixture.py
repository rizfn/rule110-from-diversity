"""Core of the construction: the 36|37|5 mixture, the block code, and elementary rules.

Cell i obeys rule 36 if i % 3 == 0, rule 37 if i % 3 == 1, rule 5 if i % 3 == 2 (Wolfram numbering, periodic lattice).
Logical 0 is the block (1, 0, 0); logical 1 is (0, 0, 0). A block is a valid code word iff its cells 1 and 2 are 0.
Three steps of the mixture equal one step of rule 110 on the logical bits."""
import numpy as np

RULES = (36, 37, 5)
TABLE = {r: np.array([(r >> k) & 1 for k in range(8)], np.uint8) for r in range(256)}
CODE = {0: (1, 0, 0), 1: (0, 0, 0)}

def neighbourhood(x):
    """index 4*left + 2*centre + right along the last axis (periodic)"""
    return (np.roll(x, 1, -1) << 2) | (x << 1) | np.roll(x, -1, -1)

def eca_step(x, rule):
    return TABLE[rule][neighbourhood(x)]

def mixture_step(x, rules=RULES):
    """one step of a spatially periodic mixture; works on 1D lattices or batches (last axis = lattice)"""
    k = neighbourhood(x); y = np.empty_like(x); L = len(rules)
    for p, r in enumerate(rules): y[..., p::L] = TABLE[r][k[..., p::L]]
    return y

def mixture_run(x, steps, rules=RULES):
    for _ in range(steps): x = mixture_step(x, rules)
    return x

def encode(bits):
    return np.array([CODE[int(b)] for b in bits], np.uint8).reshape(-1)

def is_codeword(x):
    """True for every block that is a valid code word (cells 1 and 2 white); works on batches too"""
    return (x[..., 1::3] == 0) & (x[..., 2::3] == 0)

def decode(x):
    """logical bits of a lattice whose blocks are all code words"""
    return (x[..., 0::3] == 0).astype(np.uint8)

def code_state(n_blocks, rng):
    """a perfect rule-110 code state with random logical bits"""
    return encode(rng.integers(0, 2, n_blocks))

def de_bruijn(k, n):
    """cyclic sequence over k symbols containing every word of length n exactly once"""
    a = [0] * k * n; seq = []
    def db(t, p):
        if t > n:
            if n % p == 0: seq.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]; db(t + 1, p)
            for j in range(a[t - p] + 1, k): a[t] = j; db(t + 1, t)
    db(1, 1); return np.array(seq, np.uint8)
