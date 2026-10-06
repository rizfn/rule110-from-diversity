"""Check that three steps of the 36|37|5 mixture equal one step of rule 110, for every configuration.

A block's state after 3 steps depends only on its own block and the two neighbouring blocks, so checking a periodic
ring that contains every logical word of length 7 (a de Bruijn sequence) covers every possible neighbourhood."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
import numpy as np
from mixture import mixture_run, encode, decode, is_codeword, eca_step, de_bruijn

def check(rules):
    x = de_bruijn(2, 7)
    for _ in range(5):                                   # several logical steps in a row
        y = mixture_run(encode(x), 3, rules)
        if not is_codeword(y).all() or not np.array_equal(decode(y), eca_step(x, 110)): return False
        x = eca_step(x, 110)
    return True

if __name__ == "__main__":
    print("rules 36, 37, 5 reproduce rule 110:", check((36, 37, 5)))
    print("rules 54, 37, 5 reproduce rule 110:", check((54, 37, 5)))
