"""Chain one-step relations F_R^T(enc_c x) = sigma^d enc_c(f x) into infinite emulations.
Node (R, code, target); edge to (rot_d R, code, target) with rot_d R[p] = R[(p+d) mod 6].
Keep nodes that have an infinite forward path (prune dead ends).
usage: python chains.py outputs/emul_phase_glide/onestep_rule110_P3_phasecodes_glide_k1-3_T1-12.csv"""
import csv
from collections import defaultdict, Counter
from analyze import canon_raw

def load(fn):
    edges = defaultdict(list)
    for r in csv.DictReader(open(fn)):
        R = tuple(int(r[f'r{p}']) for p in range(6)); d = int(r['d'])
        R2 = tuple(R[(p + d) % 6] for p in range(6))
        node = (R, r['codes'], int(r['target']))
        edges[node].append(((R2, r['codes'], int(r['target'])), int(r['T']), d))
    return edges

def alive(edges):
    live = set(edges)
    changed = True
    while changed:
        changed = False
        for n in list(live):
            if not any(e[0] in live for e in edges[n]):
                live.discard(n); changed = True
    return live

if __name__ == '__main__':
    import sys
    E = load(sys.argv[1]); L = alive(E)
    print(f'{len(E)} nodes with a one-step relation, {len(L)} on infinite chains')
    trip = Counter(canon_raw(n[0][:3]) for n in L)
    print('triples on infinite chains:', dict(trip))
    # example cycle per triple
    shown = set()
    for n in sorted(L):
        t = canon_raw(n[0][:3])
        if t in shown: continue
        shown.add(t)
        path, cur, seen = [], n, {}
        while cur not in seen:
            seen[cur] = len(path)
            e = next(e for e in E[cur] if e[0] in L)
            path.append((cur[0][:3], e[1], e[2])); cur = e[0]
        cyc = path[seen[cur]:]
        print(f'  {t}: code {n[1]} target {n[2]}; cycle (arrangement, T, d): {cyc}')
