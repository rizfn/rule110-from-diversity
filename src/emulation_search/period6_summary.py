"""Summary of the period-6 two-rule search (emul_search2, mode 1, k = 6): chiral targets reached by achiral vs chiral
binary patterns. Reads outputs/emul_search2/hits_chiral_P6two-rule_k6_T1-12.csv.gz."""
import csv, gzip, pathlib
from collections import defaultdict
from analyze import mirror, compl_rule
rep = lambda r: min(r, mirror(r), compl_rule(r), mirror(compl_rule(r)))
HITS = pathlib.Path(__file__).resolve().parent / 'outputs' / 'emul_search2' / 'hits_chiral_P6two-rule_k6_T1-12.csv.gz'
C = {18,183,22,151,54,147,90,165,105,150,122,161,126,129,146,182}  # class-3/4 symmetric rules

def pattern(R):  # binary pattern of a two-rule arrangement, canonical up to rotation + relabel
    a = R[0]; b = [x for x in R if x != a][0]
    m = tuple(0 if x == a else 1 for x in R)
    vs = []
    for s in range(6):
        t = m[s:] + m[:s]
        vs += [t, tuple(1 - v for v in t)]
    return min(vs)

def is_chiral(R):
    rev = R[::-1]
    return all(tuple(rev[s:] + rev[:s]) != tuple(R) for s in range(6))

stats = defaultdict(lambda: defaultdict(set))
for r in csv.DictReader(gzip.open(HITS, 'rt')):
    R = [int(r[f'r{p}']) for p in range(6)]
    if r['d'] != '0':
        continue
    ch = 'chiral pattern' if is_chiral(R) else 'achiral pattern'
    t = rep(int(r['target']))
    stats[ch][t].add((tuple(sorted(set(R))), not (set(R) & C)))

for ch in ('achiral pattern', 'chiral pattern'):
    S = stats[ch]
    print(f'{ch}: {len(S)} chiral target classes (d=0): {sorted(S)}')
    for t in (30, 106, 110, 60, 45, 41):
        if t in S:
            triv = sorted({p for p, tv in S[t] if tv})
            print(f'   class {t}: {len({p for p,_ in S[t]})} rule pairs, {len(triv)} with only class-1/2 rules {triv[:5]}')
# patterns present
pats = {pattern(tuple(int(r[f"r{p}"]) for p in range(6))) for r in csv.DictReader(gzip.open(HITS, 'rt'))}
print('distinct binary patterns among hits:', len(pats))
