"""Description length of rules (Discussion): a rule is named by its family (one of 15, equal labels: log2 15 bits)
plus its member within the family (log2 of the family size). Families: which cells a rule reads, and which of those
cells it treats as interchangeable. A rule is charged for its cheapest family."""
import itertools, math
def families():
    fams = set()
    for size in range(4):
        for S in itertools.combinations(range(3), size):                 # cells read (0 = left, 1 = centre, 2 = right)
            perms = list(itertools.permutations(range(len(S))))
            groups = set()
            for n in range(len(perms) + 1):
                for gens in itertools.combinations(perms, n):            # every subgroup of permutations of those cells
                    G = {tuple(range(len(S)))} | set(gens)
                    while True:
                        new = {tuple(g[h[i]] for i in range(len(h))) for g in G for h in G} | G
                        if new == G: break
                        G = new
                    groups.add(frozenset(G))
            for G in groups:
                orbit = {}
                for k in range(8):
                    b = ((k >> 2) & 1, (k >> 1) & 1, k & 1); key = tuple(b[i] for i in S)
                    orbit[k] = min(tuple(key[g[i]] for i in range(len(S))) for g in G) if S else ()
                reps = sorted(set(orbit.values()))
                fams.add(frozenset(sum(dict(zip(reps, o))[orbit[k]] << k for k in range(8)) for o in itertools.product([0, 1], repeat=len(reps))))
    return list(fams)
F = families(); LABEL = math.log2(len(F))
cost = {r: LABEL + min(math.log2(len(m)) for m in F if r in m) for r in range(256)}
gamma = lambda n: 2 * int(math.log2(n)) + 1                             # Elias gamma code length of an integer
print(f"{len(F)} families; label {LABEL:.2f} bits")
for r in (5, 36, 37, 110): print(f"rule {r:3d}: {cost[r]:.1f} bits")
rules = cost[36] + cost[37] + cost[5]; overhead = gamma(3) + 2 * 3 + gamma(3)      # period, code (2 blocks of 3 cells), steps
print(f"mixture 36|37|5: rules {rules:.1f} + period, code and steps {overhead} = {rules + overhead:.1f} bits (rule 110: {cost[110]:.1f})")
print(f"most expensive elementary rule: {max(cost.values()):.1f} bits; cheapest possible two-rule mixture: {2 * min(cost.values()) + gamma(2) + 4 + 1:.1f} bits")
