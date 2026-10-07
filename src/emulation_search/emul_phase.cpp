// emul_phase.cpp -- block emulation with PHASE-DEPENDENT codes.
//
// Rule arrangement: period 6, site n follows rules[n mod 6].
// Block length k with k NOT a multiple of 6. Blocks then sit at different
// phases of the arrangement, so we allow m = lcm(k,6)/k code TYPES:
//   logical cell j uses codewords (u^{(t)}_0, u^{(t)}_1),  t = j mod m,
//   u^{(t)}_0 != u^{(t)}_1.
// The supercell has L = m*k = lcm(k,6) cells and carries m logical bits.
//   k = 1 -> m = 6 ;  k = 2 -> m = 3 ;  k = 3 -> m = 2 ;  k = 4 -> m = 3
// Emulation:  F^T(enc(x)) = sigma^d(enc(f(x)))  for all x, with d a multiple
// of L (so that both the rule pattern and the code pattern are preserved).
//
// Pipeline:
//  1. uniform logical states      (ring of L cells, transition table)
//  2. logical rings (01), (001), (011) repeated to a multiple of m logical
//     cells: every logical position must decode, and must give one
//     consistent induced table f                      (bit-parallel rings)
//  3. finite windows for EVERY output type t = 0..m-1: logical-0 background,
//     256 random windows, then ALL 2^W windows (exact proof).
//
// usage: emul_phase mode kmin kmax Tmax Wmax out.csv
//   mode 0: targets {110,124,137,193}; mode 1: all chiral targets
// compile with -DPOSCTRL to read arrangements from posctrl6.txt instead.

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <array>
#include <set>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned __int128 u128;
typedef uint64_t u64;

int mirror_rule(int r) {
    int m = 0;
    for (int i = 0; i < 8; i++) {
        int l = (i >> 2) & 1, c = (i >> 1) & 1, rr = i & 1;
        if ((r >> i) & 1) m |= 1 << (4 * rr + 2 * c + l);
    }
    return m;
}
static inline u64 mask64(int n) { return n >= 64 ? ~0ull : ((1ull << n) - 1); }
static inline u64 rot(u64 w, int d, int n) {
    d = ((d % n) + n) % n;
    if (d == 0) return w;
    return ((w << d) | (w >> (n - d))) & mask64(n);
}
static inline u64 apply64(int r, u64 L, u64 C, u64 R) {
    u64 o = 0;
    for (int idx = 0; idx < 8; idx++)
        if ((r >> idx) & 1) o |= ((idx & 4) ? L : ~L) & ((idx & 2) ? C : ~C) & ((idx & 1) ? R : ~R);
    return o;
}
u64 RPM[65][6];  // ring phase masks for ring size n: bits i with i%6==p
static inline u64 ring_sim(u64 x, int n, int T, const int* rules) {
    u64 m = mask64(n);
    for (int t = 0; t < T; t++) {
        u64 L = rot(x, 1, n), R = rot(x, -1, n), out = 0;
        for (int p = 0; p < 6; p++) out |= apply64(rules[p], L, x, R) & RPM[n][p];
        x = out & m;
    }
    return x;
}

u128 PM6[6];
static inline u128 lowmask(int N) { return N >= 128 ? ~(u128)0 : (((u128)1 << N) - 1); }
static inline u128 apply128(int r, u128 L, u128 C, u128 R) {
    u128 o = 0;
    for (int idx = 0; idx < 8; idx++)
        if ((r >> idx) & 1) o |= ((idx & 4) ? L : ~L) & ((idx & 2) ? C : ~C) & ((idx & 1) ? R : ~R);
    return o;
}
static inline int floordiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
static inline int pmod(int a, int b) { return ((a % b) + b) % b; }

struct Win { int j, lo, hi, W, N, base, off; u128 mask; };
Win make_win(int k, int T, int d, int j) {
    Win w; w.j = j;
    w.lo = min(floordiv(j * k + d - T, k), j - 1);
    w.hi = max(floordiv(j * k + d + k + T - 1, k), j + 1);
    w.W = w.hi - w.lo + 1; w.N = w.W * k; w.base = (j - w.lo) * k + d;
    w.off = pmod(w.lo * k, 6); w.mask = lowmask(w.N);
    return w;
}
// codes: code[t][b] = codeword for type t, bit b
struct Code { int m, k; u64 c[6][2]; };

static inline u64 eval_win(u64 x, const Win& w, const Code& C, int T, const int* rules, int* vis) {
    int k = C.k, m = C.m;
    int re[6]; for (int p = 0; p < 6; p++) re[p] = rules[(p + w.off) % 6];
    u128 s = 0;
    for (int b = 0; b < w.W; b++) s |= (u128)C.c[pmod(w.lo + b, m)][(x >> b) & 1] << (b * k);
    for (int t = 0; t < T; t++) {
        u128 L = s << 1, R = s >> 1;
        if (vis) {
            u128 reg = lowmask(w.N - t - 1) & ~lowmask(t + 1);
            for (int idx = 0; idx < 8; idx++) {
                u128 mt = ((idx & 4) ? L : ~L) & ((idx & 2) ? s : ~s) & ((idx & 1) ? R : ~R) & reg;
                if (!mt) continue;
                for (int p = 0; p < 6; p++) if (mt & PM6[p]) vis[(p + w.off) % 6] |= 1 << idx;
            }
        }
        u128 out = 0;
        for (int p = 0; p < 6; p++) out |= PM6[p] & apply128(re[p], L, s, R);
        s = out & w.mask;
    }
    return (u64)(s >> w.base) & mask64(k);
}
static inline int nb_index(u64 x, const Win& w) {
    int c0 = w.j - w.lo;
    return 4 * ((x >> (c0 - 1)) & 1) + 2 * ((x >> c0) & 1) + ((x >> (c0 + 1)) & 1);
}

struct Hit { int r[6], k, T, d, target, vis[6]; Code C; };

int main(int argc, char** argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s mode kmin kmax Tmax Wmax out.csv\n", argv[0]); return 1; }
    int mode = atoi(argv[1]), kmin = atoi(argv[2]), kmax = atoi(argv[3]), Tmax = atoi(argv[4]), Wmax = atoi(argv[5]);
    const char* outname = argv[6];
    for (int p = 0; p < 6; p++) { PM6[p] = 0; for (int i = 0; i < 128; i++) if (i % 6 == p) PM6[p] |= (u128)1 << i; }
    for (int n = 1; n <= 64; n++) for (int p = 0; p < 6; p++) { RPM[n][p] = 0; for (int i = 0; i < n; i++) if (i % 6 == p) RPM[n][p] |= 1ull << i; }

    vector<int> sym; for (int r = 0; r < 256; r++) if (mirror_rule(r) == r) sym.push_back(r);
    vector<array<int, 6>> arrs;
#ifdef POSCTRL
    { FILE* pf = fopen("posctrl6.txt", "r"); array<int, 6> R;
      while (fscanf(pf, "%d %d %d %d %d %d", &R[0], &R[1], &R[2], &R[3], &R[4], &R[5]) == 6) arrs.push_back(R);
      fclose(pf); }
#else
    { set<array<int, 6>> seen;
      auto rot6 = [](int m, int s) { return ((m >> s) | (m << (6 - s))) & 63; };
      for (int m = 0; m < 64; m++) {
          if (rot6(m, 2) == m || rot6(m, 3) == m) continue;
          for (int a : sym) for (int b : sym) {
              if (a == b) continue;
              array<int, 6> R{}; for (int p = 0; p < 6; p++) R[p] = ((m >> p) & 1) ? b : a;
              if (seen.insert(R).second) arrs.push_back(R);
          } } }
#endif
    size_t a_begin = argc > 7 ? atoll(argv[7]) : 0, a_end = argc > 8 ? min((size_t)atoll(argv[8]), arrs.size()) : arrs.size();
    fprintf(stderr, "[%zu,%zu) ", a_begin, a_end);
    fprintf(stderr, "arrangements=%zu mode=%d k=%d..%d Tmax=%d\n", arrs.size(), mode, kmin, kmax, Tmax);
    bool isT0[256] = {false}; isT0[110] = isT0[124] = isT0[137] = isT0[193] = true;
    bool chiral[256]; for (int r = 0; r < 256; r++) chiral[r] = mirror_rule(r) != r;

    vector<Hit> hits; long long ncand = 0, nring = 0, nverify = 0, nskip = 0;
    mt19937_64 rng(777);
    auto t0 = chrono::steady_clock::now();
    vector<uint32_t> G, cur;

    for (size_t a = a_begin; a < a_end; a++) {
        int rules[6]; for (int p = 0; p < 6; p++) rules[p] = arrs[a][p];
        vector<char> found(512, 0);
        for (int k = kmin; k <= kmax; k++) {
            if (k % 6 == 0) continue;               // covered by the plain search
            int Lc = k * 6 / gcd(k, 6), m = Lc / k;
            if (Lc > 12) continue;                  // transition table limit (2^12)
            int nS = 1 << Lc;
            // logical ring lengths: multiple of m and of the pattern period
            int lam2 = m * 2 / gcd(m, 2), lam3 = m * 3 / gcd(m, 3);
            if (lam2 * k > 64 || lam3 * k > 64) continue;
            G.resize(nS); cur.resize(nS);
            for (int s = 0; s < nS; s++) { G[s] = (uint32_t)ring_sim(s, Lc, 1, rules); cur[s] = s; }
            u64 km = mask64(k);
            auto blocks_differ = [&](u64 A, u64 B) {
                for (int t = 0; t < m; t++) if (((A >> (t * k)) & km) == ((B >> (t * k)) & km)) return false;
                return true;
            };
            for (int T = 1; T <= Tmax; T++) {
                for (int s = 0; s < nS; s++) cur[s] = G[cur[s]];
                // d multiple of Lc => rotation of ring Lc by d is the identity: img = cur
                struct Cand { u64 S0, S1; int a0, a1; };
                vector<Cand> cands; vector<uint32_t> fixedv;
                for (int w = 0; w < nS; w++) {
                    uint32_t v = cur[w];
                    if (v == (uint32_t)w) { fixedv.push_back(w); continue; }
                    if (!blocks_differ(w, v)) continue;
                    if (cur[v] == v) { cands.push_back({v, (u64)w, 0, 0}); cands.push_back({(u64)w, v, 1, 1}); }
                    else if (mode == 1 && cur[v] == (uint32_t)w && (uint32_t)w < v) {
                        cands.push_back({(u64)w, v, 1, 0}); cands.push_back({v, (u64)w, 1, 0}); }
                }
                if (mode == 1)
                    for (uint32_t x : fixedv) for (uint32_t y : fixedv) if (x != y && blocks_differ(x, y)) cands.push_back({x, y, 0, 1});
                ncand += cands.size();

                for (auto& c : cands) {
                    Code C; C.m = m; C.k = k;
                    for (int t = 0; t < m; t++) { C.c[t][0] = (c.S0 >> (t * k)) & km; C.c[t][1] = (c.S1 >> (t * k)) & km; }
                    // ring sims are independent of d up to rotation: compute once
                    const int pats[3][3] = {{0, 1, -1}, {0, 0, 1}, {0, 1, 1}};
                    const int plen[3] = {2, 3, 3};
                    u64 ry[3]; int rn[3];
                    for (int pi = 0; pi < 3; pi++) {
                        int lam = (pi == 0) ? lam2 : lam3; rn[pi] = lam * k;
                        u64 x = 0;
                        for (int j = 0; j < lam; j++) x |= C.c[j % m][pats[pi][j % plen[pi]]] << (j * k);
                        ry[pi] = ring_sim(x, rn[pi], T, rules);
                    }
                    for (int d = -((T + Lc) / Lc) * Lc; d <= T + Lc; d += Lc) {
                        // ---- logical rings: decode every position, build table ----
                        int tbl = c.a0 | (c.a1 << 7), known = 1 | (1 << 7); bool ok = true;
                        for (int pi = 0; pi < 3 && ok; pi++) {
                            int lam = rn[pi] / k;
                            u64 y = rot(ry[pi], -d, rn[pi]);
                            for (int j = 0; j < lam && ok; j++) {
                                u64 o = (y >> (j * k)) & km; int t = j % m;
                                int bit = (o == C.c[t][0]) ? 0 : (o == C.c[t][1]) ? 1 : -1;
                                if (bit < 0) { ok = false; break; }
                                int P = plen[pi];
                                int l = pats[pi][(j + P - 1) % P], cc = pats[pi][j % P], r = pats[pi][(j + 1) % P];
                                int idx = 4 * l + 2 * cc + r;
                                if ((known >> idx) & 1) { if (((tbl >> idx) & 1) != bit) ok = false; }
                                else { known |= 1 << idx; tbl |= bit << idx; }
                            }
                        }
                        if (!ok || known != 255) continue;
                        if (mode == 0 && !isT0[tbl]) continue;
                        if (mode == 1 && (!chiral[tbl] || found[tbl * 2 + (d == 0)])) continue;
                        nring++;
                        // ---- windows, for every output type ----
                        vector<Win> wins;
                        for (int t = 0; t < m; t++) wins.push_back(make_win(k, T, d, t));
                        bool skip = false;
                        for (auto& w : wins) if (w.W > Wmax || w.N > 128) skip = true;
                        if (skip) { nskip++; continue; }
                        for (auto& w : wins) {
                            int t = w.j;
                            for (int idx = 0; idx < 8 && ok; idx++) {
                                int c0 = w.j - w.lo;
                                u64 x = ((u64)((idx >> 2) & 1) << (c0 - 1)) | ((u64)((idx >> 1) & 1) << c0) | ((u64)(idx & 1) << (c0 + 1));
                                if (eval_win(x, w, C, T, rules, nullptr) != C.c[t][(tbl >> idx) & 1]) ok = false;
                            }
                            u64 full = mask64(w.W);
                            for (int rr = 0; rr < 256 && ok; rr++) {
                                u64 x = rng() & full;
                                if (eval_win(x, w, C, T, rules, nullptr) != C.c[t][(tbl >> nb_index(x, w)) & 1]) ok = false;
                            }
                            if (!ok) break;
                        }
                        if (!ok) continue;
                        nverify++;
                        int vis[6] = {0};
                        for (auto& w : wins) {
                            u64 full = mask64(w.W);
                            for (u64 x = 0;; x++) {
                                if (eval_win(x, w, C, T, rules, vis) != C.c[w.j][(tbl >> nb_index(x, w)) & 1]) { ok = false; break; }
                                if (x == full) break;
                            }
                            if (!ok) break;
                        }
                        if (!ok) continue;
                        Hit h; for (int p = 0; p < 6; p++) { h.r[p] = rules[p]; h.vis[p] = vis[p]; }
                        h.k = k; h.T = T; h.d = d; h.target = tbl; h.C = C;
                        hits.push_back(h); found[tbl * 2 + (d == 0)] = 1;
                    }
                }
            }
        }
        if ((a + 1) % 8192 == 0) {
            double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
            fprintf(stderr, "  %zu/%zu hits=%zu cands=%lld ring-pass=%lld %.0fs\n", a + 1, arrs.size(), hits.size(), ncand, nring, el);
        }
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "done: hits=%zu cands=%lld ring-pass=%lld verified=%lld skipped=%lld time=%.1fs\n",
            hits.size(), ncand, nring, nverify, nskip, el);
    FILE* f = fopen(outname, "w");
    fprintf(f, "r0,r1,r2,r3,r4,r5,k,m,T,d,target,codes,vis0,vis1,vis2,vis3,vis4,vis5\n");
    for (auto& h : hits) {
        fprintf(f, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,", h.r[0], h.r[1], h.r[2], h.r[3], h.r[4], h.r[5], h.k, h.C.m, h.T, h.d, h.target);
        for (int t = 0; t < h.C.m; t++) {   // codes written cell-by-cell, left to right: type:w0/w1
            for (int b = 0; b < 2; b++) {
                for (int i = 0; i < h.k; i++) fputc('0' + ((h.C.c[t][b] >> i) & 1), f);
                fputc(b == 0 ? '/' : (t + 1 < h.C.m ? '|' : ','), f);
            }
        }
        for (int p = 0; p < 6; p++) fprintf(f, "%d%c", h.vis[p], p < 5 ? ',' : '\n');
    }
    fclose(f);
    return 0;
}
