// emul_search2.cpp  --  exhaustive block-emulation search, version 2
//
// Arrangements of mirror-symmetric ECA rules with spatial period P:
//   P = 1, 2, 3 : every arrangement (64^P)
//   P = 6       : every TWO-rule arrangement of EXACT period 6
//                 (54 binary patterns x ordered pairs a != b, deduplicated)
// Site n follows rules[n mod P].
//
// Emulation (block code, global frame shift d, d multiple of P):
//   logical bit x_j  ->  word w_{x_j} of k cells (k multiple of P)
//   F^T(enc(x)) = sigma^d(enc(f(x)))   for ALL logical configurations x.
//
// Pipeline (filters are necessary conditions, last step is an exact proof):
//   1. uniform states 0^inf, 1^inf           (ring k,   via transition table)
//   2. logical (01)^inf                      (ring 2k,  bit-parallel)
//   3. logical (001)^inf, (011)^inf          (ring 3k,  bit-parallel) -> full table
//   4. table in logical-0 background         (finite window)
//   5. 256 random windows
//   6. ALL 2^W windows that influence one output block (proof), recording
//      which neighbourhoods every phase visits.
//
// mode 0: targets {110,124,137,193} (rule 110 and its equivalents)
// mode 1: every chiral ECA target
//
// usage: emul_search2 P mode kmax Tmax Wmax out.csv [a_begin a_end]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <array>
#include <set>
#include <random>
#include <chrono>
#include <algorithm>
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

// ---------------- rings (n <= 64 cells) ----------------
static inline u64 mask64(int n) { return n >= 64 ? ~0ull : ((1ull << n) - 1); }
static inline u64 rot(u64 w, int d, int n) {  // (rot_d w)_i = w_{i-d}
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
static inline u64 ring_sim(u64 x, int n, int T, const int* rules, int P, const u64* rpm) {
    u64 m = mask64(n);
    for (int t = 0; t < T; t++) {
        u64 L = rot(x, 1, n), R = rot(x, -1, n), out = 0;
        for (int p = 0; p < P; p++) out |= apply64(rules[p], L, x, R) & rpm[p];
        x = out & m;
    }
    return x;
}

// ---------------- finite windows (<= 128 cells) ----------------
u128 PM[7][6];
static inline u128 lowmask(int N) { return N >= 128 ? ~(u128)0 : (((u128)1 << N) - 1); }
static inline u128 apply128(int r, u128 L, u128 C, u128 R) {
    u128 o = 0;
    for (int idx = 0; idx < 8; idx++)
        if ((r >> idx) & 1) o |= ((idx & 4) ? L : ~L) & ((idx & 2) ? C : ~C) & ((idx & 1) ? R : ~R);
    return o;
}
static inline u128 step_win(u128 x, const int* rules, int P, u128 mask) {
    u128 L = x << 1, R = x >> 1, out = 0;
    for (int p = 0; p < P; p++) out |= PM[P][p] & apply128(rules[p], L, x, R);
    return out & mask;
}
static inline int floordiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
struct Win { int k, T, d, lo, hi, W, N, base; u128 mask; };
Win make_win(int k, int T, int d) {
    Win w; w.k = k; w.T = T; w.d = d;
    w.lo = min(floordiv(d - T, k), -1);
    w.hi = max(floordiv(d + k + T - 1, k), 1);
    w.W = w.hi - w.lo + 1; w.N = w.W * k; w.base = (-w.lo) * k + d; w.mask = lowmask(w.N);
    return w;
}
static inline u128 build(u64 x, const Win& w, u64 w0, u64 w1) {
    u128 s = 0;
    for (int b = 0; b < w.W; b++) s |= (u128)(((x >> b) & 1) ? w1 : w0) << (b * w.k);
    return s;
}
static inline u64 eval_win(u64 x, const Win& w, u64 w0, u64 w1, const int* rules, int P, int* vis) {
    u128 s = build(x, w, w0, w1);
    for (int t = 0; t < w.T; t++) {
        if (vis) {
            u128 reg = lowmask(w.N - t - 1) & ~lowmask(t + 1);
            u128 L = s << 1, R = s >> 1;
            for (int idx = 0; idx < 8; idx++) {
                u128 mt = ((idx & 4) ? L : ~L) & ((idx & 2) ? s : ~s) & ((idx & 1) ? R : ~R) & reg;
                if (!mt) continue;
                for (int p = 0; p < P; p++) if (mt & PM[P][p]) vis[p] |= 1 << idx;
            }
        }
        s = step_win(s, rules, P, w.mask);
    }
    return (u64)(s >> w.base) & mask64(w.k);
}
static inline u64 ctx_bits(const Win& w, int idx) {
    int c0 = -w.lo;
    return ((u64)((idx >> 2) & 1) << (c0 - 1)) | ((u64)((idx >> 1) & 1) << c0) | ((u64)(idx & 1) << (c0 + 1));
}
static inline int nb_index(u64 x, const Win& w) {
    int c0 = -w.lo;
    return 4 * ((x >> (c0 - 1)) & 1) + 2 * ((x >> c0) & 1) + ((x >> (c0 + 1)) & 1);
}

struct Hit { int P, r[6], k, T, d; u64 w0, w1; int target, vis[6]; };

int main(int argc, char** argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s P mode kmax Tmax Wmax out.csv [a_begin a_end]\n", argv[0]); return 1; }
    int P = atoi(argv[1]), mode = atoi(argv[2]), kmax = atoi(argv[3]), Tmax = atoi(argv[4]), Wmax = atoi(argv[5]);
    const char* outname = argv[6];

    for (int q = 1; q <= 6; q++) for (int p = 0; p < q; p++) {
        PM[q][p] = 0;
        for (int i = 0; i < 128; i++) if (i % q == p) PM[q][p] |= (u128)1 << i;
    }
    vector<int> sym;
    for (int r = 0; r < 256; r++) if (mirror_rule(r) == r) sym.push_back(r);
    int nsym = sym.size();

    // ---------------- arrangements ----------------
    vector<array<int, 6>> arrs;
    if (P >= 1 && P <= 3) {
        long long na = 1; for (int i = 0; i < P; i++) na *= nsym;
        for (long long a = 0; a < na; a++) {
            array<int, 6> R{}; long long t = a;
            for (int p = 0; p < P; p++) { R[p] = sym[t % nsym]; t /= nsym; }
            arrs.push_back(R);
        }
    } else if (P == 6) {
        set<array<int, 6>> seen;
        auto rot6 = [](int m, int s) { return ((m >> s) | (m << (6 - s))) & 63; };
        for (int m = 0; m < 64; m++) {
            if (rot6(m, 2) == m || rot6(m, 3) == m) continue;   // smaller period
            for (int ia = 0; ia < nsym; ia++) for (int ib = 0; ib < nsym; ib++) {
                if (ia == ib) continue;
                array<int, 6> R{};
                for (int p = 0; p < 6; p++) R[p] = ((m >> p) & 1) ? sym[ib] : sym[ia];
                if (seen.insert(R).second) arrs.push_back(R);
            }
        }
    } else { fprintf(stderr, "P must be 1, 2, 3 or 6\n"); return 1; }
#ifdef POSCTRL
    // positive control: replace list by period-6 lifts read from posctrl.txt
    if (P == 6) {
        arrs.clear();
        FILE* pf = fopen("posctrl.txt", "r"); array<int, 6> R;
        while (fscanf(pf, "%d %d %d %d %d %d", &R[0], &R[1], &R[2], &R[3], &R[4], &R[5]) == 6) arrs.push_back(R);
        fclose(pf);
    }
#endif
    long long a_begin = argc > 7 ? atoll(argv[7]) : 0, a_end = argc > 8 ? atoll(argv[8]) : (long long)arrs.size();
    a_end = min(a_end, (long long)arrs.size());
    fprintf(stderr, "P=%d mode=%d arrangements=%zu (running %lld..%lld)\n", P, mode, arrs.size(), a_begin, a_end);

    bool isT0[256] = {false}; isT0[110] = isT0[124] = isT0[137] = isT0[193] = true;
    bool chiral[256]; for (int r = 0; r < 256; r++) chiral[r] = mirror_rule(r) != r;

    vector<Hit> hits;
    long long nskip = 0, nverify = 0, ncand = 0;
    mt19937_64 rng(12345);
    auto t0 = chrono::steady_clock::now();
    vector<uint32_t> G, cur, img;

    for (long long a = a_begin; a < a_end; a++) {
        int rules[6]; for (int p = 0; p < 6; p++) rules[p] = arrs[a][p];
        vector<char> found(512, 0);
        for (int k = P; k <= kmax; k += P) {
            int n1 = 1 << k;
            u64 rpm2[6] = {0}, rpm3[6] = {0};
            for (int i = 0; i < 2 * k && i < 64; i++) rpm2[i % P] |= 1ull << i;
            for (int i = 0; i < 3 * k && i < 64; i++) rpm3[i % P] |= 1ull << i;
            u64 rpm1[6] = {0};
            for (int i = 0; i < k; i++) rpm1[i % P] |= 1ull << i;
            G.resize(n1); cur.resize(n1); img.resize(n1);
            for (int s = 0; s < n1; s++) { G[s] = (uint32_t)ring_sim(s, k, 1, rules, P, rpm1); cur[s] = s; }

            for (int T = 1; T <= Tmax; T++) {
                for (int s = 0; s < n1; s++) cur[s] = G[cur[s]];
                for (int sh = 0; sh < k; sh += P) {
                    for (int w = 0; w < n1; w++) img[w] = (uint32_t)rot(cur[w], -sh, k);
                    // ---- candidate pairs (a0,a1) = (f(000), f(111)) ----
                    struct Cand { u64 w0, w1; int a0, a1; };
                    vector<Cand> cands;
                    vector<uint32_t> fixedv;
                    for (int w = 0; w < n1; w++) {
                        uint32_t v = img[w];
                        if (v == (uint32_t)w) { fixedv.push_back(w); continue; }
                        if (img[v] == v) {                 // w -> v -> v
                            cands.push_back({v, (u64)w, 0, 0});   // 0->0, 1->0
                            cands.push_back({(u64)w, v, 1, 1});   // 0->1, 1->1
                        } else if (mode == 1 && img[v] == (uint32_t)w && (uint32_t)w < v) {  // 2-cycle
                            cands.push_back({(u64)w, v, 1, 0});
                            cands.push_back({v, (u64)w, 1, 0});
                        }
                    }
                    if (mode == 1)
                        for (uint32_t x : fixedv) for (uint32_t y : fixedv) if (x != y) cands.push_back({x, y, 0, 1});
                    ncand += cands.size();

                    for (auto& c : cands) {
                        // ring 2k: logical (01)^inf -> (f(101), f(010))
                        u64 y2 = ring_sim(c.w0 | (c.w1 << k), 2 * k, T, rules, P, rpm2);
                        bool have3 = false; u64 y3[2] = {0, 0};
                        for (int d = sh - k * ((T + k) / k + 1); d <= T + k; d += k) {
                            if (d < -(T + k)) continue;
                            u64 dec = rot(y2, -d, 2 * k);
                            u64 lo_ = dec & mask64(k), hi_ = dec >> k;
                            if (mode == 0) {
                                u64 e = (c.a0 == 0) ? c.w1 : c.w0;
                                if (lo_ != e || hi_ != e) continue;
                            } else if ((lo_ != c.w0 && lo_ != c.w1) || (hi_ != c.w0 && hi_ != c.w1)) continue;
                            // ring 3k: (001)^inf and (011)^inf -> complete induced table
                            if (!have3) {
                                const int pats[2][3] = {{0, 0, 1}, {0, 1, 1}};
                                for (int pi = 0; pi < 2; pi++) {
                                    u64 x = 0;
                                    for (int b = 0; b < 3; b++) x |= (pats[pi][b] ? c.w1 : c.w0) << (b * k);
                                    y3[pi] = ring_sim(x, 3 * k, T, rules, P, rpm3);
                                }
                                have3 = true;
                            }
                            int tbl = c.a0 | (c.a1 << 7); bool ok = true;
                            const int pats[2][3] = {{0, 0, 1}, {0, 1, 1}};
                            for (int pi = 0; pi < 2 && ok; pi++) {
                                u64 y = rot(y3[pi], -d, 3 * k);
                                for (int b = 0; b < 3; b++) {
                                    u64 o = (y >> (b * k)) & mask64(k);
                                    int bit = (o == c.w0) ? 0 : (o == c.w1) ? 1 : -1;
                                    if (bit < 0) { ok = false; break; }
                                    int l = pats[pi][(b + 2) % 3], cc = pats[pi][b], r = pats[pi][(b + 1) % 3];
                                    tbl |= bit << (4 * l + 2 * cc + r);
                                }
                            }
                            if (!ok) continue;
                            if (mode == 0 && !isT0[tbl]) continue;
                            if (mode == 1 && (!chiral[tbl] || found[tbl * 2 + (d == 0)])) continue;

                            Win win = make_win(k, T, d);
                            if (win.W > Wmax || win.N > 128) { nskip++; continue; }
                            // background-0 windows
                            for (int idx = 0; idx < 8 && ok; idx++) {
                                u64 o = eval_win(ctx_bits(win, idx), win, c.w0, c.w1, rules, P, nullptr);
                                if (o != (((tbl >> idx) & 1) ? c.w1 : c.w0)) ok = false;
                            }
                            if (!ok) continue;
                            nverify++;
                            u64 full = mask64(win.W);
                            for (int rr = 0; rr < 256 && ok; rr++) {
                                u64 x = rng() & full;
                                if (eval_win(x, win, c.w0, c.w1, rules, P, nullptr) != (((tbl >> nb_index(x, win)) & 1) ? c.w1 : c.w0)) ok = false;
                            }
                            if (!ok) continue;
                            int vis[6] = {0};
                            for (u64 x = 0;; x++) {
                                if (eval_win(x, win, c.w0, c.w1, rules, P, vis) != (((tbl >> nb_index(x, win)) & 1) ? c.w1 : c.w0)) { ok = false; break; }
                                if (x == full) break;
                            }
                            if (!ok) continue;
                            Hit h; h.P = P; h.k = k; h.T = T; h.d = d; h.w0 = c.w0; h.w1 = c.w1; h.target = tbl;
                            for (int p = 0; p < 6; p++) { h.r[p] = rules[p]; h.vis[p] = vis[p]; }
                            hits.push_back(h);
                            found[tbl * 2 + (d == 0)] = 1;
                        }
                    }
                }
            }
        }
        if ((a + 1 - a_begin) % 2048 == 0) {
            double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
            fprintf(stderr, "  %lld/%lld  hits=%zu  cands=%lld  verify=%lld  %.0fs\n", a + 1 - a_begin, a_end - a_begin, hits.size(), ncand, nverify, el);
        }
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "done: hits=%zu cands=%lld verify=%lld skipped(W>Wmax)=%lld time=%.1fs\n", hits.size(), ncand, nverify, nskip, el);

    FILE* f = fopen(outname, "w");
    int Pc = max(P, 3);
    fprintf(f, "P"); for (int p = 0; p < Pc; p++) fprintf(f, ",r%d", p);
    fprintf(f, ",k,T,d,w0,w1,target"); for (int p = 0; p < Pc; p++) fprintf(f, ",vis%d", p); fprintf(f, "\n");
    for (auto& h : hits) {
        fprintf(f, "%d", h.P); for (int p = 0; p < Pc; p++) fprintf(f, ",%d", p < P ? h.r[p] : 0);
        fprintf(f, ",%d,%d,%d,%llu,%llu,%d", h.k, h.T, h.d, (unsigned long long)h.w0, (unsigned long long)h.w1, h.target);
        for (int p = 0; p < Pc; p++) fprintf(f, ",%d", p < P ? h.vis[p] : 0); fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}
