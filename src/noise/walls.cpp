// Walls created by single flips and by noise (64 independent runs per machine word, one per bit).
// compile: g++ -O3 -march=native -o walls walls.cpp
//  track NB WORDS T SEED      : one flipped non-bit cell per replica in a perfect run; print positions of the
//                               right- and left-moving walls every 100 steps (relative to the flip)
//  pairf NB WORDS SEED        : escape probability of a flipped random cell, for each of the 3 time phases
//  noisec NB WORDS T SEED EPS : code start + noise; number of walls (wall clusters) per block vs time
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <random>
#include <string>
#include <cmath>
using namespace std; typedef uint64_t u64;
struct Sys {
    int NB, N; vector<u64> x, y; mt19937_64 rng;
    Sys(int nb, u64 s) : NB(nb), N(3 * nb), x(3 * nb), y(3 * nb), rng(s) {}
    void step() {
        for (int j = 0; j < NB; j++) { int i = 3 * j;
            u64 a = x[(i + N - 1) % N], b = x[i], c = x[i + 1], d = x[i + 2], e = x[(i + 3) % N];
            y[i] = (~a & b & ~c) | (a & ~b & c); y[i + 1] = (~b & ~c & ~d) | (~b & c & ~d) | (b & ~c & d); y[i + 2] = ~c & ~e; }
        x.swap(y);
    }
    u64 cw(int j) const { return ~(x[3 * j + 1] | x[3 * j + 2]); }
    u64 allcw() const { u64 a = ~0ULL; for (int j = 0; j < NB; j++) a &= cw(j); return a; }
    void code() { for (int j = 0; j < NB; j++) { u64 b = rng(); x[3 * j] = ~b; x[3 * j + 1] = 0; x[3 * j + 2] = 0; } }
    void noise(double eps) { long total = (long)N * 64; geometric_distribution<long> G(eps);
        for (long k = G(rng); k < total; k += 1 + G(rng)) x[k / 64] ^= 1ULL << (k % 64); }
};
struct Walls { vector<u64> m0, m1, m2, w; int NB;
    Walls(int nb) : m0(nb), m1(nb), m2(nb), w(nb), NB(nb) {}
    void push(const Sys& S) { m2.swap(m1); m1.swap(m0); for (int j = 0; j < NB; j++) m0[j] = S.cw(j); }
    void compute() { for (int j = 0; j < NB; j++) { int k = (j + 1) % NB;
        w[j] = ~(m0[j] | m1[j] | m2[j]) | ~((m0[j] & m0[k]) | (m1[j] & m1[k]) | (m2[j] & m2[k])); } }
    double clusters() const { long s = 0; for (int j = 0; j < NB; j++) s += __builtin_popcountll(w[j] & ~w[(j + NB - 1) % NB]); return (double)s / NB / 64; }
};
int main(int argc, char** argv) {
    string mode = argv[1];
    if (mode == "track") {
        int NB = atoi(argv[2]), W = atoi(argv[3]); long T = atol(argv[4]); u64 seed = strtoull(argv[5], 0, 10);
        for (int w = 0; w < W; w++) {
            Sys S(NB, seed + 31 * w); S.code(); for (int t = 0; t < 300; t++) S.step();
            int c = NB / 2;
            for (int b = 0; b < 64; b++) S.x[3 * c + 1 + (S.rng() & 1)] ^= 1ULL << b;   // non-bit cell of the middle block
            Walls Wl(NB);
            for (long t = 1; t <= T; t++) {
                S.step(); Wl.push(S);
                if (t % 100 == 0) { Wl.compute();
                    for (int b = 0; b < 64; b++) {
                        int R = -1, L = 1;
                        for (int d = min(NB / 2 - 1, (int)(0.30 * t) + 400); d > 0; d--) if (Wl.w[(c + d) % NB] >> b & 1) { R = d; break; }
                        for (int d = min(NB / 2 - 1, (int)(0.30 * t) + 400); d > 0; d--) if (Wl.w[(c - d + NB) % NB] >> b & 1) { L = -d; break; }
                        printf("%d %ld %d %d\n", w * 64 + b, t, R, L);
                    } }
            }
        }
    } else if (mode == "pairf") {
        int NB = atoi(argv[2]), W = atoi(argv[3]); u64 seed = strtoull(argv[4], 0, 10);
        for (int ph = 0; ph < 3; ph++) { long esc = 0, tot = 0;
            for (int w = 0; w < W; w++) {
                Sys S(NB, seed + 97 * w + ph); S.code(); for (int t = 0; t < 300 + ph; t++) S.step();
                for (int b = 0; b < 64; b++) { long i = S.rng() % S.N; S.x[i] ^= 1ULL << b; }
                u64 healed = 0; for (int t = 0; t < 60; t++) { S.step(); healed |= S.allcw(); }
                esc += 64 - __builtin_popcountll(healed); tot += 64;
            }
            printf("phase %d escape %.4f (%ld flips)\n", ph, (double)esc / tot, tot);
        }
    } else if (mode == "noisec") {
        int NB = atoi(argv[2]), W = atoi(argv[3]); long T = atol(argv[4]); u64 seed = strtoull(argv[5], 0, 10); double eps = atof(argv[6]);
        vector<long> times; for (double t = 3; t <= T; t *= 1.08) { long v = (long)t; if (times.empty() || v != times.back()) times.push_back(v); }
        vector<double> acc(times.size(), 0);
        for (int w = 0; w < W; w++) {
            Sys S(NB, seed + 7 * w); S.code(); Walls Wl(NB); size_t ti = 0;
            for (long t = 1; t <= T + 2; t++) { S.step(); S.noise(eps); Wl.push(S);
                while (ti < times.size() && times[ti] + 2 == t) { Wl.compute(); acc[ti] += Wl.clusters() / W; ti++; } }
        }
        for (size_t i = 0; i < times.size(); i++) printf("%ld %.8f\n", times[i], acc[i]);
    }
}
