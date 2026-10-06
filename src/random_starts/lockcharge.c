/* Lock time AND net wall charge from random starts of the 36|37|5 mixture.
   64 independent runs at once (one bit each). Output: one line per run:  charge  lock_time
     charge = net phase charge of the walls at time T1 (+1 per jump-1 wall, -1 per jump-2 wall);
              0 = balanced, < 0 = excess left-moving walls, > 0 = excess right-moving walls,
              99 = no clean domains at T1 (rare)
     lock_time = first step at which every block is a code word (-1 = not locked by TMAX)
   Runs with charge > 0 never lock, so the job stops once all other runs have locked (their lock_time is then -1).
   Checkpoints every ~10 min; re-run the same command to resume.
   compile:  gcc -O3 -march=native -o lockcharge lockcharge.c
   run:      ./lockcharge NBLOCKS TMAX SEED [T1]  > outputs/lock_charge/charge_N<NBLOCKS>_tmax<TMAX>_s<SEED>.txt    (T1 defaults to 3*NBLOCKS) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
typedef uint64_t u64;
static u64 s[4];
static u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
static u64 rnd(void) { u64 r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
static long NB, N;
static void step(u64 *x, u64 *y) {
    for (long j = 0; j < NB; j++) {
        long i = 3 * j; u64 a = x[(i + N - 1) % N], b = x[i], c = x[i + 1], d = x[i + 2], e = x[(i + 3) % N];
        y[i] = (~a & b & ~c) | (a & ~b & c);                     /* rule 36 */
        y[i + 1] = (~b & ~c & ~d) | (~b & c & ~d) | (b & ~c & d); /* rule 37 */
        y[i + 2] = ~c & ~e;                                       /* rule 5  */
    }
}
/* net charge of replica bit b from code-word words at three consecutive times */
static int charge(u64 **cw, int b, int *ph) {
    for (long j = 0; j < NB; j++) {
        int m = 0, cnt = 0;
        for (int k = 0; k < 3; k++) if ((cw[k][j] >> b) & 1) { m = k; cnt++; }
        ph[j] = (cnt == 1) ? m : -1;
    }
    long start = -1;
    for (long j = 0; j < NB; j++) if (ph[j] >= 0 && ph[j] == ph[(j + 1) % NB]) { start = j; break; }
    if (start < 0) return 99;
    int cur = ph[start], q = 0;
    for (long k = 1; k <= NB; k++) {
        long j = (start + k) % NB;
        if (ph[j] >= 0 && ph[j] == ph[(j + 1) % NB] && ph[j] != cur) { q += ((ph[j] - cur + 3) % 3 == 1) ? 1 : -1; cur = ph[j]; }
    }
    return q;
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s NBLOCKS TMAX SEED [T1]\n", argv[0]); return 1; }
    NB = atol(argv[1]); N = 3 * NB; long TMAX = atol(argv[2]); u64 seed = strtoull(argv[3], 0, 10);
    long T1 = argc > 4 ? atol(argv[4]) : 3 * NB;
    char ck[256]; snprintf(ck, sizeof ck, "ckptq_N%ld_s%llu.bin", NB, (unsigned long long)seed);
    u64 *x = malloc(N * 8), *y = malloc(N * 8), locked = 0; long t = 0, lt[64]; int Q[64];
    u64 *cwbuf[3]; for (int k = 0; k < 3; k++) cwbuf[k] = malloc(NB * 8);
    int *ph = malloc(NB * sizeof(int));
    for (int b = 0; b < 64; b++) { lt[b] = -1; Q[b] = 0; }
    FILE *f = fopen(ck, "rb");
    if (f) {
        if (fread(&t, sizeof t, 1, f) != 1 || fread(&locked, sizeof locked, 1, f) != 1 || fread(lt, sizeof(long), 64, f) != 64 ||
            fread(Q, sizeof(int), 64, f) != 64 || fread(x, 8, N, f) != (size_t)N) { fprintf(stderr, "bad checkpoint\n"); return 1; }
        fclose(f); fprintf(stderr, "resumed at t=%ld\n", t);
    } else {
        s[0] = seed ^ 0x9E3779B97F4A7C15ULL; s[1] = seed * 0xBF58476D1CE4E5B9ULL + 1; s[2] = ~seed; s[3] = seed << 7 | 1;
        for (int i = 0; i < 20; i++) rnd();
        for (long i = 0; i < N; i++) x[i] = rnd();
    }
    time_t last = time(0);
    u64 need = ~0ULL;            /* runs still worth waiting for; excess right-moving runs (charge > 0) never lock */
    if (t > T1 + 2) { need = 0; for (int b = 0; b < 64; b++) if (Q[b] <= 0 || Q[b] == 99) need |= 1ULL << b; }
    while (t < TMAX && ((locked & need) != need || t < T1 + 2)) {
        step(x, y); u64 *tmp = x; x = y; y = tmp; t++;
        u64 all = ~0ULL;
        for (long j = 0; j < NB && all; j++) all &= ~(x[3 * j + 1] | x[3 * j + 2]);
        for (u64 nw = all & ~locked; nw; nw &= nw - 1) lt[__builtin_ctzll(nw)] = t;
        locked |= all;
        if (t >= T1 && t <= T1 + 2) {                         /* record the walls' net charge */
            for (long j = 0; j < NB; j++) cwbuf[t - T1][j] = ~(x[3 * j + 1] | x[3 * j + 2]);
            if (t == T1 + 2) { need = 0;
                for (int b = 0; b < 64; b++) { Q[b] = (lt[b] > 0 && lt[b] <= T1) ? 0 : charge(cwbuf, b, ph);
                                               if (Q[b] <= 0 || Q[b] == 99) need |= 1ULL << b; } }
        }
        if ((t & 65535) == 0 && time(0) - last > 600) {
            f = fopen(ck, "wb"); fwrite(&t, sizeof t, 1, f); fwrite(&locked, sizeof locked, 1, f); fwrite(lt, sizeof(long), 64, f);
            fwrite(Q, sizeof(int), 64, f); fwrite(x, 8, N, f); fclose(f); last = time(0);
            fprintf(stderr, "t=%ld locked %d/64\n", t, __builtin_popcountll(locked));
        }
    }
    for (int b = 0; b < 64; b++) printf("%d %ld\n", Q[b], lt[b]);
    remove(ck); return 0;
}
