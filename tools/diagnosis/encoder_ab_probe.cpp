/* encoder_ab_probe.cpp — re-measure the NN (mode 0) vs containment (mode 1)
 * encoder disagreement. The "4.7% of points" figure is from an old 50k sample;
 * this re-establishes it on a fresh equal-area global sample, and also asks the
 * question that matters in practice: do the two encoders disagree only in the
 * deep leaf tail, or do they land in DIFFERENT BINS at usable layers?
 *
 * For each point: encode under both modes, compare the full L29 UUID and the
 * bins at L8/L12/L16/L20.
 */
#include "hex9_c.h"

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <vector>

/* Reproducible equal-area point: lon uniform, lat = asin(2u-1). */
static uint64_t s = 0x9e3779b97f4a7c15ULL;
static double urand() {                 /* splitmix64 → [0,1) */
    s += 0x9e3779b97f4a7c15ULL;
    uint64_t z = s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return (double)(z >> 11) / 9007199254740992.0;
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    const long N = (argc > 1) ? atol(argv[1]) : 200000;
    const int  bin_layers[] = { 8, 12, 16, 20 };
    const int  NB = (int)(sizeof(bin_layers)/sizeof(bin_layers[0]));

    long full_diff = 0;
    long bin_diff[4] = {0,0,0,0};

    for (long i = 0; i < N; ++i) {
        const double lon = (urand() * 2.0 - 1.0) * 180.0;
        const double lat = std::asin(urand() * 2.0 - 1.0) * 180.0 / M_PI;

        uint8_t u0[16], u1[16];
        hex9_set_encoder(0); if (hex9_encode(lon, lat, u0) != 0) continue;
        hex9_set_encoder(1); if (hex9_encode(lon, lat, u1) != 0) continue;

        if (std::memcmp(u0, u1, 16) != 0) full_diff++;

        for (int b = 0; b < NB; ++b) {
            uint8_t b0[16], b1[16];
            hex9_bin(u0, bin_layers[b], b0);
            hex9_bin(u1, bin_layers[b], b1);
            if (std::memcmp(b0, b1, 16) != 0) bin_diff[b]++;
        }
    }

    printf("N = %ld equal-area global points\n", N);
    printf("full L29 UUID disagree : %ld  (%.3f%%)\n",
           full_diff, 100.0 * (double)full_diff / (double)N);
    for (int b = 0; b < NB; ++b)
        printf("bin L%-2d disagree       : %ld  (%.3f%%)\n",
               bin_layers[b], bin_diff[b], 100.0 * (double)bin_diff[b] / (double)N);

    /* ── Timing: each mode over the SAME fixed point set, separate loops ──── */
    const long T = (N < 40000) ? N : 40000;
    std::vector<double> lons(T), lats(T);
    for (long i = 0; i < T; ++i) {
        lons[i] = (urand() * 2.0 - 1.0) * 180.0;
        lats[i] = std::asin(urand() * 2.0 - 1.0) * 180.0 / M_PI;
    }
    using clk = std::chrono::steady_clock;
    uint8_t u[16];
    for (int mode = 0; mode <= 1; ++mode) {
        hex9_set_encoder(mode);
        hex9_encode(lons[0], lats[0], u);          /* warm */
        const auto t0 = clk::now();
        for (long i = 0; i < T; ++i) hex9_encode(lons[i], lats[i], u);
        const auto t1 = clk::now();
        const double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / (double)T;
        printf("encode %-12s: %.2f µs/point  (%.0f pts/s)\n",
               mode == 1 ? "containment" : "NN", us, 1e6 / us);
    }
    return 0;
}
