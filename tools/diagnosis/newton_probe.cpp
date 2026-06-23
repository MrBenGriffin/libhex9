/* newton_probe.cpp — PROTOTYPE of the 2D-Newton encoder inversion.
 *
 * Target: h9_lonlat_to_boct (core/h9_math.h) — the geographic→barycentric
 * INVERSION that both encoder modes share (containment's
 * lonlatdeg_to_cxcy_oid calls it too, h9_grid.h:314). Today it is a
 * BEAM=6 × DEPTH=34 search: ~1,800 trig forward projections per point,
 * the encode "pig". The forward map forward(fx,fy)→lon/lat is analytic
 * and cheap (~0.6 µs). Newton borrows the easy direction to do the hard
 * one: solve forward(fx,fy)=(lon,lat).
 *
 * This is the UNWARPED AK projection inversion — the continuous forward
 * map only (no warp, no cell corners, no seams; per the design memo the
 * warp/seam hazards live in the separate warp_inv Newton, not here). So
 * the map is smooth except at octahedron vertices/edges, and FD Jacobian
 * is clean (no seam straddling).
 *
 * The triad implemented here:
 *   - forward map  = ORACLE (cheap+exact); per-iter residual check.
 *   - FD Jacobian  = cheap PREDICTOR (provisional; vetoed by the oracle).
 *   - beam search  = fallback oracle when Newton can't make progress.
 * Seed = gnomonic L1-normalised |ECEF direction| barycentric (analytic,
 * ~0 cost): the octant face vertices are the axis directions, so the
 * linear/gnomonic barycentric of the target is just |eX|,|eY|,|eZ|
 * renormalised — an excellent seed near the AK map's near-affine core.
 *
 * The probe does NOT modify the encoder. It runs Newton ALONGSIDE the
 * beam and measures, over an equal-area global sample:
 *   - boct agreement   : max |Δ(u,v,w)| Newton-vs-beam, and #exceeding 1e-12
 *   - UUID agreement   : full L29 + bins L8/12/16/20 (the output-preserving
 *                        test, routed through the real addressing/warp)
 *   - fallback rate    : how often the guard would punt to the beam
 *   - Newton iters     : histogram
 *   - timing           : beam µs/pt vs Newton µs/pt (inversion only)
 *
 * Build:  cmake --build build --target newton_probe && ./build/newton_probe [N]
 */
#include "hex9_c.h"
#include "h9_grid.h"      /* pulls h9_math.h + h9_addressing.h */
#include "newton_invert.h"

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <vector>

static double haversine_m(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371008.8, d2r = M_PI/180.0;
    const double dlat = (lat2-lat1)*d2r, dlon = (lon2-lon1)*d2r;
    const double a = std::sin(dlat/2)*std::sin(dlat/2) +
                     std::cos(lat1*d2r)*std::cos(lat2*d2r)*std::sin(dlon/2)*std::sin(dlon/2);
    return 2.0*R*std::asin(std::min(1.0, std::sqrt(a)));
}

/* ── reproducible equal-area sample (splitmix64), as encoder_ab_probe ──────── */
static uint64_t g_s = 0x9e3779b97f4a7c15ULL;
static double urand() {
    g_s += 0x9e3779b97f4a7c15ULL;
    uint64_t z = g_s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return (double)(z >> 11) / 9007199254740992.0;
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    const long N = (argc > 1) ? atol(argv[1]) : 200000;

    long   boct_bad   = 0;          /* max|Δuvw| > 1e-12 */
    double max_duvw   = 0.0;
    long   fellback   = 0;
    long   full_diff  = 0;
    const int bin_layers[] = { 8, 12, 16, 20 };
    long   bin_diff[4]    = {0,0,0,0};
    long   iters_hist[14] = {0};
    double sum_iters = 0.0;
    double rt_beam_max = 0.0, rt_beam_sum = 0.0;   /* roundtrip geodesic, m */
    double rt_newt_max = 0.0, rt_newt_sum = 0.0;

    for (long i = 0; i < N; ++i) {
        const double lon_d = (urand() * 2.0 - 1.0) * 180.0;
        const double lat_d = std::asin(urand() * 2.0 - 1.0) * 180.0 / M_PI;
        const double lon_r = lon_d * M_PI / 180.0, lat_r = lat_d * M_PI / 180.0;

        NStat st;
        const H9BOct nb = lonlat_to_boct_newton(lon_r, lat_r, &st);
        const H9BOct bb = h9_lonlat_to_boct_beam(lon_r, lat_r);

        const double du = fabs(nb.u - bb.u), dv = fabs(nb.v - bb.v), dw = fabs(nb.w - bb.w);
        const double dm = std::max(du, std::max(dv, dw));
        if (dm > max_duvw) max_duvw = dm;
        if (dm > 1e-12) boct_bad++;

        if (st.fellback) fellback++;
        sum_iters += st.iters;
        iters_hist[std::min(st.iters, 13)]++;

        /* output-preserving test: route both bocts through the real addressing */
        uint8_t un[16], ub[16];
        h9_boct_to_uuid_beam(nb, 6, 6, un);
        h9_boct_to_uuid_beam(bb, 6, 6, ub);
        if (std::memcmp(un, ub, 16) != 0) full_diff++;

        /* roundtrip fidelity: encode→decode geodesic from the original point */
        double dlon, dlat;
        hex9_decode(un, &dlon, &dlat);
        const double rn = haversine_m(lon_d, lat_d, dlon, dlat);
        hex9_decode(ub, &dlon, &dlat);
        const double rb = haversine_m(lon_d, lat_d, dlon, dlat);
        rt_newt_sum += rn; if (rn > rt_newt_max) rt_newt_max = rn;
        rt_beam_sum += rb; if (rb > rt_beam_max) rt_beam_max = rb;
        for (int b = 0; b < 4; ++b) {
            uint8_t bn[16], bbn[16];
            hex9_bin(un, bin_layers[b], bn);
            hex9_bin(ub, bin_layers[b], bbn);
            if (std::memcmp(bn, bbn, 16) != 0) bin_diff[b]++;
        }
    }

    printf("N = %ld equal-area global points\n\n", N);
    printf("── Newton vs beam: boct agreement ──\n");
    printf("max |Δ(u,v,w)|          : %.3e\n", max_duvw);
    printf("boct |Δuvw| > 1e-12     : %ld  (%.4f%%)\n", boct_bad, 100.0*boct_bad/N);
    printf("\n── output-preserving (routed through addressing+warp) ──\n");
    printf("full L29 UUID disagree  : %ld  (%.4f%%)\n", full_diff, 100.0*full_diff/N);
    for (int b = 0; b < 4; ++b)
        printf("bin L%-2d disagree        : %ld  (%.4f%%)\n",
               bin_layers[b], bin_diff[b], 100.0*bin_diff[b]/N);
    printf("\n── roundtrip fidelity (encode→decode geodesic, warp ON) ──\n");
    printf("beam   max %.2f nm   mean %.2f nm\n", rt_beam_max*1e9, rt_beam_sum/N*1e9);
    printf("newton max %.2f nm   mean %.2f nm\n", rt_newt_max*1e9, rt_newt_sum/N*1e9);
    printf("\n── guard / convergence ──\n");
    printf("would fall back to beam : %ld  (%.4f%%)\n", fellback, 100.0*fellback/N);
    printf("mean Newton iters       : %.3f\n", sum_iters/N);
    printf("iters histogram (0..13+):");
    for (int k = 0; k < 14; ++k) printf(" %ld", iters_hist[k]);
    printf("\n");

    /* ── timing: inversion only, fixed point set, separate loops ──────────── */
    const long T = (N < 40000) ? N : 40000;
    std::vector<double> lr(T), ar(T);
    for (long i = 0; i < T; ++i) {
        lr[i] = ((urand()*2.0-1.0)*180.0) * M_PI/180.0;
        ar[i] = (std::asin(urand()*2.0-1.0)*180.0/M_PI) * M_PI/180.0;
    }
    using clk = std::chrono::steady_clock;
    volatile double sink = 0.0;

    h9_lonlat_to_boct_beam(lr[0], ar[0]);
    auto t0 = clk::now();
    for (long i = 0; i < T; ++i) { H9BOct b = h9_lonlat_to_boct_beam(lr[i], ar[i]); sink += b.u; }
    auto t1 = clk::now();
    const double beam_us = std::chrono::duration<double,std::micro>(t1-t0).count()/T;

    NStat st;
    lonlat_to_boct_newton(lr[0], ar[0], &st);
    t0 = clk::now();
    for (long i = 0; i < T; ++i) { H9BOct b = lonlat_to_boct_newton(lr[i], ar[i], &st); sink += b.u; }
    t1 = clk::now();
    const double newt_us = std::chrono::duration<double,std::micro>(t1-t0).count()/T;

    printf("\n── timing (inversion only) ──\n");
    printf("beam   : %8.3f µs/pt  (%.0f pts/s)\n", beam_us, 1e6/beam_us);
    printf("newton : %8.3f µs/pt  (%.0f pts/s)\n", newt_us, 1e6/newt_us);
    printf("speedup: %.1f×   (sink %.1f)\n", beam_us/newt_us, (double)sink);
    return 0;
}
