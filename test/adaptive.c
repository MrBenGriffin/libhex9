/* adaptive.c — hex9_adaptive_* (population-digest multi-layer grid).
 *
 * Checks the digestion invariants: conservation (emitted values sum to the
 * input weight total), ceiling respected (except single overweight points
 * and the terminal layer), floor respected for unit weights, per-point
 * assignment consistency (each point's cell is the rebin of its max_layer
 * bin at the cell's layer), and the saturated-column behaviour (one cell
 * per layer digesting ceiling-worth, remainder landing at min_layer).
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

/* deterministic LCG so the test is reproducible */
static unsigned long long rng_state = 42;
static double rng_uniform(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    /* ── uniform points, unit weights ────────────────────────────────── */
    enum { N = 20000 };
    static double lon[N], lat[N];
    for (int i = 0; i < N; ++i) {
        lon[i] = rng_uniform() * 360.0 - 180.0;
        lat[i] = rng_uniform() * 170.0 - 85.0;
    }
    const int    LMIN = 3, LMAX = 8;
    const double CEIL = 50.0, FLOOR = 5.0;
    /* adaptive takes FULL uuids (addresses, not coordinates) — encode first */
    static uint8_t uu[N * 16];
    hex9_encode_many(lon, lat, N, uu);
    hex9_adaptive *a = hex9_adaptive_create(uu, NULL, N, LMIN, LMAX, CEIL, FLOOR,
                                            err, sizeof err);
    CHECK(a != NULL, "create: %s", err);
    if (a) {
        const int m = hex9_adaptive_count(a);
        int64_t *assign = (int64_t *)malloc(N * sizeof(int64_t));
        hex9_adaptive_assign(a, assign);

        double total = 0.0;
        int64_t npts_total = 0;
        int layer_hist[30] = {0};
        for (int i = 0; i < m; ++i) {
            uint8_t u[16];
            int layer;
            double v;
            int64_t np;
            hex9_adaptive_cell(a, i, u, &layer, &v, &np);
            total += v;
            npts_total += np;
            layer_hist[layer]++;
            CHECK(layer >= LMIN && layer <= LMAX, "cell %d layer %d", i, layer);
            CHECK(np >= 1, "cell %d empty", i);
            CHECK(v <= CEIL + 1e-9 || layer == LMIN, "cell %d value %g > ceiling", i, v);
            CHECK(v >= FLOOR - 1e-9 || layer == LMIN, "cell %d value %g < floor", i, v);
        }
        CHECK(fabs(total - (double)N) < 1e-6, "conservation: %f != %d", total, N);
        CHECK(npts_total == N, "npoints: %lld != %d", (long long)npts_total, N);

        /* assignment: every point assigned; its cell is the CANONICAL bin of
         * the point AT THE CELL'S LAYER. NOT a label-prefix of the max_layer
         * bin: a bin is not a prefix of its address. At a split-hex (6/7/8)
         * leaf the canonical bin is the mode-0 parent, so a point on the
         * mode-1 side has bin nibbles that differ from its address nibbles —
         * by design (the uuid meta records the mode so the switch is
         * recoverable). The right check is direct equality against the
         * point's own canonical bin at the layer (identity path = h9_grid). */
        for (int i = 0; i < N; ++i) {
            CHECK(assign[i] >= 0 && assign[i] < m, "point %d unassigned", i);
            if (assign[i] < 0 || assign[i] >= m) break;
            uint8_t full[16], canon[16], cu[16];
            int layer;
            hex9_adaptive_cell(a, (int)assign[i], cu, &layer, NULL, NULL);
            hex9_encode(lon[i], lat[i], full);
            CHECK(hex9_k_disk(full, layer, 0, canon, 1) == 1, "point %d canon", i);
            CHECK(memcmp(cu, canon, 16) == 0,
                  "point %d: cell != canonical L%d bin of the point", i, layer);
            if (failures > 10) break;
        }
        printf("uniform: %d cells over layers", m);
        for (int l = LMIN; l <= LMAX; ++l) if (layer_hist[l]) printf(" L%d:%d", l, layer_hist[l]);
        printf("\n");
        free(assign);
        hex9_adaptive_destroy(a);
    }

    /* ── saturated column: all points in one cell ────────────────────── */
    {
        enum { M = 5000 };
        static double slon[M], slat[M];
        for (int i = 0; i < M; ++i) { slon[i] = -3.19; slat[i] = 55.95; }
        static uint8_t suu[M * 16];
        hex9_encode_many(slon, slat, M, suu);
        a = hex9_adaptive_create(suu, NULL, M, 3, 8, 50.0, 5.0, err, sizeof err);
        CHECK(a != NULL, "column create: %s", err);
        if (a) {
            const int m = hex9_adaptive_count(a);
            CHECK(m == 6, "column: %d cells != 6 (one per layer 8..3)", m);
            double total = 0.0;
            for (int i = 0; i < m; ++i) {
                int layer;
                double v;
                hex9_adaptive_cell(a, i, NULL, &layer, &v, NULL);
                total += v;
                CHECK(layer == 8 - i, "column cell %d at layer %d", i, layer);
                if (layer > 3) CHECK(v == 50.0, "column L%d value %g != 50", layer, v);
                else           CHECK(v == 5000.0 - 50.0 * 5, "column L3 value %g", v);
            }
            CHECK(total == 5000.0, "column conservation %g", total);
            hex9_adaptive_destroy(a);
        }
    }

    /* ── single overweight point lands at max_layer ──────────────────── */
    {
        double wlon[2] = { -3.19, 100.0 }, wlat[2] = { 55.95, 30.0 };
        double w[2] = { 1000.0, 1.0 };
        uint8_t wuu[2 * 16];
        hex9_encode_many(wlon, wlat, 2, wuu);
        a = hex9_adaptive_create(wuu, w, 2, 3, 8, 50.0, 5.0, err, sizeof err);
        CHECK(a != NULL, "weighted create: %s", err);
        if (a) {
            int found_heavy = 0;
            for (int i = 0; i < hex9_adaptive_count(a); ++i) {
                int layer;
                double v;
                hex9_adaptive_cell(a, i, NULL, &layer, &v, NULL);
                if (v == 1000.0) { found_heavy = 1; CHECK(layer == 8, "heavy at L%d", layer); }
            }
            CHECK(found_heavy, "overweight point not digested whole");
            hex9_adaptive_destroy(a);
        }
    }

    /* ── floor above everything: all mass reaches min_layer ─────────── */
    {
        double flon[10], flat[10];
        for (int i = 0; i < 10; ++i) { flon[i] = i * 30.0 - 150.0; flat[i] = 20.0; }
        uint8_t fuu[10 * 16];
        hex9_encode_many(flon, flat, 10, fuu);
        a = hex9_adaptive_create(fuu, NULL, 10, 2, 8, 1e9, 1e9, err, sizeof err);
        CHECK(a != NULL, "floor create: %s", err);
        if (a) {
            for (int i = 0; i < hex9_adaptive_count(a); ++i) {
                int layer;
                hex9_adaptive_cell(a, i, NULL, &layer, NULL, NULL);
                CHECK(layer == 2, "high-floor cell at L%d != 2", layer);
            }
            hex9_adaptive_destroy(a);
        }
    }

    /* ── error paths ─────────────────────────────────────────────────── */
    CHECK(hex9_adaptive_create(uu, NULL, 4, 9, 8, 10, 1, err, sizeof err) == NULL,
          "min>max accepted");
    CHECK(hex9_adaptive_create(uu, NULL, 4, 1, 8, 10, 20, err, sizeof err) == NULL,
          "floor>ceiling accepted");
    CHECK(hex9_adaptive_create(uu, NULL, 4, 1, 8, 0, 0, err, sizeof err) == NULL,
          "zero ceiling accepted");
    {   /* bin input is rejected — full uuids only (the digest re-bins
         * across layers, guaranteed only from the full uuid) */
        uint8_t bin8[16];
        hex9_bin(uu, 8, bin8);
        CHECK(hex9_adaptive_create(bin8, NULL, 1, 3, 8, 10, 1, err, sizeof err) == NULL,
              "bin input accepted");
    }

    if (failures) { printf("ADAPTIVE FAILED (%d)\n", failures); return 1; }
    printf("ADAPTIVE OK\n");
    return 0;
}
