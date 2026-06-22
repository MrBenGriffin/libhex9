/* batch.c — verify the OMP batch path is bit-identical to scalar (proves
 * reentrancy/no races) and round-trips, and report the parallel speedup. */
#include "hex9_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double wall(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("init FAILED: %s\n", err); return 1; }

    /* encode is heavy (warp-inverse Newton ~hundreds of µs/pt), so default
     * small; pass a larger N on the command line to stress / benchmark. */
    const size_t N = (argc > 1) ? (size_t)strtoul(argv[1], NULL, 10) : 20000;
    double  *lon = malloc(N * sizeof(double));
    double  *lat = malloc(N * sizeof(double));
    uint8_t *ref = malloc(N * 16);
    uint8_t *bat = malloc(N * 16);
    double  *dlon = malloc(N * sizeof(double));
    double  *dlat = malloc(N * sizeof(double));
    if (!lon || !lat || !ref || !bat || !dlon || !dlat) { printf("OOM\n"); return 1; }

    for (size_t i = 0; i < N; ++i) {
        lon[i] = -180.0 + 360.0 * ((double)i / (double)N);
        lat[i] =  -85.0 + 170.0 * ((double)((i * 2654435761u) % N) / (double)N);
    }

    /* scalar reference (single thread) */
    double t0 = wall();
    for (size_t i = 0; i < N; ++i) hex9_encode(lon[i], lat[i], ref + i * 16);
    double t_scalar = wall() - t0;

    /* batch (OMP) */
    t0 = wall();
    hex9_encode_many(lon, lat, N, bat);
    double t_batch = wall() - t0;

    /* batch must equal scalar bit-for-bit */
    size_t mism = 0;
    for (size_t i = 0; i < N * 16; ++i) if (ref[i] != bat[i]) ++mism;

    /* decode_many round-trip */
    hex9_decode_many(bat, N, dlon, dlat);
    double maxe = 0.0;
    for (size_t i = 0; i < N; ++i) {
        /* Longitude wraps: a point on the antimeridian may decode to the other
         * sign (e.g. the L30 canonical centroid returns +180 for a -180 input).
         * That is zero geodesic error, so fold the diff into [0,180]. */
        double dlo = fabs(dlon[i] - lon[i]);
        if (dlo > 180.0) dlo = 360.0 - dlo;
        double e = dlo + fabs(dlat[i] - lat[i]);
        if (e > maxe) maxe = e;
    }

    printf("N=%zu  scalar=%.3fs  batch=%.3fs  speedup=%.1fx  mismatched_bytes=%zu  rt_max=%.2e deg\n",
           N, t_scalar, t_batch, t_scalar / t_batch, mism, maxe);

    int ok = (mism == 0) && (maxe < 1e-4);
    printf(ok ? "BATCH OK\n" : "BATCH FAILED\n");
    return ok ? 0 : 1;
}
