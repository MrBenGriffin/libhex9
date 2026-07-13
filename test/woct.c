/* woct.c — w_oct 3D storage CRS: seamless round-trip, unit-octahedron
 * invariant, oid == sign(xyz), and b_oct<->xyz rotation round-trip. */
#include "hex9_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int sign_oid(const double xyz[3]) {
    return ((xyz[2] < 0) << 2) | ((xyz[1] < 0) << 1) | (xyz[0] < 0);
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("init FAILED: %s\n", err); return 1; }

    /* the historically nasty spots: both poles + all four seam meridians
     * ({lat, lon}); equator crossings at lon 0/90/180/270 and lat ±90 poles. */
    const double hard[][2] = {
        {0, 0}, {0, 90}, {0, 180}, {0, -90}, {90, 0}, {-90, 0},
        {45, 135}, {-30, -170}, {60, -45}
    };
    const size_t H = sizeof(hard) / sizeof(hard[0]);

    const size_t R = (argc > 1) ? (size_t)strtoul(argv[1], NULL, 10) : 200000;
    const size_t N = R + H;
    double *lon = malloc(N * sizeof(double));
    double *lat = malloc(N * sizeof(double));
    double *xyz = malloc(N * 3 * sizeof(double));
    double *rlon = malloc(N * sizeof(double));
    double *rlat = malloc(N * sizeof(double));
    if (!lon || !lat || !xyz || !rlon || !rlat) { printf("OOM\n"); return 1; }

    for (size_t i = 0; i < R; ++i) {
        double u = (double)((i * 2654435761u) % N) / (double)N;   /* [0,1) */
        double v = (double)((i * 40503u + 12345u) % N) / (double)N;
        lon[i] = -180.0 + 360.0 * u;
        lat[i] = asin(2.0 * v - 1.0) * 180.0 / M_PI;              /* area-uniform */
    }
    for (size_t i = 0; i < H; ++i) { lat[R + i] = hard[i][0]; lon[R + i] = hard[i][1]; }

    /* lon/lat -> c_oct (batch) */
    hex9_to_woct_many(lon, lat, N, xyz);

    /* invariants: on unit octahedron, oid == sign(xyz) */
    double max_l1 = 0.0; size_t sign_bad = 0;
    for (size_t i = 0; i < N; ++i) {
        const double *p = xyz + i * 3;
        double l1 = fabs(p[0]) + fabs(p[1]) + fabs(p[2]);
        if (fabs(l1 - 1.0) > max_l1) max_l1 = fabs(l1 - 1.0);
        /* skip on-plane (a component ~0): sign is ambiguous there by design */
        int on_plane = fabs(p[0]) < 1e-12 || fabs(p[1]) < 1e-12 || fabs(p[2]) < 1e-12;
        int oid;
        hex9_woct_to_boct(p, &(double){0}, &(double){0}, &oid);
        if (!on_plane && sign_oid(p) != oid) ++sign_bad;
    }

    /* c_oct -> lon/lat round-trip (batch) */
    hex9_from_woct_many(xyz, N, rlon, rlat);
    double max_gc = 0.0;
    for (size_t i = 0; i < N; ++i) {
        if (fabs(fabs(lat[i]) - 90.0) < 1e-9) continue;   /* lon undefined at poles */
        double dlo = fabs(rlon[i] - lon[i]); if (dlo > 180.0) dlo = 360.0 - dlo;
        double dla = fabs(rlat[i] - lat[i]);
        double gc = hypot(dlo * cos(lat[i] * M_PI / 180.0), dla);
        if (gc > max_gc) max_gc = gc;
    }

    /* b_oct <-> xyz is a pure rotation: project -> lift -> unlift == project.
     * Only defined for INTERIOR points: on a seam/pole the (cx,cy,oid) rep is
     * multivalued (canonicalised to a mode-0 sibling), so skip on-plane xyz —
     * lon/lat still round-trips there (rt_max above), which is the real test. */
    double max_bo = 0.0; size_t interior = 0;
    for (size_t i = 0; i < N; ++i) {
        double cx, cy; int oid;
        hex9_project(lon[i], lat[i], &cx, &cy, &oid);
        double p[3]; hex9_boct_to_woct(cx, cy, oid, p);
        if (fabs(p[0]) < 1e-9 || fabs(p[1]) < 1e-9 || fabs(p[2]) < 1e-9) continue;
        double cx2, cy2; int oid2;
        hex9_woct_to_boct(p, &cx2, &cy2, &oid2);
        double e = fabs(cx2 - cx) + fabs(cy2 - cy) + (oid2 != oid ? 1.0 : 0.0);
        if (e > max_bo) max_bo = e;
        ++interior;
    }

    printf("N=%zu  max|L1-1|=%.2e  sign!=oid=%zu  rt_max=%.2e deg  boct_rt_max=%.2e\n",
           N, max_l1, sign_bad, max_gc, max_bo);
    int ok = (max_l1 < 1e-12) && (sign_bad == 0) && (max_gc < 1e-4) && (max_bo < 1e-12);
    printf(ok ? "WOCT OK\n" : "WOCT FAILED\n");
    return ok ? 0 : 1;
}
