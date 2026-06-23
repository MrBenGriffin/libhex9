/* seam_transect.c — walk a physical transect across the lon-0 seam at
 * lat 89.1 and classify every rendered cell: regular hexagon, sliver
 * (vertex-radius spread), and where it sits relative to the meridian.
 * Explains the banded sliver structure seen in QGIS at L16.
 * Part of the Hex9 (H9) Project. Apache-2.0. */
#include "hex9_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RING_N 7

static void ll_to_xyz(double lon, double lat, double v[3]) {
    const double D = M_PI / 180.0;
    v[0] = cos(lat * D) * cos(lon * D);
    v[1] = cos(lat * D) * sin(lon * D);
    v[2] = sin(lat * D);
}
/* chordal angle — robust for tiny separations where acos(dot) cancels */
static double ang(const double a[3], const double b[3]) {
    const double an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const double bn = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    const double dx = a[0]/an - b[0]/bn, dy = a[1]/an - b[1]/bn,
                 dz = a[2]/an - b[2]/bn;
    return 2.0 * asin(0.5 * sqrt(dx*dx + dy*dy + dz*dz));
}

int main(int argc, char **argv) {
    /* usage: seam_transect [L] [lat] [centre_m] [half_m] [step_m] [-v] */
    const int L = (argc > 1) ? atoi(argv[1]) : 16;
    const double lat = (argc > 2) ? atof(argv[2]) : 89.1;
    const double ctr_m  = (argc > 3) ? atof(argv[3]) : 0.0;
    const double half_m = (argc > 4) ? atof(argv[4]) : 5.0;
    const double step_m = (argc > 5) ? atof(argv[5]) : 0.02;
    const int verbose = (argc > 6);
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    const double cell_area = 4.0 * M_PI / (12.0 * pow(9.0, (double)L));
    const double r = sqrt(cell_area * 2.0 / (3.0 * sqrt(3.0)));   /* rad */
    const double R_E = 6371008.8;                                  /* m  */
    const double m_per_deg_lon = (M_PI / 180.0) * R_E * cos(lat * M_PI / 180.0);
    printf("L%d lat %g: circumradius %.4g m; transect lon %+g .. %+g m, step %g m\n",
           L, lat, r * R_E, ctr_m - half_m, ctr_m + half_m, step_m);

    const double m_per_deg_lat = (M_PI / 180.0) * R_E;
    uint8_t prev[16] = {0};
    for (double x = ctr_m - half_m; x <= ctr_m + half_m; x += step_m) {
        const double lon = x / m_per_deg_lon;
        uint8_t full[16], bin[16];
        if (hex9_encode(lon, lat, full) != 0) continue;
        if (hex9_bin(full, L, bin) != 0) continue;
        if (memcmp(bin, prev, 16) == 0) continue;
        memcpy(prev, bin, 16);

        double ring[2 * RING_N];
        if (hex9_cell_ring(full, L, 0, ring, RING_N) != RING_N) {
            printf("%+8.3f m: ring FAILED\n", x);
            continue;
        }
        double v3[6][3], cen[3] = {0, 0, 0}, p3[3];
        ll_to_xyz(lon, lat, p3);
        for (int k = 0; k < 6; k++) {
            ll_to_xyz(ring[2*k], ring[2*k + 1], v3[k]);
            cen[0] += v3[k][0]; cen[1] += v3[k][1]; cen[2] += v3[k][2];
        }
        double vmin = 1e30, vmax = 0.0;
        for (int k = 0; k < 6; k++) {
            const double d = ang(cen, v3[k]);
            if (d < vmin) vmin = d;
            if (d > vmax) vmax = d;
        }
        /* lon extent of the rendered polygon, in metres */
        double lo = 1e30, hi = -1e30;
        for (int k = 0; k < 6; k++) {
            const double xm = ring[2*k] * m_per_deg_lon;
            if (xm < lo) lo = xm;
            if (xm > hi) hi = xm;
        }
        /* centroid offset from the sample point, metres (lon, lat) */
        double clon = 0.0, clat = 0.0;
        for (int k = 0; k < 6; k++) { clon += ring[2*k]; clat += ring[2*k + 1]; }
        const double offx = (clon / 6.0 - lon) * m_per_deg_lon;
        const double offy = (clat / 6.0 - lat) * m_per_deg_lat;
        const char *cls = (vmax / vmin < 1.35) ? "hex   " :
                          (vmax / vmin < 2.0)  ? "warped" : "SLIVER";
        printf("%+10.4f m: %s  vr %.2f..%.2f r  poly-lon [%+9.4f, %+9.4f] m  "
               "span %.4f m  cen-off (%+.4f, %+.4f) m\n",
               x, cls, vmin / r, vmax / r, lo, hi, hi - lo, offx, offy);
        if (verbose) {
            for (int k = 0; k < 6; k++)
                printf("      v%d (%.9f, %.9f)\n", k, ring[2*k], ring[2*k + 1]);
        }
    }
    return 0;
}
