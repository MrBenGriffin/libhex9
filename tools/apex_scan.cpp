/* apex_scan.cpp — localize the uv→lonlat discontinuity behind the four polar
 * gc_kring residuals. polar_probe showed the rendered ring is garbled: at
 * (180, 89.1) L15 the vertices (433365, 28264445) and (433364, 28264446)
 * project ~11 circumradii away while their lattice neighbours project true.
 * corner_trace cleared the lattice algebra (all six vertices in-face, one
 * frame), so the fault is inside cxcy_to_lonlat: either h9_warp_fwd or the
 * boct→lonlat projection. Scan a lattice line through the bad vertex, print
 * the warp output and the final lonlat step sizes — whichever stage jumps
 * is the culprit.
 * Part of the Hex9 (H9) Project. Apache-2.0. */
#define H9_WARP_ENABLE 1
#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"
#include "h9_cell_geom.h"
#include "h9_warp_runtime.h"
#include <cstdio>
#include <cmath>

static void ll_to_xyz(double lon, double lat, double v[3]) {
    const double D = M_PI / 180.0;
    v[0] = cos(lat * D) * cos(lon * D);
    v[1] = cos(lat * D) * sin(lon * D);
    v[2] = sin(lat * D);
}
static double ang3(const double a[3], const double b[3]) {
    const double d = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return acos(fmax(-1.0, fmin(1.0, d)));
}

int main() {
    std::string err;
    if (!h9::h9_warp_init_embedded(&err)) { printf("warp init failed: %s\n", err.c_str()); return 1; }

    const int L = 15, oid = 1;            /* pt(180, 89.1) cell frame */
    const int64_t ib = 28264445;
    const double div_f = std::pow(3.0, (double)L);
    const double cell_area = 4.0 * M_PI / (12.0 * pow(9.0, (double)L));
    const double r = sqrt(cell_area * 2.0 / (3.0 * sqrt(3.0)));
    const int mo = (int)H9_OID_MO[oid];

    /* u-scan at v = ib through bad vertex u=433365 (good at 433367) */
    for (int pass = 0; pass < 2; pass++) {
        const int64_t s3 = h9kring::pow3(L);
        printf("%s-scan (oid=%d mo=%d, 2s=%lld):\n", pass ? "v" : "u", oid, mo,
               (long long)(2 * s3));
        double prev3[3], prevb[2];
        int have = 0;
        for (int dd = -8; dd <= 8; dd++) {
            const int64_t u = pass ? 433365            : 433365 + dd;
            const int64_t v = pass ? ib + dd           : ib;
            const double cx = (double)u * H9_UV_U1 / div_f;
            const double cy = (double)v * H9_UV_V3 / div_f;
            double bx, by;
            h9_warp_fwd(cx, cy, mo, &bx, &by);
            double lon, lat;
            h9grid::cxcy_to_lonlat(cx, cy, oid, &lon, &lat, nullptr, false);
            double p3[3];
            ll_to_xyz(lon, lat, p3);
            const double dw = have ? hypot(bx - prevb[0], by - prevb[1]) : 0.0;
            const double dl = have ? ang3(prev3, p3) : 0.0;
            /* u+v vs the 2s apex edge: the suspect boundary */
            printf("  (%lld,%lld) u+v-2s=%lld  warp(%+.12e,%+.12e) dwarp=%.3e  "
                   "ll(%.7f,%.7f) step=%.2f r\n",
                   (long long)u, (long long)v, (long long)(u + v - 2 * s3),
                   bx, by, dw, lon, lat, dl / r);
            prevb[0] = bx; prevb[1] = by;
            prev3[0] = p3[0]; prev3[1] = p3[1]; prev3[2] = p3[2];
            have = 1;
        }
        printf("\n");
    }
    return 0;
}
