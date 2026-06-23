/* seam_profile.cpp — magnitude of the lateral-edge bypass discontinuity
 * along the whole seam (equator → pole).
 *
 * The bypass band (h9_warp.h edge_tol=1e-7 on the line value y−√3x+W_eq,
 * i.e. ~0.5e-7 Euclidean) is identity; just outside, the CT warp applies
 * its full local delta — so the discontinuity at the band edge IS the
 * local warp delta on the seam. This tool walks the right lateral edge of
 * the mode-0 frame, apex (0,−√6/3) → equator corner (√2/2, √6/6), offsets
 * inward to a few line-values, and prints |warp_do(p)−p|. Units: 1 b_raw
 * unit ≈ 7,077 km (the edge, pole→equator ≈ 10,007 km, has length √2).
 * Part of the Hex9 (H9) Project. Apache-2.0. */
#define H9_WARP_ENABLE 1
#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_warp_runtime.h"
#include <cstdio>
#include <cmath>

int main() {
    std::string err;
    if (!h9::h9_warp_init_embedded(&err)) { printf("warp init failed: %s\n", err.c_str()); return 1; }

    const double W  = std::sqrt(6.0) / 3.0;
    const double ax = 0.0,                  ay = -W;                    /* apex = pole  */
    const double ex = std::sqrt(2.0) / 2.0, ey = std::sqrt(6.0) / 6.0;  /* equator end  */
    const double M_PER_UNIT = 10007.543e3 / std::sqrt(2.0);            /* ≈ 7,076 km    */

    /* line values to sample: in-band, then just outside, then far field */
    const double fvals[] = { 5e-8, 2e-7, 1e-6, 1e-5, 1e-4 };
    printf("%-7s %-9s |delta| at line-value f (b_raw units; band edge at f=1e-7)\n",
           "t", "lat(deg)");
    printf("%-7s %-9s %-10s %-10s %-10s %-10s %-10s\n",
           "", "", "5e-8", "2e-7", "1e-6", "1e-5", "1e-4");
    for (int i = 1; i < 40; i++) {
        const double t  = (double)i / 40.0;
        const double px = ax + (ex - ax) * t;
        const double py = ay + (ey - ay) * t;
        double lon, lat;
        h9grid::cxcy_to_lonlat(px, py, 0, &lon, &lat, nullptr, false);
        printf("%-7.3f %-9.3f", t, lat);
        for (double f : fvals) {
            /* f(x,y) = y − √3x + W; ∇f=(−√3,1); f(P+α∇f) = 4α ⇒ α=f/4 */
            const double qx = px - std::sqrt(3.0) * (f / 4.0);
            const double qy = py + (f / 4.0);
            double wx, wy;
            h9_warp_fwd(qx, qy, 0, &wx, &wy);
            const double d = std::hypot(wx - qx, wy - qy);
            printf(" %-10.3e", d);
            if (f == 2e-7) printf("(%6.2f m)", d * M_PER_UNIT);
        }
        printf("\n");
    }
    return 0;
}
