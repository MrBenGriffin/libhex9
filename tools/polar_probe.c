/* polar_probe.c — encoder bin-transition probe for the four polar gc_kring
 * residuals (lat ±89.1, lon 0/180, L15, all "3.8 circumradii").
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2026, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * corner_trace proved these cells are NOT seam-straddlers (all vertices
 * in-face; two exactly on the seam edge; cell abuts the pole apex). Two
 * hypotheses remain:
 *   (a) gc_kring's 3-circumradii vertex-mean grounding is too tight under
 *       near-pole projection distortion — the rendered ring is truthful and
 *       the cell genuinely spans the gap;
 *   (b) the rendering mis-places the cell — the encoder's actual boundary
 *       is somewhere else entirely.
 *
 * Decisive test: walk the encoder. (1) march the geodesic from the sample
 * point to the claimed centroid, recording every bin transition — if there
 * is none, the encoder itself says point and centroid share a cell; (2) scan
 * across each rendered edge midpoint and locate the encoder's transition —
 * if it sits at the rendered edge (t ≈ 1.0 on the centroid→midpoint ray),
 * the ring matches the encoder boundary. Both holding ⇒ (a): relax the
 * grounding threshold near the poles.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define RING_N 7

static void ll_to_xyz(double lon, double lat, double v[3]) {
    const double D = M_PI / 180.0;
    v[0] = cos(lat * D) * cos(lon * D);
    v[1] = cos(lat * D) * sin(lon * D);
    v[2] = sin(lat * D);
}
static void xyz_to_ll(const double v[3], double *lon, double *lat) {
    const double n = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    *lat = asin(fmax(-1.0, fmin(1.0, v[2] / n))) * 180.0 / M_PI;
    *lon = atan2(v[1], v[0]) * 180.0 / M_PI;
}
static double ang(const double a[3], const double b[3]) {
    const double an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const double bn = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    const double d = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (an * bn);
    return acos(fmax(-1.0, fmin(1.0, d)));
}
/* spherical interpolation a→b at fraction t */
static void slerp(const double a[3], const double b[3], double t, double out[3]) {
    const double th = ang(a, b);
    if (th < 1e-15) { memcpy(out, a, 3 * sizeof(double)); return; }
    const double s = sin(th);
    const double wa = sin((1.0 - t) * th) / s, wb = sin(t * th) / s;
    for (int k = 0; k < 3; k++) out[k] = wa * a[k] + wb * b[k];
}
/* bin of the point at layer L; -1 on encode failure */
static int bin_at(double lon, double lat, int L, uint8_t out[16]) {
    uint8_t full[16];
    if (hex9_encode(lon, lat, full) != 0) return -1;
    return hex9_bin(full, L, out);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    const int L = 15;
    const double pts[4][2] = { {180, 89.1}, {-0.0, 89.1}, {0, -89.1}, {180, -89.1} };

    /* nominal circumradius, same formula as gc_kring */
    const double cell_area = 4.0 * M_PI / (12.0 * pow(9.0, (double)L));
    const double r = sqrt(cell_area * 2.0 / (3.0 * sqrt(3.0)));
    printf("L%d nominal circumradius %.6e rad\n\n", L, r);

    for (int pi = 0; pi < 4; pi++) {
        const double lon = pts[pi][0], lat = pts[pi][1];
        double p3[3];
        ll_to_xyz(lon, lat, p3);

        uint8_t full[16], bin[16];
        if (hex9_encode(lon, lat, full) != 0 || hex9_bin(full, L, bin) != 0) {
            printf("pt(%g,%g): encode/bin failed\n", lon, lat);
            continue;
        }
        double ring[2 * RING_N];
        if (hex9_cell_ring(full, L, 0, ring, RING_N) != RING_N) {
            printf("pt(%g,%g): ring failed\n", lon, lat);
            continue;
        }

        /* claimed centroid + rendered shape stats */
        double v3[6][3], cen[3] = {0, 0, 0};
        for (int k = 0; k < 6; k++) {
            ll_to_xyz(ring[2*k], ring[2*k + 1], v3[k]);
            for (int x = 0; x < 3; x++) cen[x] += v3[k][x];
        }
        double vmin = 1e30, vmax = 0.0;
        for (int k = 0; k < 6; k++) {
            const double d = ang(cen, v3[k]);
            if (d < vmin) vmin = d;
            if (d > vmax) vmax = d;
        }
        double clon, clat;
        xyz_to_ll(cen, &clon, &clat);
        printf("pt(%g,%g)  claimed centroid (%.5f,%.5f)  dist %.2f r\n",
               lon, lat, clon, clat, ang(p3, cen) / r);
        printf("  rendered vertex radii: %.2f r .. %.2f r (regular hex = 1.0)\n",
               vmin / r, vmax / r);
        for (int k = 0; k < 6; k++)
            printf("    v%d (%.6f,%.6f)  from-point %.2f r  from-centroid %.2f r\n",
                   k, ring[2*k], ring[2*k + 1], ang(p3, v3[k]) / r, ang(cen, v3[k]) / r);
        {
            double dlon, dlat, d3[3];
            if (hex9_decode(full, &dlon, &dlat) == 0) {
                ll_to_xyz(dlon, dlat, d3);
                printf("  decode(full) = (%.6f,%.6f)  from-point %.2f r\n",
                       dlon, dlat, ang(p3, d3) / r);
            }
        }

        /* pole distances — how close is this cell to the apex? */
        const double pole[3] = { 0, 0, (lat > 0) ? 1.0 : -1.0 };
        printf("  point→pole %.2f r   centroid→pole %.2f r\n",
               ang(p3, pole) / r, ang(cen, pole) / r);

        /* (1) geodesic march p → centroid: every encoder bin transition */
        const int NSTEP = 4000;
        uint8_t cur[16], nxt[16];
        memcpy(cur, bin, 16);
        int ntrans = 0;
        printf("  march p→centroid (%.2f r total):\n", ang(p3, cen) / r);
        for (int i = 1; i <= NSTEP; i++) {
            double q[3], qlon, qlat;
            slerp(p3, cen, (double)i / NSTEP, q);
            xyz_to_ll(q, &qlon, &qlat);
            if (bin_at(qlon, qlat, L, nxt) != 0) continue;
            if (memcmp(nxt, cur, 16) != 0) {
                printf("    transition %d at %.3f r (t=%.4f)\n",
                       ++ntrans, ang(p3, cen) * i / NSTEP / r, (double)i / NSTEP);
                memcpy(cur, nxt, 16);
            }
        }
        printf("    %d transition(s); centroid in %s bin as sample point\n",
               ntrans, memcmp(cur, bin, 16) == 0 ? "SAME" : "a DIFFERENT");

        /* (2) per-edge boundary scan: centroid→edge-midpoint ray, t∈[0.3,2];
         * encoder transition at t≈1.0 ⇔ rendered edge = encoder boundary */
        printf("  edge-boundary scan (transition t, rendered edge at 1.0):\n   ");
        for (int e = 0; e < 6; e++) {
            double mid[3];
            for (int k = 0; k < 3; k++)
                mid[k] = 0.5 * (v3[e][k] + v3[(e + 1) % 6][k]);
            uint8_t prev[16];
            int have = 0;
            double tstar = -1.0;
            for (int i = 0; i <= 1700; i++) {
                const double t = 0.3 + (2.0 - 0.3) * i / 1700.0;
                double q[3] = { cen[0] + (mid[0] - cen[0]) * t,
                                cen[1] + (mid[1] - cen[1]) * t,
                                cen[2] + (mid[2] - cen[2]) * t };
                double qlon, qlat;
                xyz_to_ll(q, &qlon, &qlat);
                uint8_t b[16];
                if (bin_at(qlon, qlat, L, b) != 0) continue;
                if (have && memcmp(b, prev, 16) != 0 && tstar < 0.0) tstar = t;
                memcpy(prev, b, 16);
                have = 1;
            }
            if (tstar < 0.0) printf(" e%d:none", e);
            else             printf(" e%d:%.3f", e, tstar);
        }
        printf("\n\n");
    }
    return 0;
}
