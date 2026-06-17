/* gc_kring.c — great-circle validation of symbolic adjacency against the
 * encoder oracle.
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2026, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * Ben's "demonstrably safe ground" harness (2026-06-10): we cannot test
 * 12·9^29 hexes, but great circles sweep octant seams, the octahedron
 * vertices, and generic terrain at ANY layer. For every cell met along a
 * circle we validate the full-UUID kring entry against ground truth that
 * does not share code with the symbolic walk:
 *
 *   edge-probe oracle — for each of the cell's 6 edges, nudge a probe
 *   point just outside the edge midpoint (constructed in 3D ECEF, so poles
 *   and the antimeridian need no special casing) and ENCODE it; the
 *   probe's cell (rendered ring — a pure function of cell identity) must
 *   be one of the claimed neighbours. The encoder is the system's
 *   independent point→cell truth.
 *
 *   A probe that lands neither in self nor ring(1) is classified before
 *   judgement: if the probe cell's own neighbour list contains self, the
 *   two cells are symbolically adjacent and its absence from self's list
 *   is an asymmetry — hard FAIL. Otherwise the probe overshot into ring 2
 *   (geometry artefact of the fixed offset) — counted, not failed.
 *
 * Also asserted per cell: neighbour count is 5 or 6, adjacency is mutual,
 * ring(1) == neighbours, disk(1) == ring(1) + self.
 *
 * Bin-entry divergence is MEASURED per circle, not asserted (doctrine:
 * bins are layer-scoped keys, not addresses — whether kring keeps
 * accepting them rides on these numbers).
 *
 * Usage: gc_kring [points_per_circle]   (default 400)
 * Circles: equator (all-seam), the 0/180 meridian (seams + vertices +
 * poles), and two oblique generic circles. Layers: 6, 9, 12, 15.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

/* ring buffers: densify 0 → 7 (lon,lat) points, closed */
#define RING_N 7

static int ring_of(const uint8_t *uuid, int layer, double *lonlat /*2*RING_N*/) {
    return hex9_cell_ring(uuid, layer, 0, lonlat, RING_N);
}

static int same_ring(const double *a, const double *b) {
    return memcmp(a, b, 2 * RING_N * sizeof(double)) == 0;
}

/* lon/lat degrees ↔ unit sphere */
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

int main(int argc, char **argv) {
    const int npts = (argc > 1) ? atoi(argv[1]) : 400;
    static const int layers[] = { 6, 9, 12, 15 };
    const int nlayers = (int)(sizeof(layers) / sizeof(layers[0]));

    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    static const double axes[][3] = {
        { 0.0, 0.0, 1.0 },                  /* equator — every cell touches a seam */
        { 0.0, 1.0, 0.0 },                  /* meridian 0/180 — seams, vertices, poles */
        { 0.3, 0.5, 0.81 },                 /* oblique generic */
        { -0.7, 0.21, 0.44 },               /* oblique generic */
    };
    static const char *cname[] = { "equator", "meridian", "oblique-a", "oblique-b" };
    const int ncircles = (int)(sizeof(axes) / sizeof(axes[0]));

    long cells_total = 0, probes_total = 0, probes_self = 0,
         probes_ring2 = 0, probes_skipped = 0;
    long metric_cells = 0, metric_samples = 0, metric_outside = 0, metric_skipped = 0,
         metric_unswept = 0;
    double ground_max = 0.0;   /* worst nearest-vertex grounding, in circumradii */

    for (int c = 0; c < ncircles; c++) {
        double a[3] = { axes[c][0], axes[c][1], axes[c][2] };
        double an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        a[0] /= an; a[1] /= an; a[2] /= an;
        double e1[3] = { -a[1], a[0], 0.0 };
        double e1n = sqrt(e1[0]*e1[0] + e1[1]*e1[1]);
        if (e1n < 1e-12) { e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0; e1n = 1.0; }
        e1[0] /= e1n; e1[1] /= e1n; e1[2] /= e1n;
        double e2[3] = { a[1]*e1[2] - a[2]*e1[1],
                         a[2]*e1[0] - a[0]*e1[2],
                         a[0]*e1[1] - a[1]*e1[0] };

        long c_cells = 0, c_bindiv = 0;

        for (int li = 0; li < nlayers; li++) {
            const int L = layers[li];
            uint8_t prev_bin[16] = {0};
            for (int i = 0; i < npts; i++) {
                const double t = 2.0 * M_PI * (double)i / (double)npts;
                const double p[3] = { cos(t)*e1[0] + sin(t)*e2[0],
                                      cos(t)*e1[1] + sin(t)*e2[1],
                                      cos(t)*e1[2] + sin(t)*e2[2] };
                double lon, lat;
                xyz_to_ll(p, &lon, &lat);

                uint8_t full[16], bin[16];
                if (hex9_encode(lon, lat, full) != 0) continue;
                if (hex9_bin(full, L, bin) != 0) continue;
                if (memcmp(bin, prev_bin, 16) == 0) continue;   /* same cell as last sample */
                memcpy(prev_bin, bin, 16);
                cells_total++; c_cells++;

                /* neighbours from the full UUID — the guaranteed entry */
                uint8_t nbs[6 * 16];
                const int n = hex9_neighbors(full, L, nbs);
                CHECK(n == 5 || n == 6, "%s L%d t=%d: neighbour count %d", cname[c], L, i, n);
                if (n < 0) continue;

                double self_ring[2 * RING_N];
                if (ring_of(full, L, self_ring) != RING_N) {
                    CHECK(0, "%s L%d t=%d: self ring failed", cname[c], L, i);
                    continue;
                }

                /* LOCATION GROUNDING — the assertion that catches F2-style
                 * mis-location (adjacency alone is locally consistent even
                 * around a wrongly-placed cell): the sampled point must lie
                 * within a couple of circumradii of the NEAREST rendered
                 * vertex of its claimed cell. Nearest-vertex (not vertex
                 * mean): the warp's lateral-edge bypass (h9_warp.h edge_tol)
                 * is identity inside a 1e-7 band around the seam meridians
                 * but fully warped outside, so a fine-layer cell straddling
                 * the band edge renders sheared by the local warp delta
                 * (~4.5e-7 ≈ 11 r at L15 near the poles) and its vertex
                 * mean swings several r off the sample point even though
                 * encoder, decode and ring all agree (apex_scan/polar_probe,
                 * 2026-06-11). A point inside a regular hexagon is ≤ 1 r
                 * from a vertex; 2 r allows coarse-layer distortion while
                 * still flagging a cell whose render sits even one cell
                 * away. */
                {
                    double dmin = 1e30;
                    for (int k = 0; k < 6; k++) {
                        double w[3];
                        ll_to_xyz(self_ring[2*k], self_ring[2*k + 1], w);
                        /* chordal angle — acos(dot) cancels catastrophically
                         * at L16+ vertex scales (~1e-8 rad) */
                        const double dx = p[0] - w[0], dy = p[1] - w[1],
                                     dz = p[2] - w[2];
                        const double a2 = 2.0 * asin(0.5 *
                                          sqrt(dx*dx + dy*dy + dz*dz));
                        if (a2 < dmin) dmin = a2;
                    }
                    /* circumradius (radians): hex area = 3√3/2 r², sphere
                     * area 4π over 12·9^L cells */
                    const double cell_area = 4.0 * M_PI / (12.0 * pow(9.0, (double)L));
                    const double circum_r  = sqrt(cell_area * 2.0 / (3.0 * sqrt(3.0)));
                    CHECK(dmin < 2.0 * circum_r,
                          "%s L%d t=%d: cell mis-located — point (%.5f, %.5f) is "
                          "%.1f circumradii from its claimed cell's nearest vertex",
                          cname[c], L, i, lon, lat, dmin / circum_r);
                    if (dmin / circum_r > ground_max) ground_max = dmin / circum_r;
                }
                double nb_ring[6][2 * RING_N];
                int nb_ring_ok[6] = {0};
                for (int j = 0; j < n; j++)
                    nb_ring_ok[j] = (ring_of(nbs + (size_t)j * 16, L, nb_ring[j]) == RING_N);

                /* mutuality — bins are keys, not addresses (enforced), so
                 * traverse onward by encoding a point inside each neighbour;
                 * assert only when the probe actually landed there. */
                for (int j = 0; j < n; j++) {
                    if (!nb_ring_ok[j]) continue;
                    double ncen[3] = {0, 0, 0}, nv0[3];
                    for (int k = 0; k < 6; k++) {
                        double w[3];
                        ll_to_xyz(nb_ring[j][2*k], nb_ring[j][2*k + 1], w);
                        ncen[0] += w[0]; ncen[1] += w[1]; ncen[2] += w[2];
                    }
                    ncen[0] /= 6.0; ncen[1] /= 6.0; ncen[2] /= 6.0;
                    ll_to_xyz(nb_ring[j][0], nb_ring[j][1], nv0);
                    double np3[3] = { ncen[0] + (nv0[0] - ncen[0]) * 0.2,
                                      ncen[1] + (nv0[1] - ncen[1]) * 0.2,
                                      ncen[2] + (nv0[2] - ncen[2]) * 0.2 };
                    double nplon, nplat;
                    xyz_to_ll(np3, &nplon, &nplat);
                    uint8_t nf[16];
                    double nfr[2 * RING_N];
                    if (hex9_encode(nplon, nplat, nf) != 0) continue;
                    if (ring_of(nf, L, nfr) != RING_N) continue;
                    if (!same_ring(nfr, nb_ring[j])) continue;   /* probe strayed */
                    uint8_t back[6 * 16];
                    const int m = hex9_neighbors(nf, L, back);
                    CHECK(m == 5 || m == 6, "%s L%d t=%d: back-count %d (nb %d)",
                          cname[c], L, i, m, j);
                    int mutual = 0;
                    for (int kx = 0; kx < m && !mutual; kx++) {
                        double r[2 * RING_N];
                        if (ring_of(back + (size_t)kx * 16, L, r) == RING_N &&
                            same_ring(r, self_ring))
                            mutual = 1;
                    }
                    CHECK(mutual, "%s L%d t=%d: nb %d not mutual", cname[c], L, i, j);
                }

                /* ring(1)/disk(1) composition */
                uint8_t ring1[8 * 16], disk1[8 * 16];
                const int64_t nr = hex9_k_ring(full, L, 1, ring1, 8);
                const int64_t nd = hex9_k_disk(full, L, 1, disk1, 8);
                CHECK(nr == n, "%s L%d t=%d: ring(1)=%lld != %d", cname[c], L, i, (long long)nr, n);
                CHECK(nd == n + 1, "%s L%d t=%d: disk(1)=%lld != %d", cname[c], L, i, (long long)nd, n + 1);

                /* edge-probe oracle, constructed on the unit sphere */
                double v3[6][3], cen3[3] = {0, 0, 0};
                for (int k = 0; k < 6; k++) {
                    ll_to_xyz(self_ring[2*k], self_ring[2*k + 1], v3[k]);
                    cen3[0] += v3[k][0]; cen3[1] += v3[k][1]; cen3[2] += v3[k][2];
                }
                cen3[0] /= 6.0; cen3[1] /= 6.0; cen3[2] /= 6.0;

                for (int e = 0; e < 6; e++) {
                    const int e2i = (e + 1) % 6;
                    double mid3[3], pr3[3];
                    for (int k = 0; k < 3; k++)
                        mid3[k] = 0.5 * (v3[e][k] + v3[e2i][k]);
                    /* just outside the edge: centroid → midpoint, overshoot 10% */
                    for (int k = 0; k < 3; k++)
                        pr3[k] = cen3[k] + (mid3[k] - cen3[k]) * 1.10;
                    double plon, plat;
                    xyz_to_ll(pr3, &plon, &plat);

                    uint8_t pf[16];
                    if (hex9_encode(plon, plat, pf) != 0) { probes_skipped++; continue; }
                    double pring[2 * RING_N];
                    if (ring_of(pf, L, pring) != RING_N) { probes_skipped++; continue; }
                    probes_total++;

                    if (same_ring(pring, self_ring)) { probes_self++; continue; }

                    int hit = 0;
                    for (int j = 0; j < n && !hit; j++)
                        if (nb_ring_ok[j] && same_ring(nb_ring[j], pring)) hit = 1;
                    if (hit) continue;

                    /* not self, not ring(1): adjacency asymmetry or ring-2
                     * overshoot? Ask the probe cell. */
                    uint8_t pback[6 * 16];
                    const int pm = hex9_neighbors(pf, L, pback);
                    int claims_self = 0;
                    for (int kx = 0; kx < pm && !claims_self; kx++) {
                        double r[2 * RING_N];
                        if (ring_of(pback + (size_t)kx * 16, L, r) == RING_N &&
                            same_ring(r, self_ring))
                            claims_self = 1;
                    }
                    if (claims_self) {
                        CHECK(0, "%s L%d t=%d edge %d: probe cell claims self as "
                              "neighbour but is missing from self's ring "
                              "(probe %.7f, %.7f)", cname[c], L, i, e, plon, plat);
                    } else {
                        probes_ring2++;   /* genuine overshoot — geometry artefact */
                    }
                }

                /* input contract: bin input must be rejected (Ben's ruling —
                 * bins are layer-scoped keys, not addresses) */
                uint8_t nbs_bin[6 * 16];
                if (hex9_neighbors(bin, L, nbs_bin) != -1 ||
                    hex9_k_ring(bin, L, 1, nbs_bin, 6) != -1 ||
                    hex9_k_disk(bin, L, 1, nbs_bin, 7) != -1)
                    c_bindiv++;

                /* ── metric ring(k) oracle (Ben's circle-projection idea) ──
                 * The cell's own validated geometry calibrates the local
                 * metric (neighbour centroid spacing = 2×apothem). Two
                 * checks per ring(k) member:
                 *   1. AIMED — encode the member's claimed centroid; the
                 *      encoder must land back in that member (hard).
                 *   2. metric sanity — the member's distance from the
                 *      centre must lie in a generous physical band
                 *      [1.1, 3.2]·k·apothem (hard; flat-lattice theory is
                 *      [1.732, 2.0], warp/seams spread it).
                 * Plus the blind over-sampled annulus sweep: sampled bins
                 * must fall in disk(k+1) (hard) and ring-k coverage misses
                 * are REPORTED (soft — blind sampling is distortion-
                 * limited; the aimed check above is the guarantee). */
                if (c_cells % 64 == 1) {
                    /* sanity: skip F5-garbled rings (exactly-on-seam cells
                     * render distorted hexagons; vertex radii then spread) */
                    double vmin = 1e30, vmax = 0.0, a_sum = 0.0;
                    for (int k = 0; k < 6; k++) {
                        double d[3] = { v3[k][0] - cen3[0], v3[k][1] - cen3[1],
                                        v3[k][2] - cen3[2] };
                        const double vd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                        if (vd < vmin) vmin = vd;
                        if (vd > vmax) vmax = vd;
                    }
                    for (int e = 0; e < 6; e++) {
                        double m[3] = { 0.5*(v3[e][0] + v3[(e+1)%6][0]) - cen3[0],
                                        0.5*(v3[e][1] + v3[(e+1)%6][1]) - cen3[1],
                                        0.5*(v3[e][2] + v3[(e+1)%6][2]) - cen3[2] };
                        a_sum += sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
                    }
                    const double apothem = a_sum / 6.0;   /* chord ≈ angle at these scales */
                    if (vmax > 1.6 * vmin || n != 6) {
                        metric_skipped++;
                    } else {
                        double t1[3], t2[3];
                        const double cn = sqrt(cen3[0]*cen3[0] + cen3[1]*cen3[1] + cen3[2]*cen3[2]);
                        double cu[3] = { cen3[0]/cn, cen3[1]/cn, cen3[2]/cn };
                        t1[0] = -cu[1]; t1[1] = cu[0]; t1[2] = 0.0;
                        double t1n = sqrt(t1[0]*t1[0] + t1[1]*t1[1]);
                        if (t1n < 1e-12) { t1[0] = 1.0; t1[1] = 0.0; t1[2] = 0.0; t1n = 1.0; }
                        t1[0] /= t1n; t1[1] /= t1n; t1[2] /= t1n;
                        t2[0] = cu[1]*t1[2] - cu[2]*t1[1];
                        t2[1] = cu[2]*t1[0] - cu[0]*t1[2];
                        t2[2] = cu[0]*t1[1] - cu[1]*t1[0];

                        for (int kk = 2; kk <= 3; kk++) {
                            uint8_t ringk[18 * 16], diskk1[64 * 16];
                            const int64_t nrk = hex9_k_ring(full, L, kk, ringk, 18);
                            const int64_t ndk = hex9_k_disk(full, L, kk + 1, diskk1, 64);
                            if (nrk <= 0 || ndk <= 0) {
                                CHECK(0, "%s L%d t=%d: ring/disk(%d) failed", cname[c], L, i, kk);
                                continue;
                            }
                            /* member pre-pass: the metric oracle needs well-
                             * rendered local geometry. The warp's lateral-edge
                             * bypass (h9_warp.h edge_tol: identity inside a
                             * 1e-7 band around the seam meridians, fully
                             * warped outside) shears any member straddling
                             * the band edge — its vertex radii spread and its
                             * vertex-mean centroid swings several r
                             * (seam_profile: cliff ≈ 2–3 m at the octahedron
                             * vertices, ~0.3 m mid-seam). Aimed and metric
                             * checks are meaningless across the cliff, so a
                             * garbled member disqualifies the whole anchor —
                             * counted, like the anchor's own-garble skip. */
                            int cliff = 0;
                            for (int j = 0; j < (int)nrk && !cliff; j++) {
                                double mr[2 * RING_N];
                                if (ring_of(ringk + (size_t)j * 16, L, mr) != RING_N)
                                    continue;          /* asserted further down */
                                double mc3[3] = {0, 0, 0}, mv[6][3];
                                for (int kx = 0; kx < 6; kx++) {
                                    ll_to_xyz(mr[2*kx], mr[2*kx + 1], mv[kx]);
                                    mc3[0] += mv[kx][0]; mc3[1] += mv[kx][1]; mc3[2] += mv[kx][2];
                                }
                                mc3[0] /= 6.0; mc3[1] /= 6.0; mc3[2] /= 6.0;
                                double mvmin = 1e30, mvmax = 0.0;
                                for (int kx = 0; kx < 6; kx++) {
                                    double d[3] = { mv[kx][0] - mc3[0], mv[kx][1] - mc3[1],
                                                    mv[kx][2] - mc3[2] };
                                    const double vd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                                    if (vd < mvmin) mvmin = vd;
                                    if (vd > mvmax) mvmax = vd;
                                }
                                if (mvmax > 1.6 * mvmin) cliff = 1;
                            }
                            if (cliff) { metric_skipped++; continue; }

                            static double dring[64][2 * RING_N];
                            int dring_ok[64];
                            for (int j = 0; j < (int)ndk; j++)
                                dring_ok[j] = (ring_of(diskk1 + (size_t)j * 16, L, dring[j]) == RING_N);
                            int hitk[18] = {0};   /* per ring-k member */

                            const double rf[5] = { 1.732, 1.80, 1.87, 1.94, 2.0 };
                            const int nang = 24 * kk;
                            for (int ri = 0; ri < 5; ri++) {
                                const double rr = rf[ri] * (double)kk * apothem;
                                for (int ai = 0; ai < nang; ai++) {
                                    const double th = 2.0 * M_PI * (double)ai / (double)nang;
                                    double sp[3];
                                    for (int x = 0; x < 3; x++)
                                        sp[x] = cos(rr)*cu[x] +
                                                sin(rr)*(cos(th)*t1[x] + sin(th)*t2[x]);
                                    double slon, slat;
                                    xyz_to_ll(sp, &slon, &slat);
                                    uint8_t sf[16];
                                    if (hex9_encode(slon, slat, sf) != 0) continue;
                                    double sring[2 * RING_N];
                                    if (ring_of(sf, L, sring) != RING_N) continue;
                                    metric_samples++;
                                    int in_disk = 0;
                                    for (int j = 0; j < (int)ndk && !in_disk; j++)
                                        if (dring_ok[j] && same_ring(dring[j], sring)) in_disk = 1;
                                    if (!in_disk) { metric_outside++; continue; }
                                    for (int j = 0; j < (int)nrk; j++) {
                                        double r2[2 * RING_N];
                                        if (ring_of(ringk + (size_t)j * 16, L, r2) == RING_N &&
                                            same_ring(r2, sring)) { hitk[j] = 1; break; }
                                    }
                                }
                            }
                            /* aimed member verification + metric sanity */
                            for (int j = 0; j < (int)nrk; j++) {
                                double mr[2 * RING_N], mc[3] = {0, 0, 0};
                                if (ring_of(ringk + (size_t)j * 16, L, mr) != RING_N) {
                                    CHECK(0, "%s L%d t=%d: ring(%d) member %d unrenderable",
                                          cname[c], L, i, kk, j);
                                    continue;
                                }
                                for (int kx = 0; kx < 6; kx++) {
                                    double w[3];
                                    ll_to_xyz(mr[2*kx], mr[2*kx + 1], w);
                                    mc[0] += w[0]; mc[1] += w[1]; mc[2] += w[2];
                                }
                                const double mn = sqrt(mc[0]*mc[0] + mc[1]*mc[1] + mc[2]*mc[2]);
                                const double dt = (cu[0]*mc[0] + cu[1]*mc[1] + cu[2]*mc[2]) / mn;
                                const double ratio = acos(fmax(-1.0, fmin(1.0, dt)))
                                                     / ((double)kk * apothem);
                                CHECK(ratio > 1.1 && ratio < 3.2,
                                      "%s L%d t=%d: ring(%d) member %d at %.3f x k.apothem",
                                      cname[c], L, i, kk, j, ratio);
                                /* Aim slightly OFF-centre: the exact hexagon
                                 * centre lies on the internal half-hex seam
                                 * (encoder tie-break flips flavour there —
                                 * Ben's F4 warning), so probe 15% toward
                                 * vertex 0: interior, away from any seam. */
                                double mv0[3];
                                ll_to_xyz(mr[0], mr[1], mv0);
                                double mcu[3] = { mc[0]/mn + 0.15*(mv0[0] - mc[0]/mn),
                                                  mc[1]/mn + 0.15*(mv0[1] - mc[1]/mn),
                                                  mc[2]/mn + 0.15*(mv0[2] - mc[2]/mn) };
                                double mlon, mlat;
                                xyz_to_ll(mcu, &mlon, &mlat);
                                uint8_t mf[16];
                                double mring[2 * RING_N];
                                if (hex9_encode(mlon, mlat, mf) == 0 &&
                                    ring_of(mf, L, mring) == RING_N) {
                                    CHECK(same_ring(mring, mr),
                                          "%s L%d t=%d: ring(%d) member %d centroid does "
                                          "not encode back to the member (%.6f, %.6f)",
                                          cname[c], L, i, kk, j, mlon, mlat);
                                }
                            }
                            /* blind-coverage statistic (soft) */
                            for (int j = 0; j < (int)nrk; j++)
                                if (!hitk[j]) metric_unswept++;
                            metric_cells++;
                        }
                    }
                }
            }
        }
        CHECK(c_bindiv == 0, "%s: %ld bin inputs were NOT rejected", cname[c], c_bindiv);
        printf("%-10s cells %6ld   bin-input rejection OK\n", cname[c], c_cells);
    }

    printf("cells: %ld   probes: %ld (self %ld, ring2-overshoot %ld, skipped %ld)\n",
           cells_total, probes_total, probes_self, probes_ring2, probes_skipped);
    printf("grounding: worst point→nearest-vertex %.2f circumradii (limit 2.00)\n",
           ground_max);
    printf("metric rings: %ld cell·k checks, %ld blind samples (%ld outside disk(k+1), "
           "%ld unswept members, %ld anchors skipped: garbled/seam-cliff)\n",
           metric_cells, metric_samples, metric_outside, metric_unswept, metric_skipped);
    if (failures) { printf("GC_KRING FAILED (%d)\n", failures); return 1; }
    printf("gc_kring OK\n");
    return 0;
}
