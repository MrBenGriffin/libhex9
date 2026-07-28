/* boct_io.c — the face-coordinate (bring-your-own-projection) surface.
 *
 * Proves the ⟨face⟩⟨cell⟩ seam is exactly the canonical chain's:
 *   1. project → encode_boct  ==  encode          (byte-equal, the AKW path
 *      re-entered through the public seam)
 *   2. to_woct → woct_to_boct → encode_boct == encode  (the researcher's 3D
 *      route: her [ox,oy,oz] through the projection-free chart bridge)
 *   3. decode_boct → unproject  ≈  decode          (centroid coherence)
 *   4. decode_boct → encode_boct round-trip        (centroid re-encodes to
 *      its own cell at max depth)
 *   5. cell_ring_boct → unproject  ≈  cell_ring    (per-vertex frames agree,
 *      densify 0 and 2)
 */
#include "hex9_c.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL "); printf(__VA_ARGS__); printf("\n"); ++fails; } \
} while (0)

/* Deterministic LCG sample (portable, same everywhere). */
static unsigned long long lcg = 0x9E3779B97F4A7C15ULL;
static double rnd(void) {
    lcg = lcg * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((lcg >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}

int main(void) {
    enum { N = 500 };
    double angdiff;

    /* The warp field MUST be live: the boct seam speaks the WARPED chart,
     * and an uninitialised run collapses BRAW == b_oct, making every
     * frame-confusion bug invisible (it did, during development). */
    char err[256] = {0};
    if (hex9_init(err, sizeof err)) {
        printf("FAIL init: %s\n", err);
        return 1;
    }

    for (int i = 0; i < N; ++i) {
        const double lon = rnd() * 360.0 - 180.0;
        /* equal-area latitude sample */
        const double lat = asin(2.0 * rnd() - 1.0) * 180.0 / M_PI;

        uint8_t canon[16], via_boct[16], via_woct[16];
        hex9_encode(lon, lat, canon);

        /* 1. through the public 2D seam */
        double cx, cy; int oid;
        hex9_project(lon, lat, &cx, &cy, &oid);
        CHECK(hex9_encode_boct(cx, cy, oid, via_boct) == 0,
              "encode_boct rc (%f,%f)", lon, lat);
        CHECK(memcmp(canon, via_boct, 16) == 0,
              "encode_boct != encode at (%f,%f)", lon, lat);

        /* 2. through the 3D w_oct route. The chart bridge is a rotation out
         * and back — ULP-exact, not bit-exact — so boundary-adjacent points
         * may flip LEAF-TAIL nibbles between the two routes (each route is
         * individually deterministic). Assert working-depth agreement (L20,
         * ~3 mm cells: a ULP flip there is ~1e-7 probable per point). */
        double xyz[3]; double cx2, cy2; int oid2;
        hex9_to_woct(lon, lat, xyz);
        hex9_woct_to_boct(xyz, &cx2, &cy2, &oid2);
        CHECK(hex9_encode_boct(cx2, cy2, oid2, via_woct) == 0,
              "encode_boct(woct) rc (%f,%f)", lon, lat);
        uint8_t bin_a[16], bin_b[16];
        hex9_bin(canon, 20, bin_a);
        hex9_bin(via_woct, 20, bin_b);
        CHECK(memcmp(bin_a, bin_b, 16) == 0,
              "woct route != encode at L20 bin, (%f,%f)", lon, lat);

        /* 3. centroid coherence: decode_boct + unproject vs decode */
        double dcx, dcy; int doid;
        CHECK(hex9_decode_boct(canon, &dcx, &dcy, &doid) == 0,
              "decode_boct rc (%f,%f)", lon, lat);
        double blon, blat, dlon, dlat;
        CHECK(hex9_unproject(dcx, dcy, doid, &blon, &blat) == 0,
              "unproject(decode_boct) rc (%f,%f)", lon, lat);
        hex9_decode(canon, &dlon, &dlat);
        angdiff = fabs(blon - dlon);
        if (angdiff > 180.0) angdiff = fabs(angdiff - 360.0);
        angdiff = fmax(angdiff * cos(dlat * M_PI / 180.0), fabs(blat - dlat));
        CHECK(angdiff < 1e-9, "centroid mismatch %.3e deg at (%f,%f)",
              angdiff, lon, lat);

        /* 4. centroid re-encodes to its own cell */
        uint8_t re[16];
        CHECK(hex9_encode_boct(dcx, dcy, doid, re) == 0,
              "re-encode rc (%f,%f)", lon, lat);
        CHECK(memcmp(canon, re, 16) == 0,
              "centroid re-encode != cell at (%f,%f)", lon, lat);
    }

    /* 5. ring coherence at a spread of layers, densify 0 and 2 */
    static const double pts[][2] = {
        { -3.19, 55.95 }, { -0.1276, 51.5072 }, { 139.7, 35.7 },
        { -62.396329, -15.365652 }, { 179.99, 0.01 }, { 45.0, 35.26439 },
        { 0.0, 89.9 }, { -111.728587, -83.344816 },
    };
    for (size_t p = 0; p < sizeof pts / sizeof pts[0]; ++p) {
        uint8_t u[16];
        hex9_encode(pts[p][0], pts[p][1], u);
        for (int layer = 2; layer <= 14; layer += 4) {
            for (int densify = 0; densify <= 2; densify += 2) {
                const int n = hex9_ring_npoints(densify);
                double rll[2 * 55], rcx[55], rcy[55];
                int roid[55];
                const int n1 = hex9_cell_ring(u, layer, densify, rll, n);
                const int n2 = hex9_cell_ring_boct(u, layer, densify,
                                                   rcx, rcy, roid, n);
                CHECK(n1 == n && n2 == n, "ring counts %d/%d layer %d", n1, n2, layer);
                if (n1 != n || n2 != n) continue;
                for (int v = 0; v < n; ++v) {
                    double vlon, vlat;
                    CHECK(hex9_unproject(rcx[v], rcy[v], roid[v], &vlon, &vlat) == 0,
                          "ring unproject rc p%zu l%d v%d", p, layer, v);
                    double dl = fabs(vlon - rll[2 * v]);
                    if (dl > 180.0) dl = fabs(dl - 360.0);
                    angdiff = fmax(dl * cos(rll[2 * v + 1] * M_PI / 180.0),
                                   fabs(vlat - rll[2 * v + 1]));
                    CHECK(angdiff < 1e-9,
                          "ring vertex mismatch %.3e deg p%zu l%d d%d v%d",
                          angdiff, p, layer, densify, v);
                }
            }
        }
    }

    printf("=== boct_io: %s (%d) ===\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}
