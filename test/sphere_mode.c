/* sphere_mode.c — the sphere-datum twins (2.1.0).
 *
 * Doctrine under test (hex9_c.h "TWO DATUMS, ONE REGIME"): the *_sphere
 * entry points run the identical chain minus the WGS84 authalic reduction.
 * So:
 *   1. PARITY: encode_sphere(lon, xi(lat)) must land in the same cell as
 *      encode(lon, lat), where xi is the WGS84 authalic latitude (computed
 *      here from h9_authalic.h, the same series the library uses). Where the
 *      reduction is exactly identity (equator, poles) the full UUID must be
 *      BIT-IDENTICAL; at mid-latitudes the deg->rad->deg round-trip of xi
 *      allows ~1 ULP of wobble, so parity is pinned at the L20 bin (~mm).
 *   2. DIVERGENCE: feeding the SAME numeric lon/lat to both datums must give
 *      different addresses (coarse bins agree, full UUIDs differ) — the
 *      hazard the datum-as-dataset-metadata doctrine exists for.
 *   3. Round-trips, projection twins, and the grid handle's datum memory.
 *
 * A failure in 1 means the sphere path is NOT the same chain (a second
 * regime has crept in). A failure in 2 means the reduction became a no-op
 * (sphere aux leaking into the WGS84 path or vice versa).
 */
#include "hex9_c.h"
#include "../core/h9_authalic.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const double WGS84_E2 = 0.0066943799901413199;
static H9Authalic g_aux;

static int g_fails = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); g_fails++; } \
} while (0)

/* WGS84 geodetic degrees -> authalic (spherical) degrees. */
static double xi_deg(double lat_deg)
{
    return h9_geodetic_to_authalic(&g_aux, lat_deg * (M_PI / 180.0))
           * (180.0 / M_PI);
}

static void hexify(const unsigned char *u, char *out)
{
    static const char *H = "0123456789abcdef";
    for (int i = 0; i < 16; ++i) { out[2*i] = H[u[i] >> 4]; out[2*i+1] = H[u[i] & 15]; }
    out[32] = '\0';
}

/* Mid-latitude probes (the authalic shift peaks near 45deg). */
static const double PTS[][2] = {
    { -3.188267,  55.953251 },   /* edinburgh */
    { -0.124625,  51.500745 },   /* westminster */
    { 151.215256, -33.856159 },  /* sydney */
    { -68.302832, -54.801912 },  /* ushuaia */
    { 139.691706,  35.689487 },  /* tokyo */
    {  18.423300, -33.918861 },  /* cape town */
};
static const int N_PTS = (int)(sizeof PTS / sizeof PTS[0]);

/* Exact-identity probes: xi(lat) == lat bit-exactly (series terms all
 * sin(k*pi) == 0), so full-UUID parity must be exact. */
static const double PTS_EXACT[][2] = {
    {   0.0,    0.0 },
    {  17.25,   0.0 },
    { -179.99,  0.0 },
    {   0.0,   90.0 },
    {  45.0,  -90.0 },
};
static const int N_EXACT = (int)(sizeof PTS_EXACT / sizeof PTS_EXACT[0]);

int main(void)
{
    char err[256];

    /* hex9_init is the 2.1.0 canonical name; hex9_warp_init the alias.
     * Both idempotent, both 0. */
    CHECK(hex9_init(err, sizeof err) == 0, "hex9_init: %s", err);
    CHECK(hex9_warp_init(err, sizeof err) == 0, "hex9_warp_init alias: %s", err);
    CHECK(h9_authalic_init(&g_aux, WGS84_E2), "test-local authalic series");

    /* ── 1a. exact parity where the reduction is identity ─────────────────── */
    for (int i = 0; i < N_EXACT; ++i) {
        const double lon = PTS_EXACT[i][0], lat = PTS_EXACT[i][1];
        unsigned char a[16], b[16];
        CHECK(hex9_encode(lon, lat, a) == 0, "encode exact[%d]", i);
        CHECK(hex9_encode_sphere(lon, lat, b) == 0, "encode_sphere exact[%d]", i);
        CHECK(memcmp(a, b, 16) == 0,
              "exact parity [%d] (%.4f, %.4f): identity reduction must be "
              "bit-identical", i, lon, lat);
    }

    /* ── 1b. mid-latitude parity at the L20 bin ───────────────────────────── */
    for (int i = 0; i < N_PTS; ++i) {
        const double lon = PTS[i][0], lat = PTS[i][1];
        unsigned char a[16], b[16], ba[16], bb[16];
        CHECK(hex9_encode(lon, lat, a) == 0, "encode [%d]", i);
        CHECK(hex9_encode_sphere(lon, xi_deg(lat), b) == 0, "encode_sphere [%d]", i);
        CHECK(hex9_bin(a, 20, ba) == 0 && hex9_bin(b, 20, bb) == 0, "bin [%d]", i);
        if (memcmp(ba, bb, 16) != 0) {
            char ha[33], hb[33];
            hexify(ba, ha); hexify(bb, hb);
            CHECK(0, "L20 parity [%d] (%.4f, %.4f):\n  wgs84  %s\n  sphere %s",
                  i, lon, lat, ha, hb);
        }
    }

    /* ── 2. divergence: same numbers, different datum, different address ──── */
    for (int i = 0; i < N_PTS; ++i) {
        const double lon = PTS[i][0], lat = PTS[i][1];
        unsigned char a[16], b[16], ba[16], bb[16];
        hex9_encode(lon, lat, a);
        hex9_encode_sphere(lon, lat, b);   /* same numerics, sphere datum */
        CHECK(memcmp(a, b, 16) != 0,
              "divergence [%d]: datums must give different addresses at "
              "|lat| > 30", i);
        /* ... but only below the coarse layers (max authalic shift ~21 km,
         * L2 cells ~300 km): the L2 bins must agree. */
        hex9_bin(a, 2, ba);
        hex9_bin(b, 2, bb);
        CHECK(memcmp(ba, bb, 16) == 0, "divergence [%d]: L2 bins must agree", i);
    }

    /* ── 3. sphere round-trips (encode/decode, project/unproject) ─────────── */
    for (int i = 0; i < N_PTS; ++i) {
        const double lon = PTS[i][0], lat_s = xi_deg(PTS[i][1]);
        unsigned char u[16];
        double lo, la;
        hex9_encode_sphere(lon, lat_s, u);
        CHECK(hex9_decode_sphere(u, &lo, &la) == 0, "decode_sphere [%d]", i);
        CHECK(fabs(lo - lon) < 1e-6 && fabs(la - lat_s) < 1e-6,
              "sphere encode/decode RT [%d]: (%.8f, %.8f) -> (%.8f, %.8f)",
              i, lon, lat_s, lo, la);

        double cx, cy; int oid;
        CHECK(hex9_project_sphere(lon, lat_s, &cx, &cy, &oid) == 0,
              "project_sphere [%d]", i);
        CHECK(hex9_unproject_sphere(cx, cy, oid, &lo, &la) == 0,
              "unproject_sphere [%d]", i);
        CHECK(fabs(lo - lon) < 1e-9 && fabs(la - lat_s) < 1e-9,
              "sphere project/unproject RT [%d]", i);

        /* project_sphere == project modulo the reduction: same cx/cy from
         * reduced input (the projection-level parity pin). */
        double wx, wy; int woid;
        hex9_project(lon, PTS[i][1], &wx, &wy, &woid);
        CHECK(woid == oid && fabs(wx - cx) < 1e-12 && fabs(wy - cy) < 1e-12,
              "projection parity [%d]: d=(%.3e, %.3e)", i, wx - cx, wy - cy);
    }

    /* ── 4. batch == scalar, sphere flavour ───────────────────────────────── */
    {
        double lons[6], lats[6], blon[6], blat[6];
        unsigned char many[6 * 16];
        for (int i = 0; i < N_PTS; ++i) {
            lons[i] = PTS[i][0];
            lats[i] = xi_deg(PTS[i][1]);
        }
        CHECK(hex9_encode_many_sphere(lons, lats, 6, many) == 0, "encode_many_sphere");
        for (int i = 0; i < N_PTS; ++i) {
            unsigned char u[16];
            hex9_encode_sphere(lons[i], lats[i], u);
            CHECK(memcmp(u, many + i * 16, 16) == 0, "many==scalar sphere [%d]", i);
        }
        CHECK(hex9_decode_many_sphere(many, 6, blon, blat) == 0, "decode_many_sphere");
        for (int i = 0; i < N_PTS; ++i)
            CHECK(fabs(blon[i] - lons[i]) < 1e-6 && fabs(blat[i] - lats[i]) < 1e-6,
                  "decode_many_sphere RT [%d]", i);
    }

    /* ── 5. the grid handle remembers its datum ───────────────────────────── */
    {
        const int L = 6;
        hex9_grid *gs = hex9_grid_create_sphere(-4.0, 55.0, -2.0, 56.5,
                                                L, 0, 100000, err, sizeof err);
        CHECK(gs != NULL, "grid_create_sphere: %s", err);
        hex9_grid *gw = hex9_grid_create(-4.0, 55.0, -2.0, 56.5,
                                         L, 0, 100000, err, sizeof err);
        CHECK(gw != NULL, "grid_create: %s", err);
        if (gs && gw) {
            CHECK(hex9_grid_count(gs) > 0, "sphere grid empty");
            /* Same numeric bbox, different datum: the cell sets differ. */
            int diff = (hex9_grid_count(gs) != hex9_grid_count(gw));
            if (!diff) {
                for (int i = 0; i < hex9_grid_count(gs); ++i) {
                    unsigned char us[16], uw[16];
                    hex9_grid_cell_uuid(gs, i, us);
                    hex9_grid_cell_uuid(gw, i, uw);
                    if (memcmp(us, uw, 16) != 0) { diff = 1; break; }
                }
            }
            CHECK(diff, "sphere and wgs84 grids over the same numeric bbox "
                        "must differ");
            /* Containment coherence IN the sphere datum: the handle's
             * centroid re-encodes (sphere) into the handle's own cell. */
            for (int i = 0; i < hex9_grid_count(gs); i += 7) {
                double clon, clat;
                unsigned char cu[16], eu[16], bu[16];
                hex9_grid_cell_uuid(gs, i, cu);
                hex9_grid_cell_centroid(gs, i, &clon, &clat);
                hex9_encode_sphere(clon, clat, eu);
                hex9_bin(eu, L, bu);
                CHECK(memcmp(cu, bu, 16) == 0,
                      "sphere grid containment [%d]: centroid does not re-bin "
                      "to its own cell", i);
            }
            /* cell_ring_sphere == the sphere handle's own ring. */
            {
                unsigned char cu[16];
                double ra[64], rb[64];
                hex9_grid_cell_uuid(gs, 0, cu);
                const int na = hex9_grid_cell_ring(gs, 0, 1, ra, 32);
                const int nb = hex9_cell_ring_sphere(cu, L, 1, rb, 32);
                CHECK(na == nb && na > 0, "ring point counts (%d vs %d)", na, nb);
                for (int k = 0; na == nb && k < 2 * na; ++k)
                    CHECK(ra[k] == rb[k], "ring vertex %d: %.12f vs %.12f",
                          k, ra[k], rb[k]);
            }
        }
        hex9_grid_destroy(gs);
        hex9_grid_destroy(gw);
    }

    if (g_fails) { printf("%d FAILURE(S)\n", g_fails); return 1; }
    printf("sphere_mode: all assertions passed\n");
    return 0;
}
