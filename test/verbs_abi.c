/* verbs_abi.c — behavioural laws for the C grid verbs (aim / walk_to /
 * vision_cone), the twins of the python wheel's hex9.verbs.
 *
 * Same scenario and laws as test/verbs.py: a field self-calibrated from the
 * measured hex pitch, straddling the equator octant seam on purpose. The
 * checks are behavioural, not pinned — bearings are advisory FP; the keys
 * come from the (separately pinned) canonical chain.
 */
#include "hex9_c.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); fails++; } \
} while (0)

#define LAYER 8
#define KRING 6

static const double D2R = 3.14159265358979323846 / 180.0;

static double gc_bearing(double lon1, double lat1, double lon2, double lat2) {
    const double p1 = lat1 * D2R, p2 = lat2 * D2R, dl = (lon2 - lon1) * D2R;
    double b = atan2(sin(dl) * cos(p2),
                     cos(p1) * sin(p2) - sin(p1) * cos(p2) * cos(dl)) / D2R;
    b = fmod(b, 360.0);
    return b < 0 ? b + 360.0 : b;
}

static double gc_angle(double lon1, double lat1, double lon2, double lat2) {
    const double p1 = lat1 * D2R, p2 = lat2 * D2R, dl = (lon2 - lon1) * D2R;
    double h = sin((p2 - p1) / 2) * sin((p2 - p1) / 2)
             + cos(p1) * cos(p2) * sin(dl / 2) * sin(dl / 2);
    if (h < 0) h = 0;
    if (h > 1) h = 1;
    return 2.0 * asin(sqrt(h)) / D2R;
}

static double wrap180(double a) {
    a = fmod(a + 180.0, 360.0);
    if (a < 0) a += 360.0;
    return a - 180.0;
}

static void key_of(double lon, double lat, uint8_t out[16]) {
    uint8_t full[16];
    hex9_encode(lon, lat, full);
    hex9_bin(full, LAYER, out);
}

static int in_set(const uint8_t *set, int64_t n, const uint8_t k[16]) {
    for (int64_t i = 0; i < n; ++i)
        if (memcmp(set + i * 16, k, 16) == 0) return 1;
    return 0;
}

int main(void)
{
    char err[256];
    if (hex9_init(err, sizeof err)) { printf("FAIL: init: %s\n", err); return 1; }

    const double FLON = 35.02, FLAT = -0.02;    /* just south of the seam */
    const double BEAR = 15.0, HALF = 32.0;
    uint8_t fox[16], full[16];
    key_of(FLON, FLAT, fox);
    double flon, flat;
    CHECK(hex9_decode(fox, &flon, &flat) == 0, "fox centroid");
    CHECK(hex9_encode(flon, flat, full) == 0, "fox full");

    /* pitch, self-calibrated like the module does */
    uint8_t nbs[6 * 16];
    const int nn = hex9_neighbors(full, LAYER, nbs);
    CHECK(nn >= 5, "neighbors");
    double pitch = 0;
    {
        double d[6];
        for (int i = 0; i < nn; ++i) {
            double nlon, nlat;
            hex9_decode(nbs + (size_t)i * 16, &nlon, &nlat);
            d[i] = gc_angle(flon, flat, nlon, nlat);
        }
        for (int i = 1; i < nn; ++i)
            for (int j = i; j > 0 && d[j] < d[j - 1]; --j) {
                double t = d[j]; d[j] = d[j - 1]; d[j - 1] = t;
            }
        pitch = (nn & 1) ? d[nn / 2] : 0.5 * (d[nn / 2 - 1] + d[nn / 2]);
    }
    CHECK(pitch > 0, "pitch");

    /* ── aim: inverse of every neighbour's bearing; guards ── */
    for (int i = 0; i < nn; ++i) {
        double nlon, nlat;
        hex9_decode(nbs + (size_t)i * 16, &nlon, &nlat);
        uint8_t got[16];
        CHECK(hex9_aim(fox, LAYER, gc_bearing(flon, flat, nlon, nlat), got) == 0,
              "aim rc (nb %d)", i);
        CHECK(memcmp(got, nbs + (size_t)i * 16, 16) == 0,
              "aim missed neighbour %d", i);
    }

    /* ── vision_cone: src included; cone ⊆ disk; full circle == disk;
     *    bearing law; seam crossing ── */
    const int64_t maxd = hex9_disk_ncells(KRING);
    uint8_t *disk = malloc((size_t)maxd * 16);
    uint8_t *cone = malloc((size_t)maxd * 16);
    uint8_t *cone2 = malloc((size_t)maxd * 16);
    const int64_t ndisk = hex9_k_disk(full, LAYER, KRING, disk, maxd);
    CHECK(ndisk > 0, "k_disk");
    const int64_t ncone = hex9_vision_cone(fox, LAYER, BEAR, HALF, KRING,
                                           NULL, 0, cone, maxd);
    CHECK(ncone > 0, "vision_cone rc");
    CHECK(in_set(cone, ncone, fox), "src not in cone");
    int south = 0, north = 0;
    for (int64_t i = 0; i < ncone; ++i) {
        CHECK(in_set(disk, ndisk, cone + i * 16), "cone cell outside disk");
        double clon, clat;
        hex9_decode(cone + i * 16, &clon, &clat);
        if (memcmp(cone + i * 16, fox, 16) != 0) {
            CHECK(fabs(wrap180(gc_bearing(flon, flat, clon, clat) - BEAR))
                      <= HALF, "cell outside bearing cone");
            if (clat > 0) north++; else south++;
        }
    }
    CHECK(north > 0 && south > 0, "cone did not cross the seam (%d/%d)",
          north, south);
    const int64_t nfull = hex9_vision_cone(fox, LAYER, 0.0, 180.0, KRING,
                                           NULL, 0, cone2, maxd);
    CHECK(nfull == ndisk, "full-circle cone != disk (%lld vs %lld)",
          (long long)nfull, (long long)ndisk);

    /* ── occlusion: hedge perpendicular to the bearing at mid-cone ── */
    uint8_t hedge[40 * 16];
    int64_t nhedge = 0;
    {
        const double b = BEAR * D2R;
        const double cx = FLON + 3.2 * pitch * sin(b);
        const double cy = FLAT + 3.2 * pitch * cos(b);
        for (int i = 0; i < 40; ++i) {
            const double t = (-1.6 + 3.2 * i / 39.0) * pitch;
            uint8_t hk[16];
            key_of(cx + t * cos(b), cy - t * sin(b), hk);
            if (memcmp(hk, fox, 16) == 0) continue;
            if (!in_set(hedge, nhedge, hk))
                memcpy(hedge + (nhedge++) * 16, hk, 16);
        }
    }
    CHECK(nhedge > 1, "hedge degenerate");
    const int64_t nseen = hex9_vision_cone(fox, LAYER, BEAR, HALF, KRING,
                                           hedge, (size_t)nhedge, cone2, maxd);
    CHECK(nseen > 0 && nseen < ncone, "hedge cast no shadow (%lld of %lld)",
          (long long)nseen, (long long)ncone);
    int front = 0;
    for (int64_t i = 0; i < nseen; ++i) {
        CHECK(in_set(cone, ncone, cone2 + i * 16), "obstacles revealed cells");
        if (in_set(hedge, nhedge, cone2 + i * 16)) front++;
    }
    CHECK(front > 0, "hedge front face invisible");

    /* ── walk_to: identity, single step, adjacency, detour, sealed dest ── */
    uint8_t path[512 * 16];
    int64_t np = hex9_walk_to(fox, fox, LAYER, NULL, 0, 0, path, 512);
    CHECK(np == 1 && memcmp(path, fox, 16) == 0, "walk to self");
    uint8_t step[16];
    CHECK(hex9_aim(fox, LAYER, BEAR, step) == 0, "aim step");
    np = hex9_walk_to(fox, step, LAYER, NULL, 0, 0, path, 512);
    CHECK(np == 2 && memcmp(path + 16, step, 16) == 0, "walk one step");

    /* farthest shadowed cell as the target */
    uint8_t far[16];
    int have_far = 0;
    double fard = -1;
    for (int64_t i = 0; i < ncone; ++i) {
        if (in_set(cone2, nseen, cone + i * 16)) continue;   /* not shadowed */
        double clon, clat;
        hex9_decode(cone + i * 16, &clon, &clat);
        const double d = gc_angle(flon, flat, clon, clat);
        if (d > fard) { fard = d; memcpy(far, cone + i * 16, 16); have_far = 1; }
    }
    CHECK(have_far, "no shadow target");
    const int64_t nd = hex9_walk_to(fox, far, LAYER, NULL, 0, 0, path, 512);
    CHECK(nd >= 2 && memcmp(path, fox, 16) == 0 &&
          memcmp(path + (nd - 1) * 16, far, 16) == 0, "walk endpoints");
    for (int64_t i = 0; i + 1 < nd; ++i) {
        uint8_t f2[16], nb2[6 * 16];
        double l1, l2;
        hex9_decode(path + i * 16, &l1, &l2);
        hex9_encode(l1, l2, f2);
        const int n2 = hex9_neighbors(f2, LAYER, nb2);
        CHECK(n2 > 0 && in_set(nb2, n2, path + (i + 1) * 16),
              "path step %lld not a neighbour", (long long)i);
    }
    const int64_t nv = hex9_walk_to(fox, far, LAYER, hedge, (size_t)nhedge,
                                    0, path, 512);
    CHECK(nv >= nd, "detour shorter than direct (%lld < %lld)",
          (long long)nv, (long long)nd);
    for (int64_t i = 0; i < nv; ++i)
        if (memcmp(path + i * 16, far, 16) != 0)
            CHECK(!in_set(hedge, nhedge, path + i * 16), "walk through hedge");

    /* sealed destination: ring of obstacles → no path */
    uint8_t farfull[16], ring[8 * 16];
    double rlon, rlat;
    hex9_decode(far, &rlon, &rlat);
    hex9_encode(rlon, rlat, farfull);
    const int64_t nring = hex9_k_ring(farfull, LAYER, 1, ring, 8);
    CHECK(nring > 0, "seal ring");
    np = hex9_walk_to(fox, far, LAYER, ring, (size_t)nring, 400, path, 512);
    CHECK(np == 0, "sealed dest reachable (%lld)", (long long)np);

    /* ── true-bearing correction: roundtrip, invariants, magnitude ── */
    {
        const double bs[] = {0.0, 30.0, 45.0, 90.0, 135.0, 180.0, 250.0, 359.0};
        for (int i = 0; i < 8; ++i) {
            const double v = hex9_bearing_from_true(bs[i], 55.95);
            CHECK(fabs(hex9_bearing_to_true(v, 55.95) - bs[i]) < 1e-9,
                  "true-bearing roundtrip (b=%g)", bs[i]);
            CHECK(fabs(hex9_bearing_from_true(bs[i], 90.0) - bs[i]) < 1e-9,
                  "pole must be identity (b=%g)", bs[i]);
        }
        CHECK(hex9_bearing_from_true(0.0, 45.0) == 0.0 &&
              fabs(hex9_bearing_from_true(90.0, 45.0) - 90.0) < 1e-12 &&
              fabs(hex9_bearing_from_true(180.0, 45.0) - 180.0) < 1e-12 &&
              fabs(hex9_bearing_from_true(270.0, 45.0) - 270.0) < 1e-12,
              "cardinals must be identity");
        const double dev = 45.0 - hex9_bearing_from_true(45.0, 0.0);
        CHECK(fabs(dev - 0.1918) < 1e-3, "flattening magnitude off: %g", dev);
    }

    /* E4H guard */
    uint8_t e4[16];
    CHECK(hex9_e4h_encode(FLON, FLAT, 6, 2, e4) == 0, "e4h encode");
    CHECK(hex9_aim(e4, LAYER, 0.0, step) != 0, "aim accepted e4h");
    CHECK(hex9_walk_to(e4, fox, LAYER, NULL, 0, 0, path, 512) < 0,
          "walk accepted e4h");
    CHECK(hex9_vision_cone(e4, LAYER, 0, 90, 2, NULL, 0, cone2, maxd) < 0,
          "cone accepted e4h");

    free(disk); free(cone); free(cone2);
    if (fails) { printf("FAILED: %d checks\n", fails); return 1; }
    printf("OK: C grid verbs — aim inverse, cone laws + seam + occlusion, "
           "walk laws + detour + seal, guards\n");
    return 0;
}
