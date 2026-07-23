/* cell_uv.c — the integer lattice identity surface (hex9_cell_uv, 2.1.x).
 *
 * Pins the three claims the header makes:
 *   1. PARITY: the centre key converted by pure arithmetic
 *      (ia·u1/3^L, ib·v3/3^L) equals project(decode(bin)) — the value the
 *      geometric chain recovers through the warp and a Newton solve — to
 *      Newton tolerance, same octant. Ring keys likewise, via unproject
 *      against hex9_cell_ring's lon/lat (order-aligned).
 *   2. CANONICAL KEYS: every edge-adjacent neighbour shares EXACTLY two
 *      vertex keys by integer equality — including across octant seams
 *      (the exact-on-seam resolution is per-point, so both sides report
 *      the same key). This is what makes shared-vertex pooling work.
 *   3. EXT CENSUS: exactly 12·3^L seam-chain cells per layer (the
 *      L0-descended octant-edge chains), flagged by `ext`.
 *
 * All of it is datum-free by construction (upstream of the reduction); the
 * geometric comparisons here use the WGS84 chain arbitrarily.
 */
#include "hex9_c.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fails = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); g_fails++; } \
} while (0)

static const double PTS[][2] = {
    { -3.188267,  55.953251 },   /* edinburgh */
    { -0.124625,  51.500745 },   /* westminster */
    { 151.215256, -33.856159 },  /* sydney */
    {  17.25,      0.0      },   /* equator */
    { -89.9863,    0.0091   },   /* near the (-90, 0) octahedron vertex */
    {   0.02,     89.97     },   /* near the north pole */
};
static const int N_PTS = (int)(sizeof PTS / sizeof PTS[0]);
static const int LAYERS[] = { 1, 3, 5, 8, 12 };
static const int N_LAY = (int)(sizeof LAYERS / sizeof LAYERS[0]);

static double pow3(int L) { double d = 1.0; while (L--) d *= 3.0; return d; }

int main(void)
{
    char err[256];
    CHECK(hex9_init(err, sizeof err) == 0, "hex9_init: %s", err);

    double u1, v3;
    hex9_uv_units(&u1, &v3);
    CHECK(u1 > 0 && v3 > 0, "uv units");

    /* ── 1. parity with the geometric chain ──────────────────────────── */
    for (int p = 0; p < N_PTS; ++p) {
        uint8_t full[16];
        hex9_encode(PTS[p][0], PTS[p][1], full);
        for (int li = 0; li < N_LAY; ++li) {
            const int L = LAYERS[li];
            const double div = pow3(L);
            uint8_t b[16];
            hex9_bin(full, L, b);

            int64_t cia, cib, via[6], vib[6];
            int coid, void_[6], ext;
            CHECK(hex9_cell_uv(b, L, &cia, &cib, &coid, via, vib, void_,
                               &ext) == 0, "cell_uv [%d L%d]", p, L);

            /* centre: arithmetic vs decode+project (skip ext cells — their
             * decode centroid is the 4-vertex mean, deliberately NOT the
             * lattice centre key; the header documents the distinction). */
            if (!ext) {
                double lon, lat, gx, gy;
                int goid;
                hex9_decode(b, &lon, &lat);
                hex9_project(lon, lat, &gx, &gy, &goid);
                CHECK(goid == coid, "centre oid [%d L%d]: %d vs %d",
                      p, L, goid, coid);
                CHECK(fabs(gx - (double)cia * u1 / div) < 1e-8 &&
                      fabs(gy - (double)cib * v3 / div) < 1e-8,
                      "centre parity [%d L%d]: d=(%.2e, %.2e)", p, L,
                      gx - (double)cia * u1 / div,
                      gy - (double)cib * v3 / div);
            }

            /* ring: unproject each key, compare to hex9_cell_ring */
            double ring[2 * 7];
            CHECK(hex9_cell_ring(b, L, 0, ring, 7) == 7, "ring [%d L%d]", p, L);
            for (int v = 0; v < 6; ++v) {
                double lon, lat;
                CHECK(hex9_unproject((double)via[v] * u1 / div,
                                     (double)vib[v] * v3 / div,
                                     void_[v], &lon, &lat) == 0,
                      "unproject vert [%d L%d v%d]", p, L, v);
                double dlon = fabs(lon - ring[2 * v]);
                if (dlon > 180.0) dlon = fabs(dlon - 360.0);
                /* poles: longitude is degenerate there */
                if (fabs(ring[2 * v + 1]) > 89.999999) dlon = 0.0;
                CHECK(dlon * cos(ring[2 * v + 1] * M_PI / 180.0) < 1e-7 &&
                      fabs(lat - ring[2 * v + 1]) < 1e-7,
                      "vert parity [%d L%d v%d]: (%.9f,%.9f) vs (%.9f,%.9f)",
                      p, L, v, lon, lat, ring[2 * v], ring[2 * v + 1]);
            }
        }
    }

    /* ── 2. canonical keys: neighbours share exactly two ─────────────── */
    for (int p = 0; p < N_PTS; ++p) {
        uint8_t full[16];
        hex9_encode(PTS[p][0], PTS[p][1], full);
        const int L = 6;
        int64_t cia, cib, via[6], vib[6];
        int coid, void_[6], ext;
        uint8_t b[16];
        hex9_bin(full, L, b);
        hex9_cell_uv(b, L, &cia, &cib, &coid, via, vib, void_, &ext);

        uint8_t nbs[6 * 16];
        const int nn = hex9_neighbors(full, L, nbs);
        CHECK(nn == 5 || nn == 6, "neighbor count [%d]", p);
        int triples = 0;
        for (int k = 0; k < nn; ++k) {
            int64_t nia, nib, wia[6], wib[6];
            int noid, woid[6], next;
            CHECK(hex9_cell_uv(nbs + k * 16, L, &nia, &nib, &noid,
                               wia, wib, woid, &next) == 0,
                  "neighbor cell_uv [%d k%d]", p, k);
            int shared = 0;
            for (int a = 0; a < 6; ++a)
                for (int w = 0; w < 6; ++w)
                    if (via[a] == wia[w] && vib[a] == wib[w] &&
                        void_[a] == woid[w])
                        shared++;
            /* edge-adjacency = 2 shared keys; a 5-neighbour half-hex at an
             * octahedron cone point meets ONE neighbour across two edges —
             * 3 shared keys, and only there. */
            CHECK(shared == 2 || (shared == 3 && nn == 5),
                  "shared keys [%d k%d]: %d (canonicality broken?)",
                  p, k, shared);
            if (shared == 3) triples++;
        }
        CHECK(nn == 6 ? triples == 0 : triples == 1,
              "cone-point doubled neighbor [%d]: %d triples (nn=%d)",
              p, triples, nn);
    }

    /* ── 3. ext census at L2 via batch (also pins _many == scalar) ───── */
    {
        enum { NC = 972 };
        /* enumerate all L2 cells by binning a dense sample is clumsy in C;
         * walk the 12 L0 cells' children twice instead. */
        uint8_t l0[16 * 16], kids1[12 * 9 * 16], kids2[12 * 81 * 16];
        (void)l0;
        /* collect L2 = children(children(L0)); L0 cells from binning the
         * six probes is incomplete, so derive L0 set from k-disk of one
         * cell is also clumsy — simplest: children from each of the 12 L0
         * bins obtained via hex9_bin of encodes on a coarse lattice. */
        int have = 0;
        uint8_t seen[12][16];
        for (double lat = -80.0; lat <= 80.0 && have < 12; lat += 10.0)
            for (double lon = -175.0; lon < 180.0 && have < 12; lon += 10.0) {
                uint8_t f[16], b0[16];
                hex9_encode(lon, lat, f);
                hex9_bin(f, 0, b0);
                int dup = 0;
                for (int i = 0; i < have; ++i)
                    if (memcmp(seen[i], b0, 16) == 0) { dup = 1; break; }
                if (!dup) memcpy(seen[have++], b0, 16);
            }
        CHECK(have == 12, "found %d L0 cells", have);
        int n2 = 0;
        for (int i = 0; i < have; ++i) {
            hex9_cell_children(seen[i], kids1);
            for (int c = 0; c < 9; ++c) {
                hex9_cell_children(kids1 + c * 16, kids2 + (size_t)n2 * 16);
                n2 += 9;
            }
        }
        CHECK(n2 == NC, "L2 count %d", n2);

        int64_t cia[NC], cib[NC], via[NC * 6], vib[NC * 6];
        int32_t coid[NC], void_[NC * 6], ext[NC], lay[NC];
        for (int i = 0; i < NC; ++i) lay[i] = 2;
        CHECK(hex9_cell_uv_many(kids2, lay, NC, cia, cib, coid,
                                via, vib, void_, ext) == 0, "cell_uv_many");
        int next = 0;
        for (int i = 0; i < NC; ++i) next += ext[i];
        CHECK(next == 108, "ext census at L2: %d (want 12*3^2 = 108)", next);

        /* batch == scalar spot check */
        int64_t s_cia, s_cib, s_via[6], s_vib[6];
        int s_coid, s_void[6], s_ext;
        hex9_cell_uv(kids2 + 500 * 16, 2, &s_cia, &s_cib, &s_coid,
                     s_via, s_vib, s_void, &s_ext);
        CHECK(s_cia == cia[500] && s_cib == cib[500] &&
              s_coid == (int)coid[500] && s_ext == (int)ext[500],
              "many == scalar");
    }

    if (g_fails) { printf("%d FAILURE(S)\n", g_fails); return 1; }
    printf("cell_uv: all assertions passed\n");
    return 0;
}
