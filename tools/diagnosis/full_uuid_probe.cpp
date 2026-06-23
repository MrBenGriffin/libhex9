/* full_uuid_probe.cpp — claim-1 measurement + full-uuid-from-identity check.
 *
 * For each grid cell at several layers, build the FULL (layer-29) UUID by
 * continuing the exact integer ORIGIN descent to depth 29 (the design-doc
 * prescription: cids_from_iauv at the scaled origin (ia·3^(29−L),
 * ib·3^(29−L)), depth 29), then measure two things:
 *
 *   1. h9_bin(full, L) == bin  — is the origin-descent full UUID a VALID
 *      address for this cell? (the acceptance test for Part 2)
 *   2. geodesic( decode(full), geometric centroid )  — how far does the
 *      full UUID's decoded point sit from the cell centroid, absolute and
 *      as a fraction of the cell circumradius? (claim 1: "corner, not µm")
 *
 * Build (from repo root, after cmake --build build):
 *   c++ -std=c++17 -O2 -Icore -I. tools/full_uuid_probe.cpp \
 *       build/libhex9.a -o build/full_uuid_probe   # or link the objects
 * Simplest: add as an EXCLUDE_FROM_ALL target (see CMakeLists) and
 *   cmake --build build --target full_uuid_probe && ./build/full_uuid_probe
 */
#include "hex9_c.h"
#include "h9_grid.h"

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <set>
#include <algorithm>

static double haversine_m(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371008.8;            /* mean Earth radius (m) */
    const double d2r = M_PI / 180.0;
    const double dlat = (lat2 - lat1) * d2r;
    const double dlon = (lon2 - lon1) * d2r;
    const double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
                     std::cos(lat1 * d2r) * std::cos(lat2 * d2r) *
                     std::sin(dlon / 2) * std::sin(dlon / 2);
    return 2.0 * R * std::asin(std::min(1.0, std::sqrt(a)));
}

static int64_t pow3(int e) { int64_t p = 1; for (int i = 0; i < e; i++) p *= 3; return p; }

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) {
        std::printf("warp init FAILED: %s\n", err);
        return 1;
    }

    /* Regions stressing the descent: a benign mid-latitude box, a box
     * straddling the lon=0 octant seam at the equator, and a near-polar box.
     * Each window shrinks 3×/layer so cell count stays ~comparable. */
    struct Region { const char *name; double lonc, latc; };
    const Region regions[] = {
        { "mid-lat   ", -0.10, 51.50 },
        { "equator/0 ",  0.00,  0.00 },
        { "near-pole ",  0.00, 88.50 },
    };
    const int layers[] = { 8, 11, 13, 16, 20 };

    std::printf("C = encode(interior)   D = encoder-free uuid_from_cxcy(interior)\n");
    std::printf("interior = centroid nudged 10%% toward mode-0 corner; "
                "bin = h9_bin(full,L)==cell bin\n");

    for (int ri = 0; ri < (int)(sizeof(regions)/sizeof(regions[0])); ri++) {
      const double lonc = regions[ri].lonc, latc = regions[ri].latc;
      std::printf("\n[%s] centre (%.2f, %.2f)\n", regions[ri].name, lonc, latc);
      for (int li = 0; li < (int)(sizeof(layers) / sizeof(layers[0])); li++) {
        const int L = layers[li];
        const double half = 0.20 * std::pow(3.0, 8 - L);   /* ~1255 cells/layer */
        std::vector<H9GridCell> cells;
        h9grid::enumerate(L, lonc - half, latc - half, lonc + half, latc + half, cells);
        if (cells.empty()) { std::printf("L%-3d (no cells)\n", L); continue; }

        const int64_t s = pow3(29 - L);
        const double div_f = std::pow(3.0, (double)L);
        const double ALPHA = 0.10;        /* proven sufficient interior nudge */
        long c_ok = 0, d_ok = 0;
        std::set<std::array<uint8_t,16>> d_uniq;
        std::vector<double> d_off;

        for (const auto &c : cells) {
            uint8_t bb[16];
            const int face_mode = (int)H9_OID_MO[c.oid];

            /* Mode-0 origin corner = decode of the origin-descent UUID. The
             * half-hex seam runs through the centroid; moving toward this
             * corner enters the mode-0 half. */
            uint8_t fa[16];
            h9grid::uuid_from_iauv(c.oid, c.c2, c.ia * s, c.ib * s, 29, c.ext, fa);
            double cornlon, cornlat; hex9_decode(fa, &cornlon, &cornlat);

            /* Candidate C — encode the interior point (lon/lat). Proven valid;
             * the slow reference. */
            const double clon = c.cen_lon + ALPHA * (cornlon - c.cen_lon);
            const double clat = c.cen_lat + ALPHA * (cornlat - c.cen_lat);
            uint8_t fc[16];
            if (hex9_encode(clon, clat, fc) == 0) {
                hex9_bin(fc, L, bb);
                if (std::memcmp(bb, c.uuid, 16) == 0) c_ok++;
            }

            /* Candidate D — ENCODER-FREE: build the same interior point in the
             * grid's own (cx,cy) frame and descend via uuid_from_cxcy. */
            double sum_cx = 0, sum_cy = 0; int nm = 0;
            for (int v = 0; v < 6; v++) {
                if ((int)H9_OID_MO[c.poid[v]] != face_mode) continue;
                sum_cx += (double)c.pu[v] * H9_UV_U1 / div_f;
                sum_cy += (double)c.pv[v] * H9_UV_V3 / div_f;
                nm++;
            }
            const double cen_cx = sum_cx / nm, cen_cy = sum_cy / nm;
            const double org_cx = (double)c.ia * H9_UV_U1 / div_f;
            const double org_cy = (double)c.ib * H9_UV_V3 / div_f;
            const double in_cx = cen_cx + ALPHA * (org_cx - cen_cx);
            const double in_cy = cen_cy + ALPHA * (org_cy - cen_cy);
            /* Raw lattice (cx,cy) → b_oct frame (warp_fwd), the frame
             * uuid_from_cxcy expects (it applies warp_inv internally). */
            double bx, by;
            h9_warp_fwd(in_cx, in_cy, face_mode, &bx, &by);
            uint8_t fd[16];
            h9grid::uuid_from_cxcy(bx, by, c.oid, 29, fd);
            hex9_bin(fd, L, bb);
            if (std::memcmp(bb, c.uuid, 16) == 0) d_ok++;
            { std::array<uint8_t,16> k; std::memcpy(k.data(), fd, 16); d_uniq.insert(k); }
            double dlon, dlat; hex9_decode(fd, &dlon, &dlat);
            d_off.push_back(haversine_m(dlon, dlat, c.cen_lon, c.cen_lat));
        }

        std::sort(d_off.begin(), d_off.end());
        const double d_off_med = d_off.empty() ? 0 : d_off[d_off.size()/2];
        std::printf("L%-3d n=%4zu | C(encode) %4ld/%zu | D(cxcy,free) %4ld/%zu "
                    "uniq=%zu off=%.3fm\n",
                    L, cells.size(), c_ok, cells.size(),
                    d_ok, cells.size(), d_uniq.size(), d_off_med);
      }
    }
    return 0;
}
