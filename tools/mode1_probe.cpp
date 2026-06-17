/* mode1_probe.cpp — empirically derive the mode-1-home slot/centroid map.
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2026, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * For grid points over a patch: NN-encode, run the deterministic
 * reversible-tail walk (the wip branch identity_from_uuid front end) to get
 * the leaf origin (ia, ib), slot c2 and leaf parent mode; ground truth is
 * the enumerated mesh cell whose centroid is nearest the point (interior
 * points only). For each (leaf_mode, c2) bucket, print the observed
 * (truth_centroid - origin) offsets — the table identity_from_uuid needs.
 *
 * Build: c++ -std=c++17 -O2 -Icore -DH9_WARP_ENABLE=1 tools/mode1_probe.cpp
 *        core/h9_warp_runtime.cpp core/h9_warp_embedded.cpp (via CMake dev target)
 */
#define H9_WARP_ENABLE 1
#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"
#include "h9_warp_runtime.h"
#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <string>

int main(int argc, char **argv) {
    const double lon0 = (argc > 1) ? atof(argv[1]) : 13.4053;
    const double lat0 = (argc > 2) ? atof(argv[2]) : 52.5209;
    const int    L    = (argc > 3) ? atoi(argv[3]) : 8;
    const int    N    = (argc > 4) ? atoi(argv[4]) : 160;

    std::string err;
    if (!h9::h9_warp_init_embedded(&err)) { fprintf(stderr, "warp: %s\n", err.c_str()); return 1; }

    /* ground truth: enumerate the mesh over the patch */
    const double span = 0.06;
    std::vector<H9GridCell> cells;
    h9grid::enumerate(L, lon0 - span, lat0 - span, lon0 + span, lat0 + span, cells);
    fprintf(stderr, "enumerated %zu cells\n", cells.size());
    if (cells.empty()) return 1;

    /* truth lattice positions: per cell, centroid in integer UV of its frame */
    /* (H9GridCell carries oid, c2, ia, ib — origin — plus geographic centroid) */

    /* observation buckets: (leaf_mode, c2_slot) -> map of (du, dv) -> count */
    std::map<std::pair<int,int>, std::map<std::pair<long,long>, long>> buckets;
    long used = 0, skipped = 0, mismatch_frame = 0;

    for (int iy = 0; iy < N; iy++) {
        for (int ix = 0; ix < N; ix++) {
            const double lon = lon0 - span * 0.8 + 2 * span * 0.8 * (ix + 0.5) / N;
            const double lat = lat0 - span * 0.8 + 2 * span * 0.8 * (iy + 0.5) / N;

            /* truth: nearest enumerated centroid, interior points only */
            int best = -1, second = -1;
            double bd = 1e30, sd = 1e30;
            for (size_t j = 0; j < cells.size(); j++) {
                const double dx = (cells[j].cen_lon - lon) * cos(lat * M_PI / 180.0);
                const double dy =  cells[j].cen_lat - lat;
                const double d2 = dx*dx + dy*dy;
                if (d2 < bd) { sd = bd; second = best; bd = d2; best = (int)j; }
                else if (d2 < sd) { sd = d2; second = (int)j; }
            }
            (void)second;
            if (best < 0 || bd > 0.25 * sd) { skipped++; continue; }  /* near an edge */

            /* encode + deterministic reversible-tail walk */
            H9BOct b = h9_lonlatdeg_to_boct(lon, lat);
            uint8_t full[16], nib[32];
            h9_boct_to_uuid(b, full);
            h9a_unpack(full, nib);
            if (nib[30] == 0x0Fu) { skipped++; continue; }
            const uint8_t r_mo = nib[31] & 1u;
            const int oid = (int)H9_L0HEX_BACK[nib[0]][r_mo][0];
            uint8_t c_mo = (nib[31] >> 3) & 1u;
            uint8_t c2w  = (nib[31] >> 1) & 3u;
            uint8_t rids[32] = {}; rids[0] = r_mo;
            int slot = (int)c2w, smode = (int)c_mo;
            bool ok = true;
            for (int l = 29; l >= 1; --l) {
                if (l == L) { slot = (int)c2w; smode = (int)c_mo; }
                if (nib[l] > 8) { ok = false; break; }
                const uint8_t *e = H9_HEX_REG[nib[l]][c_mo][c2w];
                if (e[0] == 0xFFu) { ok = false; break; }
                if (l <= L) rids[l] = e[0];
                c_mo = e[1]; c2w = e[2];
            }
            if (!ok) { skipped++; continue; }

            const int64_t s = h9kring::pow3(L);
            int64_t ia = 0, ib = 0, w = s / 3;
            for (int m = 1; m <= L; m++, w /= 3) {
                int ou, ov;
                if (!h9kring::cid_offset_int(H9_RID2CELL[rids[m]], &ou, &ov)) { ok = false; break; }
                ia += (int64_t)ou * w;
                ib += (int64_t)ov * w;
            }
            if (!ok) { skipped++; continue; }

            const H9GridCell &t = cells[best];
            if (t.oid != oid) { mismatch_frame++; continue; }
            const int64_t cu0 = t.ia + h9kring::H9KR_C2_DU[t.c2];
            const int64_t cv0 = t.ib + h9kring::H9KR_C2_DV[t.c2];
            buckets[{smode, slot}][{(long)(cu0 - ia), (long)(cv0 - ib)}]++;
            used++;
        }
    }

    printf("used %ld, skipped %ld, other-frame %ld\n", used, skipped, mismatch_frame);
    printf("(mode-0 truth offsets are DU={1,1,-2}, DV={1,-1,0})\n");
    for (auto &bk : buckets) {
        printf("leaf_mode=%d slot_c2=%d:\n", bk.first.first, bk.first.second);
        for (auto &ofs : bk.second)
            printf("    (du=%ld, dv=%ld) x%ld\n", ofs.first.first, ofs.first.second, ofs.second);
    }
    return 0;
}
