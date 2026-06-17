/* corner_trace.cpp — branch trace for the four polar F5 residuals.
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

int main() {
    std::string err;
    h9::h9_warp_init_embedded(&err);
    const int L = 15;
    const double pts[4][2] = { {180, 89.1}, {0, 89.1}, {0, -89.1}, {180, -89.1} };
    for (int p = 0; p < 4; p++) {
        H9BOct b = h9_lonlatdeg_to_boct(pts[p][0], pts[p][1]);
        uint8_t full[16];
        h9_boct_to_uuid(b, full);
        h9kring::H9CellId id;
        if (!h9kring::identity_from_uuid(full, L, &id)) { printf("no id\n"); continue; }
        const int64_t s = h9kring::pow3(L);
        printf("pt(%g,%g) oid=%d mo=%d ext=%d c2=%d  ia-s rel: ia=%lld ib=%lld (s=%lld)\n",
               pts[p][0], pts[p][1], id.oid, (int)H9_OID_MO[id.oid], (int)id.ext, id.c2,
               (long long)id.ia, (long long)id.ib, (long long)s);
        for (int v = 0; v < 6; v++) {
            int64_t u = (int64_t)H9_HI[0][id.c2][v][0] + id.ia;
            int64_t w = (int64_t)H9_HI[0][id.c2][v][1] + id.ib;
            const int mo = (int)H9_OID_MO[id.oid];
            const bool in0 = h9cell::uv_in_face(mo, s, u, w);
            int64_t ru = u, rv = w; int ro = id.oid;
            h9cell::resolve_uv_frame(s, &ru, &rv, &ro);
            /* boundary distances in own frame (mode-0 form: v-s, u-v-2s, u+v+2s) */
            long long d1 = (mo==0)? (long long)(w - s)     : (long long)(-s - w);
            long long d2 = (mo==0)? (long long)(u - w - 2*s): (long long)(u + w - 2*s);
            long long d3 = (mo==0)? (long long)(u + w + 2*s): (long long)(u - w + 2*s);
            printf("  v%d (%lld,%lld) in_own=%d  bdry(d1=%lld d2=%lld d3=%lld) -> oid %d (%lld,%lld)\n",
                   v, (long long)(u - id.ia), (long long)(w - id.ib), (int)in0,
                   d1, d2, d3, ro, (long long)(ru), (long long)(rv));
        }
    }
    return 0;
}
