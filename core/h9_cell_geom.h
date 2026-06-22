/* h9_cell_geom.h — identity-based cell geometry (the single source of truth).
 *
 * Renders a cell's hexagon ring and centroid from its lattice identity
 * (h9kring::H9CellId), using exactly the same per-vertex integer-UV +
 * frame-reflection construction as the grid enumerator's emit — so
 * h9_cell == h9_grid by construction, for every UUID flavour.
 *
 * WHY: the older decode path (cell_unpack's float backward walk) follows the
 * centroid-descent convention, while uuid_from_iauv (grid, encoder, k-ring,
 * labels) uses exact origin descent. For cells where the two descents
 * disagree, the old path rendered the WRONG hexagon for a canonical bin —
 * caught by the postgis_hex9 regression suite (grid_cell_agree) when the
 * extension moved into this repository.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"

namespace h9cell {

/* Integer-UV face membership for an octant frame at scale s. Mode-0 face
 * triangle: (3s, s), (-3s, s), (0, -2s); mode-1 is the v-mirror. */
static inline bool uv_in_face(int mo, int64_t s, int64_t u, int64_t v) {
    if (mo == 0) return !(v >  s || u - v >  2*s || u + v < -2*s);
    else         return !(v < -s || u + v >  2*s || u - v < -2*s);
}

/* Resolve an out-of-face vertex into the frame that owns it, via the
 * probe-validated seam maps (single edge, then two-map corner composition).
 * Ben's exact-on-seam ruling (2026-06-11): the mode-0 octant GOVERNS a
 * seam-straddling cell — its identity lives in the mode-0 frame and the
 * across-seam vertices are carried into the neighbouring frame here. The
 * index case is L0: all 12 L0 cells are centred on octahedron edges, and
 * 3x subdivision keeps a seam-centred chain at every layer. Without this,
 * out-of-face vertices fold at the seam (signed-bary projection uses
 * magnitudes), rendering garbled slivers (doctrine F5). */
static inline void resolve_uv_frame(int64_t s, int64_t *u, int64_t *v, int *oid) {
    const int mo = (int)H9_OID_MO[*oid];
    if (uv_in_face(mo, s, *u, *v)) return;
    for (int e = 0; e < 3; e++) {
        int64_t ru, rv;
        h9kring::seam_apply(e, mo, s, *u, *v, &ru, &rv);
        const int nb = H9_OID_NB[*oid][e];
        if (uv_in_face((int)H9_OID_MO[nb], s, ru, rv)) {
            *u = ru; *v = rv; *oid = nb;
            return;
        }
    }
    for (int e1 = 0; e1 < 3; e1++) {
        int64_t r1u, r1v;
        h9kring::seam_apply(e1, mo, s, *u, *v, &r1u, &r1v);
        const int n1 = H9_OID_NB[*oid][e1];
        for (int e2 = 0; e2 < 3; e2++) {
            int64_t r2u, r2v;
            h9kring::seam_apply(e2, (int)H9_OID_MO[n1], s, r1u, r1v, &r2u, &r2v);
            const int n2 = H9_OID_NB[n1][e2];
            if (uv_in_face((int)H9_OID_MO[n2], s, r2u, r2v)) {
                *u = r2u; *v = r2v; *oid = n2;
                return;
            }
        }
    }
    /* unresolved (should not happen within one cell of a face) — leave in
     * the own frame; the legacy fold applies. */
}

/* Per-vertex integer UV + owning frame for a cell identity. Mirrors
 * emit_hexes_for_oid: vertices live in the cell's own frame except v4/v5 of
 * ext (half-hex) cells, which reflect across the seam into
 * H9_OID_NB[oid][c2]'s frame; any other out-of-face vertex (seam-straddling
 * cells — the L0-descended chain) resolves via the seam maps. */
static inline void identity_vertices(const h9kring::H9CellId &id, int layer,
                                     int64_t pu[6], int64_t pv[6], int poid[6]) {
    for (int v = 0; v < 6; v++) {
        pu[v]   = (int64_t)H9_HI[0][id.c2][v][0] + id.ia;
        pv[v]   = (int64_t)H9_HI[0][id.c2][v][1] + id.ib;
        poid[v] = id.oid;
    }
    if (id.ext) {
        const int nb = H9_OID_NB[id.oid][id.c2];
        /* v4 ← reflect(v2);  v5 ← reflect(v1)  across v = 0 */
        pu[4] = (int64_t)H9_HI[0][id.c2][2][0] + id.ia;
        pv[4] = -((int64_t)H9_HI[0][id.c2][2][1] + id.ib);
        poid[4] = nb;
        pu[5] = (int64_t)H9_HI[0][id.c2][1][0] + id.ia;
        pv[5] = -((int64_t)H9_HI[0][id.c2][1][1] + id.ib);
        poid[5] = nb;
        return;
    }
    const int64_t s = h9kring::pow3(layer);
    for (int v = 0; v < 6; v++)
        resolve_uv_frame(s, &pu[v], &pv[v], &poid[v]);
}

/* Geographic centroid, grid convention: the lattice centroid (== the
 * mode-matching 6-vertex mean) for whole cells; the 4-own-vertex mean for
 * ext half-hexes. Matches H9GridCell::cen_lon/cen_lat exactly. */
static inline bool identity_centroid(const h9kring::H9CellId &id, int layer,
                                     double *lon, double *lat) {
    const double div_f = std::pow(3.0, (double)layer);
    double cx, cy;
    if (!id.ext) {
        cx = (double)(id.ia + h9kring::H9KR_C2_DU[id.c2]) * H9_UV_U1 / div_f;
        cy = (double)(id.ib + h9kring::H9KR_C2_DV[id.c2]) * H9_UV_V3 / div_f;
    } else {
        int64_t su = 0, sv = 0;
        for (int v = 0; v < 4; v++) {
            su += (int64_t)H9_HI[0][id.c2][v][0] + id.ia;
            sv += (int64_t)H9_HI[0][id.c2][v][1] + id.ib;
        }
        cx = (double)su * H9_UV_U1 / (4.0 * div_f);
        cy = (double)sv * H9_UV_V3 / (4.0 * div_f);
    }
    return h9grid::cxcy_to_lonlat(cx, cy, id.oid, lon, lat);
}

/* Antimeridian/seam normalisation + exact close, in place. */
static inline void normalize_ring(double *lons, double *lats, int n_ring) {
    for (int i = 1; i < n_ring - 1; ++i) {
        while (lons[i] - lons[i-1] >  180.0) lons[i] -= 360.0;
        while (lons[i] - lons[i-1] < -180.0) lons[i] += 360.0;
    }
    double mean = 0.0;
    for (int i = 0; i < n_ring - 1; ++i) mean += lons[i];
    mean /= (double)(n_ring - 1);
    if      (mean >  180.0) for (int i = 0; i < n_ring; ++i) lons[i] -= 360.0;
    else if (mean < -180.0) for (int i = 0; i < n_ring; ++i) lons[i] += 360.0;
    lons[n_ring - 1] = lons[0]; lats[n_ring - 1] = lats[0];
}

/* Closed lon/lat hexagon ring at `densify` (6·3^densify + 1 points) for a
 * cell identity. densify > 0 interpolates each edge in integer-consistent
 * frames with per-vertex octant-seam reflection — the same construction as
 * hex9_grid_cell_ring. Returns the point count, or -1 on error. */
static inline int identity_ring(const h9kring::H9CellId &id, int layer, int densify,
                                double *lons, double *lats) {
    if (layer < 0 || layer > H9_LMAX) return -1;
    if (densify < 0 || densify > 9 || layer + densify > H9_LMAX) return -1;
    int n_ring = 1;
    for (int i = 0; i < densify; ++i) n_ring *= 3;
    const int v_per_edge = n_ring;
    n_ring = 6 * n_ring + 1;

    int64_t pu[6], pv[6];
    int     poid[6];
    identity_vertices(id, layer, pu, pv, poid);

    const double div_f = std::pow(3.0, (double)layer);
    if (densify == 0) {
        /* project each vertex in its own frame — bit-identical to the grid
         * enumerator's stored vlon/vlat */
        for (int v = 0; v < 6; ++v) {
            const double cx = (double)pu[v] * H9_UV_U1 / div_f;
            const double cy = (double)pv[v] * H9_UV_V3 / div_f;
            h9grid::cxcy_to_lonlat(cx, cy, poid[v], &lons[v], &lats[v], nullptr, false);
        }
        lons[6] = lons[0]; lats[6] = lats[0];
        normalize_ring(lons, lats, 7);
        return 7;
    }
    int n = 0;
    for (int e = 0; e < 6; ++e) {
        const int en = (e + 1) % 6;
        const int oa = poid[e], ob = poid[en];
        const int mo_a = (int)H9_OID_MO[oa], mo_b = (int)H9_OID_MO[ob];
        int64_t ua = pu[e], va = pv[e], ub = pu[en], vb = pv[en];
        int frame_oid;
        if (oa == ob || mo_a == mo_b) frame_oid = oa;
        else if (mo_a == 0) { frame_oid = ob; va = -va; }
        else                { frame_oid = oa; vb = -vb; }
        const double cxa = (double)ua * H9_UV_U1 / div_f, cya = (double)va * H9_UV_V3 / div_f;
        const double cxb = (double)ub * H9_UV_U1 / div_f, cyb = (double)vb * H9_UV_V3 / div_f;
        for (int s = 0; s < v_per_edge; ++s) {
            const double t = (double)s / (double)v_per_edge;
            h9grid::cxcy_to_lonlat(cxa + (cxb - cxa) * t, cya + (cyb - cya) * t,
                                   frame_oid, &lons[n], &lats[n], nullptr, false);
            ++n;
        }
    }
    lons[n_ring - 1] = lons[0]; lats[n_ring - 1] = lats[0];
    normalize_ring(lons, lats, n_ring);
    return n_ring;
}

}  /* namespace h9cell */
