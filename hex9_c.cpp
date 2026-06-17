/* hex9_c.cpp — C-ABI implementation for libhex9.
 *
 * Thin extern "C" wrappers over the C++ core (core/h9_*.h). Geometry-free:
 * lon/lat degrees + 16-byte UUIDs in/out. See hex9_c.h for the contract.
 *
 * Include order mirrors the proven order from the PostGIS glue
 * (lwgeom_hex9.cpp) so the warp shim resolves for h9_grid.h's cxcy path.
 */
#define H9_WARP_ENABLE 1
#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"
#include "h9_cell_geom.h"
#include "h9_warp_runtime.h"

#include "hex9_c.h"

#include <cstring>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

/* Encoder mode: 1 = containment (grid-canonical, DEFAULT), 0 = NN (legacy).
 * Containment is canonical: encode(point) then h9_bin matches h9_grid / h9_id /
 * Python at EVERY layer including the L29 leaf; the NN beam encoder agrees only
 * down to the leaf tail (they are bin-identical through ~L20, differ in ~3.7%
 * of full UUIDs below — measured, tools/encoder_ab_probe). The PostGIS
 * extension fixes this at the default (no GUC) so h9_encode stays truly
 * IMMUTABLE; hex9_set_encoder exists only for embedders / A-B research (CLI,
 * Python, probes). Read-only during a batch loop, so OpenMP-thread-safe. */
static int g_encoder_mode = 1;

/* Shared scalar kernels — the single source of truth for one item, called by
 * both the scalar API and the batch loops (so batch == scalar by construction).
 * Must stay reentrant: locals only, read-only warp state + g_encoder_mode. */
static inline void encode_one(double lon, double lat, uint8_t out[16]) {
    if (g_encoder_mode == 1) {
        double cx, cy; int oid;
        h9grid::lonlatdeg_to_cxcy_oid(lon, lat, &cx, &cy, &oid);
        h9grid::uuid_from_cxcy_full(cx, cy, oid, out);
    } else {
        H9BOct b = h9_lonlatdeg_to_boct(lon, lat);
        h9_boct_to_uuid(b, out);
    }
}
static inline void decode_one(const uint8_t uuid[16], double *lon, double *lat) {
    /* Bin UUIDs decode to the cell's geographic centroid via the exact
     * identity path (grid convention — coherent with h9_grid/h9_cell/labels
     * for every encoding flavour). Full UUIDs keep the beam backward walk
     * (the representative point of the L29 address). */
    uint8_t nib[32];
    h9a_unpack(uuid, nib);
    if (nib[30] == 0x0Fu) {
        int bl = 30;
        while (bl >= 0 && nib[bl] == 0x0Fu) bl--;
        h9kring::H9CellId id;
        if (bl >= 0 && h9kring::identity_from_uuid(uuid, bl, &id) &&
            h9cell::identity_centroid(id, bl, lon, lat))
            return;
        /* fall through to the boct walk for malformed bins */
    }
    H9BOct b = h9_uuid_to_boct(uuid);
    h9_boct_to_lonlatdeg(b, lon, lat);
}

extern "C" const char *hex9_version(void) {
    return "libhex9 0.1.0 (" __DATE__ " " __TIME__ ")";
}

extern "C" int hex9_warp_init(char *errbuf, size_t errlen) {
    std::string err;
    const bool ok = h9::h9_warp_init_embedded(&err);
    if (!ok && errbuf && errlen) {
        std::strncpy(errbuf, err.c_str(), errlen - 1);
        errbuf[errlen - 1] = '\0';
    }
    return ok ? 0 : 1;
}

extern "C" void hex9_set_use_warp(int on) { h9::g_warp_use = (on != 0); }
extern "C" void hex9_set_encoder(int mode) { g_encoder_mode = mode; }

extern "C" int hex9_encode(double lon, double lat, uint8_t out_uuid[16]) {
    encode_one(lon, lat, out_uuid);
    return 0;
}

extern "C" int hex9_decode(const uint8_t uuid[16], double *lon, double *lat) {
    decode_one(uuid, lon, lat);
    return 0;
}

/* Canonical (grid-matching) bin at `layer` — the identity coarsening h9_grid
 * enumerates with, which resolves a split-hex (6/7/8) leaf to its canonical
 * MODE-0 parent (the cell that geometrically contains the point). h9_bin_uuid
 * alone trusts the address's reversible MODE-1 nibble ancestry, so it bins a
 * split-hex point into a cell that does NOT contain it (~20% of split-hex
 * cells; grid-verified, e.g. L3 5038 vs the containing 5058). Going through
 * identity makes hex9_bin == hex9_k_disk(.,0) == h9_grid everywhere — the
 * bin's defining property (every point inside the cell maps to it uniquely)
 * holds geometrically. The full UUID stays the reversible load-bearer; the bin
 * is derived from it. L0 has no split ambiguity; an identity failure
 * (malformed input, or a bin fed below its own layer) falls back to the
 * legacy nibble bin. */
static void canonical_bin(const uint8_t uuid[16], int layer, uint8_t out[16]) {
    if (layer >= 1) {
        h9kring::H9CellId id;
        if (h9kring::identity_from_uuid(uuid, layer, &id)) {
            h9kring::identity_to_uuid(id, layer, out);
            return;
        }
    }
    h9_bin_uuid(uuid, layer, out);
}

extern "C" int hex9_bin(const uint8_t uuid[16], int layer, uint8_t out_uuid[16]) {
    if (layer < 0 || layer > 29) return 1;
    canonical_bin(uuid, layer, out_uuid);
    return 0;
}

/* ── Batch (OpenMP) ─────────────────────────────────────────────────────────
 * Signed loop index for OpenMP portability (older specs require it). */
extern "C" int hex9_encode_many(const double *lon, const double *lat, size_t n,
                                uint8_t *out_uuid) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        encode_one(lon[i], lat[i], out_uuid + (size_t)i * 16);
    return 0;
}

extern "C" int hex9_decode_many(const uint8_t *uuid, size_t n,
                                double *lon, double *lat) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        decode_one(uuid + (size_t)i * 16, &lon[i], &lat[i]);
    return 0;
}

extern "C" int hex9_bin_many(const uint8_t *uuid, int layer, size_t n,
                             uint8_t *out_uuid) {
    if (layer < 0 || layer > 29) return 1;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        canonical_bin(uuid + (size_t)i * 16, layer, out_uuid + (size_t)i * 16);
    return 0;
}

/* ── Continuous projection (b_oct backend) ──────────────────────────────────
 * Thin surface over the grid-canonical forward kernel (the same one encode_one
 * uses in mode 1). No descent, no L29 cap — see hex9_c.h / boct_backend_design. */
extern "C" int hex9_project(double lon, double lat,
                            double *cx, double *cy, int *oid) {
    h9grid::lonlatdeg_to_boct_oid(lon, lat, cx, cy, oid);
    return 0;
}

extern "C" int hex9_project_many(const double *lon, const double *lat, size_t n,
                                 double *cx, double *cy, int *oid) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::lonlatdeg_to_boct_oid(lon[i], lat[i], &cx[i], &cy[i], &oid[i]);
    return 0;
}

extern "C" int hex9_unproject(double cx, double cy, int oid,
                              double *lon, double *lat) {
    return h9grid::boct_oid_to_lonlat(cx, cy, oid, lon, lat) ? 0 : 1;
}

extern "C" int hex9_unproject_many(const double *cx, const double *cy,
                                   const int *oid, size_t n,
                                   double *lon, double *lat) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::boct_oid_to_lonlat(cx[i], cy[i], oid[i], &lon[i], &lat[i]);
    return 0;
}

extern "C" int hex9_ring_npoints(int densify) {
    if (densify < 0 || densify > 9) return -1;
    int n = 1;
    for (int i = 0; i < densify; ++i) n *= 3;
    return 6 * n + 1;
}

/* ── Labels ─────────────────────────────────────────────────────────────────
 * h9_write_label ported from the PostGIS glue (lwgeom_hex9.cpp): body nibbles
 * 0..layer mapped to digits 0-9,a,b. Writes layer+1 chars (no NUL); the public
 * wrappers NUL-terminate. */
static int write_label(const uint8_t nibbles[32], int layer, char *buf) {
    static const char H9D[] = "0123456789ab";
    for (int i = 0; i <= layer; ++i) buf[i] = H9D[nibbles[i] & 0x0Fu];
    return layer + 1;
}

/* Labels name the CANONICAL bin (== hex9_bin == h9_grid), spelled for humans.
 * They are cosmetic, not load-bearing: the body shows the full-hex digit and
 * drops the half-hex/mode, so canonical labels are NOT cross-layer
 * prefix-nested at split-hex (6/7/8) — common-prefix locality lives in the
 * full UUID (the load-bearer) and is a downstream concern. The full UUID is
 * decoded to its canonical bin first so label(full) == label(h9_bin(full)). */
extern "C" int hex9_label(const uint8_t uuid[16], int layer, char *buf, size_t buflen) {
    if (layer < 0 || layer > 29 || !buf) return -1;
    uint8_t cb[16], nibbles[32];
    canonical_bin(uuid, layer, cb);
    h9a_unpack(cb, nibbles);
    int len = write_label(nibbles, layer, buf);
    if (len < 0 || (size_t)len + 1 > buflen) return -1;
    buf[len] = '\0';
    return len;
}

extern "C" int hex9_label_key(const uint8_t uuid[16], int layer, char *buf, size_t buflen) {
    if (layer < 0 || layer > 29 || !buf) return -1;
    /* Canonical bin: it is a bin (nibble[30]==0xF), so nibble[31] already holds
     * the layer's key tail (c2<<1 | oct_mode) — no backward walk, and no c_mo
     * display bit (the canonical mode-0 home is unambiguous, so GIS exports no
     * longer need it to pick a half-hex template). parse(label_key) rebuilds
     * this bin exactly, so parse(label_key(u,L)) == hex9_bin(u,L). */
    uint8_t cb[16], nibbles[32];
    canonical_bin(uuid, layer, cb);
    h9a_unpack(cb, nibbles);
    int len = write_label(nibbles, layer, buf);
    if (len < 0) return -1;
    const uint8_t key_tail = nibbles[31];
    static const char H9D[] = "0123456789abcdef";
    if ((size_t)len + 3 > buflen) return -1;   /* '.' + nibble + NUL */
    buf[len++] = '.';
    buf[len++] = H9D[key_tail & 0x0Fu];
    buf[len]   = '\0';
    return len;
}

/* ── Label parsing / mesh prefix ops ────────────────────────────────────────
 * See hex9_c.h. Bare labels are recovered by candidate-tail search verified
 * against the canonical re-encode (exact for canonical labels — the only
 * kind hex9_label emits for grid/k-ring/containment-encoder UUIDs); keyed
 * labels (".<k>") decode directly and are flavour-blind. */

/* Parse the body chars into nibbles; returns layer or -1. */
static int label_body_to_nibbles(const char *label, uint8_t nibbles[32],
                                 const char **tail_out) {
    static const char H9D[] = "0123456789ab";
    int len = 0;
    while (label[len] && label[len] != '.') len++;
    if (len < 1 || len > 30) return -1;
    for (int i = 0; i < len; i++) {
        const char *p = (const char *)std::memchr(H9D, label[i], 12);
        if (!p) return -1;
        const int d = (int)(p - H9D);
        if (i > 0 && d > 8) return -1;        /* layers 1+ use digits 0..8 */
        nibbles[i] = (uint8_t)d;
    }
    for (int i = len; i <= 30; i++) nibbles[i] = 0x0Fu;
    *tail_out = (label[len] == '.') ? label + len + 1 : nullptr;
    return len - 1;
}

extern "C" int hex9_parse_label(const char *label, uint8_t out_uuid[16]) {
    if (!label || !out_uuid) return -1;
    uint8_t nib[32];
    const char *tail;
    const int layer = label_body_to_nibbles(label, nib, &tail);
    if (layer < 0) return -1;

    if (tail) {                                /* keyed: ".<hex nibble>" */
        /* A keyed label IS the bin in text form: parse is the exact
         * SYNTACTIC inverse of hex9_label_key — rebuild body + key nibble,
         * NO resolution (Python h9_from_label semantics). This makes
         * parse(label_key(u, L)) == hex9_bin(u, L) an identity. (The former
         * canonicalising parse routed through the bin resolution machinery,
         * which is approximate at meta-bearing cells — the doctrine
         * fossil.)
         *
         * The key carries one bit MORE than the bin tail: key =
         * (c_mo<<3)|(c2<<1)|r_mo (c_mo lets GIS exports pick the hexagon
         * template without resolution), while bin uuids store
         * (c2<<1)|r_mo. Strip bit 3 here or the identity silently fails
         * for every c_mo=1 cell — first seen 2026-06-12 when the F6 field
         * shift landed Westminster L8 on key '.b' (bin tail 3). */
        static const char H9D[] = "0123456789abcdef";
        const char *p = (const char *)std::memchr(H9D, tail[0], 16);
        if (!p || tail[1]) return -1;
        nib[31] = (uint8_t)((p - H9D) & 0x7);
        h9a_pack(nib, out_uuid);
        return layer;
    }

    /* bare: search the 6 possible tails, require the canonical re-encode to
     * reproduce the body (rejects corrupt and non-canonical labels). */
    int found = 0;
    uint8_t result[16];
    for (int r_mo = 0; r_mo < 2; r_mo++) {
        for (int c2 = 0; c2 < 3; c2++) {
            nib[31] = (uint8_t)((c2 << 1) | r_mo);
            uint8_t bin[16];
            h9a_pack(nib, bin);
            h9kring::H9CellId id;
            if (!h9kring::identity_from_uuid(bin, layer, &id)) continue;
            uint8_t canon[16], cnib[32];
            h9kring::identity_to_uuid(id, layer, canon);
            h9a_unpack(canon, cnib);
            if (std::memcmp(cnib, nib, (size_t)layer + 1) != 0) continue;
            if (found && std::memcmp(result, canon, 16) != 0) return -1;  /* ambiguous */
            std::memcpy(result, canon, 16);
            found = 1;
        }
    }
    if (!found) return -1;
    std::memcpy(out_uuid, result, 16);
    return layer;
}

extern "C" int hex9_label_centroid(const char *label, double *lon, double *lat) {
    if (!lon || !lat) return 1;
    uint8_t bin[16];
    const int layer = hex9_parse_label(label, bin);
    if (layer < 0) return 1;
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(bin, layer, &id)) return 1;
    return h9cell::identity_centroid(id, layer, lon, lat) ? 0 : 1;
}

extern "C" int hex9_common_ancestor(const uint8_t *uuids, size_t n, int layer,
                                    char *buf, size_t buflen, uint8_t *out_uuid) {
    if (!uuids || n == 0 || layer < 0 || layer > 29 || !buf) return -1;
    uint8_t first[16], fnib[32];
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuids, layer, &id)) return -1;
    h9kring::identity_to_uuid(id, layer, first);
    h9a_unpack(first, fnib);

    int shared = layer + 1;                    /* nibbles 0..layer */
    for (size_t i = 1; i < n && shared > 0; i++) {
        uint8_t canon[16], cnib[32];
        if (!h9kring::identity_from_uuid(uuids + i * 16, layer, &id)) return -1;
        h9kring::identity_to_uuid(id, layer, canon);
        h9a_unpack(canon, cnib);
        while (shared > 0 && std::memcmp(cnib, fnib, (size_t)shared) != 0) shared--;
    }
    if (shared == 0) return -1;                /* cells span L0 hexes */

    const int anc_layer = shared - 1;
    /* Canonical labels are prefix-hierarchical, so the ancestor's label is
     * the shared prefix of the (canonical) first body, and its canonical
     * bin is recovered exactly by the label parser's verified tail search.
     * (Do NOT go through normalize_bin/h9_bin_uuid here: their key tail is
     * the deep backward-walk context, which at coarse layers is not the
     * identity-decodable canonical tail.) */
    if (hex9_label(first, anc_layer, buf, buflen) < 0) return -1;
    uint8_t anc[16];
    if (hex9_parse_label(buf, anc) != anc_layer) return -1;
    if (out_uuid) std::memcpy(out_uuid, anc, 16);
    return anc_layer;
}

/* ── Single-cell geometry ───────────────────────────────────────────────────
 * Identity-based (core/h9_cell_geom.h): the UUID resolves to its exact
 * lattice identity and renders through the same construction the grid
 * enumerator uses, so h9_cell == h9_grid for every encoding flavour. (The
 * former cell_unpack float backward-walk followed the centroid-descent
 * convention and could render the wrong hexagon for canonical bins.) */

extern "C" int hex9_cell_ring(const uint8_t uuid[16], int layer, int densify,
                              double *out_lonlat, int max_points) {
    if (layer < 1 || layer > 29) return -1;
    if (densify < 0 || densify > 9 || layer + densify > 29) return -1;
    const int n_ring = hex9_ring_npoints(densify);
    if (!out_lonlat || max_points < n_ring) return -1;

    /* identity-based rendering (h9_cell_geom.h): same construction as the
     * grid enumerator, so the cell polygon for a UUID is always the polygon
     * h9_grid emits for that UUID — regardless of encoding flavour. */
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return -1;
    std::vector<double> lons(n_ring), lats(n_ring);
    if (h9cell::identity_ring(id, layer, densify, lons.data(), lats.data()) != n_ring)
        return -1;
    for (int i = 0; i < n_ring; ++i) { out_lonlat[2*i] = lons[i]; out_lonlat[2*i+1] = lats[i]; }
    return n_ring;
}

/* ── Grid enumeration (SRF) ──────────────────────────────────────────────────
 * Wraps the core BFS (h9grid::enumerate, which dedups internally) and renders
 * per-cell rings on demand. densify=0 uses the BFS-stored vertices; densify>0
 * uses the per-vertex frame switching (ext/octant-seam reflection), ported from
 * the glue's from_cell path. */
struct hex9_grid {
    int layer;
    std::vector<H9GridCell> cells;
};

extern "C" hex9_grid *hex9_grid_create(double lon_min, double lat_min,
                                       double lon_max, double lat_max,
                                       int layer, int densify,
                                       int64_t max_cells,
                                       char *errbuf, size_t errlen) {
    auto fail = [&](const char *m) -> hex9_grid * {
        if (errbuf && errlen) { std::strncpy(errbuf, m, errlen - 1); errbuf[errlen - 1] = '\0'; }
        return nullptr;
    };
    if (layer < 1 || layer > 29) return fail("layer must be 1..29");
    if (densify < 0 || densify > 9 || layer + densify > 29) return fail("invalid densify");

    hex9_grid *g = new hex9_grid;
    g->layer = layer;
    h9grid::enumerate(layer, lon_min, lat_min, lon_max, lat_max, g->cells);
    if (max_cells > 0 && (int64_t)g->cells.size() > max_cells) {
        delete g;
        return fail("cell count exceeds max_cells");
    }
    return g;
}

extern "C" int  hex9_grid_count(const hex9_grid *g) { return g ? (int)g->cells.size() : 0; }
extern "C" void hex9_grid_cell_uuid(const hex9_grid *g, int i, uint8_t out_uuid[16]) {
    std::memcpy(out_uuid, g->cells[(size_t)i].uuid, 16);
}
extern "C" void hex9_grid_cell_id(const hex9_grid *g, int i, uint8_t out_uuid[16]) {
    std::memcpy(out_uuid, g->cells[(size_t)i].id_uuid, 16);
}
extern "C" void hex9_grid_cell_centroid(const hex9_grid *g, int i, double *lon, double *lat) {
    *lon = g->cells[(size_t)i].cen_lon;
    *lat = g->cells[(size_t)i].cen_lat;
}

extern "C" int hex9_grid_cell_ring(const hex9_grid *g, int i, int densify,
                                   double *out_lonlat, int max_points) {
    if (densify < 0 || densify > 9 || g->layer + densify > 29) return -1;
    const int n_ring = hex9_ring_npoints(densify);
    if (!out_lonlat || max_points < n_ring) return -1;
    const H9GridCell &c = g->cells[(size_t)i];
    /* same identity-based construction as hex9_cell_ring (h9_cell_geom.h) */
    h9kring::H9CellId id{c.oid, c.c2, c.ia, c.ib, (bool)c.ext};
    std::vector<double> lons(n_ring), lats(n_ring);
    if (h9cell::identity_ring(id, g->layer, densify, lons.data(), lats.data()) != n_ring)
        return -1;
    for (int k = 0; k < n_ring; ++k) { out_lonlat[2*k] = lons[k]; out_lonlat[2*k+1] = lats[k]; }
    return n_ring;
}

extern "C" void hex9_grid_destroy(hex9_grid *g) { delete g; }

/* ── Adaptive multi-layer grid (population digest) ──────────────────────────
 * See hex9_c.h for the model. Symbolic throughout: canonical binning via
 * h9kring::identity_from_uuid, parent steps via h9kring::normalize_bin. */

struct hex9_adaptive {
    /* `uuid` is the layer-scoped bin KEY (the hex9 column); `full` is a
     * representative FULL uuid of a point the cell digested, kept so geometry
     * renders through identity_from_uuid's guaranteed full path rather than the
     * F2 bin-decode fossil. */
    struct Cell { uint8_t uuid[16]; uint8_t full[16]; int layer; double value; int64_t npoints; };
    std::vector<Cell>    cells;
    std::vector<int64_t> assign;     /* per input point -> cell index */
};

extern "C" hex9_adaptive *hex9_adaptive_create(const uint8_t *uuids,
                                               const double *weight, size_t n,
                                               int min_layer, int max_layer,
                                               double ceiling, double floor_,
                                               char *errbuf, size_t errlen) {
    auto fail = [&](const char *m) -> hex9_adaptive * {
        if (errbuf && errlen) { std::strncpy(errbuf, m, errlen - 1); errbuf[errlen - 1] = '\0'; }
        return nullptr;
    };
    if (!uuids) return fail("uuids must not be NULL");
    if (min_layer < 0 || max_layer > 29 || min_layer > max_layer)
        return fail("need 0 <= min_layer <= max_layer <= 29");
    if (!(ceiling > 0.0)) return fail("ceiling must be > 0");
    if (!(floor_ >= 0.0) || floor_ > ceiling) return fail("need 0 <= floor <= ceiling");

    /* Full-uuid-native digest. Each item keeps its ORIGINAL full uuid; the
     * layer-L grouping key is recomputed fresh from that full uuid via the
     * IDENTITY path (identity_from_uuid -> identity_to_uuid) — the same
     * coarsening h9_grid enumerates with, so each emitted cell's key, geometry
     * and the grid bin it JOINs against all agree (grid-verified at split-hex
     * cells). NOT h9_bin_uuid: that trusts the full uuid's nibble ancestry and
     * mis-coarsens split-hex (6/7/8) cells to the wrong place (the F3 fossil,
     * docs/addressing-doctrine.md). The former path stored a max_layer bin and
     * re-binned the bin upward with normalize_bin — the same fossil, one level
     * worse. Input is FULL uuids only — bins are layer-scoped keys, not
     * addresses, and have already lost the upward path. */
    struct Item { uint8_t full[16]; uint8_t key[16]; double w; int64_t idx; };
    std::vector<Item> work(n);
    std::vector<char> bad(n, 0);
    {
        const ptrdiff_t N = (ptrdiff_t)n;
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < N; ++i) {
            const uint8_t *full = uuids + (size_t)i * 16;
            if ((full[15] >> 4) == 0x0Fu) { bad[i] = 2; continue; }   /* bin input */
            h9kring::H9CellId id;
            if (!h9kring::identity_from_uuid(full, max_layer, &id)) { bad[i] = 1; continue; }
            std::memcpy(work[i].full, full, 16);
            work[i].w   = weight ? weight[i] : 1.0;
            work[i].idx = (int64_t)i;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (bad[i] == 2) return fail("bin uuid input: adaptive takes FULL uuids only "
                                     "(bins are layer-scoped keys, not addresses)");
        if (bad[i])      return fail("uuid could not be binned (malformed?)");
    }

    hex9_adaptive *a = new hex9_adaptive;
    a->assign.assign(n, -1);

    auto by_key_idx = [](const Item &x, const Item &y) {
        const int c = std::memcmp(x.key, y.key, 16);
        return c ? c < 0 : x.idx < y.idx;
    };

    std::vector<Item> carry;
    for (int L = max_layer; L >= min_layer; --L) {
        /* re-key the working set at layer L straight from each FULL uuid, via
         * the grid-canonical identity path (NOT h9_bin_uuid — see above) */
        {
            const ptrdiff_t M = (ptrdiff_t)work.size();
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i = 0; i < M; ++i) {
                h9kring::H9CellId id;
                if (h9kring::identity_from_uuid(work[i].full, L, &id))
                    h9kring::identity_to_uuid(id, L, work[i].key);
                else  /* unreachable: validated decodable at init */
                    std::memset(work[i].key, 0xFF, 16);
            }
        }
        std::sort(work.begin(), work.end(), by_key_idx);
        carry.clear();
        for (size_t i = 0; i < work.size();) {
            size_t j = i;
            double V = 0.0;
            while (j < work.size() && std::memcmp(work[j].key, work[i].key, 16) == 0)
                V += work[j++].w;

            if (V >= floor_ || L == min_layer) {
                /* digest: emit the cell, first-fit whole points to the ceiling
                 * (everything when this is the terminal layer) */
                hex9_adaptive::Cell c;
                std::memcpy(c.uuid, work[i].key, 16);   /* bin key -> hex9 column */
                c.layer = L; c.value = 0.0; c.npoints = 0;
                bool have_full = false;
                const int64_t ci = (int64_t)a->cells.size();
                for (size_t k = i; k < j; ++k) {
                    if (L == min_layer || c.npoints == 0 || c.value + work[k].w <= ceiling) {
                        a->assign[(size_t)work[k].idx] = ci;
                        c.value += work[k].w;
                        c.npoints++;
                        if (!have_full) {           /* representative full uuid for geometry */
                            std::memcpy(c.full, work[k].full, 16);
                            have_full = true;
                        }
                    } else {
                        carry.push_back(work[k]);   /* excess -> parent layer */
                    }
                }
                a->cells.push_back(c);
            } else {
                /* under floor: pass the whole run to the parent layer (it will
                 * be re-keyed from each full uuid at the top of the next pass) */
                for (size_t k = i; k < j; ++k)
                    carry.push_back(work[k]);
            }
            i = j;
        }
        work.swap(carry);
        if (work.empty()) break;
    }
    return a;
}

extern "C" int hex9_adaptive_count(const hex9_adaptive *a) {
    return a ? (int)a->cells.size() : 0;
}
extern "C" void hex9_adaptive_cell(const hex9_adaptive *a, int i, uint8_t out_uuid[16],
                                   int *layer, double *value, int64_t *npoints) {
    const hex9_adaptive::Cell &c = a->cells[(size_t)i];
    if (out_uuid) std::memcpy(out_uuid, c.uuid, 16);
    if (layer)    *layer    = c.layer;
    if (value)    *value    = c.value;
    if (npoints)  *npoints  = c.npoints;
}
extern "C" void hex9_adaptive_cell_full(const hex9_adaptive *a, int i, uint8_t out_full[16]) {
    if (out_full) std::memcpy(out_full, a->cells[(size_t)i].full, 16);
}
extern "C" void hex9_adaptive_assign(const hex9_adaptive *a, int64_t *out) {
    if (out) std::memcpy(out, a->assign.data(), a->assign.size() * sizeof(int64_t));
}
extern "C" void hex9_adaptive_destroy(hex9_adaptive *a) { delete a; }

/* ── Neighbours / k-ring / k-disk ───────────────────────────────────────────
 * Thin wrappers over h9kring (core/h9_kring.h). Output bins are sorted by
 * UUID for deterministic results (matches the grid enumeration convention). */

extern "C" int hex9_neighbors(const uint8_t uuid[16], int layer, uint8_t *out_uuids) {
    if (!out_uuids || layer < 1 || layer > 29) return -1;
    if ((uuid[15] >> 4) == 0x0Fu) return -1;   /* bin input: keys, not addresses */
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return -1;
    h9kring::H9CellId nbs[6];
    const int n = h9kring::neighbors(id, layer, nbs);
    for (int i = 0; i < n; i++)
        h9kring::identity_to_uuid(nbs[i], layer, out_uuids + (size_t)i * 16);
    /* insertion sort by UUID (n <= 6) */
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0; j--) {
            uint8_t *a = out_uuids + (size_t)(j - 1) * 16, *b = out_uuids + (size_t)j * 16;
            if (std::memcmp(a, b, 16) <= 0) break;
            uint8_t t[16];
            std::memcpy(t, a, 16); std::memcpy(a, b, 16); std::memcpy(b, t, 16);
        }
    return n;
}

extern "C" int64_t hex9_disk_ncells(int k) {
    if (k < 0) return -1;
    return 1 + 3 * (int64_t)k * ((int64_t)k + 1);
}

/* Shared k_ring/k_disk body: BFS, filter, sort by UUID. */
static int64_t kring_common(const uint8_t uuid[16], int layer, int k,
                            uint8_t *out_uuids, int64_t max_cells, bool ring_only) {
    if (k < 0 || !out_uuids || max_cells <= 0 || layer < 1 || layer > 29) return -1;
    if ((uuid[15] >> 4) == 0x0Fu) return -1;   /* bin input: keys, not addresses */
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return -1;
    std::vector<std::pair<h9kring::H9CellId, int>> disk;
    if (!h9kring::k_disk(id, layer, k, disk, ring_only ? 0 : max_cells))
        return -1;                                /* disk exceeded max_cells */
    std::vector<std::array<uint8_t, 16>> uu;
    uu.reserve(disk.size());
    for (auto &cd : disk) {
        if (ring_only && cd.second != k) continue;
        std::array<uint8_t, 16> b;
        h9kring::identity_to_uuid(cd.first, layer, b.data());
        uu.push_back(b);
    }
    if ((int64_t)uu.size() > max_cells) return -1;
    std::sort(uu.begin(), uu.end(),
              [](const std::array<uint8_t,16> &a, const std::array<uint8_t,16> &b) {
                  return std::memcmp(a.data(), b.data(), 16) < 0;
              });
    for (size_t i = 0; i < uu.size(); i++)
        std::memcpy(out_uuids + i * 16, uu[i].data(), 16);
    return (int64_t)uu.size();
}

extern "C" int64_t hex9_k_ring(const uint8_t uuid[16], int layer, int k,
                               uint8_t *out_uuids, int64_t max_cells) {
    return kring_common(uuid, layer, k, out_uuids, max_cells, /*ring_only=*/true);
}

extern "C" int64_t hex9_k_disk(const uint8_t uuid[16], int layer, int k,
                               uint8_t *out_uuids, int64_t max_cells) {
    return kring_common(uuid, layer, k, out_uuids, max_cells, /*ring_only=*/false);
}

/* ── Diagnostic (internal, not stable) ──────────────────────────────────────
 * Mirrors the PostGIS h9_diag: BRAW (pre-warp) and descent (post-warp-inverse)
 * (cx,cy) the encoder feeds into descent, plus oid/mode. */
extern "C" int hex9_diag(double lon, double lat, char *buf, size_t buflen) {
    if (!buf) return -1;
    double cx, cy; int oid;
    h9grid::lonlatdeg_to_cxcy_oid(lon, lat, &cx, &cy, &oid);
    const int oct_mode = (int)H9_OID_MO[oid];
    double bx = cx, by = cy;
#if H9_WARP_ENABLE
    h9_warp_inv(cx, cy, oct_mode, &bx, &by);   /* descent runs on warp_inv(BRAW) */
#endif
    int n = std::snprintf(buf, buflen,
        "BRAW=(%.17g, %.17g) descent=(%.17g, %.17g) oid=%d mode=%d",
        cx, cy, bx, by, oid, oct_mode);
    return (n < 0 || (size_t)n >= buflen) ? -1 : n;
}
