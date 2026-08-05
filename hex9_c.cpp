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
#include "h9_e4h.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"
#include "h9_curve.h"
#include "h9_cell_geom.h"
#include "h9_warp_runtime.h"

#include "hex9_c.h"

#include <cstring>
#include <cstddef>
#include <cstdio>
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
/* The ⟨face⟩→⟨cell⟩ seam: (cx, cy, oid) → full uuid via the grid-canonical
 * containment descent, then (reclaimed layout) canonicalise to the mode-0
 * home with the 3-bit tail via the identity round-trip, so the emitted uuid
 * is exactly what bin(.,H9_LMAX) would give. L30 cells are sub-micron, so
 * the round-trip is lossless in practice (stays well under 1 um). Shared by
 * encode_one (canonical mode) and the bring-your-own-projection entry
 * hex9_encode_boct — descent is projection-blind by construction. */
static inline void encode_from_cxcy(double cx, double cy, int oid,
                                    uint8_t out[16]) {
    h9grid::uuid_from_cxcy_full(cx, cy, oid, out);
#if !H9_HAS_HTERM
    {
        h9kring::H9CellId id;
        if (h9kring::identity_from_uuid(out, H9_LMAX, &id))
            h9kring::identity_to_uuid(id, H9_LMAX, out);
    }
#endif
}

/* Same seam entered at the WARPED chart — the coordinates hex9_project /
 * hex9_woct_to_boct emit (the public b_oct). Used by hex9_encode_boct:
 * project → encode_boct is bit-identical to encode because both descend the
 * identical h9_warp_inv output. */
static inline void encode_from_boct(double cx, double cy, int oid,
                                    uint8_t out[16]) {
    h9grid::uuid_from_boct_full(cx, cy, oid, out);
#if !H9_HAS_HTERM
    {
        h9kring::H9CellId id;
        if (h9kring::identity_from_uuid(out, H9_LMAX, &id))
            h9kring::identity_to_uuid(id, H9_LMAX, out);
    }
#endif
}

static inline void encode_one(double lon, double lat, uint8_t out[16],
                              const H9Authalic *aux = &h9_g_authalic) {
    if (g_encoder_mode == 1) {
        double cx, cy; int oid;
        h9grid::lonlatdeg_to_cxcy_oid(lon, lat, &cx, &cy, &oid, aux);
        encode_from_cxcy(cx, cy, oid, out);
    } else {
        H9BOct b = h9_lonlatdeg_to_boct(lon, lat, aux);
        h9_boct_to_uuid(b, out);
#if !H9_HAS_HTERM
        /* Same canonicalisation for the legacy NN path (see encode_from_cxcy). */
        {
            h9kring::H9CellId id;
            if (h9kring::identity_from_uuid(out, H9_LMAX, &id))
                h9kring::identity_to_uuid(id, H9_LMAX, out);
        }
#endif
    }
}
static inline void decode_one(const uint8_t uuid[16], double *lon, double *lat,
                              const H9Authalic *aux = &h9_g_authalic) {
    /* Bin UUIDs decode to the cell's geographic centroid via the exact
     * identity path (grid convention — coherent with h9_grid/h9_cell/labels
     * for every encoding flavour), which already realises the canonical
     * terminal (p_mo=0, centroid t_cell). Legacy full UUIDs keep the beam
     * backward walk (the representative point of the L29 address).
     *
     * Reclaimed layout: there is no h_term, so a max-depth UUID IS a bin — its
     * deepest body nibble (H9_NIB_BODYTOP) is a real digit, not a sentinel.
     * Route it through the SAME identity/centroid path so L30 decodes exactly
     * like every shallower bin (no boct-walk special-case). */
    uint8_t nib[32];
    h9a_unpack(uuid, nib);
#if H9_HAS_HTERM
    const bool is_bin = (nib[H9_NIB_TAIL - 1] == 0x0Fu);
#else
    const bool is_bin = true;
#endif
    if (is_bin) {
        int bl = H9_NIB_BODYTOP;
        while (bl >= 0 && nib[bl] == 0x0Fu) bl--;
        h9kring::H9CellId id;
        if (bl >= 0 && h9kring::identity_from_uuid(uuid, bl, &id) &&
            h9cell::identity_centroid(id, bl, lon, lat, aux))
            return;
        /* fall through to the boct walk for malformed bins */
    }
    H9BOct b = h9_uuid_to_boct(uuid);
    h9_boct_to_lonlatdeg(b, lon, lat, aux);
}

extern "C" const char *hex9_version(void) {
    return "libhex9 " HEX9_VERSION " (" __DATE__ " " __TIME__ ")";
}

extern "C" int hex9_lmax(void) {
    return H9_LMAX;   /* 29 on HEX9_USE_L29 (legacy), 30 when reclaimed */
}

extern "C" int hex9_warp_init(char *errbuf, size_t errlen) {
    std::string err;
    /* The authalic series is the ellipsoid-specific half of the chain; the
     * warp field is the other. Neither is optional — an address cannot be
     * formed without both. */
    if (!h9_authalic_ensure()) {
        if (errbuf && errlen) {
            std::strncpy(errbuf, "authalic series init failed", errlen - 1);
            errbuf[errlen - 1] = '\0';
        }
        return 1;
    }
    const bool ok = h9::h9_warp_init_embedded(&err);
    if (!ok && errbuf && errlen) {
        std::strncpy(errbuf, err.c_str(), errlen - 1);
        errbuf[errlen - 1] = '\0';
    }
    return ok ? 0 : 1;
}

/* Canonical name since 2.1.0: the function builds the whole addressing chain
 * (authalic series + warp field), not just the warp. hex9_warp_init above is
 * the kept-forever historic alias. */
extern "C" int hex9_init(char *errbuf, size_t errlen) {
    return hex9_warp_init(errbuf, errlen);
}

extern "C" void hex9_set_use_warp(int on) { h9::g_warp_use = (on != 0); }
extern "C" void hex9_set_encoder(int mode) { g_encoder_mode = mode; }

extern "C" int hex9_encode(double lon, double lat, uint8_t out_uuid[16]) {
    encode_one(lon, lat, out_uuid);
    return 0;
}

extern "C" int hex9_decode(const uint8_t uuid[16], double *lon, double *lat) {
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_decode */
    decode_one(uuid, lon, lat);
    return 0;
}

/* ── Sphere-datum twins ──────────────────────────────────────────────────────
 * Identical chain minus the WGS84 authalic reduction: (lon, lat) are
 * already-spherical degrees (h9_g_sphere, identity both ways). Same kernels,
 * different aux — see hex9_c.h for the datum doctrine. */
extern "C" int hex9_encode_sphere(double lon, double lat, uint8_t out_uuid[16]) {
    encode_one(lon, lat, out_uuid, &h9_g_sphere);
    return 0;
}

extern "C" int hex9_decode_sphere(const uint8_t uuid[16], double *lon, double *lat) {
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_decode */
    decode_one(uuid, lon, lat, &h9_g_sphere);
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
    if (layer == 0) {
        /* L0 bin = the root hex itself (no split ambiguity), canonicalised to
         * the mode-0 octant rep: keep nibble[0], sentinel the body, and write
         * the canonical key_tail (c2_canon<<1, r_mo=p_mo=0). Matches Python
         * h9_bin_pts' L0 branch; the general h9_bin_uuid path corrupts the root
         * nibble at layer 0. Malformed root (>11) falls through to legacy. */
        uint8_t nib[32];
        h9a_unpack(uuid, nib);
        if (nib[0] <= 11) {
            uint8_t obn[32];
            obn[0] = nib[0];
            for (int i = 1; i <= H9_NIB_TAIL - 1; ++i) obn[i] = 0x0Fu;
            obn[H9_NIB_TAIL] = (uint8_t)((H9_L0HEX_BACK[nib[0]][0][1] & 3u) << 1);
            h9a_pack(obn, out);
            return;
        }
    }
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
    if (layer < 0 || layer > H9_LMAX) return 1;
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_bin */
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
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuid + i * 16)) return 1;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        decode_one(uuid + (size_t)i * 16, &lon[i], &lat[i]);
    return 0;
}

extern "C" int hex9_encode_many_sphere(const double *lon, const double *lat,
                                       size_t n, uint8_t *out_uuid) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        encode_one(lon[i], lat[i], out_uuid + (size_t)i * 16, &h9_g_sphere);
    return 0;
}

extern "C" int hex9_decode_many_sphere(const uint8_t *uuid, size_t n,
                                       double *lon, double *lat) {
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuid + i * 16)) return 1;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        decode_one(uuid + (size_t)i * 16, &lon[i], &lat[i], &h9_g_sphere);
    return 0;
}

extern "C" int hex9_bin_many(const uint8_t *uuid, int layer, size_t n,
                             uint8_t *out_uuid) {
    if (layer < 0 || layer > H9_LMAX) return 1;
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuid + i * 16)) return 1;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        canonical_bin(uuid + (size_t)i * 16, layer, out_uuid + (size_t)i * 16);
    return 0;
}

/* ── Canonical cell parent / ancestor (mode-0 convention) ─────────────────
 * Thin surface over h9kring::h9_cell_parent_uuid / h9_cell_ancestor_uuid —
 * see the doctrine comment there. */
extern "C" int hex9_cell_parent(const uint8_t uuid[16], uint8_t out_uuid[16]) {
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_bin/split */
    return h9kring::h9_cell_parent_uuid(uuid, out_uuid) ? 0 : 1;
}

extern "C" int hex9_cell_parent_many(const uint8_t *uuid, size_t n,
                                     uint8_t *out_uuid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= h9e4h_is_marked(uuid + (size_t)i * 16) ? 1
            : (h9kring::h9_cell_parent_uuid(uuid + (size_t)i * 16,
                                            out_uuid + (size_t)i * 16) ? 0 : 1);
    return rc;
}

extern "C" int hex9_cell_ancestor(const uint8_t uuid[16], int layer,
                                  uint8_t out_uuid[16]) {
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_bin/split */
    return h9kring::h9_cell_ancestor_uuid(uuid, layer, out_uuid) ? 0 : 1;
}

extern "C" int hex9_cell_ancestor_many(const uint8_t *uuid, int layer, size_t n,
                                       uint8_t *out_uuid) {
    if (layer < 0 || layer > H9_LMAX) return 1;
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= h9e4h_is_marked(uuid + (size_t)i * 16) ? 1
            : (h9kring::h9_cell_ancestor_uuid(uuid + (size_t)i * 16, layer,
                                              out_uuid + (size_t)i * 16) ? 0 : 1);
    return rc;
}

/* ── Hamiltonian curve addressing ───────────────────────────────────────────
 * Surface over core/h9_curve.h (36-state transducer; tables verbatim from
 * hhg9). Input normalisation policy: a curve-uuid passes through; a
 * canonical bin is trusted as-is (like the cell parent/ancestor family); a
 * full/reversible uuid is re-binned canonically at its own depth through
 * canonical_bin — the SAME path hex9_bin takes, so curve(encode(pt)) and
 * curve(bin(encode(pt), L)) agree with h9_grid / Python at every layer. */

/* Shared kernel: uuid (h9 or curve) -> curve rows + layer. */
static int curve_rows_of(const uint8_t uuid[16], uint8_t *rows) {
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H: no curve on tails */
    if (h9curve::is_curve_uuid(uuid))
        return h9curve::curve_unpack(uuid, rows);
    bool rebin = false;
    const int L = h9curve::curve_input_layer(uuid, &rebin);
    if (L < 0) return -1;
    uint8_t bin[16];
    if (rebin)
        canonical_bin(uuid, L, bin);
    else
        std::memcpy(bin, uuid, 16);
    if (!h9curve::curve_rows_from_bin(bin, L, rows, nullptr)) return -1;
    return L;
}

extern "C" int hex9_curve(const uint8_t uuid[16], uint8_t out_curve[16]) {
    uint8_t rows[H9_NLEVELS];
    const int L = curve_rows_of(uuid, rows);
    if (L < 0) return 1;
    h9curve::curve_pack_rows(rows, L, out_curve);
    return 0;
}

extern "C" int hex9_curve_many(const uint8_t *uuid, size_t n, uint8_t *out_curve) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= hex9_curve(uuid + (size_t)i * 16, out_curve + (size_t)i * 16);
    return rc;
}

extern "C" int hex9_curve_decode(const uint8_t curve[16], uint8_t out_uuid[16]) {
    if (h9e4h_is_marked(curve)) return 1;   /* E4H is neither curve nor h9 */
    if (!h9curve::is_curve_uuid(curve)) {           /* h9-uuid passes through */
        if (h9curve::curve_input_layer(curve, nullptr) < 0) return 1;
        std::memcpy(out_uuid, curve, 16);
        return 0;
    }
    return h9curve::curve_decode_to_bin(curve, out_uuid) ? 0 : 1;
}

extern "C" int hex9_curve_decode_many(const uint8_t *curve, size_t n,
                                      uint8_t *out_uuid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= hex9_curve_decode(curve + (size_t)i * 16, out_uuid + (size_t)i * 16);
    return rc;
}

extern "C" int hex9_is_curve(const uint8_t uuid[16]) {
    return h9curve::is_curve_uuid(uuid) ? 1 : 0;
}

extern "C" int hex9_curve_layer(const uint8_t curve[16]) {
    uint8_t rows[H9_NLEVELS];
    if (!h9curve::is_curve_uuid(curve)) return -1;
    return h9curve::curve_unpack(curve, rows);
}

extern "C" int hex9_curve_bin(const uint8_t curve[16], int layer,
                              uint8_t out_curve[16]) {
    uint8_t rows[H9_NLEVELS];
    if (!h9curve::is_curve_uuid(curve)) return 1;
    const int L = h9curve::curve_unpack(curve, rows);
    if (L < 0 || layer < 0 || layer > L) return 1;
    h9curve::curve_pack_rows(rows, layer, out_curve);
    return 0;
}

extern "C" int hex9_curve_index(const uint8_t uuid[16], char *buf, size_t buflen) {
    uint8_t rows[H9_NLEVELS];
    const int L = curve_rows_of(uuid, rows);
    if (L < 0) return -1;
    return h9curve::curve_index_decimal(rows, L, buf, buflen);
}

extern "C" int hex9_curve_pack(const char *index_dec, int layer,
                               uint8_t out_curve[16]) {
    if (layer < 0 || layer > H9_LMAX) return 1;
    uint8_t rows[H9_NLEVELS];
    if (!h9curve::curve_rows_from_decimal(index_dec, layer, rows)) return 1;
    h9curve::curve_pack_rows(rows, layer, out_curve);
    return 0;
}

extern "C" int hex9_curve_label(const uint8_t uuid[16], char *buf, size_t buflen) {
    uint8_t rows[H9_NLEVELS];
    const int L = curve_rows_of(uuid, rows);
    if (L < 0) return -1;
    if (buflen < (size_t)L + 3) return -1;          /* 'c' + slot + L ranks + NUL */
    buf[0] = 'c';
    buf[1] = (char)(rows[0] < 10 ? '0' + rows[0] : 'a' + rows[0] - 10);
    for (int k = 1; k <= L; ++k) buf[1 + k] = (char)('0' + rows[k]);
    buf[L + 2] = '\0';
    return L + 2;
}

extern "C" int hex9_curve_parse_label(const char *label, uint8_t out_curve[16]) {
    if (!label || label[0] != 'c' || !label[1]) return 1;
    const char sc = label[1];
    uint8_t rows[H9_NLEVELS];
    if (sc >= '0' && sc <= '9') rows[0] = (uint8_t)(sc - '0');
    else if (sc == 'a' || sc == 'b') rows[0] = (uint8_t)(sc - 'a' + 10);
    else return 1;
    int L = 0;
    for (const char *p = label + 2; *p; ++p) {
        if (*p < '0' || *p > '8' || L >= H9_LMAX) return 1;
        rows[++L] = (uint8_t)(*p - '0');
    }
    h9curve::curve_pack_rows(rows, L, out_curve);
    return 0;
}

extern "C" int64_t hex9_curve_ncells(int from_layer, int to_layer) {
    if (from_layer < 0 || to_layer < from_layer || to_layer > H9_LMAX) return -1;
    const int g = to_layer - from_layer;
    if (g > 19) return -1;                           /* 9^20 overflows int64 */
    int64_t n = 1;
    for (int i = 0; i < g; ++i) n *= 9;
    return n;
}

/* Shared input resolution for the enumeration family: h9 bin trusted,
 * full uuid re-binned canonically, curve-uuid decoded. Returns the cell's
 * own layer, or -1. */
static int curve_cell_input(const uint8_t uuid[16], uint8_t anc[16]) {
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H: no lineage on tails */
    if (h9curve::is_curve_uuid(uuid)) {
        if (!h9curve::curve_decode_to_bin(uuid, anc)) return -1;
    } else {
        bool rebin = false;
        const int L = h9curve::curve_input_layer(uuid, &rebin);
        if (L < 0) return -1;
        if (rebin) canonical_bin(uuid, L, anc);
        else       std::memcpy(anc, uuid, 16);
    }
    uint8_t nib[32];
    h9a_unpack(anc, nib);
    return h9curve::curve_body_layer(nib);
}

extern "C" int64_t hex9_curve_cells(const uint8_t uuid[16], int layer,
                                    uint8_t *out_bins, uint8_t *out_curves,
                                    int64_t max_cells) {
    if (layer < 0 || layer > H9_LMAX || !out_bins) return -1;
    uint8_t anc[16];
    const int from = curve_cell_input(uuid, anc);
    if (from < 0) return -1;
    const int64_t want = hex9_curve_ncells(from, layer);
    if (want < 0 || want > max_cells) return -1;
    uint8_t rows[H9_NLEVELS];
    int state = 0;
    if (!h9curve::curve_rows_from_bin(anc, from, rows, &state)) return -1;
    int64_t count = 0;
    if (!h9curve::curve_cells_rec(anc, from, state, layer, rows,
                                  out_bins, out_curves, &count))
        return -1;
    return count;
}

extern "C" int hex9_cell_children(const uint8_t uuid[16], uint8_t *out_uuids) {
    if (!out_uuids) return 1;
    uint8_t cell[16];
    const int L = curve_cell_input(uuid, cell);
    if (L < 0 || L >= H9_LMAX) return 1;
    uint8_t rows[H9_NLEVELS];
    int state = 0;
    if (!h9curve::curve_rows_from_bin(cell, L, rows, &state)) return 1;
    uint8_t kids[9][16];
    if (!h9curve::curve_children_ranked(cell, L, state, kids)) return 1;
    std::memcpy(out_uuids, kids, 9 * 16);
    return 0;
}

extern "C" int64_t hex9_owned_cells(const uint8_t uuid[16], int layer,
                                    uint8_t *out_bins, uint8_t *out_curves,
                                    int64_t max_cells) {
    if (layer < 0 || layer > H9_LMAX || !out_bins) return -1;
    uint8_t zone[16];
    const int from = curve_cell_input(uuid, zone);
    if (from < 0 || layer < from) return -1;
    return h9curve::curve_owned_cells(zone, from, layer,
                                      out_bins, out_curves, max_cells);
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

/* Sphere-datum twins of the projection surface (see hex9_encode_sphere). */
extern "C" int hex9_project_sphere(double lon, double lat,
                                   double *cx, double *cy, int *oid) {
    h9grid::lonlatdeg_to_boct_oid(lon, lat, cx, cy, oid, &h9_g_sphere);
    return 0;
}

extern "C" int hex9_project_many_sphere(const double *lon, const double *lat,
                                        size_t n,
                                        double *cx, double *cy, int *oid) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::lonlatdeg_to_boct_oid(lon[i], lat[i], &cx[i], &cy[i], &oid[i],
                                      &h9_g_sphere);
    return 0;
}

extern "C" int hex9_unproject_sphere(double cx, double cy, int oid,
                                     double *lon, double *lat) {
    return h9grid::boct_oid_to_lonlat(cx, cy, oid, lon, lat, &h9_g_sphere) ? 0 : 1;
}

extern "C" int hex9_unproject_many_sphere(const double *cx, const double *cy,
                                          const int *oid, size_t n,
                                          double *lon, double *lat) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::boct_oid_to_lonlat(cx[i], cy[i], oid[i], &lon[i], &lat[i],
                                   &h9_g_sphere);
    return 0;
}

/* ── Warped octahedral cartesian CRS (w_oct; the 3D storage baseline) ────────
 * Seamless 3D coordinate for the sphere — see hex9_c.h. xyz is on the unit
 * octahedron; oid == sign(xyz). b_oct <-> xyz is pure rotation; lon/lat <-> xyz
 * crosses the authalic warp (via hex9_project/unproject). */
extern "C" int hex9_boct_to_woct(double cx, double cy, int oid, double xyz[3]) {
    if (oid < 0 || oid > 7) return 1;
    h9grid::boct_to_woct(cx, cy, oid, xyz);
    return 0;
}

extern "C" int hex9_woct_to_boct(const double xyz[3], double *cx, double *cy,
                                 int *oid) {
    h9grid::woct_to_boct(xyz, cx, cy, oid);
    return 0;
}

extern "C" int hex9_to_woct(double lon, double lat, double xyz[3]) {
    double cx, cy; int oid;
    h9grid::lonlatdeg_to_boct_oid(lon, lat, &cx, &cy, &oid);
    h9grid::boct_to_woct(cx, cy, oid, xyz);
    return 0;
}

extern "C" int hex9_from_woct(const double xyz[3], double *lon, double *lat) {
    double cx, cy; int oid;
    h9grid::woct_to_boct(xyz, &cx, &cy, &oid);
    return h9grid::boct_oid_to_lonlat(cx, cy, oid, lon, lat) ? 0 : 1;
}

extern "C" int hex9_boct_to_woct_many(const double *cx, const double *cy,
                                      const int *oid, size_t n, double *xyz) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::boct_to_woct(cx[i], cy[i], oid[i], xyz + (size_t)i * 3);
    return 0;
}

extern "C" int hex9_woct_to_boct_many(const double *xyz, size_t n,
                                      double *cx, double *cy, int *oid) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i)
        h9grid::woct_to_boct(xyz + (size_t)i * 3, &cx[i], &cy[i], &oid[i]);
    return 0;
}

extern "C" int hex9_to_woct_many(const double *lon, const double *lat, size_t n,
                                 double *xyz) {
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < N; ++i) {
        double cx, cy; int oid;
        h9grid::lonlatdeg_to_boct_oid(lon[i], lat[i], &cx, &cy, &oid);
        h9grid::boct_to_woct(cx, cy, oid, xyz + (size_t)i * 3);
    }
    return 0;
}

extern "C" int hex9_from_woct_many(const double *xyz, size_t n,
                                   double *lon, double *lat) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i) {
        double cx, cy; int oid;
        h9grid::woct_to_boct(xyz + (size_t)i * 3, &cx, &cy, &oid);
        rc |= h9grid::boct_oid_to_lonlat(cx, cy, oid, &lon[i], &lat[i]) ? 0 : 1;
    }
    return rc;
}

extern "C" int hex9_ring_npoints(int densify) {
    if (densify < 0 || densify > 9) return -1;
    int n = 1;
    for (int i = 0; i < densify; ++i) n *= 3;
    return 6 * n + 1;
}

/* ── Face-coordinate addressing (bring your own ⟨polyhedron⟩) ────────────────
 * The ⟨face⟩⟨cell⟩ layers without AKW: mint and read addresses directly at
 * the (cx, cy, oid) chart. Descent is projection-blind, so a caller with her
 * own sphere→octahedron map gets the full cell machinery — and the resulting
 * addresses are NON-CANONICAL: same 128-bit space, different meaning. The
 * projection identity is dataset metadata; never mix projections within one
 * dataset (see docs/projection.md). Always the grid-canonical containment
 * descent — hex9_set_encoder does not apply here. */
extern "C" int hex9_encode_boct(double cx, double cy, int oid,
                                uint8_t out_uuid[16]) {
    if (!out_uuid || oid < 0 || oid > 7 ||
        !std::isfinite(cx) || !std::isfinite(cy)) return 1;
    encode_from_boct(cx, cy, oid, out_uuid);
    return 0;
}

extern "C" int hex9_encode_boct_many(const double *cx, const double *cy,
                                     const int *oid, size_t n,
                                     uint8_t *out_uuid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= hex9_encode_boct(cx[i], cy[i], oid[i], out_uuid + (size_t)i * 16);
    return rc;
}

/* Cell centroid at the uuid's bin layer, in face coordinates — the same
 * lattice centroid hex9_decode unprojects, emitted before the unprojection
 * (projection-free; the chart is the returned oid's). */
extern "C" int hex9_decode_boct(const uint8_t uuid[16],
                                double *cx, double *cy, int *oid) {
    if (!uuid || !cx || !cy || !oid) return 1;
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: hex9_e4h_decode */
    uint8_t nib[32];
    h9a_unpack(uuid, nib);
#if H9_HAS_HTERM
    /* Legacy layout: only bin uuids carry the identity path; full L29
     * addresses would need the beam walk, which is lon/lat-coupled. */
    if (nib[H9_NIB_TAIL - 1] != 0x0Fu) return 1;
#endif
    int bl = H9_NIB_BODYTOP;
    while (bl >= 0 && nib[bl] == 0x0Fu) bl--;
    h9kring::H9CellId id;
    if (bl >= 0 && h9kring::identity_from_uuid(uuid, bl, &id) &&
        h9cell::identity_centroid_boct(id, bl, cx, cy, oid))
        return 0;
    return 1;
}

extern "C" int hex9_decode_boct_many(const uint8_t *uuid, size_t n,
                                     double *cx, double *cy, int *oid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= hex9_decode_boct(uuid + (size_t)i * 16, &cx[i], &cy[i], &oid[i]);
    return rc;
}

/* Cell boundary ring in face coordinates: (cx, cy, oid) per vertex, each in
 * the chart of ITS oid (per-vertex frames, same seam-reflection choice the
 * lon/lat ring uses before unprojection). Closed ring — hex9_ring_npoints
 * points, last repeats first. Returns point count, or -1 on error. */
extern "C" int hex9_cell_ring_boct(const uint8_t uuid[16], int layer,
                                   int densify,
                                   double *out_cx, double *out_cy,
                                   int *out_oid, int max_points) {
    if (!out_cx || !out_cy || !out_oid) return -1;
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: no h9 cell ring */
    if (layer < 0 || layer > H9_LMAX) return -1;
    if (densify < 0 || densify > 9 || layer + densify > H9_LMAX) return -1;
    const int n_ring = hex9_ring_npoints(densify);
    if (max_points < n_ring) return -1;
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return -1;
    return h9cell::identity_ring_boct(id, layer, densify, out_cx, out_cy, out_oid);
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
    if (layer < 0 || layer > H9_LMAX || !buf) return -1;
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: hex9_e4h_label */
    uint8_t cb[16], nibbles[32];
    canonical_bin(uuid, layer, cb);
    h9a_unpack(cb, nibbles);
    int len = write_label(nibbles, layer, buf);
    if (len < 0 || (size_t)len + 1 > buflen) return -1;
    buf[len] = '\0';
    return len;
}

extern "C" int hex9_label_key(const uint8_t uuid[16], int layer, char *buf, size_t buflen) {
    if (layer < 0 || layer > H9_LMAX || !buf) return -1;
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: hex9_e4h_label */
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
    if (len < 1 || len > H9_NIB_BODYTOP + 1) return -1;   /* body L0..L{BODYTOP} */
    for (int i = 0; i < len; i++) {
        const char *p = (const char *)std::memchr(H9D, label[i], 12);
        if (!p) return -1;
        const int d = (int)(p - H9D);
        if (i > 0 && d > 8) return -1;        /* layers 1+ use digits 0..8 */
        nibbles[i] = (uint8_t)d;
    }
    for (int i = len; i <= H9_NIB_TAIL - 1; i++) nibbles[i] = 0x0Fu;
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
    if (!uuids || n == 0 || layer < 0 || layer > H9_LMAX || !buf) return -1;
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuids + i * 16)) return -1;
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
    /* The ancestor's label IS the shared nibble prefix (write_label is a
     * pure nibble->char map); its canonical bin is recovered by the label
     * parser's verified tail search. Do NOT label via hex9_label(first,
     * anc_layer): that re-coarsens the bin through the identity/context
     * walk — the F3 fossil path, unguaranteed for bins below their own
     * layer (it mislabelled GB's L1 ancestor '43' as '65' on the L30
     * layout). The prefix needs no walk at all. */
    if (buflen < (size_t)shared + 1) return -1;
    write_label(fnib, anc_layer, buf);
    buf[shared] = '\0';
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

/* ── Cell lattice identity (integer UV) — see hex9_c.h ─────────────────── */

/* Canonical pool key for a vertex. resolve_uv_frame canonicalises only
 * STRICTLY out-of-face vertices; a vertex ON a face boundary is in-face for
 * every frame that shares the edge, so each cell would otherwise report it
 * in its own frame and pooling would miss the match (found by test: the
 * L0-descended seam chains). Rule: enumerate the equivalent representations
 * across every edge the vertex lies on (seam maps; corners compose) and
 * return the lexicographically smallest (oid, ia, ib). Deterministic, frame
 * independent, exact integers throughout. */
static void uv_canonical(int64_t s, int64_t *u, int64_t *v, int *oid) {
    struct Rep { int64_t u, v; int oid; };
    Rep reps[8];
    int n = 0;
    reps[n++] = { *u, *v, *oid };
    for (int i = 0; i < n; ++i) {
        const int mo = (int)H9_OID_MO[reps[i].oid];
        for (int e = 0; e < 3; ++e) {
            int64_t ru, rv;
            h9kring::seam_apply(e, mo, s, reps[i].u, reps[i].v, &ru, &rv);
            const int nb = H9_OID_NB[reps[i].oid][e];
            if (!h9cell::uv_in_face((int)H9_OID_MO[nb], s, ru, rv)) continue;
            bool dup = false;
            for (int k = 0; k < n; ++k)
                if (reps[k].u == ru && reps[k].v == rv && reps[k].oid == nb) {
                    dup = true;
                    break;
                }
            if (!dup && n < 8) reps[n++] = { ru, rv, nb };
        }
    }
    /* only boundary vertices have orbit > 1; interior verts fall through */
    int best = 0;
    for (int i = 1; i < n; ++i)
        if (reps[i].oid < reps[best].oid ||
            (reps[i].oid == reps[best].oid &&
             (reps[i].u < reps[best].u ||
              (reps[i].u == reps[best].u && reps[i].v < reps[best].v))))
            best = i;
    *u = reps[best].u;
    *v = reps[best].v;
    *oid = reps[best].oid;
}

extern "C" int hex9_cell_uv(const uint8_t uuid[16], int layer,
                            int64_t *c_ia, int64_t *c_ib, int *c_oid,
                            int64_t v_ia[6], int64_t v_ib[6], int v_oid[6],
                            int *ext) {
    if (layer < 0 || layer > H9_LMAX) return 1;
    if (h9e4h_is_marked(uuid)) return 1;    /* E4H input: no h9 lattice id */
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return 1;
    /* centre key stays in the cell's GOVERNING frame (mode-0 side for seam
     * cells — the identity's own oid), so c_oid keeps its lineage meaning;
     * vertex keys are pool keys and get the canonical representative. */
    *c_ia  = id.ia + h9kring::H9KR_C2_DU[id.c2];
    *c_ib  = id.ib + h9kring::H9KR_C2_DV[id.c2];
    *c_oid = id.oid;
    *ext   = id.ext ? 1 : 0;
    h9cell::identity_vertices(id, layer, v_ia, v_ib, v_oid);
    const int64_t s = h9kring::pow3(layer);
    for (int v = 0; v < 6; ++v)
        uv_canonical(s, &v_ia[v], &v_ib[v], &v_oid[v]);
    return 0;
}

extern "C" int hex9_cell_uv_many(const uint8_t *uuid, const int32_t *layer,
                                 size_t n,
                                 int64_t *c_ia, int64_t *c_ib, int32_t *c_oid,
                                 int64_t *v_ia, int64_t *v_ib, int32_t *v_oid,
                                 int32_t *ext) {
    int rc = 0;
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuid + i * 16)) return 1;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i) {
        int co, ex;
        int vo[6];
        rc |= hex9_cell_uv(uuid + (size_t)i * 16, (int)layer[i],
                           &c_ia[i], &c_ib[i], &co,
                           v_ia + (size_t)i * 6, v_ib + (size_t)i * 6, vo,
                           &ex);
        c_oid[i] = (int32_t)co;
        ext[i]   = (int32_t)ex;
        for (int v = 0; v < 6; ++v) v_oid[(size_t)i * 6 + v] = (int32_t)vo[v];
    }
    return rc;
}

extern "C" void hex9_uv_units(double *u1, double *v3) {
    *u1 = H9_UV_U1;
    *v3 = H9_UV_V3;
}

static int cell_ring_impl(const uint8_t uuid[16], int layer, int densify,
                          double *out_lonlat, int max_points,
                          const H9Authalic *aux) {
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: no h9 cell ring */
    if (layer < 0 || layer > H9_LMAX) return -1;
    if (densify < 0 || densify > 9 || layer + densify > H9_LMAX) return -1;
    const int n_ring = hex9_ring_npoints(densify);
    if (!out_lonlat || max_points < n_ring) return -1;

    /* identity-based rendering (h9_cell_geom.h): same construction as the
     * grid enumerator, so the cell polygon for a UUID is always the polygon
     * h9_grid emits for that UUID — regardless of encoding flavour. */
    h9kring::H9CellId id;
    if (!h9kring::identity_from_uuid(uuid, layer, &id)) return -1;
    std::vector<double> lons(n_ring), lats(n_ring);
    if (h9cell::identity_ring(id, layer, densify, lons.data(), lats.data(), aux) != n_ring)
        return -1;
    for (int i = 0; i < n_ring; ++i) { out_lonlat[2*i] = lons[i]; out_lonlat[2*i+1] = lats[i]; }
    return n_ring;
}

extern "C" int hex9_cell_ring(const uint8_t uuid[16], int layer, int densify,
                              double *out_lonlat, int max_points) {
    return cell_ring_impl(uuid, layer, densify, out_lonlat, max_points,
                          &h9_g_authalic);
}

extern "C" int hex9_cell_ring_sphere(const uint8_t uuid[16], int layer,
                                     int densify,
                                     double *out_lonlat, int max_points) {
    return cell_ring_impl(uuid, layer, densify, out_lonlat, max_points,
                          &h9_g_sphere);
}

/* ── Grid enumeration (SRF) ──────────────────────────────────────────────────
 * Wraps the core BFS (h9grid::enumerate, which dedups internally) and renders
 * per-cell rings on demand. densify=0 uses the BFS-stored vertices; densify>0
 * uses the per-vertex frame switching (ext/octant-seam reflection), ported from
 * the glue's from_cell path. */
struct hex9_grid {
    int layer;
    /* Datum of every lon/lat this handle emits (and of the bbox it was built
     * from). Set at create, immutable after — the handle knows what it is, so
     * the accessors need no _sphere twins. */
    const H9Authalic *aux;
    std::vector<H9GridCell> cells;
};

static hex9_grid *grid_create_impl(double lon_min, double lat_min,
                                   double lon_max, double lat_max,
                                   int layer, int densify,
                                   int64_t max_cells,
                                   char *errbuf, size_t errlen,
                                   const H9Authalic *aux) {
    auto fail = [&](const char *m) -> hex9_grid * {
        if (errbuf && errlen) { std::strncpy(errbuf, m, errlen - 1); errbuf[errlen - 1] = '\0'; }
        return nullptr;
    };
    if (layer < 0 || layer > H9_LMAX) {
        char m[32];
        std::snprintf(m, sizeof m, "layer must be 0..%d", H9_LMAX);
        return fail(m);
    }
    if (densify < 0 || densify > 9 || layer + densify > H9_LMAX) return fail("invalid densify");

    hex9_grid *g = new hex9_grid;
    g->layer = layer;
    g->aux   = aux;
    // enumerate bails mid-BFS once the budget is exceeded, so a deep whole-world
    // request fails fast here instead of materialising 12·9^layer nodes.
    if (!h9grid::enumerate(layer, lon_min, lat_min, lon_max, lat_max, g->cells,
                           max_cells, aux)) {
        delete g;
        return fail("cell count exceeds max_cells");
    }
    return g;
}

extern "C" hex9_grid *hex9_grid_create(double lon_min, double lat_min,
                                       double lon_max, double lat_max,
                                       int layer, int densify,
                                       int64_t max_cells,
                                       char *errbuf, size_t errlen) {
    return grid_create_impl(lon_min, lat_min, lon_max, lat_max, layer, densify,
                            max_cells, errbuf, errlen, &h9_g_authalic);
}

extern "C" hex9_grid *hex9_grid_create_sphere(double lon_min, double lat_min,
                                              double lon_max, double lat_max,
                                              int layer, int densify,
                                              int64_t max_cells,
                                              char *errbuf, size_t errlen) {
    return grid_create_impl(lon_min, lat_min, lon_max, lat_max, layer, densify,
                            max_cells, errbuf, errlen, &h9_g_sphere);
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
    if (densify < 0 || densify > 9 || g->layer + densify > H9_LMAX) return -1;
    const int n_ring = hex9_ring_npoints(densify);
    if (!out_lonlat || max_points < n_ring) return -1;
    const H9GridCell &c = g->cells[(size_t)i];
    /* same identity-based construction as hex9_cell_ring (h9_cell_geom.h) */
    h9kring::H9CellId id{c.oid, c.c2, c.ia, c.ib, (bool)c.ext};
    std::vector<double> lons(n_ring), lats(n_ring);
    if (h9cell::identity_ring(id, g->layer, densify, lons.data(), lats.data(),
                              g->aux) != n_ring)
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
    for (size_t i = 0; i < n; ++i)
        if (h9e4h_is_marked(uuids + i * 16))
            return fail("E4H uuid input — the digest takes h9 full uuids");
    if (min_layer < 0 || max_layer > H9_LMAX || min_layer > max_layer) {
        char msg[80];
        std::snprintf(msg, sizeof msg,
                      "need 0 <= min_layer <= max_layer <= %d", H9_LMAX);
        return fail(msg);
    }
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
    if (!out_uuids || layer < 0 || layer > H9_LMAX) return -1;
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: h9 adjacency only */
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
    if (k < 0 || !out_uuids || max_cells <= 0 || layer < 0 || layer > H9_LMAX) return -1;
    if (h9e4h_is_marked(uuid)) return -1;   /* E4H input: h9 adjacency only */
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

/* ── E4H: aperture-4 structural tail ────────────────────────────────────────
 * ABI over core/h9_e4h.h (the exact classifier; tables frozen from the hhg9
 * reference by tools/gen_e4h_tables.py). See hex9_c.h for the model and
 * docs/universality.md for the contract: the seed FP program + the 2^-46
 * snap + the integer descent ARE the E4H address definition. */

/* host-bin context: key_tail state + the LATTICE HEX CENTRE in the host's
 * governing chart. The centre convention matters: hhg9's h9_dec (the frame
 * the reference builds on) uses the whole-hexagon lattice centre — for a
 * seam-straddling (ext) host that point lies ON the octant seam, NOT at the
 * mode-0 half's representative centroid that hex9_decode_boct emits (the F2
 * distinction hex9_cell_uv documents). Using the representative here shifted
 * every ext-host frame by half a cell and broke corpus parity on all of
 * them; the lattice centre agrees with the reference to <= 1 ulp
 * (cell-relative <= 7e-12 — far under the corpus margin filter). */
static int e4h_host_info(const uint8_t host[16], int layer,
                         int *hoid, int *c2, int *p_mo,
                         double *hcx, double *hcy) {
    uint8_t nib[32];
    h9a_unpack(host, nib);
    const uint8_t kt = nib[H9_NIB_TAIL];
    *c2   = (kt >> 1) & 3;
    *p_mo = (kt >> 3) & 1;
    if (*c2 > 2) return E4H_EBAD;
    int64_t c_ia, c_ib, v_ia[6], v_ib[6];
    int c_oid, v_oid[6], ext;
    if (hex9_cell_uv(host, layer, &c_ia, &c_ib, &c_oid,
                     v_ia, v_ib, v_oid, &ext)) return E4H_EBAD;
    double u1, v3, p3 = 1.0;
    hex9_uv_units(&u1, &v3);
    for (int i = 0; i < layer; ++i) p3 *= 3.0;
    *hcx = (double)c_ia * u1 / p3;
    *hcy = (double)c_ib * v3 / p3;
    *hoid = c_oid;
    return 0;
}

static int e4h_encode_one(double lon, double lat, int layer, int depth,
                          int sphere, uint8_t out[16]) {
    if (!out || layer < 0 || depth < 0 ||
        layer + 2 + depth > H9_NIB_BODYTOP) return E4H_EBAD;
    double px, py;
    int g;
    int rc = sphere ? hex9_project_sphere(lon, lat, &px, &py, &g)
                    : hex9_project(lon, lat, &px, &py, &g);
    if (rc) return E4H_EBAD;
    uint8_t full[16], host[16];
    rc = sphere ? hex9_encode_sphere(lon, lat, full)
                : hex9_encode(lon, lat, full);
    if (rc) return E4H_EBAD;
    if (hex9_bin(full, layer, host)) return E4H_EBAD;
    int hoid, c2, p_mo;
    double hcx, hcy;
    if (e4h_host_info(host, layer, &hoid, &c2, &p_mo, &hcx, &hcy))
        return E4H_EBAD;
    rc = h9e4h_unfold_point(&px, &py, g, hoid);
    if (rc) return rc;
    double wx, wy;
    h9e4h_to_frame(px, py, hcx, hcy, p_mo, c2, layer, &wx, &wy);
    int half;
    uint8_t digits[28];
    rc = h9e4h_descend(wx, wy, hoid, c2, depth, &half, digits);
    if (rc) return rc;
    uint8_t nib[32];
    h9a_unpack(host, nib);
    nib[layer + 1] = E4H_MARK;
    nib[layer + 2] = (uint8_t)half;
    for (int i = 0; i < depth; ++i) nib[layer + 3 + i] = digits[i];
    h9a_pack(nib, out);
    return 0;
}

static int e4h_decode_one(const uint8_t u[16], int sphere, int partner,
                          double *lon, double *lat) {
    if (!u || !lon || !lat) return E4H_EBAD;
    uint8_t nib[32], hostnib[32], digits[28], host[16];
    h9a_unpack(u, nib);
    int layer, half, nd;
    int rc = h9e4h_split_nibs(nib, hostnib, &layer, &half, digits, &nd);
    if (rc) return rc;
    h9a_pack(hostnib, host);
    int hoid, c2, p_mo;
    double hcx, hcy;
    if (e4h_host_info(host, layer, &hoid, &c2, &p_mo, &hcx, &hcy))
        return E4H_EBAD;
    /* partner probe: the mirror across the trapezoid's long side */
    const double pr = partner ? -E4H_CEN_RE : E4H_CEN_RE;
    double rr, ri;
    rc = h9e4h_compose(half, digits, nd, hoid, c2, pr, E4H_CEN_IM, &rr, &ri);
    if (rc) return rc;
    double zx, zy;
    h9e4h_from_frame(rr, ri, hcx, hcy, p_mo, c2, layer, &zx, &zy);
    int oid = hoid;
    h9e4h_fold(&zx, &zy, &oid);
    rc = sphere ? hex9_unproject_sphere(zx, zy, oid, lon, lat)
                : hex9_unproject(zx, zy, oid, lon, lat);
    return rc ? E4H_EBAD : 0;
}

extern "C" int hex9_e4h_encode(double lon, double lat, int layer, int depth,
                               uint8_t out_uuid[16]) {
    return e4h_encode_one(lon, lat, layer, depth, 0, out_uuid);
}
extern "C" int hex9_e4h_encode_sphere(double lon, double lat, int layer,
                                      int depth, uint8_t out_uuid[16]) {
    return e4h_encode_one(lon, lat, layer, depth, 1, out_uuid);
}
extern "C" int hex9_e4h_encode_many(const double *lon, const double *lat,
                                    size_t n, int layer, int depth,
                                    uint8_t *out_uuid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= e4h_encode_one(lon[i], lat[i], layer, depth, 0,
                             out_uuid + (size_t)i * 16);
    return rc;
}
extern "C" int hex9_e4h_encode_many_sphere(const double *lon, const double *lat,
                                           size_t n, int layer, int depth,
                                           uint8_t *out_uuid) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= e4h_encode_one(lon[i], lat[i], layer, depth, 1,
                             out_uuid + (size_t)i * 16);
    return rc;
}
extern "C" int hex9_e4h_decode(const uint8_t uuid[16],
                               double *lon, double *lat) {
    return e4h_decode_one(uuid, 0, 0, lon, lat);
}
extern "C" int hex9_e4h_decode_sphere(const uint8_t uuid[16],
                                      double *lon, double *lat) {
    return e4h_decode_one(uuid, 1, 0, lon, lat);
}
extern "C" int hex9_e4h_decode_many(const uint8_t *uuid, size_t n,
                                    double *lon, double *lat) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= e4h_decode_one(uuid + (size_t)i * 16, 0, 0, &lon[i], &lat[i]);
    return rc;
}
extern "C" int hex9_e4h_decode_many_sphere(const uint8_t *uuid, size_t n,
                                           double *lon, double *lat) {
    int rc = 0;
    const ptrdiff_t N = (ptrdiff_t)n;
    #pragma omp parallel for schedule(static) reduction(|:rc)
    for (ptrdiff_t i = 0; i < N; ++i)
        rc |= e4h_decode_one(uuid + (size_t)i * 16, 1, 0, &lon[i], &lat[i]);
    return rc;
}
extern "C" int hex9_e4h_partner(const uint8_t uuid[16],
                                double *lon, double *lat) {
    return e4h_decode_one(uuid, 0, 1, lon, lat);
}
extern "C" int hex9_e4h_partner_sphere(const uint8_t uuid[16],
                                       double *lon, double *lat) {
    return e4h_decode_one(uuid, 1, 1, lon, lat);
}

extern "C" int hex9_e4h_split(const uint8_t uuid[16], uint8_t out_host[16],
                              int *half, uint8_t digits[28], int *ndigits) {
    if (!uuid || !out_host || !half || !digits || !ndigits) return E4H_EBAD;
    uint8_t nib[32], hostnib[32];
    int layer;
    h9a_unpack(uuid, nib);
    const int rc = h9e4h_split_nibs(nib, hostnib, &layer, half, digits,
                                    ndigits);
    if (rc) return rc;
    h9a_pack(hostnib, out_host);
    return 0;
}

extern "C" int hex9_e4h_bin(const uint8_t uuid[16], int depth,
                            uint8_t out_uuid[16]) {
    if (!uuid || !out_uuid || depth < 0) return E4H_EBAD;
    uint8_t nib[32], hostnib[32], digits[28];
    int layer, half, nd;
    h9a_unpack(uuid, nib);
    const int rc = h9e4h_split_nibs(nib, hostnib, &layer, &half, digits, &nd);
    if (rc) return rc;
    if (depth > nd) return E4H_EBAD;      /* deepening is unreachable */
    for (int i = layer + 3 + depth; i <= H9_NIB_TAIL - 1; ++i) nib[i] = 0x0Fu;
    h9a_pack(nib, out_uuid);
    return 0;
}

extern "C" int hex9_e4h_depth(const uint8_t uuid[16]) {
    if (!uuid) return -1;
    uint8_t nib[32], hostnib[32], digits[28];
    int layer, half, nd;
    h9a_unpack(uuid, nib);
    if (h9e4h_split_nibs(nib, hostnib, &layer, &half, digits, &nd)) return -1;
    return nd;
}

extern "C" int hex9_is_e4h(const uint8_t uuid[16]) {
    return uuid ? h9e4h_is_marked(uuid) : 0;
}

extern "C" int hex9_e4h_label(const uint8_t uuid[16], char *buf,
                              size_t buflen) {
    if (!uuid || !buf) return -1;
    uint8_t nib[32], hostnib[32], digits[28];
    int layer, half, nd;
    h9a_unpack(uuid, nib);
    if (h9e4h_split_nibs(nib, hostnib, &layer, &half, digits, &nd)) return -1;
    const int len = (layer + 1) + 2 + nd;      /* body + 'E' + half + digits */
    if ((size_t)len + 1 > buflen) return -1;
    write_label(hostnib, layer, buf);
    buf[layer + 1] = 'E';
    buf[layer + 2] = (char)('0' + half);
    for (int i = 0; i < nd; ++i) buf[layer + 3 + i] = (char)('0' + digits[i]);
    buf[len] = '\0';
    return len;
}

extern "C" int hex9_e4h_parse_label(const char *label, uint8_t out_uuid[16]) {
    if (!label || !out_uuid) return E4H_EBAD;
    const char *e = std::strchr(label, 'E');
    if (!e || e == label) return E4H_EGRAM;
    const size_t blen = (size_t)(e - label);
    if (blen > (size_t)H9_NIB_BODYTOP + 1) return E4H_EGRAM;
    char body[40];
    std::memcpy(body, label, blen);
    body[blen] = '\0';
    uint8_t host[16];
    const int layer = hex9_parse_label(body, host);
    if (layer < 0 || (size_t)layer + 1 != blen) return E4H_EGRAM;
    const char *p = e + 1;
    if (*p != '0' && *p != '1') return E4H_EGRAM;
    const int half = *p - '0';
    ++p;
    uint8_t digits[28];
    int nd = 0;
    for (; *p; ++p) {
        if (*p < '0' || *p > '5' || nd >= 28) return E4H_EGRAM;
        digits[nd++] = (uint8_t)(*p - '0');
    }
    if (layer + 2 + nd > H9_NIB_BODYTOP) return E4H_EGRAM;
    uint8_t nib[32];
    h9a_unpack(host, nib);
    if (nib[layer] == 0x0Fu || nib[layer + 1] != 0x0Fu) return E4H_EGRAM;
    nib[layer + 1] = E4H_MARK;
    nib[layer + 2] = (uint8_t)half;
    for (int i = 0; i < nd; ++i) nib[layer + 3 + i] = digits[i];
    h9a_pack(nib, out_uuid);
    return 0;
}

/* ── Grid verbs: aim / walk_to / vision_cone ────────────────────────────────
 * C twins of the python wheel's hex9.verbs (hhg9 fox-and-rabbits PoC).
 * Cell-first: bin KEYS in and out; full-uuid re-derivation internal (the
 * guaranteed path). Bearings/distances are advisory FP over centroids —
 * they steer, they never mint. The behavioural reference is the python
 * module; the primitives (decode/encode/bin/neighbors/k_disk) are the very
 * same C functions the wheel calls, so the twins differ only in scalar
 * bearing arithmetic. */

#include <array>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

namespace {

constexpr double H9V_DEG = 3.14159265358979323846 / 180.0;

double h9v_wrap180(double a) {
    a = std::fmod(a + 180.0, 360.0);
    if (a < 0) a += 360.0;
    return a - 180.0;
}

/* initial great-circle bearing, degrees clockwise from north */
double h9v_bearing(double lon1, double lat1, double lon2, double lat2) {
    const double p1 = lat1 * H9V_DEG, p2 = lat2 * H9V_DEG;
    const double dl = (lon2 - lon1) * H9V_DEG;
    const double y = std::sin(dl) * std::cos(p2);
    const double x = std::cos(p1) * std::sin(p2)
                   - std::sin(p1) * std::cos(p2) * std::cos(dl);
    double b = std::atan2(y, x) / H9V_DEG;
    b = std::fmod(b, 360.0);
    if (b < 0) b += 360.0;
    return b;
}

/* central angle, degrees */
double h9v_angle(double lon1, double lat1, double lon2, double lat2) {
    const double p1 = lat1 * H9V_DEG, p2 = lat2 * H9V_DEG;
    const double dl = (lon2 - lon1) * H9V_DEG;
    double h = std::sin((p2 - p1) / 2) * std::sin((p2 - p1) / 2)
             + std::cos(p1) * std::cos(p2)
               * std::sin(dl / 2) * std::sin(dl / 2);
    if (h < 0) h = 0;
    if (h > 1) h = 1;
    return 2.0 * std::asin(std::sqrt(h)) / H9V_DEG;
}

using H9Key = std::array<uint8_t, 16>;
struct H9KeyHash {
    size_t operator()(const H9Key &k) const {
        uint64_t a, b;
        std::memcpy(&a, k.data(), 8);
        std::memcpy(&b, k.data() + 8, 8);
        return (size_t)(a * 0x9e3779b97f4a7c15ull ^ b);
    }
};

H9Key h9v_key(const uint8_t *p) { H9Key k; std::memcpy(k.data(), p, 16); return k; }

/* bin key -> centroid (representative point) */
int h9v_centroid(const uint8_t key[16], double *lon, double *lat) {
    return hex9_decode(key, lon, lat);
}

/* bin key -> full uuid (guaranteed re-derivation: bins are not addresses) */
int h9v_full(const uint8_t key[16], uint8_t full[16]) {
    double lon, lat;
    if (hex9_decode(key, &lon, &lat)) return 1;
    return hex9_encode(lon, lat, full);
}

/* neighbours of a key: decode->encode->neighbors; returns count or -1 */
int h9v_neighbors(const uint8_t key[16], int layer, uint8_t out[6 * 16]) {
    uint8_t full[16];
    if (h9v_full(key, full)) return -1;
    return hex9_neighbors(full, layer, out);
}

/* median centroid step to the neighbours — the local hex pitch (numpy
 * median convention: mean of the two middle values when count is even) */
double h9v_pitch(const uint8_t key[16], int layer) {
    uint8_t nbs[6 * 16];
    const int n = h9v_neighbors(key, layer, nbs);
    if (n <= 0) return -1.0;
    double lon0, lat0, d[6];
    if (h9v_centroid(key, &lon0, &lat0)) return -1.0;
    for (int i = 0; i < n; ++i) {
        double lon1, lat1;
        if (h9v_centroid(nbs + (size_t)i * 16, &lon1, &lat1)) return -1.0;
        d[i] = h9v_angle(lon0, lat0, lon1, lat1);
    }
    std::sort(d, d + n);
    return (n & 1) ? d[n / 2] : 0.5 * (d[n / 2 - 1] + d[n / 2]);
}

}  /* namespace */

/* WGS84 first eccentricity squared (2f - f^2, f = 1/298.257223563) */
static const double H9V_E2 = 0.0066943799901413165;

static double h9v_true_conv(double bearing, double lat, bool from_true) {
    const double b = bearing * H9V_DEG;
    const double s = std::sin(lat * H9V_DEG);
    const double mn = (1.0 - H9V_E2) / (1.0 - H9V_E2 * s * s);   /* M/N */
    const double q = from_true ? mn : 1.0 / mn;
    double out = std::atan2(q * std::sin(b), std::cos(b)) / H9V_DEG;
    out = std::fmod(out, 360.0);
    return out < 0 ? out + 360.0 : out;
}

extern "C" double hex9_bearing_from_true(double bearing, double lat) {
    return h9v_true_conv(bearing, lat, true);
}

extern "C" double hex9_bearing_to_true(double bearing, double lat) {
    return h9v_true_conv(bearing, lat, false);
}

extern "C" int hex9_aim(const uint8_t key[16], int layer, double bearing,
                        uint8_t out_key[16]) {
    if (!key || !out_key || layer < 0 || layer > H9_LMAX ||
        !std::isfinite(bearing)) return 1;
    if (h9e4h_is_marked(key)) return 1;
    uint8_t nbs[6 * 16];
    const int n = h9v_neighbors(key, layer, nbs);
    if (n <= 0) return 1;
    double lon0, lat0;
    if (h9v_centroid(key, &lon0, &lat0)) return 1;
    int best = -1;
    double bd = 0.0;
    for (int i = 0; i < n; ++i) {
        double lon1, lat1;
        if (h9v_centroid(nbs + (size_t)i * 16, &lon1, &lat1)) return 1;
        const double diff = std::fabs(h9v_wrap180(
            h9v_bearing(lon0, lat0, lon1, lat1) - bearing));
        if (best < 0 || diff < bd) { bd = diff; best = i; }
    }
    std::memcpy(out_key, nbs + (size_t)best * 16, 16);
    return 0;
}

extern "C" int64_t hex9_walk_to(const uint8_t src_key[16],
                                const uint8_t dest_key[16], int layer,
                                const uint8_t *obstacles, size_t n_obstacles,
                                int64_t max_expand,
                                uint8_t *out_keys, int64_t max_cells) {
    if (!src_key || !dest_key || !out_keys || max_cells <= 0 ||
        layer < 0 || layer > H9_LMAX) return -1;
    if (h9e4h_is_marked(src_key) || h9e4h_is_marked(dest_key)) return -1;
    if (max_expand <= 0) max_expand = 5000;

    std::unordered_set<H9Key, H9KeyHash> block;
    for (size_t i = 0; i < n_obstacles; ++i)
        block.insert(h9v_key(obstacles + i * 16));

    double dlon, dlat, pitch;
    if (h9v_centroid(dest_key, &dlon, &dlat)) return -1;
    pitch = h9v_pitch(src_key, layer);
    if (pitch <= 0) return -1;
    const H9Key src = h9v_key(src_key), dst = h9v_key(dest_key);

    /* the python module's heap order exactly: (f, then g, then key bytes) */
    struct Node { double f; int64_t g; H9Key key; };
    struct NodeGt {
        bool operator()(const Node &a, const Node &b) const {
            if (a.f != b.f) return a.f > b.f;
            if (a.g != b.g) return a.g > b.g;
            return std::memcmp(a.key.data(), b.key.data(), 16) > 0;
        }
    };
    auto h = [&](const H9Key &k) -> double {
        double lon, lat;
        if (h9v_centroid(k.data(), &lon, &lat)) return 1e300;
        return h9v_angle(lon, lat, dlon, dlat) / pitch * 0.95;
    };

    std::priority_queue<Node, std::vector<Node>, NodeGt> open;
    std::unordered_map<H9Key, H9Key, H9KeyHash> came;
    std::unordered_map<H9Key, int64_t, H9KeyHash> g;
    open.push({h(src), 0, src});
    came.emplace(src, src);
    g.emplace(src, 0);

    for (int64_t it = 0; it < max_expand && !open.empty(); ++it) {
        const Node cur = open.top();
        open.pop();
        if (cur.key == dst) {
            std::vector<H9Key> path{cur.key};
            while (!(path.back() == src))
                path.push_back(came.at(path.back()));
            if ((int64_t)path.size() > max_cells) return -1;
            for (size_t i = 0; i < path.size(); ++i)
                std::memcpy(out_keys + i * 16,
                            path[path.size() - 1 - i].data(), 16);
            return (int64_t)path.size();
        }
        if (cur.g > g.at(cur.key)) continue;
        uint8_t nbs[6 * 16];
        const int n = h9v_neighbors(cur.key.data(), layer, nbs);
        if (n <= 0) continue;
        for (int i = 0; i < n; ++i) {
            const H9Key kk = h9v_key(nbs + (size_t)i * 16);
            if (block.count(kk) && !(kk == dst)) continue;
            auto it2 = g.find(kk);
            if (it2 == g.end() || cur.g + 1 < it2->second) {
                g[kk] = cur.g + 1;
                came[kk] = cur.key;
                open.push({(double)(cur.g + 1) + h(kk), cur.g + 1, kk});
            }
        }
    }
    return 0;    /* no path within max_expand */
}

extern "C" int64_t hex9_vision_cone(const uint8_t src_key[16], int layer,
                                    double bearing, double half_angle, int k,
                                    const uint8_t *obstacles,
                                    size_t n_obstacles,
                                    uint8_t *out_keys, int64_t max_cells) {
    if (!src_key || !out_keys || max_cells <= 0 || k < 0 ||
        layer < 0 || layer > H9_LMAX ||
        !std::isfinite(bearing) || !std::isfinite(half_angle)) return -1;
    if (h9e4h_is_marked(src_key)) return -1;

    double slon, slat;
    if (h9v_centroid(src_key, &slon, &slat)) return -1;
    uint8_t full[16];
    if (hex9_encode(slon, slat, full)) return -1;
    const int64_t maxd = hex9_disk_ncells(k);
    std::vector<uint8_t> disk((size_t)maxd * 16);
    const int64_t nd = hex9_k_disk(full, layer, k, disk.data(), maxd);
    if (nd < 0) return -1;
    const H9Key src = h9v_key(src_key);

    std::unordered_set<H9Key, H9KeyHash> block;
    for (size_t i = 0; i < n_obstacles; ++i)
        block.insert(h9v_key(obstacles + i * 16));

    const double pitch = n_obstacles ? h9v_pitch(src_key, layer) : 1.0;
    if (pitch <= 0) return -1;
    const double ax = std::cos(slat * H9V_DEG) * std::cos(slon * H9V_DEG);
    const double ay = std::cos(slat * H9V_DEG) * std::sin(slon * H9V_DEG);
    const double az = std::sin(slat * H9V_DEG);

    std::vector<H9Key> seen{src};
    for (int64_t i = 0; i < nd; ++i) {
        const uint8_t *cell = disk.data() + (size_t)i * 16;
        const H9Key ck = h9v_key(cell);
        if (ck == src) continue;
        double clon, clat;
        if (h9v_centroid(cell, &clon, &clat)) return -1;
        const double br = h9v_bearing(slon, slat, clon, clat);
        if (std::fabs(h9v_wrap180(br - bearing)) > half_angle) continue;
        bool blocked = false;
        if (n_obstacles) {
            /* sight line sampled at a third of the hex pitch, each sample
             * binned; an obstacle blocks unless it IS the source/target */
            const double dist = h9v_angle(slon, slat, clon, clat);
            const double bx = std::cos(clat * H9V_DEG) * std::cos(clon * H9V_DEG);
            const double by = std::cos(clat * H9V_DEG) * std::sin(clon * H9V_DEG);
            const double bz = std::sin(clat * H9V_DEG);
            long ns = (long)std::ceil(dist / (pitch / 3.0));
            if (ns < 2) ns = 2;
            for (long j = 0; j < ns && !blocked; ++j) {
                const double t = ((double)j + 0.5) / (double)ns;
                double px = ax * (1 - t) + bx * t;
                double py = ay * (1 - t) + by * t;
                double pz = az * (1 - t) + bz * t;
                const double nrm = std::sqrt(px * px + py * py + pz * pz);
                px /= nrm; py /= nrm; pz /= nrm;
                double c = pz;
                if (c < -1) c = -1;
                if (c > 1) c = 1;
                const double plat = std::asin(c) / H9V_DEG;
                const double plon = std::atan2(py, px) / H9V_DEG;
                uint8_t sf[16], sb[16];
                if (hex9_encode(plon, plat, sf) ||
                    hex9_bin(sf, layer, sb)) return -1;
                const H9Key sk = h9v_key(sb);
                if (block.count(sk) && !(sk == src) && !(sk == ck))
                    blocked = true;
            }
        }
        if (!blocked) seen.push_back(ck);
    }
    if ((int64_t)seen.size() > max_cells) return -1;
    std::sort(seen.begin(), seen.end(),
              [](const H9Key &a, const H9Key &b) {
                  return std::memcmp(a.data(), b.data(), 16) < 0;
              });
    for (size_t i = 0; i < seen.size(); ++i)
        std::memcpy(out_keys + i * 16, seen[i].data(), 16);
    return (int64_t)seen.size();
}
