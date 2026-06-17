/* h9_kring.h — symbolic neighbour / k-ring / k-disk for the H9 mesh.
 *
 * Built on the seam algebra derived and exhaustively validated by
 * tools/kring_probe.cpp against h9grid::enumerate_global geometry
 * (L1..L4, 472k steps: every neighbour of every cell resolves uniquely
 * and matches geometric edge-adjacency; zero ambiguity).
 *
 * Cell identity & lattice
 * ───────────────────────
 *   Every cell at `layer` has exactly one lattice identity (oid, c2, ia, ib):
 *   a mode-0 leaf supercell origin (ia, ib) in oid's integer-UV frame at
 *   scale s = 3^layer, plus the hex slot c2 (0..2).  Centroids sit at
 *   origin + {(1,1),(1,-1),(-2,0)}[c2]; the 6 hex neighbours are the
 *   centroid steps (0,±2),(±3,±1).  Tree membership of a candidate origin
 *   is decided by exact integer descent convergence (cids_from_iauv's walk
 *   landing precisely on the origin at a mode-0 leaf) — no floating point.
 *
 * Octant seams
 * ────────────
 *   A step that leaves oid's supercell tree resolves in an edge-adjacent
 *   octant frame.  Crossing edge e (0..2 — the seam's c2 index) of `oid`:
 *     p' = R[e]·p + s·T[mo][e],     target octant = H9_OID_NB[oid][e]
 *   with R[0] = I and R[1], R[2] = ∓120° lattice rotations (all det +1);
 *   T depends only on the source octant's mode.  Octahedron-vertex corner
 *   regions resolve via the composition of two edge maps.
 *
 *   (hhg9 prior art: grid.py get_canonical_gp pins the on-seam v-flip and
 *   uvc.py hex_step does within-octant stepping; these maps extend both
 *   across the seam — on the seam lines they reduce exactly to the v-flip.)
 *
 * Vertex hexagons (ext)
 * ─────────────────────
 *   The 6 octahedron vertices each carry one hexagon split into two
 *   half-hex cells (ext = true) owned by *diagonal* octants with identical
 *   (c2, ia, ib).  A half-hex has 4 side neighbours plus its partner half:
 *   degree 5.  All other cells have degree 6.
 */
#pragma once

#include <cstdint>
#include <cstring>
#include <set>
#include <array>
#include <vector>

#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"

namespace h9kring {

struct H9CellId {
    int     oid;     /* octant frame 0..7 */
    int     c2;      /* hex slot within the supercell, 0..2 */
    int64_t ia, ib;  /* mode-0 leaf supercell origin, integer UV at 3^layer */
    bool    ext;     /* half-hex of an octahedron-vertex hexagon */
    bool operator==(const H9CellId &o) const {
        return oid == o.oid && c2 == o.c2 && ia == o.ia && ib == o.ib;
    }
};

/* c2 hex centroid offset from supercell origin (= template vertex 6). */
static const int64_t H9KR_C2_DU[3] = { 1, 1, -2 };
static const int64_t H9KR_C2_DV[3] = { 1, -1, 0 };

/* The 6 neighbour steps on the centroid lattice. */
static const int64_t H9KR_STEP[6][2] = {
    {0, 2}, {3, 1}, {3, -1}, {0, -2}, {-3, -1}, {-3, 1}
};

/* Seam maps (probe-derived): crossing edge e of an octant with mode mo,
 *   p' = (R[e]·p)/2 + s·T[mo][e]  into octant H9_OID_NB[oid][e]. */
static const int H9KR_SEAM_R[3][4] = {   /* numerators over /2 */
    { 2,  0,  0,  2 },
    {-1,  3, -1, -1 },
    {-1, -3,  1, -1 },
};
static const int H9KR_SEAM_T[2][3][2] = {
    { {0, -2}, { 3,  1}, {-3,  1} },     /* source octant mode 0 */
    { {0,  2}, {-3, -1}, { 3, -1} },     /* source octant mode 1 */
};

static inline void seam_apply(int e, int mo, int64_t s, int64_t u, int64_t v,
                              int64_t *ou, int64_t *ov) {
    *ou = (H9KR_SEAM_R[e][0]*u + H9KR_SEAM_R[e][1]*v) / 2 + s*H9KR_SEAM_T[mo][e][0];
    *ov = (H9KR_SEAM_R[e][2]*u + H9KR_SEAM_R[e][3]*v) / 2 + s*H9KR_SEAM_T[mo][e][1];
}

static inline int64_t pow3(int layer) {
    int64_t s = 1;
    for (int i = 0; i < layer; i++) s *= 3;
    return s;
}

/* ext flag for a candidate identity — the gated outside_oct test the BFS
 * emit uses (mode-0 octants own the seam halves; mode-1 never set ext). */
static inline bool ext_for(int oid, int c2, int64_t ia, int64_t ib, int64_t s) {
    if (H9_OID_MO[oid] != 0) return false;
    auto outside = [s](int64_t u, int64_t v) -> bool {
        return (v > s) || (u - v > 2*s) || (u + v < -2*s);
    };
    const int64_t u4 = (int64_t)H9_HI[0][c2][4][0] + ia;
    const int64_t v4 = (int64_t)H9_HI[0][c2][4][1] + ib;
    const int64_t u5 = (int64_t)H9_HI[0][c2][5][0] + ia;
    const int64_t v5 = (int64_t)H9_HI[0][c2][5][1] + ib;
    return outside(u4, v4) && outside(u5, v5);
}

/* Exact integer descent from oid's root towards supercell origin (ia, ib).
 * Mirrors cids_from_iauv's walk (classification + NN fallback) but only
 * tracks convergence.  True iff (oid, ia, ib) is a mode-0 leaf supercell of
 * oid's BFS tree — the membership test everything here rests on. */
static inline bool origin_in_tree(int oid, int64_t ia, int64_t ib, int layer) {
    int64_t scale = pow3(layer);
    int     p_mo = (int)H9_OID_MO[oid];
    int64_t ou = 0, ov = 0;
    for (int k = 0; k < layer; k++) {
        const int (*ofs)[2] = (p_mo == 1) ? H9_MODE1_OFS        : H9_MODE0_OFS;
        const int *cm       = (p_mo == 1) ? H9_MODE1_CHILD_MODE : H9_MODE0_CHILD_MODE;
        const uint8_t *ctab = (p_mo == 1) ? H9_UP_CIDS          : H9_DN_CIDS;
        const uint8_t cid   = h9grid::classify_band_int(ia - ou, ib - ov, scale, p_mo);
        const int64_t cs    = scale / 3;
        int j = -1;
        for (int kk = 0; kk < 9; kk++) if (ctab[kk] == cid) { j = kk; break; }
        if (j < 0) {
            int64_t best = INT64_MAX; int bj = 0;
            for (int kk = 0; kk < 9; kk++) {
                const int64_t du = (ia - ou) - (int64_t)ofs[kk][0] * cs;
                const int64_t dv = (ib - ov) - (int64_t)ofs[kk][1] * cs;
                const int64_t d2 = du*du + 3*dv*dv;
                if (d2 < best) { best = d2; bj = kk; }
            }
            j = bj;
        }
        ou += (int64_t)ofs[j][0] * cs;
        ov += (int64_t)ofs[j][1] * cs;
        p_mo = cm[j];
        scale = cs;
    }
    return ou == ia && ov == ib && p_mo == 0;
}

/* Resolve a centroid (cu, cv) in oid's frame to a cell identity.
 * (The probe verified at most one c2 candidate ever converges.) */
static inline bool resolve_centroid(int oid, int64_t cu, int64_t cv, int layer,
                                    int64_t s, H9CellId *out) {
    for (int c2 = 0; c2 < 3; c2++) {
        const int64_t ia = cu - H9KR_C2_DU[c2], ib = cv - H9KR_C2_DV[c2];
        if (origin_in_tree(oid, ia, ib, layer)) {
            out->oid = oid; out->c2 = c2; out->ia = ia; out->ib = ib;
            out->ext = ext_for(oid, c2, ia, ib, s);
            return true;
        }
    }
    return false;
}

/* Resolve the cell whose centroid lies at (tu, tv) in `frame_oid`'s
 * coordinates: first in the frame's own tree (when include_own), then
 * across the 3 edge seams, then around the corners (two-map compositions).
 * The probe showed resolution is unique for every target within one step
 * of an existing cell.  When include_own is false the result can never be
 * the cell at (tu, tv) itself, since all other frames carry a different
 * oid — used for the half-hex partner probe. */
static inline bool resolve_frames(int frame_oid, int layer, int64_t s,
                                  int64_t tu, int64_t tv, bool include_own,
                                  H9CellId *out) {
    if (include_own && resolve_centroid(frame_oid, tu, tv, layer, s, out))
        return true;
    const int moA = (int)H9_OID_MO[frame_oid];
    for (int e = 0; e < 3; e++) {
        const int nb  = H9_OID_NB[frame_oid][e];
        const int moB = (int)H9_OID_MO[nb];
        int64_t u1, v1;
        seam_apply(e, moA, s, tu, tv, &u1, &v1);
        if (resolve_centroid(nb, u1, v1, layer, s, out))
            return true;
        for (int e2 = 0; e2 < 3; e2++) {
            const int nc = H9_OID_NB[nb][e2];
            if (nc == frame_oid) continue;
            int64_t u2, v2;
            seam_apply(e2, moB, s, u1, v1, &u2, &v2);
            if (resolve_centroid(nc, u2, v2, layer, s, out))
                return true;
        }
    }
    return false;
}

/* All neighbours of `id` at `layer`: 6 for ordinary cells, 5 for half-hexes
 * (4 sides + the partner half).  Returns the count. */
static inline int neighbors(const H9CellId &id, int layer, H9CellId out[6]) {
    const int64_t s  = pow3(layer);
    const int64_t cu = id.ia + H9KR_C2_DU[id.c2];
    const int64_t cv = id.ib + H9KR_C2_DV[id.c2];
    int n = 0;
    auto add = [&](const H9CellId &nb) {
        for (int i = 0; i < n; i++) if (out[i] == nb) return;
        out[n++] = nb;
    };
    H9CellId nb;
    /* the partner half of a vertex hexagon: zero step, foreign frames only */
    if (id.ext && resolve_frames(id.oid, layer, s, cu, cv, /*include_own=*/false, &nb))
        add(nb);
    for (int st = 0; st < 6; st++) {
        if (resolve_frames(id.oid, layer, s, cu + H9KR_STEP[st][0], cv + H9KR_STEP[st][1],
                           /*include_own=*/true, &nb))
            add(nb);
    }
    return n;
}

/* ── UUID interface ───────────────────────────────────────────────────────── */

/* Normalise a (full or bin) UUID to the canonical bin at `layer`.
 * Full UUIDs go through h9_bin_uuid.  Bin UUIDs at exactly `layer` pass
 * through; deeper bins are re-binned by walking the registration context
 * down to `layer` (the hex9_label_key pattern).  Returns false when the
 * input is a bin shallower than `layer`.
 *
 * FOSSIL (docs/addressing-doctrine.md, F3): "canonical" is aspirational —
 * the bin→coarser walk shares h9_bin_uuid's split-hex (6/7/8) ancestry
 * failure, and the emitted key tail is the deep-walk context, which is not
 * identity-decodable for meta-bearing cells (F2). Do not build core
 * functionality on this; derive coarse bins from the FULL uuid. */
static inline bool normalize_bin(const uint8_t uuid[16], int layer, uint8_t out[16]) {
    uint8_t nib[32];
    h9a_unpack(uuid, nib);
    if (nib[30] != 0x0Fu) {                 /* full UUID */
        h9_bin_uuid(uuid, layer, out);
        return true;
    }
    int bl = 30;                            /* bin layer = last valid nibble */
    while (bl >= 0 && nib[bl] == 0x0Fu) bl--;
    if (bl < layer) return false;
    if (bl == layer) { std::memcpy(out, uuid, 16); return true; }
    /* walk the key_tail context from bl down to `layer` */
    const uint8_t r_mo = nib[31] & 1u;
    const uint8_t c2_bl = (nib[31] >> 1) & 3u;
    /* bin key_tail stores no c_mo; recover it the way cell_unpack does —
     * try both, prefer the one whose backward pass closes on the L0 hex. */
    const uint8_t expected_c2 = H9_L0HEX_BACK[nib[0]][r_mo][1];
    uint8_t kt_c2 = c2_bl;
    for (int try_cmo = 0; try_cmo < 2; ++try_cmo) {
        uint8_t tc_mo = (uint8_t)try_cmo, tc2 = c2_bl, at_layer_c2 = c2_bl;
        for (int l = bl; l >= 1; --l) {
            if (l == layer) at_layer_c2 = tc2;
            const uint8_t *e = H9_HEX_REG[nib[l]][tc_mo][tc2];
            tc_mo = e[1]; tc2 = e[2];
        }
        if (tc2 == expected_c2 || try_cmo == 1) { kt_c2 = at_layer_c2; break; }
    }
    for (int i = layer + 1; i <= 30; ++i) nib[i] = 0x0Fu;
    nib[31] = (uint8_t)((kt_c2 << 1) | r_mo);
    h9a_pack(nib, out);
    return true;
}

/* Integer offset of a child cid within its parent supercell — the integer
 * mirror of h9a_cid_offset.  cid → offset is well-defined without knowing
 * the parent mode: the 6 cids shared by both child tables carry the same
 * offset in each. */
static inline bool cid_offset_int(uint8_t cid, int *u, int *v) {
    for (int j = 0; j < 9; j++)
        if (H9_UP_CIDS[j] == cid) { *u = H9_MODE1_OFS[j][0]; *v = H9_MODE1_OFS[j][1]; return true; }
    for (int j = 0; j < 9; j++)
        if (H9_DN_CIDS[j] == cid) { *u = H9_MODE0_OFS[j][0]; *v = H9_MODE0_OFS[j][1]; return true; }
    return false;
}

/* Recover the lattice identity from a UUID (full, or bin at >= layer).
 *
 * Every cell has TWO supercell homes — a mode-0 leaf and a mode-1 leaf
 * (possibly in different octant frames) — and encoders may produce either
 * flavour (uuid_from_iauv emits the canonical mode-0 form; the default
 * beam/NN point encoder can land on the mode-1 form for seam cells).
 * Decoding is flavour-blind: the backward registration pass plus per-cid
 * integer offsets (the integer mirror of cell_unpack) give the leaf
 * supercell origin exactly; centroid = origin + the key_tail c2 slot
 * offset (slot identity confirmed across L1..L4); the centroid resolves
 * to the canonical mode-0 identity through the frame cascade.
 *
 * Vertex hexagons: the two half-hexes share one centroid, so resolution
 * alone is ambiguous there.  Each of the 4 octants around the vertex
 * belongs to exactly one half (its own, or the one whose ext-reflection
 * it hosts: H9_OID_NB[half.oid][half.c2]); the encoding octant picks the
 * half, and we swap to the partner when it points the other way. */
static inline bool identity_from_uuid(const uint8_t uuid[16], int layer,
                                      H9CellId *out) {
    if (layer < 0 || layer > 29) return false;   /* layer 0 = the 12 L0 cells */

    uint8_t nib[32];
    h9a_unpack(uuid, nib);
    const uint8_t r_mo = nib[31] & 1u;
    if (nib[0] > 11) return false;
    const int oid = (int)H9_L0HEX_BACK[nib[0]][r_mo][0];

    uint8_t rids[32] = {};
    rids[0] = r_mo;
    int c2;                                       /* leaf slot at `layer` */
    int leaf_mode = 0;                            /* leaf parent mode at `layer` */

    if (nib[30] != 0x0Fu) {
        /* Full UUID — the reversible tail carries the complete meta
         * (tail.py: p_mo in bit 3, p_c2 in bits 1-2, r_mo in bit 0; h_term
         * in nibble[30]), so the registration walk is DETERMINISTIC, the
         * mirror of Python h9_dec/hex_digits_reg: seed (c_mo, c2) from the
         * tail and walk L29 -> 1.  The context (c2) on arriving at `layer`
         * is the leaf slot — the same value h9_bin_uuid emits as the bin
         * key c2.  (The former path normalised to a bin first, discarding
         * p_mo, then guessed it back via L0-closure — the F2 mis-location
         * of meta-bearing cells in docs/addressing-doctrine.md.) */
        uint8_t c_mo = (nib[31] >> 3) & 1u;
        uint8_t c2w  = (nib[31] >> 1) & 3u;
        int slot = (int)c2w;                      /* layer == 29 case */
        int smode = (int)c_mo;
        for (int l = 29; l >= 1; --l) {
            if (l == layer) { slot = (int)c2w; smode = (int)c_mo; }
            if (nib[l] > 8) return false;
            const uint8_t *e = H9_HEX_REG[nib[l]][c_mo][c2w];
            if (e[0] == 0xFFu) return false;
            if (l <= layer) rids[l] = e[0];
            c_mo = e[1]; c2w = e[2];
        }
        if (layer == 0) { slot = (int)c2w; smode = (int)c_mo; }
        c2 = slot;
        leaf_mode = smode;
    } else {
        /* Bin UUID — FOSSIL path (docs/addressing-doctrine.md, F3): the key
         * tail has no p_mo, so it is recovered with cell_unpack's two-
         * candidate L0-closure trick — unguaranteed at split-hex ancestry. */
        uint8_t bin[16];
        if (!normalize_bin(uuid, layer, bin)) return false;
        h9a_unpack(bin, nib);
        c2 = (int)((nib[31] >> 1) & 3u);          /* key_tail c2 == leaf slot */

        const uint8_t expected_c2 = H9_L0HEX_BACK[nib[0]][r_mo][1];
        for (int try_cmo = 0; try_cmo < 2; ++try_cmo) {
            uint8_t tc_mo = (uint8_t)try_cmo, tc2 = (uint8_t)c2;
            uint8_t tmp[32] = {};
            tmp[0] = r_mo;
            bool valid = true;
            for (int l = layer; l >= 1; --l) {
                if (nib[l] > 8) { valid = false; break; }
                const uint8_t *e = H9_HEX_REG[nib[l]][tc_mo][tc2];
                if (e[0] == 0xFFu) { valid = false; break; }
                tmp[l] = e[0]; tc_mo = e[1]; tc2 = e[2];
            }
            if (!valid) { if (try_cmo == 1) return false; continue; }
            if (tc2 == expected_c2 || try_cmo == 1) {
                std::memcpy(rids, tmp, (size_t)(layer + 1));
                break;
            }
        }
    }

    /* leaf supercell origin: Σ cid_offset · 3^(layer-m) over the rid path */
    const int64_t s = pow3(layer);
    int64_t ia = 0, ib = 0, w = s / 3;
    for (int m = 1; m <= layer; m++, w /= 3) {
        int ou, ov;
        if (!cid_offset_int(H9_RID2CELL[rids[m]], &ou, &ov)) return false;
        ia += (int64_t)ou * w;
        ib += (int64_t)ov * w;
    }

    /* centroid → canonical identity through the frame cascade.
     * The slot offset lives in the leaf parent's frame: under a MODE-1 leaf
     * parent the slot c2 takes the mode-0 offset of (c2+1) mod 3 — derived
     * empirically by tools/mode1_probe.cpp (clean single-valued table,
     * ~1,700 obs/bucket, consistent across octants and layers 7..10).
     * Without it, addresses whose leaf parent is mode-1 — one half of every
     * hexagon's interior — resolve one cell off. */
    const int sc = (c2 + leaf_mode) % 3;
    const int64_t cu = ia + H9KR_C2_DU[sc];
    const int64_t cv = ib + H9KR_C2_DV[sc];
    if (!resolve_frames(oid, layer, s, cu, cv, /*include_own=*/true, out))
        return false;

    /* vertex-hexagon half disambiguation by encoding octant */
    if (out->ext && out->oid != oid && H9_OID_NB[out->oid][out->c2] != oid) {
        H9CellId partner;
        const int64_t pu = out->ia + H9KR_C2_DU[out->c2];
        const int64_t pv = out->ib + H9KR_C2_DV[out->c2];
        if (resolve_frames(out->oid, layer, s, pu, pv, /*include_own=*/false, &partner)
            && partner.ext)
            *out = partner;
    }
    return true;
}

static inline void identity_to_uuid(const H9CellId &id, int layer, uint8_t out[16]) {
    h9grid::uuid_from_iauv(id.oid, id.c2, id.ia, id.ib, layer, id.ext, out);
}

/* ── k-disk / k-ring (BFS on the symbolic neighbour step) ─────────────────── */

/* All cells within graph distance k of `start` (start itself at distance 0).
 * Appends (id, distance) pairs to `out` in BFS order.  `max_cells` > 0 caps
 * the result; returns false when the cap would be exceeded. */
static inline bool k_disk(const H9CellId &start, int layer, int k,
                          std::vector<std::pair<H9CellId, int>> &out,
                          int64_t max_cells = 0) {
    std::set<std::array<int64_t, 4>> seen;
    auto key = [](const H9CellId &c) {
        return std::array<int64_t, 4>{ (int64_t)c.oid, (int64_t)c.c2, c.ia, c.ib };
    };
    out.clear();
    out.push_back({start, 0});
    seen.insert(key(start));
    size_t head = 0;
    while (head < out.size()) {
        const H9CellId cur = out[head].first;
        const int      d   = out[head].second;
        head++;
        if (d == k) continue;
        H9CellId nbs[6];
        const int n = neighbors(cur, layer, nbs);
        for (int i = 0; i < n; i++) {
            if (!seen.insert(key(nbs[i])).second) continue;
            if (max_cells > 0 && (int64_t)out.size() >= max_cells) return false;
            out.push_back({nbs[i], d + 1});
        }
    }
    return true;
}

}  /* namespace h9kring */
