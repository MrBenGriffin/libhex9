/* h9_curve.h — the H9 Hamiltonian space-filling curve (36-state transducer).
 *
 * Port of hhg9/h9/curve_tables.py + addressing.x_adr_curve +
 * uuid_address.h9_curve_* (tables verbatim; provenance: the hexsphere
 * machine-closure derivation, h9curve_rowtables_L5c — see the Python
 * module docstring).
 *
 * The curve's tree is the LINEAGE tree — iterated one-generation canonical
 * parents (h9_cell_fold_uuid one layer up; a single deep fold is WRONG on
 * hexagon-band cells) — and a cell's curve address is emitted by a row
 * transducer walked down that chain:
 *
 *     S     = ROOT_STATE[root_digit]        # 12 L0 hexagons
 *     slot  = AXIOM_POS[root_digit]         # curve order of the roots
 *     per layer:  r = T1[S][foster][digit]  # base-9 rank digit
 *                 S = T2[S][r]
 *
 * where `digit` is the cell's own body nibble at that layer and `foster`
 * marks a lineage parent whose address is not the cell's body prefix.
 * T1 is PARTIAL: -1 entries cannot occur on canonical lineage chains, so a
 * -1 lookup is a totality violation (non-canonical input), never a table
 * gap.  T2 is total.  Each state's T1 row is rank-bijective, which makes
 * the curve address symbolically invertible down to the (foster, digit)
 * string; only foster address-SPELLINGS need the geometric child oracle.
 *
 * Packed curve-uuid layout (32 nibbles, mirrors Python h9_curve_uuid):
 *   nibble 0    = 0xC  (type marker — an h9-uuid's nibble 0 is the root
 *                 digit 0..11, so this is positionally unambiguous)
 *   nibble 1    = axiom slot 0..11
 *   nibble 1+k  = base-9 rank at layer k (k = 1..L)
 *   the rest    = 0xF sentinels (full depth L30 lands exactly on nibble 31)
 * Byte order at a fixed layer IS curve order, and unlike an h9-uuid body a
 * curve-uuid truncates EXACTLY: dropping rank digits gives the lineage
 * ancestor's curve address (index/9 == parent index).
 */
#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

#include "h9_kring.h"

namespace h9curve {

static const uint8_t H9CURVE_MARK = 0x0C;

/* ── Tables (verbatim string literals from hhg9/h9/curve_tables.py) ───────── */

static const char *H9C_AXIOM      = "0457362ab918";
static const char *H9C_ROOT_STATE = "141420535320";

static const char *H9C_T1[36][2] = {
    {"285136074", "........."},
    {"803652714", "........."},
    {"825174036", "........."},
    {"714263805", "........."},
    {"063714852", "........."},
    {"075236184", "........."},
    {"085236174", "........."},
    {"0.371.8.2", ".6...4.5."},
    {"8.365.7.4", ".0...2.1."},
    {"8.517.0.6", ".2...4.3."},
    {"174625083", "........."},
    {"306147285", "........."},
    {"741562830", "........."},
    {"138265047", "........."},
    {"6.375.8.4", ".0...2.1."},
    {"471632580", "........."},
    {"528361407", "........."},
    {"813652704", "........."},
    {"147326058", "........."},
    {"2.513.0.4", ".8...6.7."},
    {"750623841", "........."},
    {"0.523.1.4", ".8...6.7."},
    {"582741603", "........."},
    {"7.426.8.5", ".1...3.0."},
    {"7.156.8.0", ".4...2.3."},
    {"1.462.0.3", ".7...5.8."},
    {"360527481", "........."},
    {"417256308", "........."},
    {"603752814", "........."},
    {"1.732.0.8", ".4...6.5."},
    {"5.836.4.7", ".2...1.0."},
    {"4.725.3.8", ".1...6.0."},
    {"3.052.4.1", ".6...7.8."},
    {"4.163.5.0", ".7...2.8."},
    {"0.523.1.4", ".7...6.8."},
    {"8.365.7.4", ".1...2.0."},
};

static const char *H9C_T2[36] = {
    "820376375", "h9a19dc75", "8a19a19ab", "e6gf19a1j", "m3763763l",
    "h9ik76h2j", "h9ik76375", "m37l376nl", "h9a89do75", "8p19a89ab",
    "e6376rq1j", "m3763763l", "e6376rq1j", "e4519dc75", "h9a89as7l",
    "e6gf19a1j", "820376375", "e4519dc75", "e6gf19a1j", "89037l375",
    "h9ik76h2j", "h9tk7l375", "8a19a19ab", "e6uf89a1j", "el376vq1j",
    "e637lrw1j", "h9a19as4l", "e6376rq1j", "h9a19as4l", "e6gx19a8j",
    "820n7637y", "e6n76rq8j", "z9a19ps4l", "elgf19p1j", "h9tk7lh2j",
    "e4589do75",
};

struct H9CurveTables {
    uint8_t axiom[12];       /* curve position -> root digit */
    uint8_t axiom_pos[12];   /* root digit -> curve position */
    uint8_t root_state[12];  /* root digit -> start state */
    int8_t  t1[36][2][9];    /* [state][foster][digit] -> rank, -1 unreachable */
    uint8_t t2[36][9];       /* [state][rank] -> next state (total) */
    uint8_t inv_f[36][9];    /* [state][rank] -> foster (T1 rows are rank-bijective) */
    uint8_t inv_d[36][9];    /* [state][rank] -> digit */
};

static inline const H9CurveTables &curve_tables(void) {
    static const H9CurveTables T = []() {
        H9CurveTables t;
        auto b36 = [](char c) -> uint8_t {
            return (uint8_t)(c <= '9' ? c - '0' : c - 'a' + 10);
        };
        for (int i = 0; i < 12; ++i) {
            t.axiom[i] = b36(H9C_AXIOM[i]);          /* hex digit, b36 == b16 here */
            t.axiom_pos[t.axiom[i]] = (uint8_t)i;
            t.root_state[i] = b36(H9C_ROOT_STATE[i]);
        }
        for (int s = 0; s < 36; ++s)
            for (int r = 0; r < 9; ++r) {
                t.t2[s][r] = b36(H9C_T2[s][r]);
                t.inv_f[s][r] = 0xFFu;
                t.inv_d[s][r] = 0xFFu;
            }
        for (int s = 0; s < 36; ++s)
            for (int f = 0; f < 2; ++f)
                for (int d = 0; d < 9; ++d) {
                    const char c = H9C_T1[s][f][d];
                    const int8_t r = (c == '.') ? (int8_t)-1 : (int8_t)(c - '0');
                    t.t1[s][f][d] = r;
                    if (r >= 0) {
                        t.inv_f[s][r] = (uint8_t)f;
                        t.inv_d[s][r] = (uint8_t)d;
                    }
                }
        return t;
    }();
    return T;
}

/* ── Small helpers ────────────────────────────────────────────────────────── */

/* Deepest non-sentinel body nibble (the bin layer). */
static inline int curve_body_layer(const uint8_t nib[32]) {
    int L = 0;
    for (int l = 1; l <= H9_NIB_BODYTOP; ++l) {
        if (nib[l] == 0x0Fu) break;
        L = l;
    }
    return L;
}

static inline bool is_curve_uuid(const uint8_t uuid[16]) {
    return (uuid[0] >> 4) == H9CURVE_MARK;
}

/* Layer + validity scan of an h9 input; also reports whether it carries a
 * reversible (full-uuid) tail and must be canonically re-binned by the
 * caller (the ABI owns the guaranteed full->bin path, canonical_bin).
 * Returns the layer, or -1 on a malformed root/body nibble. */
static inline int curve_input_layer(const uint8_t in[16], bool *needs_rebin) {
    uint8_t nib[32];
    h9a_unpack(in, nib);
    if (nib[0] > 11) return -1;
    const int L = curve_body_layer(nib);
    for (int l = 1; l <= L; ++l)
        if (nib[l] > 8) return -1;
    bool reversible = ((nib[H9_NIB_TAIL] >> 3) & 1u) != 0;   /* p_mo tail bit */
#if H9_HAS_HTERM
    if (nib[H9_NIB_HTERM] != 0x0Fu) reversible = true;       /* legacy full uuid */
#endif
    if (needs_rebin) *needs_rebin = reversible;
    return L;
}

/* ── Forward encode: canonical bin -> curve rows ──────────────────────────── */

/* rows[0] = axiom slot, rows[k] = base-9 rank at layer k (k = 1..L).
 * When end_state is non-NULL the terminal transducer state is written (the
 * state a descent from this cell continues in).  Returns false on a
 * malformed chain or a T1 miss (non-canonical input — the totality test). */
static inline bool curve_rows_from_bin(const uint8_t bin[16], int L,
                                       uint8_t *rows, int *end_state) {
    const H9CurveTables &T = curve_tables();
    /* Lineage chain: iterated ONE-generation canonical folds (deep ownership
     * differs on hexagon-band cells, so one deep fold would be wrong). */
    uint8_t chain[H9_NLEVELS][32];
    uint8_t cur[16], par[16];
    std::memcpy(cur, bin, 16);
    h9a_unpack(cur, chain[L]);
    for (int k = L; k >= 1; --k) {
        if (!h9kring::h9_cell_fold_uuid(cur, k - 1, par)) return false;
        h9a_unpack(par, chain[k - 1]);
        std::memcpy(cur, par, 16);
    }
    const uint8_t root = chain[0][0];
    if (root > 11) return false;
    rows[0] = T.axiom_pos[root];
    int state = T.root_state[root];
    for (int k = 1; k <= L; ++k) {
        const int foster = std::memcmp(chain[k], chain[k - 1], (size_t)k) != 0;
        const uint8_t digit = chain[k][k];
        if (digit > 8) return false;
        const int8_t r = T.t1[state][foster][digit];
        if (r < 0) return false;
        rows[k] = (uint8_t)r;
        state = T.t2[state][r];
    }
    if (end_state) *end_state = state;
    return true;
}

/* rows (slot + L ranks) -> packed curve-uuid. */
static inline void curve_pack_rows(const uint8_t *rows, int L, uint8_t out[16]) {
    uint8_t nib[32];
    for (int i = 0; i < 32; ++i) nib[i] = 0x0Fu;
    nib[0] = H9CURVE_MARK;
    for (int k = 0; k <= L; ++k) nib[1 + k] = rows[k];
    h9a_pack(nib, out);
}

/* Packed curve-uuid -> rows.  Returns the layer, or -1 when the input is
 * not a well-formed curve uuid (marker, slot range, rank range, clean
 * sentinel run). */
static inline int curve_unpack(const uint8_t cu[16], uint8_t *rows) {
    uint8_t nib[32];
    h9a_unpack(cu, nib);
    if (nib[0] != H9CURVE_MARK || nib[1] > 11) return -1;
    rows[0] = nib[1];
    int L = 0;
    for (int l = 1; l <= H9_LMAX; ++l) {
        const uint8_t v = nib[1 + l];
        if (v == 0x0Fu) {
            for (int m = l; m <= 30; ++m)
                if (nib[1 + m] != 0x0Fu) return -1;
            return L;
        }
        if (v > 8) return -1;
        rows[l] = v;
        L = l;
    }
    return L;
}

/* ── Curve index as a decimal string (12·9^L exceeds uint64 above L18) ────── */

/* rows -> decimal index string.  buf gets the NUL-terminated decimal numeral
 * (<= 31 digits at L30).  Returns the string length, or -1 if buflen is too
 * small (40 bytes always suffice). */
static inline int curve_index_decimal(const uint8_t *rows, int L,
                                      char *buf, size_t buflen) {
    uint8_t dec[40];                     /* little-endian decimal digits */
    int nd = 0;
    unsigned seed = rows[0];
    do { dec[nd++] = (uint8_t)(seed % 10u); seed /= 10u; } while (seed);
    for (int k = 1; k <= L; ++k) {
        unsigned carry = rows[k];
        for (int i = 0; i < nd; ++i) {
            const unsigned v = (unsigned)dec[i] * 9u + carry;
            dec[i] = (uint8_t)(v % 10u);
            carry = v / 10u;
        }
        while (carry) { dec[nd++] = (uint8_t)(carry % 10u); carry /= 10u; }
    }
    if (buflen < (size_t)nd + 1) return -1;
    for (int i = 0; i < nd; ++i) buf[i] = (char)('0' + dec[nd - 1 - i]);
    buf[nd] = '\0';
    return nd;
}

/* Decimal index string + layer -> rows.  Returns false on a malformed
 * numeral or an index out of range for the layer (slot residue > 11). */
static inline bool curve_rows_from_decimal(const char *s, int L, uint8_t *rows) {
    uint8_t dec[64];                     /* big-endian decimal digits */
    int nd = 0;
    if (!s || !*s) return false;
    for (const char *p = s; *p; ++p) {
        if (*p < '0' || *p > '9' || nd >= (int)sizeof(dec)) return false;
        dec[nd++] = (uint8_t)(*p - '0');
    }
    /* Peel base-9 digits off the low end: repeated big-decimal divmod 9. */
    for (int k = L; k >= 1; --k) {
        unsigned rem = 0;
        for (int i = 0; i < nd; ++i) {
            const unsigned v = rem * 10u + dec[i];
            dec[i] = (uint8_t)(v / 9u);
            rem = v % 9u;
        }
        rows[k] = (uint8_t)rem;
        while (nd > 1 && dec[0] == 0) { std::memmove(dec, dec + 1, (size_t)--nd); }
    }
    unsigned slot = 0;
    for (int i = 0; i < nd; ++i) {
        slot = slot * 10u + dec[i];
        if (slot > 11u) return false;
    }
    rows[0] = (uint8_t)slot;
    return true;
}

/* ── The child oracle: 9 canonical children of a cell ─────────────────────
 *
 * Seed: the layer-(L+1) cell containing the parent's mode-0-interior nudged
 * centroid — the SAME decode + containment descent as h9_cell_bin_at_uuid
 * (h9_kring.h), run one layer DOWN instead of up, returning the lattice
 * identity.  Candidates: k_disk(seed, 3) — 37 cells, a strict superset of
 * the 9 children (the parent spans 3 child widths).  Filter: the EXACT
 * address-space fold h9_cell_fold_uuid == parent, so correctness never
 * rests on the float seed (it only has to land within the disk, which
 * holds with huge margin against the nudge's ~3^-L noise). */

/* Lattice identity of the layer-(top_l + gens) cell containing the input
 * cell's interior point.  Body of h9_cell_bin_at_uuid descending `gens`
 * generations DOWN and returning the identity; see that function for the
 * doctrine commentary.  The seed is a geometric candidate only — every
 * caller re-adjudicates with the exact address-space fold. */
static inline bool curve_child_seed(const uint8_t in[16], int gens,
                                    h9kring::H9CellId *out, int *child_layer) {
    uint8_t nib[32];
    h9a_unpack(in, nib);
    if (nib[0] > 11 || gens < 1) return false;
    const uint8_t r_mo = nib[H9_NIB_TAIL] & 1u;
    const uint8_t p_c2 = (nib[H9_NIB_TAIL] >> 1) & 3u;
    int top_l = 0;
    for (int l = 1; l <= H9_NIB_BODYTOP; ++l) {
        if (nib[l] == 0x0Fu) break;
        top_l = l;
    }
    if (top_l + gens > H9_LMAX) return false;   /* would exceed lmax */
    const int layer = top_l + gens;
    if (child_layer) *child_layer = layer;

    uint8_t cids[H9_NLEVELS];
    {
        uint8_t c_mo = 0, c2 = p_c2;
        uint8_t rids[H9_NLEVELS] = {};
        rids[0] = r_mo;
        for (int l = top_l; l >= 1; --l) {
            const uint8_t *e = H9_HEX_REG[nib[l]][c_mo][c2];
            rids[l] = e[0]; c_mo = e[1]; c2 = e[2];
        }
        for (int i = 0; i <= top_l; ++i) cids[i] = H9_RID2CELL[rids[i]];
    }

    static const double cen_sx[3] = {  1.0 * H9_UV_U1,  1.0 * H9_UV_U1, -2.0 * H9_UV_U1 };
    static const double cen_sy[3] = {  1.0 * H9_UV_V3, -1.0 * H9_UV_V3,  0.0 };
    double org_x = 0.0, org_y = 0.0;
    double cen_x = cen_sx[p_c2], cen_y = cen_sy[p_c2];
    for (int l = top_l - 1; l >= 0; --l) {
        double ox, oy;
        h9a_cid_offset(cids[l + 1], &ox, &oy);
        org_x = org_x / 3.0 + ox;  org_y = org_y / 3.0 + oy;
        cen_x = cen_x / 3.0 + ox;  cen_y = cen_y / 3.0 + oy;
    }

    const double NUDGE = 0.10;
    double cx = cen_x + NUDGE * (org_x - cen_x);
    double cy = cen_y + NUDGE * (org_y - cen_y);
    const int oid = (int)H9_L0HEX_BACK[nib[0]][r_mo][0];
    h9grid::h9_preamble_nudge(&cx, &cy, (int)H9_OID_MO[oid]);

    const double R3 = std::sqrt(3.0);
    int     p_mo = (int)H9_OID_MO[oid];
    int     leaf_mode = p_mo;
    int64_t ia = 0, ib = 0;
    int64_t iscale = h9kring::pow3(layer);
    uint8_t drids[H9_NLEVELS + 1];
    drids[0] = (uint8_t)p_mo;
    for (int l = 0; l <= layer; ++l) {
        const double  xdot = cx * R3;
        const uint8_t cid  = h9grid::classify_band(xdot, cy, p_mo);
        const double  *offx = (p_mo == 1) ? H9.UP_X    : H9.DN_X;
        const double  *offy = (p_mo == 1) ? H9.UP_Y    : H9.DN_Y;
        const uint8_t *ctab = (p_mo == 1) ? H9_UP_CIDS : H9_DN_CIDS;
        const int     *mtab = (p_mo == 1) ? H9_UP_MODE : H9_DN_MODE;
        const int   (*ofs)[2] = (p_mo == 1) ? H9_MODE1_OFS : H9_MODE0_OFS;
        int j = -1;
        for (int k = 0; k < 9; ++k) if (ctab[k] == cid) { j = k; break; }
        if (j < 0) {
            j = h9grid::h9_recover_cell(&cx, &cy, p_mo, ctab);
            if (j < 0) {
                double best = 1e300; int bj = 0;
                for (int k = 0; k < 9; ++k) {
                    const double dx = cx - offx[k], dy = cy - offy[k];
                    const double d2 = dx*dx + dy*dy;
                    if (d2 < best) { best = d2; bj = k; }
                }
                j = bj;
            }
        }
        drids[l + 1] = H9_CELL2RID[ctab[j]];
        if (l < layer) {
            iscale /= 3;
            ia += (int64_t)ofs[j][0] * iscale;
            ib += (int64_t)ofs[j][1] * iscale;
        } else {
            leaf_mode = p_mo;
        }
        cx = (cx - offx[j]) * 3.0;
        cy = (cy - offy[j]) * 3.0;
        p_mo = mtab[j];
    }

    const uint8_t *se = H9_REG_HEX[drids[layer - 1] & 1u][drids[layer]][drids[layer + 1]];
    if (se[0] == 0xFFu) return false;
    const int sc = ((int)se[1] + leaf_mode) % 3;
    const int64_t s = h9kring::pow3(layer);
    h9kring::H9CellId id;
    if (!h9kring::resolve_frames(oid, layer, s,
                                 ia + h9kring::H9KR_C2_DU[sc],
                                 ib + h9kring::H9KR_C2_DV[sc],
                                 /*include_own=*/true, &id))
        return false;
    if (id.ext && id.oid != oid && H9_OID_NB[id.oid][id.c2] != oid) {
        h9kring::H9CellId partner;
        const int64_t pu = id.ia + h9kring::H9KR_C2_DU[id.c2];
        const int64_t pv = id.ib + h9kring::H9KR_C2_DV[id.c2];
        if (h9kring::resolve_frames(id.oid, layer, s, pu, pv,
                                    /*include_own=*/false, &partner)
            && partner.ext)
            id = partner;
    }
    *out = id;
    return true;
}

/* The 9 canonical children of a layer-L bin, as bins at L+1, UNORDERED.
 * Returns true and fills out_child[0..8]; false on error. */
static inline bool curve_children(const uint8_t parent[16], int parent_layer,
                                  uint8_t out_child[9][16]) {
    h9kring::H9CellId seed;
    int k = 0;
    if (!curve_child_seed(parent, 1, &seed, &k)) return false;
    std::vector<std::pair<h9kring::H9CellId, int>> disk;
    for (int radius = 3; radius <= 4; ++radius) {
        disk.clear();
        if (!h9kring::k_disk(seed, k, radius, disk)) return false;
        int n = 0;
        for (const auto &pr : disk) {
            uint8_t cu[16], pu[16];
            h9kring::identity_to_uuid(pr.first, k, cu);
            if (parent_layer == 0
                    ? h9kring::h9_cell_fold_uuid(cu, 0, pu)
                    : h9kring::h9_cell_fold_uuid(cu, parent_layer, pu)) {
                if (std::memcmp(pu, parent, 16) == 0) {
                    if (n < 9) std::memcpy(out_child[n], cu, 16);
                    ++n;
                }
            }
        }
        if (n == 9) return true;
        if (n > 9) return false;         /* impossible by doctrine; bail */
    }
    return false;
}

/* ── Constructive inverse: curve-uuid -> canonical bin ────────────────────
 *
 * Forward-fit from the root (mirror of Python h9_curve_decode): the T2
 * state thread is known from the ranks alone; each rank inverts to
 * (foster, digit) through the rank-bijective T1 row; the matching child is
 * selected from the parent's 9 canonical children by that (foster, digit)
 * — unique by the same bijection.  O(L) with one child-oracle call per
 * level. */
static inline bool curve_decode_to_bin(const uint8_t cu[16], uint8_t out[16]) {
    const H9CurveTables &T = curve_tables();
    uint8_t rows[H9_NLEVELS];
    const int L = curve_unpack(cu, rows);
    if (L < 0) return false;
    const uint8_t root = T.axiom[rows[0]];
    uint8_t nib[32];
    for (int i = 0; i < 32; ++i) nib[i] = 0x0Fu;
    nib[0] = root;
    nib[H9_NIB_TAIL] = (uint8_t)((H9_L0HEX_BACK[root][0][1] & 3u) << 1);
    uint8_t cur[16];
    h9a_pack(nib, cur);
    int state = T.root_state[root];
    uint8_t cur_nib[32];
    h9a_unpack(cur, cur_nib);
    for (int k = 1; k <= L; ++k) {
        const uint8_t r = rows[k];
        const uint8_t want_f = T.inv_f[state][r];
        const uint8_t want_d = T.inv_d[state][r];
        if (want_f == 0xFFu) return false;
        uint8_t kids[9][16];
        if (!curve_children(cur, k - 1, kids)) return false;
        bool found = false;
        for (int j = 0; j < 9; ++j) {
            uint8_t kn[32];
            h9a_unpack(kids[j], kn);
            const uint8_t foster = std::memcmp(kn, cur_nib, (size_t)k) != 0;
            if (foster == want_f && kn[k] == want_d) {
                std::memcpy(cur, kids[j], 16);
                std::memcpy(cur_nib, kn, 32);
                found = true;
                break;
            }
        }
        if (!found) return false;
        state = T.t2[state][r];
    }
    std::memcpy(out, cur, 16);
    return true;
}

/* The 9 canonical children of a layer-L bin IN CURVE-RANK ORDER (the
 * deterministic order the curve induces; slot j holds the rank-j child).
 * `state` is the parent's terminal transducer state (curve_rows_from_bin).
 * Returns false on an oracle failure or a rank-bijection violation. */
static inline bool curve_children_ranked(const uint8_t cell[16], int layer,
                                         int state, uint8_t out[9][16]) {
    const H9CurveTables &T = curve_tables();
    uint8_t kids[9][16];
    if (!curve_children(cell, layer, kids)) return false;
    uint8_t pn[32], kn[32];
    h9a_unpack(cell, pn);
    int by_rank[9];
    for (int r = 0; r < 9; ++r) by_rank[r] = -1;
    for (int j = 0; j < 9; ++j) {
        h9a_unpack(kids[j], kn);
        const int foster = std::memcmp(kn, pn, (size_t)(layer + 1)) != 0;
        const uint8_t digit = kn[layer + 1];
        if (digit > 8) return false;
        const int8_t r = T.t1[state][foster][digit];
        if (r < 0 || by_rank[r] >= 0) return false;
        by_rank[r] = j;
    }
    for (int r = 0; r < 9; ++r) {
        if (by_rank[r] < 0) return false;
        std::memcpy(out[r], kids[by_rank[r]], 16);
    }
    return true;
}

/* ── Curve-ordered LINEAGE descendant generator ───────────────────────────
 *
 * Emit every layer-`to_layer` LINEAGE descendant of a bin (iterated
 * one-generation children — the curve tree), in curve order — depth-first
 * by rank, the state threaded down, so the emit order IS ascending curve
 * index (and out_curves, when given, is ascending as a byte string).
 * Cost: one child-oracle call per internal node. */
static inline bool curve_cells_rec(const uint8_t cell[16], int layer, int state,
                                   int to_layer, uint8_t *rows,
                                   uint8_t *out_bins, uint8_t *out_curves,
                                   int64_t *count) {
    if (layer == to_layer) {
        std::memcpy(out_bins + *count * 16, cell, 16);
        if (out_curves)
            curve_pack_rows(rows, layer, out_curves + *count * 16);
        ++*count;
        return true;
    }
    const H9CurveTables &T = curve_tables();
    uint8_t kids[9][16];
    if (!curve_children_ranked(cell, layer, state, kids)) return false;
    for (int r = 0; r < 9; ++r) {
        rows[layer + 1] = (uint8_t)r;
        if (!curve_cells_rec(kids[r], layer + 1, T.t2[state][r],
                             to_layer, rows, out_bins, out_curves, count))
            return false;
    }
    return true;
}

/* ── OWNED sub-zone enumeration (the OWNERSHIP relation, downward) ────────
 *
 * Every layer-`to_layer` cell whose canonical deep-fold ancestor
 * (h9_cell_fold_uuid — the OWNERSHIP relation, exact at any depth) is
 * `zone`: exactly 9^g cells, and the owned sets over all zones at the
 * zone's layer PARTITION the layer (cf. OGC API-DGGS issue #108's "owned
 * sub-zone": counts exactly aperture^d, bounded excursion — only rim
 * splits protrude, by their far half).  Distinct from the LINEAGE set
 * (curve_cells_rec) beyond one generation: lineage is transitive with
 * ~1/6 displaced area; ownership is geometrically bounded but not
 * transitive.  The two coincide for a single generation.
 *
 * Candidates: k_disk of radius 3^g + 2 around the interior-point seed
 * (covers the zone's footprint with margin); the filter is the exact
 * fold, so — as everywhere in this module — correctness never rests on
 * the float seed.  Output is CURVE-SORTED (the curve induces the owned
 * subset's order), bins + optional curve uuids.  Returns the count
 * (9^g), or -1 on error / a count violation. */
static inline int64_t curve_owned_cells(const uint8_t zone[16], int zone_layer,
                                        int to_layer,
                                        uint8_t *out_bins, uint8_t *out_curves,
                                        int64_t max_cells) {
    const int g = to_layer - zone_layer;
    if (g < 0 || to_layer > H9_LMAX || !out_bins) return -1;
    int64_t want = 1;
    for (int i = 0; i < g; ++i) want *= 9;
    if (want > max_cells) return -1;
    if (g == 0) {
        std::memcpy(out_bins, zone, 16);
        if (out_curves) {
            uint8_t rows[H9_NLEVELS];
            if (!curve_rows_from_bin(zone, zone_layer, rows, nullptr)) return -1;
            curve_pack_rows(rows, zone_layer, out_curves);
        }
        return 1;
    }

    h9kring::H9CellId seed;
    int k = 0;
    if (!curve_child_seed(zone, g, &seed, &k) || k != to_layer) return -1;
    int64_t radius = 2;
    for (int i = 0; i < g; ++i) radius *= 3;
    radius = radius / 2 + 2;                       /* 3^g + 2 */
    std::vector<std::pair<h9kring::H9CellId, int>> disk;
    if (!h9kring::k_disk(seed, to_layer, (int)radius, disk,
                         /*max_cells=*/want * 8 + 64))
        return -1;

    std::vector<std::array<uint8_t, 16>> owned;
    owned.reserve((size_t)want);
    for (const auto &pr : disk) {
        uint8_t cu[16], pu[16];
        h9kring::identity_to_uuid(pr.first, to_layer, cu);
        if (h9kring::h9_cell_fold_uuid(cu, zone_layer, pu)
                && std::memcmp(pu, zone, 16) == 0) {
            std::array<uint8_t, 16> a;
            std::memcpy(a.data(), cu, 16);
            owned.push_back(a);
        }
    }
    if ((int64_t)owned.size() != want) return -1;  /* 9^g by doctrine */

    /* Curve-sort: compute each owned cell's curve address; byte order at
     * one layer IS curve order.  (Owned cells are not one lineage subtree,
     * so each needs its own chain walk — O(9^g · L) table lookups.) */
    std::vector<std::array<uint8_t, 16>> curves((size_t)want);
    std::vector<int> order((size_t)want);
    for (int64_t i = 0; i < want; ++i) {
        uint8_t rows[H9_NLEVELS];
        if (!curve_rows_from_bin(owned[(size_t)i].data(), to_layer, rows, nullptr))
            return -1;
        curve_pack_rows(rows, to_layer, curves[(size_t)i].data());
        order[(size_t)i] = (int)i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return std::memcmp(curves[(size_t)a].data(), curves[(size_t)b].data(), 16) < 0;
    });
    for (int64_t i = 0; i < want; ++i) {
        std::memcpy(out_bins + i * 16, owned[(size_t)order[(size_t)i]].data(), 16);
        if (out_curves)
            std::memcpy(out_curves + i * 16, curves[(size_t)order[(size_t)i]].data(), 16);
    }
    return want;
}

}  /* namespace h9curve */
