/* e4h.c — E4H structural laws, ABI-level (no corpus).
 *
 * Asserts the tail algebra the docs promise, using only the public surface:
 *   - count law: exactly 2*4^d addresses decode over one host at depth d,
 *     and every one round-trips decode -> encode byte-exactly (the leaf
 *     centroid lies in its own leaf — exact nesting);
 *   - truncation = binning (suffix-local, unlike the a9 digits);
 *   - the nibble budget boundary (layer + 2 + depth <= body top);
 *   - grammar rejects (double marker, bad half, bad digit, digit after pad);
 *   - matched pairs: the partner shares the final class digit and the
 *     partner-of-partner is the original address (involution);
 *   - the 0xE marker guard on every h9/curve entry point.
 */
#include "hex9_c.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); fails++; } \
} while (0)

static void unpack(const uint8_t u[16], uint8_t nib[32])
{
    for (int i = 0; i < 16; ++i) { nib[2*i] = u[i] >> 4; nib[2*i+1] = u[i] & 0x0F; }
}

static void pack(const uint8_t nib[32], uint8_t u[16])
{
    for (int i = 0; i < 16; ++i) u[i] = (uint8_t)((nib[2*i] << 4) | nib[2*i+1]);
}

/* enumerate all (half, digit^d) candidates over base's host; count the ones
 * that decode, and round-trip each through encode */
static void count_law(const uint8_t base[16], int layer, int depth)
{
    uint8_t nib[32];
    unpack(base, nib);
    long count = 0, expect = 2;
    for (int i = 0; i < depth; ++i) expect *= 4;

    long total = 2;
    for (int i = 0; i < depth; ++i) total *= 6;
    for (long c = 0; c < total; ++c) {
        long v = c;
        nib[layer + 2] = (uint8_t)(v % 2); v /= 2;
        for (int i = 0; i < depth; ++i) {
            nib[layer + 3 + i] = (uint8_t)(v % 6);
            v /= 6;
        }
        uint8_t u[16];
        pack(nib, u);
        double lon, lat;
        if (hex9_e4h_decode(u, &lon, &lat)) continue;
        count++;
        uint8_t r[16];
        CHECK(hex9_e4h_encode(lon, lat, layer, depth, r) == 0,
              "re-encode rc (c=%ld)", c);
        CHECK(memcmp(u, r, 16) == 0, "decode->encode moved (c=%ld)", c);
    }
    CHECK(count == expect, "count law L%d D%d: %ld != %ld",
          layer, depth, count, expect);
}

int main(void)
{
    char err[256];
    if (hex9_init(err, sizeof err)) { printf("FAIL: init: %s\n", err); return 1; }

    const double LON = -3.19, LAT = 55.95;      /* edinburgh host */

    /* ── count law + roundtrip, depths 1..3 ── */
    for (int d = 1; d <= 3; ++d) {
        uint8_t u[16];
        CHECK(hex9_e4h_encode(LON, LAT, 6, d, u) == 0, "encode L6 D%d", d);
        count_law(u, 6, d);
    }

    /* ── truncation is binning ── */
    {
        uint8_t u8[16];
        CHECK(hex9_e4h_encode(LON, LAT, 6, 8, u8) == 0, "encode L6 D8");
        CHECK(hex9_e4h_depth(u8) == 8, "depth() != 8");
        for (int d = 0; d <= 8; ++d) {
            uint8_t t[16], e[16];
            CHECK(hex9_e4h_bin(u8, d, t) == 0, "e4h_bin d=%d", d);
            CHECK(hex9_e4h_encode(LON, LAT, 6, d, e) == 0, "encode d=%d", d);
            CHECK(memcmp(t, e, 16) == 0, "truncation != binning at d=%d", d);
        }
        uint8_t t[16];
        CHECK(hex9_e4h_bin(u8, 9, t) != 0, "bin deeper than tail must fail");
    }

    /* ── nibble budget boundary ── */
    {
        uint8_t u[16];
        const int top = hex9_lmax();            /* body top nibble index */
        CHECK(hex9_e4h_encode(LON, LAT, 6, top - 8, u) == 0,
              "budget: L6 full depth");
        CHECK(hex9_e4h_encode(LON, LAT, 6, top - 7, u) != 0,
              "budget: L6 overdepth must fail");
        CHECK(hex9_e4h_encode(LON, LAT, 0, top - 2, u) == 0,
              "budget: L0 full depth");
        CHECK(hex9_e4h_encode(LON, LAT, 0, top - 1, u) != 0,
              "budget: L0 overdepth must fail");
        CHECK(hex9_e4h_encode(LON, LAT, -1, 1, u) != 0, "negative layer");
        CHECK(hex9_e4h_encode(LON, LAT, 6, -1, u) != 0, "negative depth");
    }

    /* ── grammar rejects ── */
    {
        uint8_t u[16], nib[32], m[16];
        double lon, lat;
        CHECK(hex9_e4h_encode(LON, LAT, 6, 4, u) == 0, "encode L6 D4");

        unpack(u, nib); nib[20] = 0xE; pack(nib, m);       /* second marker */
        CHECK(hex9_e4h_decode(m, &lon, &lat) != 0, "double marker accepted");

        unpack(u, nib); nib[8] = 2; pack(nib, m);          /* half > 1 */
        CHECK(hex9_e4h_decode(m, &lon, &lat) != 0, "half=2 accepted");

        unpack(u, nib); nib[9] = 6; pack(nib, m);          /* digit > 5 */
        CHECK(hex9_e4h_decode(m, &lon, &lat) != 0, "digit=6 accepted");

        unpack(u, nib); nib[9] = 0xD; pack(nib, m);        /* reserved value */
        CHECK(hex9_e4h_decode(m, &lon, &lat) != 0, "digit=0xD accepted");

        unpack(u, nib); nib[20] = 3; pack(nib, m);         /* digit after pad */
        CHECK(hex9_e4h_decode(m, &lon, &lat) != 0, "digit-after-pad accepted");

        CHECK(hex9_e4h_decode(u, &lon, &lat) == 0, "control decode failed");
    }

    /* ── split / label / parse ── */
    {
        uint8_t u[16], host[16], digits[28], v[16];
        int half, nd;
        char lab[64];
        CHECK(hex9_e4h_encode(LON, LAT, 6, 4, u) == 0, "encode L6 D4");
        CHECK(hex9_e4h_split(u, host, &half, digits, &nd) == 0, "split rc");
        CHECK(nd == 4 && (half == 0 || half == 1), "split fields");
        CHECK(hex9_is_e4h(u) == 1 && hex9_is_e4h(host) == 0, "is_e4h");
        uint8_t hb[16];
        CHECK(hex9_bin(host, 6, hb) == 0 && memcmp(hb, host, 16) == 0,
              "split host is not the canonical L6 bin");
        const int len = hex9_e4h_label(u, lab, sizeof lab);
        CHECK(len == 7 + 2 + 4, "label length %d", len);
        CHECK(lab[7] == 'E' && lab[8] == (char)('0' + half), "label grammar");
        CHECK(hex9_e4h_parse_label(lab, v) == 0 && memcmp(u, v, 16) == 0,
              "label parse roundtrip");
        CHECK(hex9_e4h_parse_label("0031586E9", v) != 0, "bad digit parsed");
        CHECK(hex9_e4h_parse_label("0031586E2", v) != 0, "bad half parsed");
    }

    /* ── matched pairs: shared final digit, involution ── */
    {
        uint8_t u[16], host[16], digits[28], v[16], vh[16], vd[28];
        int half, nd, vhalf, vnd, checked = 0;
        double lon, lat;
        for (int d = 2; d <= 4 && checked < 2; ++d) {
            CHECK(hex9_e4h_encode(LON, LAT, 6, d, u) == 0, "encode pair L6");
            CHECK(hex9_e4h_split(u, host, &half, digits, &nd) == 0, "split");
            if (digits[nd - 1] == 0) continue;      /* centre child: partner
                                                     * is the host's other
                                                     * half — no digit law */
            CHECK(hex9_e4h_partner(u, &lon, &lat) == 0, "partner rc");
            CHECK(hex9_e4h_encode(lon, lat, 6, d, v) == 0, "partner encode");
            CHECK(memcmp(u, v, 16) != 0, "partner == self");
            CHECK(hex9_e4h_split(v, vh, &vhalf, vd, &vnd) == 0, "psplit");
            CHECK(vnd == nd && vd[vnd - 1] == digits[nd - 1],
                  "matched pair final digit");
            uint8_t w[16];
            CHECK(hex9_e4h_partner(v, &lon, &lat) == 0, "partner^2 rc");
            CHECK(hex9_e4h_encode(lon, lat, 6, d, w) == 0, "partner^2 enc");
            CHECK(memcmp(u, w, 16) == 0, "partner involution");
            checked++;
        }
        CHECK(checked >= 1, "no non-centre pair exercised");
    }

    /* ── marker guards: every h9/curve entry point rejects E4H input ── */
    {
        uint8_t u[16], full[16], curve[16], out[16 * 16];
        double lon, lat, cx, cy;
        int oid;
        int64_t ia6[6], ib6[6];
        int64_t ia, ib;
        int o6[6], ext;
        char buf[64];
        CHECK(hex9_e4h_encode(LON, LAT, 6, 2, u) == 0, "guard encode");
        CHECK(hex9_encode(LON, LAT, full) == 0, "plain encode");
        CHECK(hex9_curve(full, curve) == 0, "plain curve");
        CHECK(hex9_is_e4h(full) == 0 && hex9_is_e4h(curve) == 0,
              "is_e4h false positives");

        CHECK(hex9_decode(u, &lon, &lat) != 0, "guard decode");
        CHECK(hex9_decode_sphere(u, &lon, &lat) != 0, "guard decode_sphere");
        CHECK(hex9_decode_many(u, 1, &lon, &lat) != 0, "guard decode_many");
        CHECK(hex9_bin(u, 4, out) != 0, "guard bin");
        CHECK(hex9_bin_many(u, 4, 1, out) != 0, "guard bin_many");
        CHECK(hex9_cell_parent(u, out) != 0, "guard cell_parent");
        CHECK(hex9_cell_ancestor(u, 3, out) != 0, "guard cell_ancestor");
        CHECK(hex9_curve(u, out) != 0, "guard curve");
        CHECK(hex9_curve_decode(u, out) != 0, "guard curve_decode");
        CHECK(hex9_curve_index(u, buf, sizeof buf) < 0, "guard curve_index");
        CHECK(hex9_curve_label(u, buf, sizeof buf) < 0, "guard curve_label");
        CHECK(hex9_cell_children(u, out) != 0, "guard cell_children");
        CHECK(hex9_curve_cells(u, 8, out, NULL, 16) < 0, "guard curve_cells");
        CHECK(hex9_owned_cells(u, 8, out, NULL, 16) < 0, "guard owned_cells");
        CHECK(hex9_neighbors(u, 6, out) < 0, "guard neighbors");
        CHECK(hex9_k_ring(u, 6, 1, out, 16) < 0, "guard k_ring");
        CHECK(hex9_k_disk(u, 6, 1, out, 16) < 0, "guard k_disk");
        CHECK(hex9_label(u, 6, buf, sizeof buf) < 0, "guard label");
        CHECK(hex9_label_key(u, 6, buf, sizeof buf) < 0, "guard label_key");
        CHECK(hex9_common_ancestor(u, 1, 6, buf, sizeof buf, NULL) < 0,
              "guard common_ancestor");
        CHECK(hex9_decode_boct(u, &cx, &cy, &oid) != 0, "guard decode_boct");
        CHECK(hex9_cell_ring(u, 6, 0, (double *)out, 8) < 0,
              "guard cell_ring");
        CHECK(hex9_cell_uv(u, 6, &ia, &ib, &oid, ia6, ib6, o6, &ext) != 0,
              "guard cell_uv");
    }

    if (fails) { printf("FAILED: %d checks\n", fails); return 1; }
    printf("OK: e4h laws — count law, truncation=binning, budget, grammar, "
           "matched pairs, marker guards\n");
    return 0;
}
