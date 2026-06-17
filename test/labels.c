/* labels.c — hex9_parse_label / hex9_label_centroid / hex9_common_ancestor.
 *
 * Round-trips labels (bare and keyed) against canonical bins at interior,
 * seam and polar cells; checks the keyed form is flavour-blind (labels of
 * beam-encoded bins parse to the same cell); checks ancestor/prefix
 * structure over k-disks and against the grid enumerator's centroids.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

static void roundtrip(double lon, double lat, int layer, const char *tag) {
    uint8_t full[16], bin[16], canon[16], parsed[16];
    char lbl[40], lblk[40];
    hex9_encode(lon, lat, full);
    hex9_bin(full, layer, bin);
    /* canonical form via k_disk(0) on the FULL uuid (bins are keys, not
     * addresses — adjacency input is full-UUID only) */
    CHECK(hex9_k_disk(full, layer, 0, canon, 1) == 1, "%s L%d: canonicalise", tag, layer);

    /* bare label of the canonical bin parses back to it */
    CHECK(hex9_label(canon, layer, lbl, sizeof lbl) == layer + 1, "%s L%d: label", tag, layer);
    CHECK(hex9_parse_label(lbl, parsed) == layer, "%s L%d: parse('%s')", tag, layer, lbl);
    CHECK(memcmp(parsed, canon, 16) == 0, "%s L%d: parse('%s') != canon", tag, layer, lbl);

    /* keyed label: a keyed label IS the bin in text form, and parse is its
     * exact SYNTACTIC inverse — parse(label_key(bin)) == bin, an identity
     * (no resolution; bins are keys, not addresses). */
    CHECK(hex9_label_key(bin, layer, lblk, sizeof lblk) > 0, "%s L%d: label_key", tag, layer);
    CHECK(hex9_parse_label(lblk, parsed) == layer, "%s L%d: parse('%s')", tag, layer, lblk);
    CHECK(memcmp(parsed, bin, 16) == 0, "%s L%d: parse('%s') != bin", tag, layer, lblk);

    /* centroid is finite and decodes near the cell (within a few cell sizes) */
    double clon, clat, dlon, dlat;
    CHECK(hex9_label_centroid(lbl, &clon, &clat) == 0, "%s L%d: centroid", tag, layer);
    hex9_decode(canon, &dlon, &dlat);
    double dd = fabs(clon - dlon);
    if (dd > 180.0) dd = 360.0 - dd;
    const double cell_deg = 60.0 / pow(3.0, layer) + 1e-9;
    CHECK(dd * cos(clat * M_PI / 180.0) < 4 * cell_deg && fabs(clat - dlat) < 4 * cell_deg,
          "%s L%d: centroid (%f,%f) far from decode (%f,%f)", tag, layer, clon, clat, dlon, dlat);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    int layers[] = { 1, 4, 8, 15, 29 };
    for (size_t i = 0; i < sizeof layers / sizeof *layers; ++i) {
        roundtrip(-3.19, 55.95, layers[i], "edinburgh");
        roundtrip(12.0, 0.001, layers[i], "equator");
        roundtrip(0.001, 45.0, layers[i], "meridian");
        roundtrip(13.7, 89.999, layers[i], "pole");
        roundtrip(0.001, -0.001, layers[i], "vertex");
    }

    /* keyed-parse identity sweep over a world lattice — the load-bearing check:
     * parse(label_key(u,L)) == hex9_bin(u,L). Labels now name the CANONICAL
     * (grid-matching) bin, so label/label_key/parse are all the identity
     * family: the keyed tail is the canonical bin tail (c2<<1 | oct_mode, no
     * c_mo display bit), and parse reconstructs the bin exactly. The sweep must
     * exercise split-hex (6/7/8) cells (where the canonical home differs from
     * the address descent) or it has gone blind — asserted via keys_seen. */
    {
        unsigned keys_seen = 0;
        int split_hex_seen = 0;
        for (int yi = 0; yi < 20; yi++) {
            for (int xi = 0; xi < 20; xi++) {
                const double lon = -178.0 + 356.0 * xi / 19.0;
                const double lat = -85.0 + 170.0 * yi / 19.0;
                uint8_t full[16], bin[16], parsed[16];
                char lblk[40], lbl[40];
                if (hex9_encode(lon, lat, full) != 0) continue;
                for (int L = 1; L <= 11; L += 2) {
                    hex9_bin(full, L, bin);
                    const int kl = hex9_label_key(full, L, lblk, sizeof lblk);
                    CHECK(kl > 2, "sweep (%g,%g) L%d: label_key", lon, lat, L);
                    if (kl <= 2) continue;
                    const char kc = lblk[kl - 1];
                    keys_seen |= 1u << (kc <= '9' ? kc - '0' : kc - 'a' + 10);
                    CHECK(hex9_parse_label(lblk, parsed) == L,
                          "sweep (%g,%g) L%d: parse('%s')", lon, lat, L, lblk);
                    CHECK(memcmp(parsed, bin, 16) == 0,
                          "sweep (%g,%g) L%d: parse('%s') != bin", lon, lat, L, lblk);
                    if (hex9_label(full, L, lbl, sizeof lbl) == L + 1) {
                        const char leaf = lbl[L];
                        if (leaf == '6' || leaf == '7' || leaf == '8') split_hex_seen++;
                    }
                }
            }
        }
        CHECK(keys_seen != 0, "keyed sweep saw no tails (keys 0x%04x)", keys_seen);
        CHECK(split_hex_seen > 0,
              "keyed sweep exercised no split-hex (6/7/8) leaves — canonical "
              "resolution untested, guard blind");
    }

    /* grid cross-check: every grid cell's label parses back to its uuid and
     * label_centroid reproduces the grid centroid exactly */
    hex9_grid *g = hex9_grid_create(-0.3, 51.3, 0.1, 51.7, 8, 0, 0, err, sizeof err);
    CHECK(g != NULL, "grid create: %s", err);
    if (g) {
        const int n = hex9_grid_count(g);
        /* 1247 under the F6 field (2026-06-12 re-baseline; was 1255 pre-F6 —
         * the global field shift moved boundary cells across the bbox) */
        CHECK(n == 1247, "grid count %d != 1247", n);
        for (int i = 0; i < n; i++) {
            uint8_t uu[16], parsed[16];
            char lbl[40];
            hex9_grid_cell_uuid(g, i, uu);
            hex9_label(uu, 8, lbl, sizeof lbl);
            CHECK(hex9_parse_label(lbl, parsed) == 8 && memcmp(parsed, uu, 16) == 0,
                  "grid cell %d: label '%s' round-trip", i, lbl);
            double glon, glat, clon, clat;
            hex9_grid_cell_centroid(g, i, &glon, &glat);
            CHECK(hex9_label_centroid(lbl, &clon, &clat) == 0 &&
                  fabs(clon - glon) < 1e-9 && fabs(clat - glat) < 1e-9,
                  "grid cell %d: centroid mismatch (%.12f,%.12f) vs (%.12f,%.12f)",
                  i, clon, clat, glon, glat);
            if (failures > 8) break;
        }
        hex9_grid_destroy(g);
    }

    /* common ancestor: a k-disk shares a proper prefix; each member's label
     * starts with the ancestor label */
    {
        uint8_t full[16], canon[16];
        hex9_encode(-3.19, 55.95, full);
        hex9_k_disk(full, 12, 0, canon, 1);
        uint8_t disk[19 * 16];
        const int64_t nd = hex9_k_disk(full, 12, 2, disk, 19);
        CHECK(nd == 19, "disk(2) = %lld", (long long)nd);
        char anc[40], anc_lbl_check[40];
        uint8_t anc_uuid[16];
        const int al = hex9_common_ancestor(disk, (size_t)nd, 12, anc, sizeof anc, anc_uuid);
        CHECK(al >= 0 && al < 12, "ancestor layer %d", al);
        CHECK((int)strlen(anc) == al + 1, "ancestor label len");
        CHECK(hex9_label(anc_uuid, al, anc_lbl_check, sizeof anc_lbl_check) == al + 1 &&
              strcmp(anc, anc_lbl_check) == 0, "ancestor uuid/label agree");
        for (int i = 0; i < nd; i++) {
            char lbl[40];
            hex9_label(disk + (size_t)i * 16, 12, lbl, sizeof lbl);
            CHECK(strncmp(lbl, anc, (size_t)al + 1) == 0,
                  "disk cell %d label '%s' lacks ancestor prefix '%s'", i, lbl, anc);
        }
        /* ancestor of a single cell is itself */
        CHECK(hex9_common_ancestor(canon, 1, 12, anc, sizeof anc, NULL) == 12,
              "self ancestor != layer");
        /* far-apart cells have no common ancestor */
        uint8_t two[32];
        memcpy(two, canon, 16);
        hex9_encode(151.2, -33.87, full);          /* sydney */
        hex9_bin(full, 12, two + 16);
        CHECK(hex9_common_ancestor(two, 2, 12, anc, sizeof anc, NULL) == -1,
              "edinburgh+sydney should have no ancestor");
    }

    /* parse rejects garbage */
    {
        uint8_t u[16];
        CHECK(hex9_parse_label("", u) == -1, "empty accepted");
        CHECK(hex9_parse_label("c12", u) == -1, "bad L0 digit accepted");
        CHECK(hex9_parse_label("49", u) == -1, "digit 9 beyond L0 accepted");
        CHECK(hex9_parse_label("432.q", u) == -1, "bad key accepted");
        CHECK(hex9_parse_label("432.", u) == -1, "empty key accepted");
    }

    if (failures) { printf("LABELS FAILED (%d)\n", failures); return 1; }
    printf("LABELS OK\n");
    return 0;
}
