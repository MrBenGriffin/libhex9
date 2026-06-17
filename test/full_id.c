/* full_id.c — the grid's full identity UUID (h9_grid_cell_id) is a valid,
 * unique, reversible address for every cell.
 *
 * For each cell: h9_bin(id, layer) must equal the cell's own bin (so the
 * identity re-bins exactly — the property a raw bin's own coarsening does NOT
 * guarantee, doctrine F3/F4), and every cell's identity must be distinct.
 * Exercised on a benign mid-latitude box, a box straddling the lon=0 octant
 * seam, and a near-polar box — the seam/pole cells are where the centroid's
 * half-hex-seam instability would otherwise bind (see full_id_from_cell).
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

static int cmp16(const void *a, const void *b) {
    return memcmp(a, b, 16);
}

static void check_window(double lonc, double latc, int layer, const char *tag) {
    /* Shrink the window 3×/layer so the cell count stays bounded (~1–2k). */
    const double half = 0.20 * pow(3.0, 8 - layer);
    char err[256] = {0};
    hex9_grid *g = hex9_grid_create(lonc - half, latc - half,
                                    lonc + half, latc + half,
                                    layer, 0, 0, err, sizeof err);
    CHECK(g != NULL, "%s L%d: grid create: %s", tag, layer, err);
    if (!g) return;

    const int n = hex9_grid_count(g);
    CHECK(n > 0, "%s L%d: empty grid", tag, layer);

    uint8_t *ids = (uint8_t *)malloc((size_t)n * 16);
    int bad_bin = 0;
    for (int i = 0; i < n; ++i) {
        uint8_t bin[16], id[16], binback[16];
        hex9_grid_cell_uuid(g, i, bin);
        hex9_grid_cell_id(g, i, id);
        hex9_bin(id, layer, binback);
        if (memcmp(binback, bin, 16) != 0) bad_bin++;
        memcpy(ids + (size_t)i * 16, id, 16);
    }
    CHECK(bad_bin == 0, "%s L%d: %d/%d ids do not bin back to their cell",
          tag, layer, bad_bin, n);

    /* Uniqueness: sort the ids and look for adjacent duplicates. */
    qsort(ids, (size_t)n, 16, cmp16);
    int dups = 0;
    for (int i = 1; i < n; ++i)
        if (memcmp(ids + (size_t)(i-1) * 16, ids + (size_t)i * 16, 16) == 0) dups++;
    CHECK(dups == 0, "%s L%d: %d duplicate identity uuids among %d cells",
          tag, layer, dups, n);

    free(ids);
    hex9_grid_destroy(g);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    const int layers[] = { 8, 11, 13, 16, 20 };
    const double centres[][2] = { { -0.10, 51.50 }, { 0.0, 0.0 }, { 0.0, 88.5 } };
    const char  *tags[]       = { "mid-lat", "equator-seam", "near-pole" };

    for (int r = 0; r < 3; ++r)
        for (int li = 0; li < (int)(sizeof(layers)/sizeof(layers[0])); ++li)
            check_window(centres[r][0], centres[r][1], layers[li], tags[r]);

    if (failures) { printf("full_id: %d FAILURES\n", failures); return 1; }
    printf("full_id: OK\n");
    return 0;
}
