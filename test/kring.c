/* kring.c — hex9_neighbors / hex9_k_ring / hex9_k_disk checks.
 *
 * The seam algebra itself is exhaustively validated against the geometric
 * mesh by tools/kring_probe.cpp (L1..L4) and against the encoder oracle by
 * test/gc_kring.c. This test exercises the C ABI at interior, seam, polar
 * and deep layers, checking the invariants that hold everywhere: neighbour
 * counts (6, or 5 on octahedron-vertex half-hexes), mutual adjacency,
 * distinctness, ring/disk arithmetic — and the INPUT CONTRACT (Ben's
 * ruling, 2026-06-11): adjacency input is a FULL UUID only; bins are
 * layer-scoped keys, not addresses, and are rejected. Output bins cannot
 * be traversed further — onward steps go through a point + hex9_encode,
 * which is how mutuality is checked here.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

static int contains(const uint8_t *uuids, int n, const uint8_t want[16]) {
    for (int i = 0; i < n; ++i)
        if (memcmp(uuids + (size_t)i * 16, want, 16) == 0) return 1;
    return 0;
}

/* Encode a point inside the cell keyed by `bin` and return its full UUID.
 * (Output bins are keys: the only sanctioned way to traverse onward is via
 * a geographic point.) Interior point = ring vertex-mean nudged 20% toward
 * vertex 0, which stays off the internal half-hex seam. */
static int full_inside(const uint8_t bin[16], int layer, uint8_t out_full[16]) {
    double ring[14];
    if (hex9_cell_ring(bin, layer, 0, ring, 7) != 7) return 0;
    double clon = 0.0, clat = 0.0;
    for (int k = 0; k < 6; ++k) { clon += ring[2*k]; clat += ring[2*k + 1]; }
    clon /= 6.0; clat /= 6.0;
    const double plon = clon + (ring[0] - clon) * 0.2;
    const double plat = clat + (ring[1] - clat) * 0.2;
    return hex9_encode(plon, plat, out_full) == 0;
}

/* Neighbour invariants for the cell containing (lon, lat) at `layer`.
 * Returns the neighbour count. */
static int check_cell(double lon, double lat, int layer, const char *tag) {
    uint8_t full[16], bin[16], canon[16], nbs[6 * 16];
    hex9_encode(lon, lat, full);
    hex9_bin(full, layer, bin);

    /* k_disk(k=0) on the full UUID is the cell itself in canonical form —
     * what all k-ring output uses. */
    CHECK(hex9_k_disk(full, layer, 0, canon, 1) == 1, "%s L%d: canonicalise failed", tag, layer);

    /* input contract: bins are keys, not addresses — rejected */
    uint8_t rej[6 * 16];
    CHECK(hex9_neighbors(bin, layer, rej) == -1, "%s L%d: bin input not rejected", tag, layer);
    CHECK(hex9_k_ring(bin, layer, 1, rej, 6) == -1, "%s L%d: bin k_ring not rejected", tag, layer);
    CHECK(hex9_k_disk(bin, layer, 1, rej, 7) == -1, "%s L%d: bin k_disk not rejected", tag, layer);

    const int n = hex9_neighbors(full, layer, nbs);
    CHECK(n == 5 || n == 6, "%s L%d: count %d not in {5,6}", tag, layer, n);
    if (n < 0) return n;

    /* distinct, none equal to self */
    for (int i = 0; i < n; ++i) {
        CHECK(memcmp(nbs + (size_t)i * 16, canon, 16) != 0, "%s L%d: nb %d == self", tag, layer, i);
        for (int j = i + 1; j < n; ++j)
            CHECK(memcmp(nbs + (size_t)i * 16, nbs + (size_t)j * 16, 16) != 0,
                  "%s L%d: nb %d == nb %d", tag, layer, i, j);
    }

    /* mutuality: encode a point inside each neighbour and ask for its
     * neighbours; assert only when the probe actually landed in the
     * intended neighbour (pole/half-hex geometry can deflect the interior
     * construction — those probes are skipped, not failed). */
    for (int i = 0; i < n; ++i) {
        uint8_t nf[16], ncanon[16], back[6 * 16];
        if (!full_inside(nbs + (size_t)i * 16, layer, nf)) continue;
        if (hex9_k_disk(nf, layer, 0, ncanon, 1) != 1) continue;
        if (memcmp(ncanon, nbs + (size_t)i * 16, 16) != 0) continue;  /* probe strayed */
        const int m = hex9_neighbors(nf, layer, back);
        CHECK(m == 5 || m == 6, "%s L%d: back-count %d", tag, layer, m);
        CHECK(m > 0 && contains(back, m, canon), "%s L%d: nb %d not mutual", tag, layer, i);
    }

    /* k=1 ring == neighbours; k=1 disk == neighbours + self */
    uint8_t ring[8 * 16], disk[8 * 16];
    const int64_t nr = hex9_k_ring(full, layer, 1, ring, 8);
    const int64_t nd = hex9_k_disk(full, layer, 1, disk, 8);
    CHECK(nr == n, "%s L%d: ring(1)=%lld != %d", tag, layer, (long long)nr, n);
    CHECK(nd == n + 1, "%s L%d: disk(1)=%lld != %d", tag, layer, (long long)nd, n + 1);
    CHECK(nd > 0 && contains(disk, (int)nd, canon), "%s L%d: disk(1) misses centre", tag, layer);
    for (int i = 0; i < n; ++i)
        CHECK(contains(ring, (int)nr, nbs + (size_t)i * 16), "%s L%d: ring(1) misses nb %d",
              tag, layer, i);
    return n;
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    CHECK(hex9_disk_ncells(0) == 1 && hex9_disk_ncells(1) == 7 &&
          hex9_disk_ncells(2) == 19 && hex9_disk_ncells(3) == 37, "disk_ncells formula");
    CHECK(hex9_disk_ncells(-1) == -1, "disk_ncells(-1)");

    /* interior cells across shallow and deep layers */
    int layers[] = { 2, 5, 8, 15, 22, 29 };
    for (size_t i = 0; i < sizeof layers / sizeof *layers; ++i) {
        const int n = check_cell(-3.19, 55.95, layers[i], "edinburgh");
        CHECK(n == 6, "edinburgh L%d: expected 6 neighbours, got %d", layers[i], n);
    }

    /* seam straddlers: equator and the lon-0 meridian */
    check_cell(12.0, 0.001, 8, "equator");
    check_cell(12.0, -0.001, 8, "equator-south");
    check_cell(0.001, 45.0, 8, "meridian");
    check_cell(-0.001, 45.0, 8, "meridian-west");
    check_cell(90.001, 0.001, 8, "junction-90-0");

    /* poles: octahedron vertices — half-hex country (counts may be 5) */
    check_cell(0.0, 89.999, 6, "north-pole");
    check_cell(13.7, 89.9999, 12, "north-pole-12");
    check_cell(-120.0, -89.999, 6, "south-pole");

    /* the octahedron vertex at (0, 0): vertex hexagon halves */
    int five_seen = 0;
    double probes[][2] = { {0.001, 0.001}, {-0.001, 0.001}, {0.001, -0.001}, {-0.001, -0.001} };
    for (int i = 0; i < 4; ++i) {
        const int n = check_cell(probes[i][0], probes[i][1], 4, "vertex-0-0");
        if (n == 5) five_seen++;
    }
    CHECK(five_seen > 0, "no half-hex (degree 5) found at the (0,0) octahedron vertex at L4");

    /* k-disk growth at an interior cell */
    uint8_t full[16];
    hex9_encode(-3.19, 55.95, full);
    for (int k = 1; k <= 5; ++k) {
        const int64_t cap = hex9_disk_ncells(k);
        uint8_t *buf = (uint8_t *)malloc((size_t)cap * 16);
        const int64_t nd = hex9_k_disk(full, 8, k, buf, cap);
        const int64_t nr = hex9_k_ring(full, 8, k, buf, cap);
        CHECK(nd == cap, "disk(k=%d)=%lld != %lld", k, (long long)nd, (long long)cap);
        CHECK(nr == 6 * k, "ring(k=%d)=%lld != %d", k, (long long)nr, 6 * k);
        free(buf);
    }

    /* error paths */
    uint8_t buf[7 * 16];
    CHECK(hex9_neighbors(full, 0, buf) == -1, "layer 0 not rejected");
    CHECK(hex9_neighbors(full, 30, buf) == -1, "layer 30 not rejected");
    CHECK(hex9_k_disk(full, 8, 1, buf, 6) == -1, "undersized buffer not rejected");
    CHECK(hex9_k_ring(full, 8, 0, buf, 7) == 1, "ring(0) should be the cell itself");

    if (failures) { printf("KRING FAILED (%d)\n", failures); return 1; }
    printf("KRING OK\n");
    return 0;
}
