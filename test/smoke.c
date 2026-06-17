/* smoke.c — minimal standalone check: warp init + encode/decode round-trip. */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(void) {
    char err[256] = {0};
    printf("%s\n", hex9_version());

    if (hex9_warp_init(err, sizeof err)) {
        printf("warp init FAILED: %s\n", err);
        return 1;
    }
    printf("warp init OK\n");

    /* Edinburgh-ish */
    const double lon = -3.19, lat = 55.95;
    uint8_t uuid[16];
    hex9_encode(lon, lat, uuid);

    printf("uuid: ");
    for (int i = 0; i < 16; ++i) printf("%02x", uuid[i]);
    printf("\n");

    double dlon, dlat;
    hex9_decode(uuid, &dlon, &dlat);
    const double err_lon = fabs(dlon - lon), err_lat = fabs(dlat - lat);
    printf("in (%.6f, %.6f) -> out (%.6f, %.6f)  |dlon|=%.2e |dlat|=%.2e\n",
           lon, lat, dlon, dlat, err_lon, err_lat);

    /* Layer-29 cell is ~95 nm; round-trip should be far under 1e-4 deg. */
    if (err_lon > 1e-4 || err_lat > 1e-4) {
        printf("ROUND-TRIP TOO LARGE\n");
        return 1;
    }
    printf("round-trip OK\n");

    /* labels: bin to L8, label it; expect the CANONICAL (grid-matching) bin
     * 432177468 — the L8 cell that geometrically contains Edinburgh (grid- and
     * ST_Contains-verified). The address-family value 432177478 (hhg9 oracle /
     * pre-canonical h9_bin) is the mode-1 sibling and does NOT contain the
     * point; h9_bin now resolves split-hex (6/7/8) leaves to their mode-0
     * home. (Digit 7: address 7 -> canonical 6.) */
    uint8_t binned[16];
    hex9_bin(uuid, 8, binned);
    char lbl[40] = {0}, lblk[40] = {0};
    int ln  = hex9_label(binned, 8, lbl, sizeof lbl);
    int lnk = hex9_label_key(uuid, 8, lblk, sizeof lblk);
    printf("label(L8)=%s  label_key(L8)=%s  (len %d/%d)\n", lbl, lblk, ln, lnk);
    if (ln < 0 || lnk < 0 || strncmp(lbl, "432177468", 9) != 0) {
        printf("LABEL MISMATCH\n");
        return 1;
    }
    printf("label OK\n");

    /* Coarsening: ALWAYS re-bin from the FULL uuid (the load-bearer). Bins are
     * canonical layer-scoped keys, NOT nested addresses — coarsening a bin by a
     * further bin (the prefix-cut lapse) is not guaranteed and at split-hex
     * (6/7/8, which Edinburgh is) diverges from the full re-bin; see
     * test/bin_prefix_guard. Re-binning at the bin's own layer IS the identity.
     * (Regression: the backward pass once walked the 0xF sentinel nibbles of
     * bin input through H9_HEX_REG — out-of-bounds, build-dependent.) */
    uint8_t bin88[16];
    hex9_bin(binned, 8, bin88);
    if (memcmp(binned, bin88, 16) != 0) {
        printf("REBIN (own-layer identity) MISMATCH\n");
        return 1;
    }
    printf("rebin OK\n");

    /* cell_ring: densify 0 → 7 pts, densify 1 → 19; ring must be closed. */
    double ring[2 * 19];
    int np0 = hex9_cell_ring(binned, 8, 0, ring, 19);
    int closed = (np0 > 0) && ring[0] == ring[2 * (np0 - 1)] && ring[1] == ring[2 * (np0 - 1) + 1];
    printf("cell_ring(L8,d0): n=%d v0=(%.7f, %.7f) closed=%d\n", np0, ring[0], ring[1], closed);
    printf("binned uuid: ");
    for (int i = 0; i < 16; ++i) printf("%02x", binned[i]);
    printf("\n");
    int np1 = hex9_cell_ring(binned, 8, 1, ring, 19);
    printf("cell_ring(L8,d1): n=%d\n", np1);
    if (np0 != 7 || np1 != 19 || !closed) { printf("CELL_RING WRONG\n"); return 1; }
    printf("cell_ring OK\n");
    return 0;
}
