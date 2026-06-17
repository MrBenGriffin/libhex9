/* bin_prefix_guard.c — a standing guard against the "prefix-cut" lapse.
 *
 * Tempting but WRONG assumption: "to coarsen a bin/address to layer Lc, just
 * cut its label to the first Lc+1 chars (or trust the leading nibbles)."
 *
 * It is false by design. A split-hex body (6/7/8) names two half-hexes that
 * share one parent; when such a cell is the LEAF we must pick a parent, and by
 * design that is always the MODE-0 parent. A point on the mode-1 side descends
 * through mode-1 nibbles, but its canonical bin is the mode-0 parent — so the
 * bin's body differs from the address's body (recoverable only via the uuid
 * meta). Bins are layer-scoped KEYS, not nested addresses.
 *
 * The reliable coarsening is to re-bin from the FULL uuid via the identity path
 * — hex9_bin / hex9_k_disk(full, L, 0) — which is exactly what h9_grid
 * enumerates. (hex9_bin resolves split-hex leaves to their mode-0 home, so it
 * IS canonical; the lapse that remains is cutting a cosmetic LABEL string.)
 *
 * Two guards on a global corpus:
 *   1. the label prefix-cut must keep diverging (> 0) — else the corpus stopped
 *      exercising split-hex cells and the guard has gone BLIND;
 *   2. hex9_bin must equal the canonical bin everywhere (== 0 divergences) —
 *      pins hex9_bin to the identity/grid coarsening against regression.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

/* deterministic LCG (same generator as test/adaptive.c) */
static unsigned long long rng_state = 42;
static double rng_uniform(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11) / (double)(1ULL << 53);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    enum { N = 20000 };
    const int L_FINE = 12;
    const int L_COARSE[] = { 1, 2, 3, 4, 6, 8 };
    const int NC = (int)(sizeof L_COARSE / sizeof L_COARSE[0]);

    long checked = 0;
    long prefix_div = 0;   /* label prefix-cut  != canonical bin (the named lapse) */
    long nibble_div = 0;   /* hex9_bin nibble-trust != canonical bin (same family) */
    int  shown = 0;

    for (int i = 0; i < N; ++i) {
        const double lon = rng_uniform() * 360.0 - 180.0;
        const double lat = rng_uniform() * 170.0 - 85.0;
        uint8_t full[16];
        hex9_encode(lon, lat, full);

        /* canonical fine bin + its label (identity path = h9_grid) */
        uint8_t fine[16];
        if (hex9_k_disk(full, L_FINE, 0, fine, 1) != 1) continue;
        char fine_lbl[40];
        if (hex9_label(fine, L_FINE, fine_lbl, sizeof fine_lbl) < 0) continue;

        for (int c = 0; c < NC; ++c) {
            const int Lc = L_COARSE[c];
            if (Lc >= L_FINE) continue;

            /* the RIGHT answer: re-bin the full uuid at the coarse layer */
            uint8_t canon[16];
            if (hex9_k_disk(full, Lc, 0, canon, 1) != 1) continue;
            char canon_lbl[40];
            if (hex9_label(canon, Lc, canon_lbl, sizeof canon_lbl) < 0) continue;
            checked++;

            /* lapse #1 — prefix-cut the fine label to Lc+1 chars and parse it */
            char cut[40];
            memcpy(cut, fine_lbl, (size_t)(Lc + 1));
            cut[Lc + 1] = '\0';
            uint8_t cut_uuid[16];
            const int parsed = hex9_parse_label(cut, cut_uuid);
            const int prefix_wrong = (parsed != Lc) || memcmp(cut_uuid, canon, 16) != 0;
            if (prefix_wrong) prefix_div++;

            /* lapse #2 — trust leading nibbles via hex9_bin (h9_bin_uuid) */
            uint8_t nb[16];
            int nibble_wrong = 0;
            if (hex9_bin(full, Lc, nb) == 0)
                nibble_wrong = memcmp(nb, canon, 16) != 0;
            if (nibble_wrong) nibble_div++;

            if ((prefix_wrong || nibble_wrong) && shown < 8) {
                char nb_lbl[40] = "?"; hex9_label(nb, Lc, nb_lbl, sizeof nb_lbl);
                printf("  split-hex divergence @L%d  pt=(%.4f,%.4f)\n", Lc, lon, lat);
                printf("    canonical (k_disk from full) = %s\n", canon_lbl);
                printf("    prefix-cut  of L%d label      = %s  %s\n",
                       L_FINE, cut, prefix_wrong ? "<-- WRONG" : "(ok)");
                printf("    hex9_bin    (nibble trust)    = %s  %s\n",
                       nb_lbl, nibble_wrong ? "<-- WRONG" : "(ok)");
                shown++;
            }
        }
    }

    printf("\nbin-prefix guard: %ld (point,layer) pairs checked\n", checked);
    printf("  label prefix-cut wrong       : %ld (%.3f%%)\n", prefix_div,
           checked ? 100.0 * (double)prefix_div / (double)checked : 0.0);
    printf("  hex9_bin != canonical bin    : %ld (%.3f%%)\n", nibble_div,
           checked ? 100.0 * (double)nibble_div / (double)checked : 0.0);

    /* Guard 1 — the LAPSE must stay catchable: prefix-cutting a (cosmetic)
     * label to coarsen is wrong, because canonical labels are NOT cross-layer
     * prefix-nested at split-hex (6/7/8). A zero count means the corpus stopped
     * exercising split-hex cells and the guard has gone blind. */
    CHECK(prefix_div > 0,
          "prefix-cut guard went BLIND (0 divergences) — corpus no longer "
          "exercises split-hex cells; the 'labels prefix-cut to coarsen' lapse "
          "can no longer be caught here");

    /* Guard 2 — pin the fix: hex9_bin MUST equal the canonical (grid-matching)
     * bin everywhere. If it ever regresses to the address-nibble bin, split-hex
     * cells diverge and this jumps off zero. */
    CHECK(nibble_div == 0,
          "hex9_bin diverged from the canonical bin in %ld cases — it must be "
          "the identity/grid coarsening (split-hex leaf -> mode-0 home)", nibble_div);

    printf("\nDIRECTIVE: to COARSEN, re-bin from the FULL uuid — hex9_bin / "
           "hex9_k_disk\n  (C) or h9_bin / h9_kdisk (SQL); all are the canonical "
           "identity coarsening\n  and == h9_grid. Do NOT prefix-cut a label string "
           "(labels are cosmetic and\n  not cross-layer prefix-nested at split-hex "
           "6/7/8: the canonical leaf is the\n  mode-0 parent, whose body differs from "
           "the address body). Bins are\n  layer-scoped keys, not nested addresses. "
           "See docs/addressing-doctrine.md.\n");

    printf("\n%s (%d failures)\n", failures ? "BIN-PREFIX GUARD FAILED" : "bin-prefix guard OK", failures);
    return failures ? 1 : 0;
}
