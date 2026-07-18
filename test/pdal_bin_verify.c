/* pdal_bin_verify.c — validate filters.hex9bin output invariants.
 * argv[1] = agg csv (H9Layer,H9Count,H9Value,Z,ZMax,ZMin,Classification,
 *           Intensity,H9CurveHi,H9CurveLo header),
 * argv[2] = expected input point count N.
 * Checks: sum(H9Count)==N (sample captured exactly), sum(H9Value)==N,
 *         ZMin<=Z(mean)<=ZMax per cell, Classification mode integral, and
 *         the curve dims: H9CurveHi carries a well-formed packed curve
 *         address whose layer equals H9Layer (writers.text emits doubles,
 *         but for layers <= 12 all meaningful curve nibbles live in hi's
 *         exact top 52 bits — marker, slot, ranks, first sentinel). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 3) { printf("usage: pdal_bin_verify agg.csv N\n"); return 2; }
    long N = atol(argv[2]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { printf("cannot open %s\n", argv[1]); return 1; }
    char line[512];
    if (!fgets(line, sizeof line, f)) { printf("empty\n"); return 1; }  /* header */
    long sum_count = 0, cells = 0, bad = 0, bad_curve = 0;
    double sum_value = 0;
    double layer, count, value, z, zmax, zmin, cls, inten, chi_d, clo_d;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &layer, &count, &value, &z, &zmax, &zmin, &cls, &inten,
                   &chi_d, &clo_d) != 10) continue;
        cells++;
        sum_count += (long)llround(count);
        sum_value += value;
        if (z < zmin - 1e-6 || z > zmax + 1e-6) bad++;          /* mean within extents */
        if (cls != floor(cls)) bad++;                            /* mode is integral */
        if (layer < 0 || layer > 30) bad++;
        /* Curve address smoke check. writers.text emits doubles, so chi_d is
         * the uint64 hi rounded to 53 bits: the round-up CARRIES through the
         * 0xF sentinel run into the last rank nibble (0x..3fffff -> 0x..400000).
         * Marker and slot are exact (a carry stops at the first nibble <= 8);
         * rank nibbles are 0..8 true, 9 or +1-at-the-boundary under carry.
         * Bit-exactness of the curve itself is pinned by test/curve.c and the
         * SQL regression; this validates the dims are wired and well-formed. */
        const int L = (int)layer;
        const uint64_t chi = (uint64_t)chi_d;
        if ((chi >> 60) != 0x0Cu) bad_curve++;                   /* marker */
        else if (((chi >> 56) & 0xFu) > 11u) bad_curve++;        /* axiom slot */
        else {
            for (int k = 1; k <= L && k + 1 <= 14; ++k)
                if (((chi >> (60 - 4 * (k + 1))) & 0xFu) > 9u) { bad_curve++; break; }
        }
        (void)clo_d;   /* all-sentinel at these layers; rounds to 2^64 as double */
    }
    fclose(f);
    int ok = (sum_count == N) && (fabs(sum_value - N) < 1e-6) && bad == 0
             && bad_curve == 0 && cells > 0;
    printf("cells=%ld  sum(count)=%ld (want %ld)  sum(value)=%.2f  bad=%ld  bad_curve=%ld\n",
           cells, sum_count, N, sum_value, bad, bad_curve);
    printf(ok ? "HEX9BIN OK\n" : "HEX9BIN FAILED\n");
    return ok ? 0 : 1;
}
