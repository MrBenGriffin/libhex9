/* pdal_bin_verify.c — validate filters.hex9bin output invariants.
 * argv[1] = agg csv (H9Layer,H9Count,H9Value,Z,ZMax,ZMin,Classification header),
 * argv[2] = expected input point count N.
 * Checks: sum(H9Count)==N (sample captured exactly), sum(H9Value)==N,
 *         ZMin<=Z(mean)<=ZMax per cell, Classification mode integral. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
    if (argc < 3) { printf("usage: pdal_bin_verify agg.csv N\n"); return 2; }
    long N = atol(argv[2]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { printf("cannot open %s\n", argv[1]); return 1; }
    char line[512];
    if (!fgets(line, sizeof line, f)) { printf("empty\n"); return 1; }  /* header */
    long sum_count = 0, cells = 0, bad = 0;
    double sum_value = 0;
    double layer, count, value, z, zmax, zmin, cls;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &layer, &count, &value, &z, &zmax, &zmin, &cls) != 7) continue;
        cells++;
        sum_count += (long)llround(count);
        sum_value += value;
        if (z < zmin - 1e-6 || z > zmax + 1e-6) bad++;          /* mean within extents */
        if (cls != floor(cls)) bad++;                            /* mode is integral */
        if (layer < 0 || layer > 29) bad++;
    }
    fclose(f);
    int ok = (sum_count == N) && (fabs(sum_value - N) < 1e-6) && bad == 0 && cells > 0;
    printf("cells=%ld  sum(count)=%ld (want %ld)  sum(value)=%.2f  bad=%ld\n",
           cells, sum_count, N, sum_value, bad);
    printf(ok ? "HEX9BIN OK\n" : "HEX9BIN FAILED\n");
    return ok ? 0 : 1;
}
