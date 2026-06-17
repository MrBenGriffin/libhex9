/* cellgeom.c — coherence of single-cell geometry with the grid enumerator.
 *
 * hex9_cell_ring must reproduce hex9_grid_cell_ring exactly for every grid
 * cell (the regression that motivated this: the old cell_unpack decode path
 * rendered the wrong hexagon for canonical bins of cells where centroid- and
 * origin-descent disagree — e.g. Westminster at L8), and hex9_decode of a
 * bin must equal the grid centroid.
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

static void check_bbox(double x0, double y0, double x1, double y1, int layer,
                       const char *tag) {
    char err[256] = {0};
    hex9_grid *g = hex9_grid_create(x0, y0, x1, y1, layer, 0, 0, err, sizeof err);
    CHECK(g != NULL, "%s: grid create: %s", tag, err);
    if (!g) return;
    const int n = hex9_grid_count(g);
    CHECK(n > 0, "%s: empty grid", tag);

    for (int d = 0; d <= 1; ++d) {
        const int np = hex9_ring_npoints(d);
        double *rg = (double *)malloc((size_t)np * 2 * sizeof(double));
        double *rc = (double *)malloc((size_t)np * 2 * sizeof(double));
        for (int i = 0; i < n; ++i) {
            uint8_t uu[16];
            hex9_grid_cell_uuid(g, i, uu);
            CHECK(hex9_grid_cell_ring(g, i, d, rg, np) == np, "%s: grid ring %d", tag, i);
            CHECK(hex9_cell_ring(uu, layer, d, rc, np) == np, "%s: cell ring %d", tag, i);
            CHECK(memcmp(rg, rc, (size_t)np * 2 * sizeof(double)) == 0,
                  "%s: cell %d d%d ring differs from grid", tag, i, d);
            if (d == 0) {
                double glon, glat, dlon, dlat;
                hex9_grid_cell_centroid(g, i, &glon, &glat);
                hex9_decode(uu, &dlon, &dlat);
                CHECK(fabs(glon - dlon) < 1e-12 && fabs(glat - dlat) < 1e-12,
                      "%s: cell %d decode (%.9f,%.9f) != centroid (%.9f,%.9f)",
                      tag, i, dlon, dlat, glon, glat);
            }
            if (failures > 8) { free(rg); free(rc); hex9_grid_destroy(g); return; }
        }
        free(rg); free(rc);
    }
    printf("%s: %d cells coherent (d0+d1 rings, centroids)\n", tag, n);
    hex9_grid_destroy(g);
}

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }

    /* Westminster — contains the convention-divergent cell that exposed the bug */
    check_bbox(-0.135, 51.500, -0.120, 51.512, 8, "westminster-L8");
    /* the original 1255-cell London bbox */
    check_bbox(-0.3, 51.3, 0.1, 51.7, 8, "london-L8");
    /* seam + pole coverage at coarser/finer layers */
    check_bbox(-5.0, -3.0, 5.0, 3.0, 5, "equator-meridian-L5");
    check_bbox(-180.0, 85.0, 180.0, 90.0, 4, "polar-L4");
    check_bbox(-0.13, 51.50, -0.125, 51.505, 12, "fine-L12");

    if (failures) { printf("CELLGEOM FAILED (%d)\n", failures); return 1; }
    printf("CELLGEOM OK\n");
    return 0;
}
