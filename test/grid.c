/* grid.c — enumerate a bbox and print count + cell-0 (uuid/centroid/ring) for
 * cross-checking against the PostGIS h9_grid over the same rectangular bbox. */
#include "hex9_c.h"
#include <stdio.h>

int main(void) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("init FAILED: %s\n", err); return 1; }

    const double lo0 = -0.3, la0 = 51.3, lo1 = 0.1, la1 = 51.7;  /* London-ish */
    const int layer = 8;

    hex9_grid *g = hex9_grid_create(lo0, la0, lo1, la1, layer, 0, 5000000, err, sizeof err);
    if (!g) { printf("grid_create FAILED: %s\n", err); return 1; }

    const int n = hex9_grid_count(g);
    printf("grid L%d bbox(%.2f,%.2f,%.2f,%.2f) count=%d\n", layer, lo0, la0, lo1, la1, n);

    if (n > 0) {
        uint8_t u[16];
        hex9_grid_cell_uuid(g, 0, u);
        double clon, clat;
        hex9_grid_cell_centroid(g, 0, &clon, &clat);
        double ring[2 * 7];
        int nr = hex9_grid_cell_ring(g, 0, 0, ring, 7);
        printf("cell0 uuid=");
        for (int i = 0; i < 16; ++i) printf("%02x", u[i]);
        printf(" centroid=(%.7f, %.7f) ring_n=%d v0=(%.7f, %.7f)\n",
               clon, clat, nr, ring[0], ring[1]);
    }
    hex9_grid_destroy(g);
    return n > 0 ? 0 : 1;
}
