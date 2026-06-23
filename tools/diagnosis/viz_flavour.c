/* viz_flavour.c — dump data to visualise the mode-1-half flavour resolution
 * defect (docs/addressing-doctrine.md, the wip/identity-reversible-walk
 * blocker).
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2026, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * For an N×N grid of points covering a cell and its 2-ring at layer L,
 * encode each point and report:
 *   - its WALK BIN (h9_bin — byte-identical to the Python reference), and
 *   - its IDENTITY cell (the ring rendered from the full uuid, hashed),
 * plus the outlines of the disk(2) cells (rendered from kring's canonical
 * uuids) for overlay.
 *
 * Usage: viz_flavour lon lat L N [encoder]  > dump.csv
 *        encoder: 0 = nn/beam (default), 1 = containment
 * Lines:  P,lon,lat,walkbin_hex,ringhash      (grid samples)
 *         R,cellidx,vtx,lon,lat               (outline vertices)
 */
#include "hex9_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RING_N 7

/* FNV-1a over the ring doubles — stable id for "which cell" */
static unsigned long ring_hash(const double *ring) {
    const unsigned char *b = (const unsigned char *)ring;
    unsigned long h = 1469598103934665603UL;
    for (size_t i = 0; i < 2 * RING_N * sizeof(double); i++) {
        h ^= b[i];
        h *= 1099511628211UL;
    }
    return h;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: viz_flavour lon lat L N\n"); return 2; }
    const double lon0 = atof(argv[1]), lat0 = atof(argv[2]);
    const int L = atoi(argv[3]), N = atoi(argv[4]);

    char err[256];
    if (hex9_warp_init(err, sizeof err)) { fprintf(stderr, "warp: %s\n", err); return 1; }
    if (argc > 5) hex9_set_encoder(atoi(argv[5]));

    uint8_t centre[16];
    if (hex9_encode(lon0, lat0, centre) != 0) { fprintf(stderr, "encode failed\n"); return 1; }

    /* disk(2) outlines from canonical kring emission */
    uint8_t disk2[19 * 16];
    const int64_t nd = hex9_k_disk(centre, L, 2, disk2, 19);
    if (nd <= 0) { fprintf(stderr, "k_disk failed\n"); return 1; }
    double bb_lon_min = 1e30, bb_lon_max = -1e30, bb_lat_min = 1e30, bb_lat_max = -1e30;
    for (int j = 0; j < (int)nd; j++) {
        double ring[2 * RING_N];
        if (hex9_cell_ring(disk2 + (size_t)j * 16, L, 0, ring, RING_N) != RING_N) continue;
        for (int k = 0; k < RING_N; k++) {
            printf("R,%d,%d,%.12f,%.12f\n", j, k, ring[2*k], ring[2*k + 1]);
            if (k < 6) {
                if (ring[2*k]   < bb_lon_min) bb_lon_min = ring[2*k];
                if (ring[2*k]   > bb_lon_max) bb_lon_max = ring[2*k];
                if (ring[2*k+1] < bb_lat_min) bb_lat_min = ring[2*k + 1];
                if (ring[2*k+1] > bb_lat_max) bb_lat_max = ring[2*k + 1];
            }
        }
    }

    /* sample grid over the disk(2) bbox (slightly shrunk to stay relevant) */
    const double pad = -0.02;
    const double dl = (bb_lon_max - bb_lon_min) * pad, dphi = (bb_lat_max - bb_lat_min) * pad;
    bb_lon_min -= dl; bb_lon_max += dl; bb_lat_min -= dphi; bb_lat_max += dphi;

    for (int iy = 0; iy < N; iy++) {
        for (int ix = 0; ix < N; ix++) {
            const double lon = bb_lon_min + (bb_lon_max - bb_lon_min) * (ix + 0.5) / N;
            const double lat = bb_lat_min + (bb_lat_max - bb_lat_min) * (iy + 0.5) / N;
            uint8_t full[16], bin[16];
            if (hex9_encode(lon, lat, full) != 0) continue;
            if (hex9_bin(full, L, bin) != 0) continue;
            double ring[2 * RING_N];
            if (hex9_cell_ring(full, L, 0, ring, RING_N) != RING_N) continue;
            char binhex[33];
            for (int b = 0; b < 16; b++) sprintf(binhex + 2*b, "%02x", bin[b]);
            printf("P,%.12f,%.12f,%s,%lu\n", lon, lat, binhex, ring_hash(ring));
        }
    }
    return 0;
}
