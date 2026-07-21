#!/usr/bin/env python3
"""
gen_pole_seam_ref.py — parity reference for test/pole_seam.c.

A DELIBERATELY ADVERSARIAL companion to gen_via_sphere_ref.py. That file
samples the sphere broadly (equal-area random + a handful of landmarks); this
one samples only the places where the chain is hard, densely:

  * the polar approach — lat = ±(90 − ε) for ε geometrically spaced from 1e-12
    to 2 degrees, plus exactly ±90. The authalic series, the atan2/hypot
    latitude recovery and the octant-vertex snap all converge here at once.
  * the octant seams — the 8 face boundaries, sampled ON the seam and at
    ±1e-9 / ±1e-13 degrees either side, so the oid tie-break and the
    seam-fold frame are exercised from both directions.
  * the octant CORNERS — lat = ±atan(1/√2) crossed with the seam longitudes.
    This is the triple point where three faces and two seams meet, and it is
    the single worst-conditioned family of points on the grid.
  * the poles crossed with every seam longitude, where longitude is formally
    undefined but the encoder must still be deterministic.

Why this exists: the hhg9 pytests cover the poles (test_octant_seams::
test_poles_mode_only, at exactly ±90) and the seams (test_seam_greenwich_bulk),
but the bulk seam sweep deliberately clamps to ±88° "to avoid poles/vertices" —
so the pole×seam INTERACTION band, 88°..90°, is the one region neither test
reaches. That band is where the seam/pole/authalic bugs historically lived.

Writes test_data/pole_seam_ref.bin, same layout as via_sphere_ref.bin:

    int64 n
    then n rows of 8 f64:
      lon, lat, cx, cy, oid, rt_lon, rt_lat, pad0

Run with `hhg9` importable:  python3 tools/support/gen_pole_seam_ref.py
"""
import struct
import sys

import numpy as np

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
sys.path.insert(0, HHG9)

from hhg9 import Registrar, Points          # noqa: E402

OUT = '/Users/ben/Documents/Projects/libhex9/test_data/pole_seam_ref.bin'

# The octant seams in longitude: the four meridians where faces meet, plus the
# antimeridian reached from both sides.
SEAM_LONS = np.array([-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0, 180.0])

# Offsets either side of a seam. 1e-13 deg is ~11 nm — below the L30 cell size,
# so these pairs must land in the same cell; 1e-9 deg is ~0.11 mm, which is
# larger than an L23 cell, so those pairs may legitimately differ deep in the
# tail. Both are useful: the first pins tie-break stability, the second pins
# that the seam itself is not a discontinuity in the projection.
SEAM_EPS = np.array([0.0, -1e-13, 1e-13, -1e-9, 1e-9])

# atan(1/sqrt(2)) — the latitude of the octahedron's equatorial vertices, i.e.
# the octant corner / triple point.
CORNER_LAT = np.degrees(np.arctan(1.0 / np.sqrt(2.0)))


def build_points():
    """Return an (n, 2) array of [lat, lon]."""
    rows = []

    # 1. Polar approach at a spread of longitudes (including off-seam ones, so
    #    a pole bug that only shows away from a seam is still caught).
    eps = np.geomspace(1e-12, 2.0, 60)
    pole_lons = np.array([0.0, 10.0, 45.0, 90.0, 123.456, 180.0, -73.21])
    for lon in pole_lons:
        for e in eps:
            rows.append([90.0 - e, lon])
            rows.append([-(90.0 - e), lon])

    # 2. The poles themselves, from every seam longitude. Longitude is formally
    #    undefined at a pole; the encoder must still be deterministic and must
    #    agree with Python about which cell it picks.
    for lon in SEAM_LONS:
        rows.append([90.0, lon])
        rows.append([-90.0, lon])

    # 3. Seams: each seam meridian, straddled, over a latitude ladder that
    #    reaches INTO the polar band the pytests stop short of.
    seam_lats = np.concatenate([
        np.linspace(-88.0, 88.0, 45),
        np.array([89.0, 89.9, 89.99, 89.999, 89.999999]),
        -np.array([89.0, 89.9, 89.99, 89.999, 89.999999]),
        np.array([0.0, 1e-9, -1e-9]),            # equator, straddled
    ])
    for lon in SEAM_LONS:
        for d in SEAM_EPS:
            for lat in seam_lats:
                rows.append([lat, lon + d])

    # 4. Octant corners — the triple points — straddled in both axes.
    for lon in SEAM_LONS:
        for dlon in SEAM_EPS:
            for dlat in SEAM_EPS:
                rows.append([CORNER_LAT + dlat, lon + dlon])
                rows.append([-CORNER_LAT + dlat, lon + dlon])

    pts = np.array(rows, dtype=np.float64)
    # Keep latitudes physical; longitudes may sit just outside ±180 by an
    # epsilon and MUST be left that way — wrapping them here would silently
    # remove the antimeridian straddle this file exists to test.
    pts[:, 0] = np.clip(pts[:, 0], -90.0, 90.0)
    return pts


def main():
    reg = Registrar()
    reg.set_ellipsoid(a=6378137.0, inv_f=298.257223563, name='WGS84',
                      via_sphere=True)
    b = reg.domain('b_oct')
    b.no_lib()
    reg.projection('gcd_brw').set_parallel(threshold=1_000_000)
    g = reg.domain('g_gcd')

    pts = build_points()

    fwd = reg.project(Points(pts.copy(), g), ['g_gcd', 'b_oct'])
    back = reg.project(fwd, ['b_oct', 'g_gcd'])

    rows = np.column_stack([
        pts[:, 1], pts[:, 0],                     # lon, lat
        fwd.coords[:, 0], fwd.coords[:, 1],
        np.asarray(fwd.oid, dtype=np.float64),
        back.coords[:, 1], back.coords[:, 0],     # rt_lon, rt_lat
        np.zeros(len(pts)),
    ])
    with open(OUT, 'wb') as fh:
        fh.write(struct.pack('<q', len(rows)))
        fh.write(rows.astype('<f8').tobytes(order='C'))
    print(f'wrote {OUT}: {len(rows)} rows')


if __name__ == '__main__':
    main()
