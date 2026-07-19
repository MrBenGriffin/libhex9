#!/usr/bin/env python3
"""
gen_via_sphere_ref.py — parity reference for test/via_sphere.c.

Runs the hhg9 via-sphere chain (Registrar(via_sphere=True), pure Python,
wedge-fold warp — the frozen reference) over landmarks + seeded random
points and writes:

    test_data/via_sphere_ref.bin
      int64 n
      then n rows of 8 f64:
        lon, lat, cx, cy, oid, rt_lon, rt_lat, pad0

The C test enables hex9_set_via_sphere(1) and checks hex9_project /
hex9_unproject against these rows.

Run with `hhg9` importable:  python3 tools/support/gen_via_sphere_ref.py
"""
import struct
import sys

import numpy as np

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
sys.path.insert(0, HHG9)

from hhg9 import Registrar, Points          # noqa: E402

OUT = '/Users/ben/Documents/Projects/libhex9/test_data/via_sphere_ref.bin'


def main():
    reg = Registrar()
    reg.set_ellipsoid(a=6378137.0, inv_f=298.257223563, name='WGS84',
                      via_sphere=True)
    b = reg.domain('b_oct')
    b.no_lib()
    reg.projection('gcd_brw').set_parallel(threshold=1_000_000)
    g = reg.domain('g_gcd')

    rng = np.random.default_rng(4007)
    n = 4000
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    lon = rng.uniform(-180.0, 180.0, n)
    extra = np.array([
        [89.9999, 10.0], [-89.9999, -120.0], [0.0, 45.0001],
        [0.0001, -45.0], [45.0, 90.0001], [-45.0, 0.0001],
        [51.4778, -0.0015], [35.6586, 139.7454], [88.0, 1.0],
        [-88.0, 179.0],
    ])
    pts = np.vstack([np.column_stack([lat, lon]), extra])   # [lat, lon]

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
