#!/usr/bin/env python3
"""
gen_fund_warp_ref.py — parity reference for test/fund_warp.cpp.

Samples deterministic points over the mode-0 face (interior + on-edge +
near-corner + near-mirror-line), runs the hhg9 FOLD warp (the frozen
Python reference, wedge-CT §3.7b) forward and inverse, and writes:

    test_data/fund_warp_ref.bin
      int64 n
      then n rows of 8 f64:
        x, y, mo, do_x, do_y, undo_of_do_x, undo_of_do_y, pad0

The C++ test builds the v4 blob state and checks warp_do/warp_undo
against these rows.

Run with `hhg9` importable:  python3 tools/support/gen_fund_warp_ref.py
"""
import struct
import sys

import numpy as np

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
sys.path.insert(0, HHG9)

from hhg9.h9 import H9K                                       # noqa: E402
from hhg9.domains.octahedral_barycentric import AuthalicWarp  # noqa: E402

NPZ = HHG9 + '/hhg9/data/Sphere_l6_fund_warp_data.npz'
OUT = '/Users/ben/Documents/Projects/libhex9/test_data/fund_warp_ref.bin'

VF, VC, TR = H9K.limits.VF, H9K.limits.VC, H9K.limits.TR
R3, W = H9K.radical.R3, H9K.Ẇ
rng = np.random.default_rng(4007)


def face_sample(n):
    a, b = rng.random(n), rng.random(n)
    flip = a + b > 1.0
    a[flip], b[flip] = 1.0 - a[flip], 1.0 - b[flip]
    v0, v1, v2 = np.array([0.0, VF]), np.array([TR, VC]), np.array([-TR, VC])
    return v0 + np.outer(a, v1 - v0) + np.outer(b, v2 - v0)


def main():
    warp = AuthalicWarp(NPZ, fold=True)

    pts = [face_sample(3000)]
    # on-edge (both laterals), near-edge, near mirror lines, near corners
    t = rng.uniform(0.0, TR, 300)
    pts.append(np.column_stack([t, R3 * t - W]))
    pts.append(np.column_stack([-t, R3 * t - W]))
    t = rng.uniform(1e-3, TR - 1e-3, 300)
    pts.append(np.column_stack([t, R3 * t - W + rng.uniform(0, 2e-4, 300)]))
    y = rng.uniform(VF * 0.98, VC * 0.98, 300)
    pts.append(np.column_stack([rng.uniform(-2e-4, 2e-4, 300), y]))
    for corner in ([0.0, VF], [TR, VC], [-TR, VC]):
        c = np.asarray(corner)
        pts.append(c + 3e-3 * (face_sample(200) - c))
    pts.append(np.column_stack([rng.uniform(-TR * 0.9, TR * 0.9, 300),
                                np.full(300, VC)]))          # equator line
    p = np.vstack(pts)

    rows = []
    for mo in (0, 1):
        q = p if mo == 0 else np.column_stack([p[:, 0], -p[:, 1]])
        fwd = warp.do(q, mo=mo)
        back = warp.undo(fwd, mo=mo)
        n = len(q)
        rows.append(np.column_stack([q, np.full(n, float(mo)), fwd, back,
                                     np.zeros(n)]))
    allr = np.vstack(rows)
    with open(OUT, 'wb') as fh:
        fh.write(struct.pack('<q', len(allr)))
        fh.write(allr.astype('<f8').tobytes(order='C'))
    print(f'wrote {OUT}: {len(allr)} rows '
          f'(RT floor: {np.abs(allr[:, 5:7] - allr[:, 0:2]).max():.3e})')


if __name__ == '__main__':
    main()
