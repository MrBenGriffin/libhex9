#!/usr/bin/env python3
"""
export_warp_fund_v4.py — write the v4 fundamental-domain ".h9warp" sidecar:
the Sphere-L6 wedge field (deltas + final CT gradients) for the wedge-fold
warp (hex9 wedge-CT, h9_proj L6_MOBIUS_WARP.md §3.7b).

The v4 payload covers ONLY the (2 3 6) fundamental wedge of the mode-0
face (corners C = pole face corner, M = lateral-edge midpoint, G = face
centroid): 399,675 vertices at L6 instead of the 2,394,766-vertex face.
The C side folds queries into the wedge by the D3 group and unfolds the
displacement (core/h9_warp_fund.h), and derives halo vertices/gradients
by exact reflection — so, as with v3, shipping the FINAL per-vertex
gradients (hhg9's fold build: halo-padded Bell-Sibson estimate + fold-C1
/ seam-tangency projection) makes the C field match hhg9 by construction.

Inputs:
  hhg9/data/Sphere_l6_fund_warp_data.npz
      (source_pts, target_pts, edge_class, corners, layer)
  the hhg9 fold build itself (AuthalicWarp(npz, fold=True)) for the
      final projected gradients — ~10 s.

----------------------------------------------------------------------------
V4 FORMAT — header layout IDENTICAL to v2/v3, with:
  version = 0x0004
  level   = 6, mode = 0
  count   = wedge vertex count (399,675 at L6)
  flags   bit1 = HAS_GRADIENTS (required), bit2 = FUND (wedge payload)
  payload = count * 6 * f64 per vertex:
      (dx, dy, ddx/dx, ddx/dy, ddy/dx, ddy/dy)
  vertex order: tri_mesh(level, mode) lex order FILTERED by the closed
      wedge test  x >= -1e-9  and  p·N2 >= -1e-9  (N2 the G-M median
      normal, sign toward C) — the C loader re-derives the same list
      from its own tri_mesh, so no coordinates are shipped.
  crc32 over the payload bytes, as in v2/v3.
----------------------------------------------------------------------------
Run from anywhere with `hhg9` importable:
    python3 tools/support/export_warp_fund_v4.py
"""

import os
import struct
import sys
import zlib

import numpy as np

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
sys.path.insert(0, HHG9)

from hhg9.h9.polygon import tri_mesh                      # noqa: E402
from hhg9.h9 import H9K                                   # noqa: E402
from hhg9.domains.octahedral_barycentric import AuthalicWarp  # noqa: E402

NPZ_FUND = os.path.join(HHG9, 'hhg9/data/Sphere_l6_fund_warp_data.npz')
OUT_PATH = '/Users/ben/Documents/Projects/libhex9/core/Sphere_l6_fund.f64g.h9warp'

MAGIC = b'H9WP'
VERSION = 4
LEVEL = 6
MODE = 0
DTYPE_F64 = 0
FLAG_HAS_GRADIENTS = 1 << 1
FLAG_FUND = 1 << 2
SPHERE_A = 6371007.1809   # authalic radius (informational)
WEDGE_TOL = 1e-9


def _qkeys(a, tol=1e-10):
    k = int(np.ceil(-np.log2(tol)))
    return np.rint(np.ldexp(a, k)).astype(np.int64)


def wedge_mask(verts):
    """Closed-wedge membership, matching stage3b_fund_l6.in_wedge and the
    C loader (h9_warp_fund.h)."""
    c = np.array([0.0, H9K.limits.VF])
    b = np.array([H9K.limits.TR, H9K.limits.VC])
    m = 0.5 * (c + b)
    u2 = m / np.linalg.norm(m)
    n2 = np.array([-u2[1], u2[0]])
    if n2 @ c < 0:
        n2 = -n2
    return (verts[:, 0] >= -WEDGE_TOL) & (verts @ n2 >= -WEDGE_TOL)


def main():
    data = np.load(NPZ_FUND, allow_pickle=True)
    src = np.asarray(data['source_pts'], dtype=np.float64)
    dst = np.asarray(data['target_pts'], dtype=np.float64)
    assert int(data['layer']) == LEVEL

    print('building hhg9 fold warp (final gradients)...', flush=True)
    warp = AuthalicWarp(NPZ_FUND, fold=True)
    n_real = len(warp.src)
    assert n_real == len(src)
    gdx = np.asarray(warp.ct_dx.grad[:n_real, 0, :], dtype=np.float64)
    gdy = np.asarray(warp.ct_dy.grad[:n_real, 0, :], dtype=np.float64)
    assert np.isfinite(gdx).all() and np.isfinite(gdy).all()
    # the fold build snaps on-line targets exactly — export ITS deltas
    deltas_npz = warp.dst - warp.src

    # canonical order: tri_mesh lex order filtered by the wedge test
    verts = tri_mesh(LEVEL, MODE)[0]
    wsel = wedge_mask(verts)
    wverts = verts[wsel]
    if len(wverts) != len(src):
        print(f'FATAL: wedge vert count {len(wverts)} != npz count {len(src)}',
              file=sys.stderr)
        sys.exit(1)

    nk = np.stack([_qkeys(src[:, 0]), _qkeys(src[:, 1])], axis=1)
    vk = np.stack([_qkeys(wverts[:, 0]), _qkeys(wverts[:, 1])], axis=1)
    no = np.lexsort((nk[:, 0], nk[:, 1]))
    vo = np.lexsort((vk[:, 0], vk[:, 1]))
    assert np.array_equal(nk[no], vk[vo]), 'point-set mismatch vs tri_mesh wedge'
    perm = np.empty(len(wverts), dtype=np.int64)
    perm[vo] = no
    err = np.max(np.abs(src[perm] - wverts))
    print(f'tri_mesh wedge order check: max coord err {err:.3e}')

    n = len(wverts)
    payload_arr = np.empty((n, 6), dtype=np.float64)
    payload_arr[:, 0:2] = deltas_npz[perm]
    payload_arr[:, 2:4] = gdx[perm]
    payload_arr[:, 4:6] = gdy[perm]
    payload = payload_arr.tobytes(order='C')
    crc = zlib.crc32(payload) & 0xFFFFFFFF

    header = bytearray(56)
    header[0:4] = MAGIC
    struct.pack_into('<H', header, 4, VERSION)
    header[6] = LEVEL
    header[7] = MODE
    struct.pack_into('<I', header, 8, n)
    header[12] = DTYPE_F64
    header[13] = FLAG_HAS_GRADIENTS | FLAG_FUND
    struct.pack_into('<d', header, 16, 1.0)
    struct.pack_into('<d', header, 24, 1.0)
    struct.pack_into('<d', header, 32, SPHERE_A)
    struct.pack_into('<d', header, 40, 0.0)
    struct.pack_into('<I', header, 48, crc)

    with open(OUT_PATH, 'wb') as fh:
        fh.write(header)
        fh.write(payload)
    print(f'wrote {OUT_PATH}: {56 + len(payload):,} bytes '
          f'({n:,} wedge vertices x 48 B + 56 B header), crc 0x{crc:08x}')


if __name__ == '__main__':
    main()
