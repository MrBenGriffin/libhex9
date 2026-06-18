#!/usr/bin/env python3
"""
export_warp_v3.py — write the v3 ".h9warp" sidecar: per-vertex deltas PLUS the
final CT gradients (F6 edge-tangent field, Option A of
tools/support/warp-port-brief-f6-cside.md).

Why gradients are shipped: tangent vertex deltas alone do not give a tangent
interpolant — hhg9's loader applies ghost-orbit padding, x-mirror gradient
symmetrisation and the edge-tangent gradient projection before CT evaluation.
Shipping the FINAL per-vertex gradients lets the C side skip all of that and
be bit-exact with hhg9 by construction (same values + same gradients + same
per-triangle CT formula => same field).

Inputs:
  hhg9/data/WGS84_l5_warp_data.npz                  (source_pts, target_pts)
  experimental/sinkhorn/output/WGS84_l5_warp_grads.npz
      (source_pts, target_pts, grad_dx, grad_dy — exported from the hhg9
       loader's self.ct_dx.grad / self.ct_dy.grad, REAL vertices only,
       source_pts order; see the port brief's exporter snippet)

----------------------------------------------------------------------------
V3 FORMAT — header layout IDENTICAL to v2 (export_warp_deltas.py), with:
  version = 0x0003
  flags  bit1 = HAS_GRADIENTS (bit0 mirror flag must be 0 in v3 for now)
  payload = count * 6 * f64, per vertex in tri_mesh lex order:
      (dx, dy, ddx/dx, ddx/dy, ddy/dx, ddy/dy)
  crc32 over the payload bytes, as in v2.
v1/v2 readers reject v3 via the version field.
----------------------------------------------------------------------------
Run from the hhg9 repo root:
    python3 /Users/ben/Documents/Projects/hex9_cli/hex9_cli/export_warp_v3.py

----------------------------------------------------------------------------
PROVENANCE: this script was authored + run inline in a prior libhex9 Claude
session (it produced the currently-shipped core/WGS84_l5_warp_f6.full.f64g.h9warp)
but was never committed as a file. Recovered verbatim from the session
transcript on 2026-06-18 and confirmed as THE producer: both the delta columns
and the gradient columns of the shipped blob are bit-identical (max diff 0.0
over 266,815 verts) to NPZ_DATA / NPZ_GRAD below. Requires a python env with
`hhg9` importable. See tools/support/README.md.
----------------------------------------------------------------------------
"""

import os
import struct
import sys
import zlib

import numpy as np

from hhg9.h9.polygon import tri_mesh

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
NPZ_DATA = os.path.join(HHG9, 'hhg9/data/WGS84_l5_warp_data.npz')
NPZ_GRAD = os.path.join(HHG9, 'experimental/sinkhorn/output/WGS84_l5_warp_grads.npz')
OUT_PATH = '/Users/ben/Documents/Projects/libhex9/core/WGS84_l5_warp_f6.full.f64g.h9warp'

MAGIC = b'H9WP'
VERSION = 3
LEVEL = 5
MODE = 0
DTYPE_F64 = 0
FLAG_HAS_GRADIENTS = 1 << 1
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563


def _qkeys(a, tol=1e-10):
    k = int(np.ceil(-np.log2(tol)))
    return np.rint(np.ldexp(a, k)).astype(np.int64)


def perm_npz_to_trimesh(npz_src):
    """perm such that npz_src[perm] == tri_mesh verts (lex order)."""
    verts, _e, _t = tri_mesh(levels=LEVEL, mode=MODE)
    assert verts.shape == npz_src.shape, (verts.shape, npz_src.shape)
    nk = np.stack([_qkeys(npz_src[:, 0]), _qkeys(npz_src[:, 1])], axis=1)
    vk = np.stack([_qkeys(verts[:, 0]), _qkeys(verts[:, 1])], axis=1)
    no = np.lexsort((nk[:, 0], nk[:, 1]))
    vo = np.lexsort((vk[:, 0], vk[:, 1]))
    assert np.array_equal(nk[no], vk[vo]), 'point-set mismatch vs tri_mesh'
    perm = np.empty(len(verts), dtype=np.int64)
    perm[vo] = no
    err = np.max(np.abs(npz_src[perm] - verts))
    print(f'tri_mesh order check: max coord err {err:.3e}')
    return perm


def main():
    data = np.load(NPZ_DATA, allow_pickle=True)
    grad = np.load(NPZ_GRAD, allow_pickle=True)
    src = np.asarray(data['source_pts'], dtype=np.float64)
    dst = np.asarray(data['target_pts'], dtype=np.float64)
    gsrc = np.asarray(grad['source_pts'], dtype=np.float64)
    gdx = np.asarray(grad['grad_dx'], dtype=np.float64)
    gdy = np.asarray(grad['grad_dy'], dtype=np.float64)
    if not (np.array_equal(src, gsrc)):
        print('FATAL: grads npz source_pts differ from data npz', file=sys.stderr)
        sys.exit(1)
    assert np.isfinite(gdx).all() and np.isfinite(gdy).all()

    perm = perm_npz_to_trimesh(src)
    deltas = (dst - src)[perm]
    gdx, gdy = gdx[perm], gdy[perm]

    n = len(deltas)
    payload_arr = np.empty((n, 6), dtype=np.float64)
    payload_arr[:, 0:2] = deltas
    payload_arr[:, 2:4] = gdx
    payload_arr[:, 4:6] = gdy
    payload = payload_arr.tobytes(order='C')
    crc = zlib.crc32(payload) & 0xFFFFFFFF

    header = bytearray(56)
    header[0:4] = MAGIC
    struct.pack_into('<H', header, 4, VERSION)
    header[6] = LEVEL
    header[7] = MODE
    struct.pack_into('<I', header, 8, n)
    header[12] = DTYPE_F64
    header[13] = FLAG_HAS_GRADIENTS
    struct.pack_into('<d', header, 16, 1.0)
    struct.pack_into('<d', header, 24, 1.0)
    struct.pack_into('<d', header, 32, WGS84_A)
    struct.pack_into('<d', header, 40, WGS84_F)
    struct.pack_into('<I', header, 48, crc)

    with open(OUT_PATH, 'wb') as fh:
        fh.write(header)
        fh.write(payload)
    print(f'wrote {OUT_PATH}: {56 + len(payload):,} bytes '
          f'({n:,} vertices x 48 B + 56 B header), crc 0x{crc:08x}')

    # spot-check (port brief): right-edge apex slide profile, tri_mesh order
    verts, _e, _t = tri_mesh(levels=LEVEL, mode=MODE)
    R3, W = np.sqrt(3.0), np.sqrt(6.0) / 3.0
    on = np.abs(R3 * verts[:, 0] - verts[:, 1] - W) < 1e-12
    t_hat = np.array([1.0, R3]) / 2.0
    M = 7.076e6
    for ytgt, expect in ((-0.81482, 2746.9), (0.40321, -5347.7)):
        i = np.flatnonzero(on)[np.argmin(np.abs(verts[on][:, 1] - ytgt))]
        got = float(deltas[i] @ t_hat) * M
        status = 'OK' if abs(got - expect) < 0.1 else 'MISMATCH'
        print(f'  spot y={verts[i, 1]:+.5f}: slide {got:+.1f} m (expect {expect:+.1f}) {status}')


if __name__ == '__main__':
    main()
