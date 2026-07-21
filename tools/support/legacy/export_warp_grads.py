#!/usr/bin/env python3
"""
export_warp_grads.py — produce WGS84_l5_warp_grads.npz, the gradient input that
export_warp_v3.py needs.

This is step 2 of the v3 warp pipeline:
    1. (sinkhorn training)            -> hhg9/data/WGS84_l5_warp_data.npz   (deltas)
    2. export_warp_grads.py           -> WGS84_l5_warp_grads.npz            (gradients)
    3. export_warp_v3.py              -> core/WGS84_l5_warp_f6.full.f64g.h9warp
    4. CMake .incbin                  -> embedded in libhex9

The FINAL per-vertex CT gradients are read off the hhg9 AuthalicWarp loader
AFTER it has applied ghost-orbit padding + x-mirror gradient symmetrisation +
edge-tangent gradient projection (Option A of the F6 C-side port brief). scipy
stores grad as (n_padded, 1, 2); the first n entries are the REAL vertices, in
source_pts order (ghosts follow).

PROVENANCE: recovered 2026-06-18 from the "Exporter" snippet in the F6 C-side
port brief (tools/support/warp-port-brief-f6-cside.md, lines 68-79) — it had only
ever existed as prose in the brief, never as a committed script. Run from the
hhg9 repo root (needs `hhg9` importable). See tools/support/README.md.
"""

import os

import numpy as np

from hhg9.domains.octahedral_barycentric import AuthalicWarp

HHG9 = '/Users/ben/Documents/Projects/PyCharm/hex9'
NPZ_DATA = os.path.join(HHG9, 'hhg9/data/WGS84_l5_warp_data.npz')
OUT_PATH = os.path.join(HHG9, 'experimental/sinkhorn/output/WGS84_l5_warp_grads.npz')


def main():
    # H9_WARP_EDGE is irrelevant here — we only read the constructed gradients.
    w = AuthalicWarp(NPZ_DATA)
    n = len(w.src)
    np.savez(OUT_PATH,
             source_pts=w.src, target_pts=w.dst,
             grad_dx=np.asarray(w.ct_dx.grad)[:n, 0, :],
             grad_dy=np.asarray(w.ct_dy.grad)[:n, 0, :])
    print(f'wrote {OUT_PATH}: {n:,} real vertices (grad_dx/grad_dy, source_pts order)')


if __name__ == '__main__':
    main()
