#!/usr/bin/env python3
"""
export_warp_deltas.py — generate the compact sidecar that replaces the 578MB
`h9_wgs84_ct_mesh.h` (plus the two smaller `h9_wgs84_warp_*.h` headers) used by
the C++ extension.

Premise (see memory/warp_distribution_plan.md):
    `source_pts` from `WGS84_l5_warp_data.npz` is analytically derivable from
    `hhg9.h9.polygon.tri_mesh(levels=5, mode=0)`. Therefore the only payload
    that must be shipped is the per-vertex displacement (`target_pts - source_pts`)
    in the same vertex order that `tri_mesh` produces.

This script:
  1. Verifies the premise — `tri_mesh(5, 0)` matches `npz['source_pts']` (modulo ordering)
     and prints the worst absolute mismatch.
  2. Reorders the npz deltas into `tri_mesh` vertex order.
  3. Writes a versioned binary sidecar `WGS84_l5_warp.h9warp`.

Run from the hex9_cli source directory:
    python export_warp_deltas.py

----------------------------------------------------------------------------
SIDECAR FORMAT  ".h9warp"   (little-endian throughout)
----------------------------------------------------------------------------
off  size  field         meaning
  0    4   magic         ASCII "H9WP"
  4    2   version       uint16, currently 0x0002
  6    1   level         uint8  (5 for the L5 lattice)
  7    1   mode          uint8  (0 = canonical "down" octant)
  8    4   count         uint32 number of stored vertex deltas
 12    1   dtype         uint8  0=f64, 1=f32, 2=int16, 3=int32
 13    1   flags         uint8  bit0: 1=mirrored (only x>=0 stored; x<0 derived
                                          by (-dx, +dy) of mirror)
                                bits 1-7: reserved (must be 0)
 14    2   reserved      zero-padded
 16    8   scale_x       float64 (multiplier applied to stored value)
 24    8   scale_y       float64 (multiplier applied to stored value)
 32    8   ellipsoid_a   float64 (informational: WGS84 semi-major)
 40    8   ellipsoid_f   float64 (informational: WGS84 flattening)
 48    4   crc32_payload uint32 CRC-32 of the payload bytes
 52    4   reserved      zero-padded
 56   ...  payload       count * 2 * sizeof(dtype) bytes

Payload ordering: tri_mesh produces all N vertices in lex-sorted quantized-key
order (see hhg9.h9.polygon._unique_rows_tol). When `flags & 1` is clear, the
payload is one (dx, dy) per vertex in that order (count == N). When set, the
payload contains one (dx, dy) per vertex with quantized_x >= 0, in tri_mesh
order; vertices with quantized_x < 0 are derived from their (-x, y) mirror as
(-dx, +dy). The seam (x == 0) is stored once with dx ~ 0 by symmetry.

For dtype=int16: stored value v must be multiplied by scale_axis to recover the
true delta in b_oct units. For f32/f64: scale_x/scale_y are 1.0 (informational).

v1 readers see flags=0 implicitly (the byte was zero-padded reserved).
----------------------------------------------------------------------------
"""

import os
import struct
import zlib
import math
import numpy as np

from hhg9.h9.polygon import tri_mesh

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH = os.path.join(SCRIPT_DIR, 'WGS84_l5_warp_data.npz')
OUT_PATH = os.path.join(SCRIPT_DIR, 'WGS84_l5_warp.h9warp')

MAGIC = b'H9WP'
VERSION = 2
LEVEL = 5
MODE = 0

DTYPE_F64 = 0
DTYPE_F32 = 1
DTYPE_I16 = 2
DTYPE_I32 = 3

FLAG_MIRRORED = 1 << 0

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563


def _quantize_match(a: np.ndarray, b: np.ndarray, tol: float = 1e-10):
    """Return per-row sort key based on power-of-two quantization (same scheme as
    `_unique_rows_tol` in polygon.py)."""
    k = int(np.ceil(-np.log2(tol)))
    return np.rint(np.ldexp(a, k)).astype(np.int64), np.rint(np.ldexp(b, k)).astype(np.int64)


def verify_lattice(npz_src: np.ndarray):
    """Return (perm, verts) such that `npz_src[perm] == verts` (within tol),
    where verts is tri_mesh's lex-sorted output. Aborts on mismatch."""
    verts, _edges, _tris = tri_mesh(levels=LEVEL, mode=MODE)
    print(f"tri_mesh(L={LEVEL}, mode={MODE}): {verts.shape[0]} verts, "
          f"{_tris.shape[0]} tris, {_edges.shape[0]} edges")
    print(f"npz['source_pts']:                {npz_src.shape[0]} pts")

    if verts.shape != npz_src.shape:
        raise RuntimeError(
            f"shape mismatch: tri_mesh -> {verts.shape}, npz -> {npz_src.shape}\n"
            f"Premise broken — tri_mesh does not reproduce source_pts at L={LEVEL}."
        )

    # Build a quantized lookup so we can find each tri_mesh vertex in npz_src.
    nq0, nq1 = _quantize_match(npz_src[:, 0], npz_src[:, 1])
    npz_key = np.stack([nq0, nq1], axis=1)
    vq0, vq1 = _quantize_match(verts[:, 0], verts[:, 1])
    verts_key = np.stack([vq0, vq1], axis=1)

    # lexsort indices for both
    npz_order = np.lexsort((npz_key[:, 0], npz_key[:, 1]))
    verts_order = np.lexsort((verts_key[:, 0], verts_key[:, 1]))
    if not np.array_equal(npz_key[npz_order], verts_key[verts_order]):
        # find first diff
        n = npz_key[npz_order]
        v = verts_key[verts_order]
        diff = np.where(np.any(n != v, axis=1))[0]
        first = diff[0] if diff.size else -1
        raise RuntimeError(
            f"point-set mismatch: first diverging sorted row at {first}\n"
            f"  npz   : {npz_src[npz_order][first]}\n"
            f"  tri_mesh: {verts[verts_order][first]}\n"
            f"Premise broken — same count but different points."
        )

    # Map tri_mesh vertex i -> index into npz_src
    perm = np.empty(verts.shape[0], dtype=np.int64)
    perm[verts_order] = npz_order
    err = np.max(np.abs(npz_src[perm] - verts))
    print(f"Max |npz_src[perm] - tri_mesh_verts| = {err:.3e}  "
          f"(expected ~0; quantization tol = 1e-10)")
    return perm, verts


def pack(out_path: str, deltas: np.ndarray, dtype_code: int, flags: int = 0):
    if dtype_code == DTYPE_F64:
        scale_x = scale_y = 1.0
        payload = deltas.astype(np.float64, copy=False).tobytes(order='C')
    elif dtype_code == DTYPE_F32:
        scale_x = scale_y = 1.0
        payload = deltas.astype(np.float32, copy=False).tobytes(order='C')
    elif dtype_code in (DTYPE_I16, DTYPE_I32):
        if dtype_code == DTYPE_I16:
            i_max = 32767
            np_dtype = np.int16
            label = 'int16'
        else:
            i_max = 2147483647
            np_dtype = np.int32
            label = 'int32'
        # Per-axis scale with ~5% headroom over max |delta|.
        max_x = float(np.max(np.abs(deltas[:, 0])))
        max_y = float(np.max(np.abs(deltas[:, 1])))
        scale_x = max_x * 1.05 / i_max
        scale_y = max_y * 1.05 / i_max
        q = np.empty_like(deltas, dtype=np.int64)
        q[:, 0] = np.rint(deltas[:, 0] / scale_x).astype(np.int64)
        q[:, 1] = np.rint(deltas[:, 1] / scale_y).astype(np.int64)
        q = q.astype(np_dtype)
        recovered = q.astype(np.float64) * np.array([scale_x, scale_y])
        rt_err = float(np.max(np.abs(recovered - deltas)))
        ground_nm = rt_err * 6.37e6 * 1e9
        print(f"  {label} round-trip max error: {rt_err:.3e} b_oct "
              f"(~{ground_nm:.1f} nm on ground)")
        payload = q.tobytes(order='C')
    else:
        raise ValueError(f"unknown dtype_code {dtype_code}")

    crc = zlib.crc32(payload) & 0xFFFFFFFF

    header = bytearray(56)
    header[0:4] = MAGIC
    struct.pack_into('<H', header, 4, VERSION)
    header[6] = LEVEL
    header[7] = MODE
    struct.pack_into('<I', header, 8, deltas.shape[0])
    header[12] = dtype_code
    header[13] = flags
    # 14..15 reserved (zero)
    struct.pack_into('<d', header, 16, scale_x)
    struct.pack_into('<d', header, 24, scale_y)
    struct.pack_into('<d', header, 32, WGS84_A)
    struct.pack_into('<d', header, 40, WGS84_F)
    struct.pack_into('<I', header, 48, crc)
    # 52..55 reserved (zero)

    with open(out_path, 'wb') as fh:
        fh.write(header)
        fh.write(payload)
    return len(header) + len(payload)


def main():
    print(f"Loading {os.path.normpath(NPZ_PATH)} ...")
    repo = np.load(NPZ_PATH, allow_pickle=True)
    src = np.asarray(repo['source_pts'], dtype=np.float64)
    dst = np.asarray(repo['target_pts'], dtype=np.float64)
    repo.close()

    print("\n--- Step 1: verify tri_mesh reproduces source_pts ---")
    perm, verts = verify_lattice(src)

    print("\n--- Step 2: reorder deltas into tri_mesh vertex order ---")
    deltas_full = (dst - src)[perm]
    print(f"  delta range: x [{deltas_full[:, 0].min():+.6e}, "
          f"{deltas_full[:, 0].max():+.6e}]  "
          f"y [{deltas_full[:, 1].min():+.6e}, {deltas_full[:, 1].max():+.6e}]")
    print(f"  max |delta|: x={np.max(np.abs(deltas_full[:, 0])):.6e}  "
          f"y={np.max(np.abs(deltas_full[:, 1])):.6e}")

    # Mirror-symmetric half: keep vertices with quantized x >= 0.
    # Use the same quantization scheme as _unique_rows_tol so the predicate is
    # exact and matches what the C++ loader will compute.
    tol = 1e-12
    k = int(np.ceil(-np.log2(tol)))
    qx = np.rint(np.ldexp(verts[:, 0], k)).astype(np.int64)
    keep_mask = qx >= 0
    deltas_half = deltas_full[keep_mask]
    n_half = int(keep_mask.sum())
    n_seam = int((qx == 0).sum())
    print(f"  mirror split: kept {n_half:,} (x>=0)  "
          f"seam {n_seam}  dropped {len(verts) - n_half:,} (x<0)")

    seam_only = qx[keep_mask] == 0
    if seam_only.any():
        seam_dx = float(np.max(np.abs(deltas_half[seam_only, 0])))
        print(f"  seam (x=0) max |dx| in stored half = {seam_dx:.3e} (expect ~0)")

    print("\n--- Step 3: write sidecars (full and mirrored, each dtype) ---")
    variants = [
        ('full',     deltas_full, 0),
        ('mirrored', deltas_half, FLAG_MIRRORED),
    ]
    for variant_label, ds, flags in variants:
        for code, label in [(DTYPE_F64, 'f64'),
                            (DTYPE_F32, 'f32'),
                            (DTYPE_I32, 'i32'),
                            (DTYPE_I16, 'i16')]:
            out = OUT_PATH.replace(
                '.h9warp', f'.{variant_label}.{label}.h9warp')
            size = pack(out, ds, code, flags=flags)
            print(f"  {variant_label:8s} {label}: {size:>10,} bytes  "
                  f"-> {os.path.basename(out)}")

    print("\nDone. For a 1 µm-on-ground target, mirrored.i32 (~1.07 MB) gives "
          "~27 nm precision at f32 size. mirrored.f64 (~2.14 MB) is exact.")


if __name__ == '__main__':
    main()
