"""
h9warp_io.py — Python reader for the .h9warp sidecar format.

Companion to export_warp_deltas.py.  Used both for testing (verifying that the
sidecar round-trips through AuthalicWarp) and as the reference implementation
that the C++ loader must match.

Format spec lives in export_warp_deltas.py's module docstring.
"""

import struct
import zlib
import numpy as np

from hhg9.h9.polygon import tri_mesh

_MAGIC = b'H9WP'
_HEADER_LEN = 56

_DTYPE_F64 = 0
_DTYPE_F32 = 1
_DTYPE_I16 = 2
_DTYPE_I32 = 3

_FLAG_MIRRORED = 1 << 0


def load_h9warp(path: str):
    """Read a .h9warp file and return (source_pts, target_pts, info).

    source_pts is reconstructed from `tri_mesh(level, mode)`; target_pts is
    source_pts + deltas. info is a dict of header fields for diagnostics.
    """
    with open(path, 'rb') as fh:
        blob = fh.read()
    if len(blob) < _HEADER_LEN or blob[:4] != _MAGIC:
        raise ValueError(f"{path}: not an H9WP file")
    version = struct.unpack_from('<H', blob, 4)[0]
    level = blob[6]
    mode = blob[7]
    count = struct.unpack_from('<I', blob, 8)[0]
    dtype_code = blob[12]
    flags = blob[13] if version >= 2 else 0
    scale_x = struct.unpack_from('<d', blob, 16)[0]
    scale_y = struct.unpack_from('<d', blob, 24)[0]
    ellipsoid_a = struct.unpack_from('<d', blob, 32)[0]
    ellipsoid_f = struct.unpack_from('<d', blob, 40)[0]
    crc_stored = struct.unpack_from('<I', blob, 48)[0]
    payload = blob[_HEADER_LEN:]
    crc_actual = zlib.crc32(payload) & 0xFFFFFFFF
    if crc_actual != crc_stored:
        raise ValueError(
            f"{path}: CRC mismatch (stored {crc_stored:#x}, actual {crc_actual:#x})")

    if dtype_code == _DTYPE_F64:
        arr = np.frombuffer(payload, dtype='<f8').reshape(count, 2).astype(np.float64)
        deltas_stored = arr
    elif dtype_code == _DTYPE_F32:
        arr = np.frombuffer(payload, dtype='<f4').reshape(count, 2).astype(np.float64)
        deltas_stored = arr
    elif dtype_code == _DTYPE_I16:
        arr = np.frombuffer(payload, dtype='<i2').reshape(count, 2).astype(np.float64)
        arr[:, 0] *= scale_x
        arr[:, 1] *= scale_y
        deltas_stored = arr
    elif dtype_code == _DTYPE_I32:
        arr = np.frombuffer(payload, dtype='<i4').reshape(count, 2).astype(np.float64)
        arr[:, 0] *= scale_x
        arr[:, 1] *= scale_y
        deltas_stored = arr
    else:
        raise ValueError(f"{path}: unknown dtype_code {dtype_code}")

    verts, _edges, _tris = tri_mesh(levels=level, mode=mode)
    n_full = verts.shape[0]

    if flags & _FLAG_MIRRORED:
        tol = 1e-12
        k = int(np.ceil(-np.log2(tol)))
        qx = np.rint(np.ldexp(verts[:, 0], k)).astype(np.int64)
        qy = np.rint(np.ldexp(verts[:, 1], k)).astype(np.int64)
        keep_mask = qx >= 0
        if int(keep_mask.sum()) != count:
            raise ValueError(
                f"{path}: mirrored count {count} != x>=0 vertex count {int(keep_mask.sum())}")
        # Payload is in tri_mesh order over the x>=0 subset.
        pos_indices = np.where(keep_mask)[0]   # tri_mesh index of each stored row
        deltas_full = np.empty((n_full, 2), dtype=np.float64)
        deltas_full[pos_indices] = deltas_stored
        # Map (qx, qy) -> subset row, then reflect to x<0 vertices.
        key_to_subset = {(int(qx[j]), int(qy[j])): subset_i
                         for subset_i, j in enumerate(pos_indices)}
        neg_idx = np.where(qx < 0)[0]
        for i in neg_idx:
            subset_i = key_to_subset.get((-int(qx[i]), int(qy[i])))
            if subset_i is None:
                raise ValueError(
                    f"{path}: no mirror counterpart for vertex {i} "
                    f"(qx={qx[i]}, qy={qy[i]})")
            deltas_full[i, 0] = -deltas_stored[subset_i, 0]
            deltas_full[i, 1] = +deltas_stored[subset_i, 1]
    else:
        if count != n_full:
            raise ValueError(
                f"{path}: full count {count} != tri_mesh vertex count {n_full}")
        deltas_full = deltas_stored

    source_pts = verts.copy()
    target_pts = source_pts + deltas_full
    info = dict(
        version=version, level=level, mode=mode, count=count,
        dtype_code=dtype_code, flags=flags,
        scale_x=scale_x, scale_y=scale_y,
        ellipsoid_a=ellipsoid_a, ellipsoid_f=ellipsoid_f,
        n_full=n_full,
        # `deltas` is the raw stored (dx, dy) per tri_mesh vertex, with the
        # mirror unfolded. Distinct from (target_pts - source_pts) which can
        # lose 1 ULP through float roundoff. Useful for bit-exact tests
        # against the C++ loader.
        deltas=deltas_full)
    return source_pts, target_pts, info


def write_npz_from_h9warp(h9warp_path: str, npz_path: str):
    """Convert .h9warp -> .npz with the same schema as WGS84_l5_warp_data.npz,
    so AuthalicWarp / b_oct.set_warp can ingest it directly."""
    src, dst, info = load_h9warp(h9warp_path)
    np.savez(npz_path, source_pts=src, target_pts=dst, layer=np.int64(info['level']))
    return info
