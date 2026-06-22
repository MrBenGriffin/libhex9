#!/usr/bin/env python3
"""Golden regression anchor for the L29 UUID layout.

A `USE_L29=ON` build must reproduce GOLDEN_L29 byte-for-byte: this is the
trip-wire that proves the H9_LMAX cap refactor (and any later plumbing) did
NOT change observable behaviour on the legacy layout.

This intentionally does NOT apply to a USE_L29=OFF (L30) build, which is a
different on-disk format; skip it there (it is expected to differ).

Run:  PYTHONPATH=build python3 test/golden_l29.py
"""
import sys, hashlib, numpy as np, hex9_ext as h

GOLDEN_L29 = "a85141ce85514a965beb6daa5c05edf038904aaa5938c0d2dd9ade9bcffd43e8"

def digest():
    h.warp_init()
    rng = np.random.default_rng(12345); N = 50000
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, N))); lon = rng.uniform(-180, 180, N)
    lon, lat = np.ascontiguousarray(lon), np.ascontiguousarray(lat)
    full = h.encode(lon, lat)
    parts = [full.tobytes()]
    for L in range(1, 30):
        parts.append(h.bin(full, L).tobytes())
    dlon, dlat = h.decode(full)
    parts += [dlon.tobytes(), dlat.tobytes()]
    return hashlib.sha256(b"".join(parts)).hexdigest()

if __name__ == "__main__":
    d = digest()
    ok = d == GOLDEN_L29
    print(f"libhex9 {h.version()}")
    print(f"digest : {d}")
    print(f"golden : {GOLDEN_L29}")
    print("MATCH" if ok else "MISMATCH (regression on L29, or an L30 build — skip there)")
    sys.exit(0 if ok else 1)
