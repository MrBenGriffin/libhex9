#!/usr/bin/env python3
"""Canonical-address invariants for the Hex9 cell id (the GeoPlegma ZoneId).

Proves, against the compiled module, that a cell's canonical key
(x_list + (c2,r_mo), i.e. what `bin` emits) is simultaneously:
  * UNIQUE      - exhaustive enumeration at L4/L5 gives 12*9^L distinct keys
  * REVERSIBLE  - decode(encode(pt)) round-trips to < 1 um (geodesic, WGS84)
  * TRAVERSABLE - decode(bin)->encode->bin lands in the SAME cell (no fossil),
                  so parent/children compose from existing primitives
  * SEAM/DEFECT CLEAN - holds across octant seams and the six valence-4 vertices.

Run:  PYTHONPATH=build python3 test/canonical_invariants.py
Env:  H9_DEEP (default 6) caps the seam/vertex idempotence depth.
"""
import os, sys, numpy as np, hex9_ext as h
try:
    from geographiclib.geodesic import Geodesic
    GEOD = Geodesic.WGS84
except ImportError:
    GEOD = None

def vview(u):  return np.ascontiguousarray(u).view(np.dtype((np.void, 16))).ravel()
def nibs(u):
    o = np.empty((len(u), 32), np.uint8); o[:, 0::2] = u >> 4; o[:, 1::2] = u & 0xF; return o

FAILED = []
def check(name, ok, extra=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}")
    if not ok: FAILED.append(name)

def main():
    h.warp_init()
    rng = np.random.default_rng(0)

    # 1+2+3: whole-world enumeration -> uniqueness, centroid round-trip, traversability
    for L in (4, 5):
        uu, cent, _ = h.grid(-180., -90., 180., 90., L, 0)
        n, exp = len(uu), 12 * 9 ** L
        check(f"L{L}: bins globally unique", len(np.unique(vview(uu))) == n == exp,
              f"(n={n}, expect={exp})")
        lon, lat = np.ascontiguousarray(cent[:, 0]), np.ascontiguousarray(cent[:, 1])
        b = h.bin(h.encode(lon, lat), L)
        check(f"L{L}: bin(encode(centroid)) == own cell", int((vview(b) == vview(uu)).sum()) == n)
        dlon, dlat = h.decode(uu)                                   # bin -> interior pt
        re = h.bin(h.encode(np.ascontiguousarray(dlon), np.ascontiguousarray(dlat)), L)
        sib = int((vview(re) != vview(uu)).sum())
        check(f"L{L}: decode(bin)->encode->bin == bin (traversable, no fossil)", sib == 0,
              f"({sib} siblings)")
        if L == 4:
            check("L4: bins carry no h_term (nibble30==0xF) yet unique",
                  bool(np.all(nibs(uu)[:, 30] == 0xF)))

    # 4: seam + six valence-4 vertices, decode->encode->bin idempotence
    deep = int(os.environ.get("H9_DEEP", "6"))
    specials = [(0., 51.48), (180., 0.), (0., 0.), (90., 0.), (-90., 0.), (0., 89.9), (0., -89.9)]
    for L in (5, deep):
        bad = tot = 0
        for clon, clat in specials:
            p = rng.normal(0, 0.3, (3000, 2)) + (clon, clat)
            p[:, 0] = (p[:, 0] + 180) % 360 - 180; p[:, 1] = np.clip(p[:, 1], -89.999, 89.999)
            f1 = h.encode(np.ascontiguousarray(p[:, 0]), np.ascontiguousarray(p[:, 1]))
            dl, da = h.decode(f1)
            f2 = h.encode(np.ascontiguousarray(dl), np.ascontiguousarray(da))
            bad += int((vview(h.bin(f1, L)) != vview(h.bin(f2, L))).sum()); tot += len(p)
        check(f"L{L}: seam/vertex decode->encode->bin idempotent", bad == 0, f"({bad}/{tot})")

    # 5: geodesic round-trip error on 20k random points
    if GEOD is not None:
        N = 20000
        lat = np.degrees(np.arcsin(rng.uniform(-1, 1, N))); lon = rng.uniform(-180, 180, N)
        dlon, dlat = h.decode(h.encode(np.ascontiguousarray(lon), np.ascontiguousarray(lat)))
        e = np.fromiter((GEOD.Inverse(lat[i], lon[i], dlat[i], dlon[i])['s12'] for i in range(N)),
                        float, N)
        print(f"  round-trip geodesic: max={e.max()*1e6:.4f} um  median={np.median(e)*1e6:.4f} um")
        check("max round-trip error < 1 um", e.max() < 1e-6, f"(max={e.max()*1e6:.4f} um)")
    else:
        print("  (geographiclib not installed - skipping geodesic round-trip)")

    print("\nRESULT:", "ALL PASS" if not FAILED else f"FAILED: {FAILED}")
    sys.exit(1 if FAILED else 0)

if __name__ == "__main__":
    main()
