#!/usr/bin/env python3
"""Deep-layer invariants — the L30-reclamation oracle.

canonical_invariants.py proves the shallow (L4/L5/L6) canonical-id properties by
exhaustive enumeration; it cannot reach the reclaimed nibble 30. This probe
exercises the DEEPEST addressable layers (lmax-3 .. lmax) by sampling, so it is
the falsifiable check for the HEX9_USE_L29=OFF (L30) branch:

  * DEEP ROUND-TRIP   - decode(encode(pt)) geodesic error < 1 um at full depth
                        (finer grid => should be <= the L29 figure).
  * DEEP TRAVERSABLE  - decode(bin_L)->encode->bin_L == bin_L for L up to lmax
                        (no fossil at depth; the canonical id stays an address).
  * DEEP UNIQUE       - well-separated random points get distinct bins at lmax
                        (no addressing collisions in the reclaimed range).
  * LAYOUT WITNESS    - on an L30 build nibble[30] of a full uuid is a body digit
                        (0..8), addressing reaches 30, and bin(full,30) is valid;
                        on L29 it is h_term and 30 is out of range.

Layout-agnostic: drives h.lmax(), so a clean run on the L29 build validates the
probe itself before the L30 branch exists. Then re-run on the USE_L29=OFF build.

Run:  PYTHONPATH=build python3 test/l30_invariants.py
Env:  H9_N (default 20000) sample size.
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
def cc(a):     return np.ascontiguousarray(a)

FAILED = []
def check(name, ok, extra=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}")
    if not ok: FAILED.append(name)

def main():
    h.warp_init()
    LMAX = h.lmax()
    print(f"libhex9 {h.version()}  —  lmax={LMAX}")
    rng = np.random.default_rng(20260621)
    N = int(os.environ.get("H9_N", "20000"))

    # Uniform-on-sphere sample + the seam/vertex specials folded in.
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, N)))
    lon = rng.uniform(-180, 180, N)
    specials = np.array([(0., 51.48), (180., 0.), (0., 0.), (90., 0.),
                         (-90., 0.), (0., 89.9), (0., -89.9)])
    lon[:len(specials)], lat[:len(specials)] = specials[:, 0], specials[:, 1]
    lon, lat = cc(lon), cc(lat)
    full = h.encode(lon, lat)

    # 1: deep round-trip at full depth.
    if GEOD is not None:
        dlon, dlat = h.decode(full)
        e = np.fromiter((GEOD.Inverse(lat[i], lon[i], dlat[i], dlon[i])['s12']
                         for i in range(N)), float, N)
        print(f"  deep round-trip geodesic: max={e.max()*1e6:.4f} um  median={np.median(e)*1e6:.4f} um")
        check("deep round-trip < 1 um", e.max() < 1e-6, f"(max={e.max()*1e6:.4f} um)")
    else:
        print("  (geographiclib not installed — skipping geodesic round-trip)")

    # 2: deep traversability — decode(bin_L)->encode->bin_L idempotent, no fossil.
    for L in range(max(1, LMAX - 3), LMAX + 1):
        b = h.bin(full, L)
        dl, da = h.decode(b)
        re = h.bin(h.encode(cc(dl), cc(da)), L)
        sib = int((vview(re) != vview(b)).sum())
        check(f"L{L}: decode(bin)->encode->bin == bin (traversable, no fossil)",
              sib == 0, f"({sib}/{N} siblings)")

    # 3: deep uniqueness — distinct, well-separated points get distinct bins.
    bmax = h.bin(full, LMAX)
    u = len(np.unique(vview(bmax)))
    # Random sphere points at L>=25 are essentially never co-cellular (cells ~cm),
    # so collisions would signal an addressing defect, not a coincidence.
    check(f"L{LMAX}: {N} sampled points yield distinct bins", u == N, f"(unique={u}/{N})")

    # 4: layout witness — what the reclaimed nibble 30 carries.
    n30 = nibs(full)[:, 30]
    if LMAX >= 30:
        check("L30 build: full-uuid nibble[30] is a body digit (0..8)",
              bool(np.all(n30 <= 8)), f"(max={int(n30.max())})")
        try:
            h.bin(full[:1], 30); ok30 = True
        except Exception:
            ok30 = False
        check("L30 build: bin(full, 30) is in range", ok30)
    else:
        check("L29 build: full-uuid nibble[30] is h_term (0..11, the look-ahead RID)",
              bool(np.all(n30 <= 11)), f"(max={int(n30.max())})")
        try:
            h.bin(full[:1], 30); rng30 = True
        except Exception:
            rng30 = False
        check("L29 build: layer 30 is out of range", not rng30)

    print("\nRESULT:", "ALL PASS" if not FAILED else f"FAILED: {FAILED}")
    sys.exit(1 if FAILED else 0)

if __name__ == "__main__":
    main()
