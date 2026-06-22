"""Python smoke test for the nanobind hex9_ext module.

Run: PYTHONPATH=build python3 test/smoke.py
Exercises numpy-array encode/decode, checks the round-trip, cross-checks one
known point, and times the batch path. Importing numpy + hex9_ext together is
the moment a duplicate-OpenMP-runtime clash would abort — if we get past it,
we're clear.
"""
import time
import numpy as np
import hex9_ext as h9   # <- if this aborts with "OMP: Error #15", we have a libomp clash

print("version:", h9.version())
print("lmax:", h9.lmax(), "(29 = legacy layout, 30 = reclaimed h_term)")
h9.warp_init()

n = 100_000
rng = np.random.default_rng(0)
lon = rng.uniform(-180.0, 180.0, n)
lat = rng.uniform(-85.0, 85.0, n)

t0 = time.perf_counter()
uu = h9.encode(lon, lat)
dt = time.perf_counter() - t0
print(f"encode {n} pts -> {uu.shape} {uu.dtype} in {dt:.3f}s ({n/dt:,.0f} pts/s, OMP)")

dlon, dlat = h9.decode(uu)
# Fold longitude wraparound: an antimeridian point may decode to the other sign
# (e.g. the L30 canonical centroid returns +180 for a -180 input) — zero geodesic
# error, so the naive |Δlon| must be reduced into [0,180].
dlo = np.abs(dlon - lon); dlo = np.minimum(dlo, 360.0 - dlo)
err = float(dlo.max() + np.abs(dlat - lat).max())
print(f"decode round-trip max err: {err:.2e} deg")

# Cross-check one known point against the hhg9 oracle (432177478...).
# (F6 field re-baseline 2026-06-12: full uuid byte-identical to hhg9;
# was 432177468 under the pre-F6 field.)
u0 = h9.encode(np.array([-3.19]), np.array([55.95]))[0]
print("edinburgh uuid:", u0.tobytes().hex())

# bin to layer 8 and confirm shape
b = h9.bin(uu[:5], 8)
print("bin(layer=8) sample shape:", b.shape, b.dtype)

assert err < 1e-4, "round-trip too large"
assert u0.tobytes().hex().startswith("432177478"), "edinburgh prefix mismatch"

# cell(): single hexagon ring
binned = h9.bin(u0.reshape(1, 16), 8)[0]
ring0 = h9.cell(binned, 8, 0)
ring1 = h9.cell(binned, 8, 1)
print(f"cell d0 ring {ring0.shape}, d1 ring {ring1.shape}, v0={tuple(round(x,7) for x in ring0[0])}")
assert ring0.shape == (7, 2) and ring1.shape == (19, 2)

# neighbors / k_ring / k_disk (batch shape check on 5 random cells)
# input is FULL uuids only — bins are layer keys, not addresses
nbs_r, counts_r = h9.neighbors(uu[:5], 8)
assert nbs_r.shape == (5, 6, 16) and counts_r.shape == (5,)
assert set(counts_r.tolist()) <= {5, 6}, f"unexpected neighbour counts {counts_r}"
# single-cell semantics on the Edinburgh cell (full uuid input)
nbs, counts = h9.neighbors(u0.reshape(1, 16), 8)
ring1 = h9.k_ring(u0, 8, 1)
disk2 = h9.k_disk(u0, 8, 2)
assert ring1.shape[0] == counts[0], "k_ring(1) != neighbour count"
assert disk2.shape == (h9.disk_ncells(2), 16), f"k_disk(2) {disk2.shape}"
# input contract: bin input is rejected (per-row count -1 / raise)
_, rej = h9.neighbors(binned.reshape(1, 16), 8)
assert rej[0] == -1, "bin input to neighbors not rejected"
try:
    h9.k_ring(binned, 8, 1)
    raise AssertionError("bin input to k_ring not rejected")
except RuntimeError:
    pass
print(f"neighbors: counts={counts.tolist()}, ring1={ring1.shape[0]}, disk2={disk2.shape[0]}")

# labels: round-trip + ancestor over a disk
canon = h9.k_disk(u0, 8, 0)[0]
lbl = h9.label(canon, 8)
uu_p, layer_p = h9.parse_label(lbl)
assert layer_p == 8 and bytes(uu_p) == bytes(canon), "label round-trip"
lon0, lat0 = h9.label_centroid(lbl)
assert abs(lon0 - -3.19) < 1.0 and abs(lat0 - 55.95) < 1.0, f"centroid ({lon0},{lat0})"
anc = h9.common_ancestor(disk2, 8)
assert anc is not None and anc[1] < 8 and lbl.startswith(anc[0]), f"ancestor {anc}"
print(f"label: {lbl} centroid=({lon0:.4f},{lat0:.4f}) ancestor={anc[0]}@L{anc[1]}")

# grid(): bbox enumeration — 1247 under the F6 field (2026-06-12
# re-baseline; was 1255 pre-F6)
uuids, cents, rings = h9.grid(-0.3, 51.3, 0.1, 51.7, 8)
print(f"grid: n={uuids.shape[0]} uuids{uuids.shape} cents{cents.shape} rings{rings.shape}")
assert uuids.shape == (1247, 16), f"grid count {uuids.shape[0]} != 1247"
assert rings.shape == (1247, 6, 2)
# densified grid -> list of (19,2) rings
_, _, rings_d1 = h9.grid(-0.3, 51.3, 0.1, 51.7, 8, densify=1)
assert isinstance(rings_d1, list) and rings_d1[0].shape == (19, 2)
print(f"grid d1: list of {len(rings_d1)} rings, first {rings_d1[0].shape}")

# adaptive digest: conservation + assignment shape (input = FULL uuids —
# addresses, not coordinates; bins are rejected)
au, al, av, an, aa = h9.adaptive(uu[:2000], 3, 8, 50.0, 5.0)
assert au.shape[1] == 16 and av.sum() == 2000.0 and an.sum() == 2000
assert aa.shape == (2000,) and aa.min() >= 0 and aa.max() < au.shape[0]
print(f"adaptive: {au.shape[0]} cells, layers {sorted(set(al.tolist()))}, sum={av.sum():.0f}")

print("PY SMOKE OK")
