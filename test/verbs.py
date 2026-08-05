"""verbs.py — behavioural checks for the cell-first grid verbs.

Runs against the installed wheel (`hex9.verbs`) or the dev tree
(PYTHONPATH=build, flat `hex9_ext` module + verbs.py loaded by path —
the module's own dev fallback resolves the core). The field straddles
the equator octant seam on purpose: the seam must change nothing.

The checks are behavioural, not pinned: bearings are advisory FP, and
the addresses the verbs return are minted by the (separately pinned)
canonical chain — so we assert grid laws (membership, adjacency,
occlusion, detour), not byte goldens.
"""
import sys

import numpy as np

try:
    from hex9 import verbs as V
    from hex9 import _core as H
except ImportError:                    # dev tree: flat module, load by path
    import importlib.util
    import pathlib
    import hex9_ext as H
    _p = pathlib.Path(__file__).resolve().parent.parent / 'python' / 'hex9' / 'verbs.py'
    _spec = importlib.util.spec_from_file_location('verbs', _p)
    V = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(V)

LAYER = 8
FOX = (35.02, -0.02)                   # (lon, lat) — just south of the seam
BEARING, HALF, K = 15.0, 32.0, 6

fox_key = V.encode_keys(*FOX, LAYER)[0]
# self-calibrate the field to the layer's cell size, like the verbs do
pitch = V._pitch_of(fox_key, LAYER)
assert 2 * pitch < abs(FOX[1]) < (K - 1) * pitch, "fox must straddle-range the seam"

# 1. bind/centroid round-trip: a key's centroid re-binds to the same key
lon = np.array([34.7, 35.0, 35.3, FOX[0]])
lat = np.array([-0.3, 0.0, 0.35, FOX[1]])
keys = V.encode_keys(lon, lat, LAYER)
clon, clat = V.centroids(keys)
assert V.encode_keys(clon, clat, LAYER) == keys, "centroid re-bind moved a key"

# 2. multi-layer bind agrees with per-layer binds
multi = V.encode_keys_multi(lon, lat, (6, LAYER))
assert multi[LAYER] == keys
assert multi[6] == V.encode_keys(lon, lat, 6)

# 3. aim is the inverse of a neighbour's bearing, for every neighbour
nbs, cnt = H.neighbors(V._full([fox_key]), LAYER)
nb = np.asarray(nbs, dtype=np.uint8)[0, :int(np.asarray(cnt).reshape(-1)[0])]
flon, flat_ = (float(v[0]) for v in V.centroids([fox_key]))
nlon, nlat = (np.asarray(a) for a in H.decode(nb))
for i in range(nb.shape[0]):
    b = float(V.gc_bearing(flon, flat_, nlon[i], nlat[i]))
    assert V.aim(fox_key, LAYER, b) == bytes(nb[i]), "aim missed its neighbour"

# 4. batched aim == scalar aim, including one shared bearing
bearings = np.array([0.0, 60.0, 120.0, 275.0])
batch = V.aim(keys, LAYER, bearings)
assert batch == [V.aim(k, LAYER, float(b)) for k, b in zip(keys, bearings)]
assert V.aim(keys, LAYER, 15.0) == [V.aim(k, LAYER, 15.0) for k in keys]

# 5. vision_cone: source always visible; cone ⊆ disk; full circle = disk;
#    every member is inside the cone by the verb's own bearing law
disk_keys = {bytes(r) for r in np.asarray(
    H.k_disk(V._full([fox_key])[0], LAYER, K), dtype=np.uint8).reshape(-1, 16)}
cone = V.vision_cone(fox_key, LAYER, BEARING, HALF, K)
assert fox_key in cone
assert cone <= disk_keys
assert V.vision_cone(fox_key, LAYER, 0.0, 180.0, K) == disk_keys
mlon, mlat = V.centroids(sorted(cone - {fox_key}))
mbr = V.gc_bearing(flon, flat_, mlon, mlat)
assert np.all(np.abs(V.wrap180(mbr - BEARING)) <= HALF), "cell outside cone"

# 6. the cone crosses the octant seam (fox south of the equator, looking north)
assert np.any(mlat > 0) and np.any(mlat < 0), "cone did not cross the seam"

# 7. occlusion: a hedge across the cone shadows cells but stays visible itself
#    (laid perpendicular to the bearing at mid-cone, sized from the pitch)
t = np.linspace(-1.6, 1.6, 40) * pitch
b = np.radians(BEARING)
cx = FOX[0] + 3.2 * pitch * np.sin(b)
cy = FOX[1] + 3.2 * pitch * np.cos(b)
hedge = set(V.encode_keys(cx + t * np.cos(b), cy - t * np.sin(b), LAYER))
hedge.discard(fox_key)
seen = V.vision_cone(fox_key, LAYER, BEARING, HALF, K, hedge)
assert seen <= cone, "obstacles revealed cells"
shadow = cone - seen
assert shadow, "hedge cast no shadow"
# an obstacle never occludes itself (you can see the hedge, not through
# it) — its front face is visible; only cells BEHIND the nearest hedge
# cell may fall in shadow (hedge cells can shadow each other end-on)
assert hedge & seen, "the hedge's front face must be visible"
hlon, hlat = V.centroids(sorted(hedge & cone))
slon2, slat2 = V.centroids(sorted(shadow))
assert (np.min(V.gc_angle(flon, flat_, slon2, slat2))
        > np.min(V.gc_angle(flon, flat_, hlon, hlat))), \
    "a shadow cell sits nearer than the whole hedge"

# 8. walk_to: identity, single step, adjacency, obstacle detour
assert V.walk_to(fox_key, fox_key, LAYER) == [fox_key]
step = V.aim(fox_key, LAYER, BEARING)
assert V.walk_to(fox_key, step, LAYER) == [fox_key, step]

far = sorted(shadow)[0] if shadow else sorted(cone)[-1]
path = V.walk_to(fox_key, far, LAYER)
assert path is not None and path[0] == fox_key and path[-1] == far
for a, b in zip(path, path[1:]):
    nbs, cnt = H.neighbors(V._full([a]), LAYER)
    row = np.asarray(nbs, dtype=np.uint8)[0, :int(np.asarray(cnt).reshape(-1)[0])]
    assert b in {bytes(r) for r in row}, "path step is not a neighbour"

target = sorted(k for k in shadow if k not in hedge)
if target:
    detour = V.walk_to(fox_key, target[0], LAYER, hedge)
    direct = V.walk_to(fox_key, target[0], LAYER)
    assert detour is not None and not (set(detour) - {target[0]}) & hedge
    assert len(detour) >= len(direct), "detour shorter than direct path"

# 9. true-bearing correction: roundtrip, invariants, known magnitude
bs = np.array([0.0, 30.0, 45.0, 90.0, 135.0, 180.0, 250.0, 359.0])
ls = np.array([-60.0, -0.02, 0.0, 33.0, 55.95, 89.0])
for l0 in ls:
    v = V.from_true_bearing(bs, l0)
    assert np.allclose(V.to_true_bearing(v, l0), bs, atol=1e-9), "true-bearing roundtrip"
assert np.allclose(V.from_true_bearing(np.array([0.0, 90.0, 180.0, 270.0]), 45.0),
                   [0.0, 90.0, 180.0, 270.0], atol=1e-12), "cardinals must be identity"
assert np.allclose(V.from_true_bearing(bs, 90.0), bs, atol=1e-9), "poles must be identity"
dev = 45.0 - float(V.from_true_bearing(45.0, 0.0))
assert abs(dev - 0.1918) < 1e-3, f"flattening magnitude off: {dev}"

# 10. a destination sealed inside obstacles is unreachable
ring = {bytes(r) for r in np.asarray(
    H.k_ring(V._full([far])[0], LAYER, 1), dtype=np.uint8).reshape(-1, 16)}
assert V.walk_to(fox_key, far, LAYER, ring, max_expand=400) is None

print(f"verbs: OK (layer {LAYER}, cone {len(cone)} cells, "
      f"{len(shadow)} shadowed, walk {len(path) - 1} steps over the seam)")
sys.exit(0)
