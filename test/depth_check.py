"""Validate a beam-search DEPTH change: encode a fixed sample (+ boundary stress
points) and compare UUIDs against a saved baseline. Zero diffs = the change is
output-preserving. Run once to save the baseline, again after the change.

  PYTHONPATH=build python3 test/depth_check.py        # save baseline / compare
"""
import os, sys, time
import numpy as np
import hex9_ext as h9

h9.warp_init()
n = 500_000
rng = np.random.default_rng(42)
lon = rng.uniform(-180.0, 180.0, n)
lat = rng.uniform(-89.9, 89.9, n)

# boundary / pole / seam stress points
b_lon = np.array([0., 90., 180., -90., -179.999, 0., 45., -3.19, 13.405, 0.0001])
b_lat = np.array([0.,  0.,  0.,   0.,    0.,    89.999, -89.999, 55.95, 52.52, -0.0001])
lon = np.concatenate([lon, b_lon])
lat = np.concatenate([lat, b_lat])

t = time.perf_counter()
uu = h9.encode(lon, lat)
dt = time.perf_counter() - t
print(f"encoded {len(lon)} pts in {dt:.3f}s ({len(lon)/dt:,.0f} pts/s)")

ref = "/tmp/h9_depth_baseline.npy"
if not os.path.exists(ref):
    np.save(ref, uu)
    print("SAVED baseline (DEPTH=40). Now change DEPTH, rebuild, rerun.")
    sys.exit(0)

base = np.load(ref)
diffs = int(np.count_nonzero(np.any(uu != base, axis=1)))
print(f"diffs vs baseline: {diffs} / {len(lon)}")
sys.exit(0 if diffs == 0 else 2)
