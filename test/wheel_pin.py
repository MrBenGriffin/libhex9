"""wheel_pin.py — universality pins for published hex9 wheels.

Runs against the INSTALLED wheel (cibuildwheel test-command), on every
platform a wheel is built for, before anything ships. The pins are frozen
regime_pin goldens (test_data/regime_pin.tsv, deterministic chain
2026-07-27): the same lon/lat MUST mint the bit-identical uuid on every
platform. A failure is a regime violation or a broken build — never a
tolerance issue. Do not loosen; re-derive deliberately or fix the build.
"""
import numpy as np
import hex9

# (lon, lat, full-depth uuid hex) — includes s005/s030, the two points that
# historically flipped under FMA-contraction variance.
PINS = [
    (-62.39632882728877, -15.365651509427153, "b6707406817251531762571256488134"),
    (-111.72858692570949, -83.34481569615366, "92850838378016070631204602706585"),
    (139.12172450959264, -4.4624304013754319, "13478410802038168037034263288802"),
]

lon = np.array([p[0] for p in PINS])
lat = np.array([p[1] for p in PINS])
want = [p[2] for p in PINS]

# 1. full-depth universality
uuids = hex9.encode(lon, lat)
got = [row.tobytes().hex() for row in uuids]
assert got == want, f"UNIVERSALITY VIOLATION:\n got {got}\nwant {want}"

# 2. the boct seam reproduces the canonical chain bit-for-bit
cx, cy, oid = hex9.project(lon, lat)
assert (hex9.encode_boct(cx, cy, oid) == uuids).all(), "boct seam != encode"

# 3. round-trip sanity (decode is the cell centroid: sub-micron at L30)
rlon, rlat = hex9.decode(uuids)
assert np.max(np.abs(rlat - lat)) < 1e-6, "decode round-trip out of tolerance"

# 4. package identity
assert hex9.__version__ == hex9.version().split()[1]

print(f"hex9 {hex9.__version__}: universality pins OK "
      f"({len(PINS)} full-depth, boct seam bit-exact)")
