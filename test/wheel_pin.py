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

# (lon, lat, full-depth uuid hex, curve uuid hex) — includes s005/s030, the
# two points that historically flipped under FMA-contraction variance.
PINS = [
    (-62.39632882728877, -15.365651509427153,
     "b6707406817251531762571256488134", "c6750051227032726616543232058462"),
    (-111.72858692570949, -83.34481569615366,
     "92850838378016070631204602706585", "c7542844518408200782438382578143"),
    (139.12172450959264, -4.4624304013754319,
     "13478410802038168037034263288802", "ca768411048602461681067147744668"),
]

lon = np.array([p[0] for p in PINS])
lat = np.array([p[1] for p in PINS])
want = [p[2] for p in PINS]
want_curve = [p[3] for p in PINS]

# 1. full-depth universality — cell uuids, curve uuids, and curve labels
#    (a full-depth curve label IS the curve uuid's hex spelling)
uuids = hex9.encode(lon, lat)
got = [row.tobytes().hex() for row in uuids]
assert got == want, f"UNIVERSALITY VIOLATION:\n got {got}\nwant {want}"
curves = hex9.curve(uuids)
got_c = [row.tobytes().hex() for row in curves]
assert got_c == want_curve, f"UNIVERSALITY VIOLATION (curve):\n got {got_c}\nwant {want_curve}"
for i, wc in enumerate(want_curve):
    assert hex9.curve_label(uuids[i]) == wc, "curve_label != pinned curve"
    assert (hex9.curve_parse_label(wc) == curves[i]).all(), "curve_parse_label roundtrip"
assert hex9.curve_label(uuids) == want_curve, "batch curve_label != pins"
assert hex9.label(uuids, 8) == [hex9.label(u, 8) for u in uuids], "batch label != per-row"

# 1b. the curve-marker guard: cell labels of curve addresses are meaningless,
#     and must raise rather than mint plausible-looking garbage.
try:
    hex9.label(curves[0], 8)
    raise SystemExit("marker guard FAILED: label() accepted a curve-uuid")
except RuntimeError:
    pass

# 1c. E4H pins — the exact aperture-4 tail classifier, verbatim rows from
#     test_data/e4h_pin.tsv (hhg9-minted), probing to the full nibble budget
#     (layer 0, depth 28). Byte-exact or it does not ship.
E4H_PINS = [
    (-106.20243356509908, -69.900329880850606, 6, 2,
     "9856124e012ffffffffffffffffffff1", "9856124E012"),
    (-44.700874355144002, -37.799880486031171, 8, 4,
     "a02518026e01333ffffffffffffffff0", "a02518026E01333"),
    (81.390053692695759, 21.713188322189765, 2, 16,
     "544e00110512512055001ffffffffff0", "544E00110512512055001"),
    (-157.43716898579854, -30.803097692654944, 0, 28,
     "3e110015540001004050555110544410", "3E11001554000100405055511054441"),
]
for elon, elat, elayer, edepth, ewant, elabel in E4H_PINS:
    eu = hex9.e4h_encode(np.array([elon]), np.array([elat]),
                         layer=elayer, depth=edepth)
    got_e = eu[0].tobytes().hex()
    assert got_e == ewant, \
        f"UNIVERSALITY VIOLATION (e4h L{elayer}D{edepth}):\n got {got_e}\nwant {ewant}"
    assert hex9.e4h_label(eu[0]) == elabel, "e4h_label != pinned label"
    assert (hex9.e4h_parse_label(elabel) == eu[0]).all(), "e4h label roundtrip"
    assert hex9.is_e4h(eu[0]) and not hex9.is_e4h(uuids[0]), "is_e4h"
    assert int(hex9.e4h_depth(eu)[0]) == edepth, "e4h_depth"
    # truncation = binning: the tail is suffix-local
    b = hex9.e4h_bin(eu, max(0, edepth - 1))
    eb = hex9.e4h_encode(np.array([elon]), np.array([elat]),
                         layer=elayer, depth=max(0, edepth - 1))
    assert (b == eb).all(), "e4h truncation != binning"

# 1c-bis. hexagon binning (the GIS surface): the canonical hexagon key is
#     the matched pair's MODE-0 member — pinned byte-exactly per E4H pin
#     (the depth-28 pin is mode-1, so its key is its partner's address).
E4H_HEX = [
    "9856124e012ffffffffffffffffffff1",
    "a02518026e01333ffffffffffffffff0",
    "544e00110512512055001ffffffffff0",
    "3e110015540001004050555110544510",
]
for (elon, elat, elayer, edepth, ewant, _), ehex in zip(E4H_PINS, E4H_HEX):
    eu = hex9.e4h_encode(np.array([elon]), np.array([elat]),
                         layer=elayer, depth=edepth)
    hx = hex9.e4h_hex(eu)
    assert hx[0].tobytes().hex() == ehex, \
        f"UNIVERSALITY VIOLATION (e4h_hex L{elayer}D{edepth})"
    assert (hex9.e4h_hex(hx) == hx).all(), "e4h_hex not idempotent"
    assert int(hex9.e4h_mode(hx)[0]) == 0, "hex key not mode-0"
    plon, plat = hex9.e4h_partner(eu)
    pv = hex9.e4h_encode(np.array([plon[0]]), np.array([plat[0]]),
                         layer=elayer, depth=edepth)
    assert (hex9.e4h_hex(pv) == hx).all(), "pair members disagree on hex key"

# 1d. the E4H marker guard: h9 machinery must reject 0xE-marked uuids.
_e4h_u = hex9.e4h_encode(np.array([0.5]), np.array([0.5]), layer=6, depth=2)
try:
    hex9.decode(_e4h_u)
    raise SystemExit("marker guard FAILED: decode() accepted an E4H uuid")
except RuntimeError:
    pass
try:
    hex9.bin(_e4h_u, 4)
    raise SystemExit("marker guard FAILED: bin() accepted an E4H uuid")
except RuntimeError:
    pass

# 2. the boct seam reproduces the canonical chain bit-for-bit
cx, cy, oid = hex9.project(lon, lat)
assert (hex9.encode_boct(cx, cy, oid) == uuids).all(), "boct seam != encode"

# 3. round-trip sanity (decode is the cell centroid: sub-micron at L30)
rlon, rlat = hex9.decode(uuids)
assert np.max(np.abs(rlat - lat)) < 1e-6, "decode round-trip out of tolerance"

# 4. package identity
assert hex9.__version__ == hex9.version().split()[1]

# 5. the shipped runtime selftest (hex9.selftest) reproduces its own pins —
#    the embedder-facing check ships working on every wheel platform.
#    (9 h9/curve pins + 1 full-budget e4h pin)
assert hex9.selftest() == 10

print(f"hex9 {hex9.__version__}: universality pins OK "
      f"({len(PINS)} full-depth, boct seam bit-exact)")
