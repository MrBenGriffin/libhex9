"""selftest — runtime universality pins for the installed hex9 wheel.

The wheel's build platforms are gated on test/wheel_pin.py before anything
ships, but the machine a wheel eventually *runs* on is not ours to test.
This module carries a subset of the frozen regime pins
(test_data/regime_pin.tsv, deterministic chain 2026-07-27) so an embedder
can verify their own environment conforms: the same lon/lat MUST mint the
bit-identical uuid here as everywhere else. A failure is a regime violation
or a nonconforming build/platform — never a tolerance issue. Doctrine:
libhex9 docs/universality.md.
"""
import numpy as np

from . import _core

# name, lon, lat, cell uuid, curve uuid — verbatim rows from regime_pin.tsv:
# landmarks, the poles, the antimeridian, an octant seam, and s005/s030 (the
# two equal-area sample points that historically flipped under FMA variance).
_PINS = (
    ("edinburgh", -3.1899999999999999, 55.950000000000003,
     "43217747878386156377686620685721", "c1142537043612886887584815842403"),
    ("null-island", 0.0, 0.0,
     "02222222222222222222222222222220", "c0333333333333333333333333333333"),
    ("north-pole", 0.0, 90.0,
     "42222222222222222222222222222222", "c1533333333333333333333333333333"),
    ("south-pole", 0.0, -90.0,
     "a2222222222222222222222222222222", "c7533333333333333333333333333333"),
    ("antimeridian-eq", 180.0, 0.0,
     "12222222222222222222222222222220", "ca355555555555555555555555555555"),
    ("prime-45n", 0.0, 45.0,
     "43333683683682222333683248342220", "c1165646104104533365284745776455"),
    ("seam-90e-45n", 90.0, 45.0,
     "53333323323342228333326282682224", "c2116560657655733416560734316533"),
    ("s005", -62.39632882728877, -15.365651509427153,
     "b6707406817251531762571256488134", "c6750051227032726616543232058462"),
    ("s030", -111.72858692570949, -83.34481569615366,
     "92850838378016070631204602706585", "c7542844518408200782438382578143"),
)


# lon, lat, layer, depth, e4h uuid — verbatim from test_data/e4h_pin.tsv
_E4H_PIN = (-157.43716898579854, -30.803097692654944, 0, 28,
            "3e110015540001004050555110544410")


def _first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i
    return -1


def selftest():
    """Re-mint the embedded regime pins on this machine.

    Returns the number of points verified. Raises RuntimeError if any minted
    address differs from the frozen goldens: this build or platform is then
    nonconforming (see libhex9 docs/universality.md) and addresses it mints
    must not be stored or exchanged.
    """
    lon = np.array([p[1] for p in _PINS])
    lat = np.array([p[2] for p in _PINS])

    uuids = _core.encode(lon, lat)
    curves = _core.curve(uuids)
    failures = []
    for i, (name, _, _, want_u, want_c) in enumerate(_PINS):
        got_u = uuids[i].tobytes().hex()
        got_c = curves[i].tobytes().hex()
        if got_u != want_u:
            failures.append(f"{name}: uuid {got_u} != {want_u} "
                            f"(first diff at nibble {_first_diff(got_u, want_u)})")
        elif got_c != want_c:
            failures.append(f"{name}: curve {got_c} != {want_c} "
                            f"(first diff at nibble {_first_diff(got_c, want_c)})")
    if failures:
        raise RuntimeError(
            "UNIVERSALITY VIOLATION — this environment does not reproduce "
            "the hex9 addressing regime; addresses minted here must not be "
            "stored (docs/universality.md):\n  " + "\n  ".join(failures))

    # The boct seam reproduces the canonical chain bit-for-bit.
    cx, cy, oid = _core.project(lon, lat)
    if not (_core.encode_boct(cx, cy, oid) == uuids).all():
        raise RuntimeError("UNIVERSALITY VIOLATION — boct seam != encode")

    # Round-trip sanity (decode is the cell centroid: sub-micron at L30);
    # poles excluded — longitude is degenerate there.
    rlon, rlat = _core.decode(uuids)
    finite = [i for i, p in enumerate(_PINS) if abs(p[2]) != 90.0]
    if max(abs(rlat[i] - lat[i]) for i in finite) > 1e-6:
        raise RuntimeError("decode round-trip out of tolerance")

    # One E4H pin (test_data/e4h_pin.tsv row L0D28-rand00): the exact
    # aperture-4 tail classifier at the FULL nibble budget — the deepest
    # address the format can hold, minted bit-identically or not at all.
    eu = _core.e4h_encode(np.array([_E4H_PIN[0]]), np.array([_E4H_PIN[1]]),
                          layer=_E4H_PIN[2], depth=_E4H_PIN[3])
    got_e = eu[0].tobytes().hex()
    if got_e != _E4H_PIN[4]:
        raise RuntimeError(
            "UNIVERSALITY VIOLATION — e4h tail classifier: "
            f"{got_e} != {_E4H_PIN[4]} "
            f"(first diff at nibble {_first_diff(got_e, _E4H_PIN[4])})")

    return len(_PINS) + 1
