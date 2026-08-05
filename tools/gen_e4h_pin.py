# Part of the Hex9 (H9) Project
# Copyright ©2026, Ben Griffin
# Licensed under the Apache License, Version 2.0

"""gen_e4h_pin.py — generate test_data/e4h_pin.tsv from the hhg9 reference.

The corpus is the E4H conformance interchange: every row's uuid was minted
by hhg9/h9/e4h.py (the normative reference), and libhex9's exact classifier
must reproduce it BYTE-IDENTICALLY (test/e4h_parity.c); decode and partner
representatives must agree to 1e-9 degrees.

Sampling: area-uniform global points across all 8 octants, the equator and
meridian seam bands, cone-point rings (the 6 octahedron vertices), and
near-pole points — per (layer, depth) regime, with depths probing far past
the CSP-verified range (Ben's ruling 2026-08-05: pin > 4), up to the full
nibble budget (layer 0, depth 28).

Knife-edge filtering: a point whose classification margin (best minus
second-best signed-distance score, in the current residual frame) at ANY
level is below MARGIN is dropped — on a cut line either answer is a valid
address, and the frozen C program (snap 2^-46) and the reference (double
descent) may legitimately differ there. Same doctrine as the reference's
own 856-point census vs the geometric PoC.

Run:  python tools/gen_e4h_pin.py [path-to-hhg9-repo]
"""
import math
import os
import subprocess
import sys

import numpy as np

HHG9 = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/Documents/Projects/PyCharm/hex9')
sys.path.insert(0, HHG9)

from hhg9 import Points, Registrar                     # noqa: E402
from hhg9.h9 import e4h                                # noqa: E402
from hhg9.h9.uuid_address import h9_bin_pts            # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'test_data', 'e4h_pin.tsv')
MARGIN = 1e-5
REGIMES = [                    # (layer, depth) — depth probes to full budget
    (6, 2), (8, 4), (5, 9), (2, 16), (10, 18), (0, 28),
]

rng = np.random.default_rng(20260805)


def sample_points():
    """(name, lat, lon) tuples covering octants, seams, cones, poles."""
    pts = []
    n = 40
    lon = rng.uniform(-180.0, 180.0, n)
    lat = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n)))
    pts += [('rand%02d' % i, lat[i], lon[i]) for i in range(n)]
    for j, (name, la0, lo0) in enumerate(
            [('eq', 0.0, None), ('m0', None, 0.0), ('m90', None, 90.0),
             ('m180', None, 180.0)]):
        for i in range(5):
            d = rng.uniform(-0.2, 0.2)
            lo = rng.uniform(-180, 180) if lo0 is None else lo0 + d
            la = rng.uniform(-60, 60) if la0 is None else la0 + d
            pts.append(('seam-%s%02d' % (name, i), la, lo))
    cones = [(90.0, 0.0), (-90.0, 0.0), (0.0, 0.0), (0.0, 90.0),
             (0.0, 180.0), (0.0, -90.0)]
    for j, (cla, clo) in enumerate(cones):
        for i in range(4):
            th = rng.uniform(0, 2 * math.pi)
            r = rng.uniform(0.05, 0.5)
            la = max(-89.99, min(89.99, cla + r * math.cos(th)))
            lo = clo + (r * math.sin(th) / max(0.05, math.cos(math.radians(la))))
            pts.append(('cone%d-%02d' % (j, i), la, lo))
    return pts


def classify_margin(w, cand):
    us = [(w - b) / a for a, b in cand]
    ss = [e4h._score(u) for u in us]
    k = int(np.argmax(ss))
    srt = sorted(ss, reverse=True)
    return k, us[k], srt[0] - srt[1]


def min_margin(lat, lon, layer, depth, reg):
    """Replicate the reference descent, tracking the classify margins.
    Returns None where the reference raises (two seams / degenerate)."""
    g_gcd, b_oct = reg.domain('g_gcd'), reg.domain('b_oct')
    bp = reg.project(Points(np.array([[lat, lon]]), g_gcd), [g_gcd, b_oct])
    p = bp.coords[0, :2].copy()
    g = int(np.asarray(bp.oid)[0])
    host = h9_bin_pts(bp, layer)[0]
    info = e4h._host_info(host, b_oct)
    fr = e4h._host_frame(info)
    _c, o, _c2, _mo, _lay = info
    if g != o:
        m = e4h._unfolds().get((o, g))
        if m is None:
            return None
        p = np.array([p[0], p[1], 1.0]) @ m
    w = e4h._fwd(fr, complex(p[0], p[1]))
    mm = math.inf
    half, w, mg = classify_margin(w, e4h._HALVES)
    mm = min(mm, mg)
    for _ in range(depth):
        _k, w, mg = classify_margin(w, e4h._MAPS)
        mm = min(mm, mg)
    return mm


def main():
    reg = Registrar()
    rows = []
    dropped = 0
    for (layer, depth) in REGIMES:
        for name, la, lo in sample_points():
            mm = min_margin(la, lo, layer, depth, reg)
            if mm is None or mm < MARGIN:
                dropped += 1
                continue
            u = e4h.h9e_encode([la], [lo], layer=layer, depth=depth, reg=reg)[0]
            lab = e4h.h9e_label(u)
            dla, dlo = e4h.h9e_decode([u], reg=reg)
            pla, plo = e4h.h9e_partner_point([u], reg=reg)
            rows.append('\t'.join([
                'L%dD%d-%s' % (layer, depth, name),
                '%.17g' % lo, '%.17g' % la, str(layer), str(depth),
                '%032x' % u.int, lab,
                '%.17g' % dlo[0], '%.17g' % dla[0],
                '%.17g' % plo[0], '%.17g' % pla[0],
            ]))
    try:
        commit = subprocess.check_output(
            ['git', '-C', HHG9, 'rev-parse', '--short', 'HEAD'],
            text=True).strip()
    except Exception:
        commit = 'unknown'
    with open(OUT, 'w') as f:
        f.write('# e4h_pin.tsv — E4H conformance corpus, generated by '
                'tools/gen_e4h_pin.py\n')
        f.write('# reference: hhg9 %s (hhg9/h9/e4h.py); margin filter %g; '
                'regimes %s\n' % (commit, MARGIN, REGIMES))
        f.write('# name\tlon\tlat\tlayer\tdepth\tuuid\tlabel\tdec_lon\t'
                'dec_lat\tpar_lon\tpar_lat\n')
        f.write('\n'.join(rows) + '\n')
    print('wrote %s: %d rows (%d dropped by margin/seam filter)'
          % (OUT, len(rows), dropped))


if __name__ == '__main__':
    main()
