# Part of the Hex9 (H9) Project
# Copyright ©2026, Ben Griffin
# Licensed under the Apache License, Version 2.0

"""gen_e4h_tables.py — freeze the E4H constant tables into core/h9_e4h_tables.h.

The E4H reference implementation is hhg9/h9/e4h.py (the symbolic/structural
form, 2026-08-05). This generator imports it and emits:

  * EXACT integer tables for the descent classifier in Z[1/2, sqrt3]:
    the rep-4 child maps as u = M_k (w - B_k) with M_k = 4*conj(a_k)
    (integer components — exactly (w - b_k)/a_k since |a_k|^2 = 1/4),
    the trapezoid edges for the signed-distance score, and the digit
    algebra (ROT, piece inversion, BASES). Doubles from the reference
    are matched to (A + B*sqrt3)/2^k form and asserted exact; the C
    kernel then runs these as integers, so classification decisions
    can never flip on any platform at any depth.
  * FP seed tables (frozen as hex-float literals): the state frames
    A[p_mo][c2], the per-(mode, corner) seam unfold/fold affines, the
    sv triangles for the in-face test, and the decode probes. These
    feed the det-math FP seed (project → host frame → snap); after the
    snap at 2^-E4H_SIGMA everything is integer.
  * A rigorous worst-case bit bound for the descent components
    (per-level infinity-norm growth), emitted as a static_assert.

Regenerating against a different hhg9 commit and getting a different
header is a REGIME CHANGE for E4H addresses — see docs/universality.md.

Run:  python tools/gen_e4h_tables.py [path-to-hhg9-repo]
"""
import math
import os
import subprocess
import sys

import numpy as np

HHG9 = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/Documents/Projects/PyCharm/hex9')
sys.path.insert(0, HHG9)

from hhg9.h9 import e4h                                     # noqa: E402
from hhg9.h9.polygon import H9P, SV_CORNER_AXIS, fold_to_octant  # noqa: E402
from hhg9.h9.constants import H9O                           # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'core', 'h9_e4h_tables.h')
SIGMA = 46          # snap scale: w -> round(w * 2^SIGMA); part of the address
S3 = math.sqrt(3.0)


def exact(v, kmax=4, bound=16):
    """Match a double v to (A + B*sqrt3)/2^k exactly (A, B, k small ints)."""
    for k in range(kmax + 1):
        s = v * (1 << k)
        for b in range(-bound, bound + 1):
            a = s - b * S3
            if abs(a - round(a)) < 1e-9:
                A = int(round(a))
                if abs((A + b * S3) / (1 << k) - v) < 1e-12:
                    return A, b, k
    raise SystemExit(f'no exact Z[1/2,sqrt3] form for {v!r}')


def cd(x):
    """C hex-float literal for a double (bit-exact)."""
    return float(x).hex()


# ── pull the reference constants ─────────────────────────────────────
A_FRAMES = e4h._state_frames()               # (2,3) complex128
UNFOLDS = e4h._unfolds()                     # {(o,g): 3x2}
MAPS, HALVES, EDGES, Z0 = e4h._MAPS, e4h._HALVES, e4h._EDGES, e4h._Z0
CEN = e4h._CEN
ROT = e4h._ROT
BASES = e4h._bases()
assert tuple(ROT) == (0, 2, 3, 4), 'child rotations moved'
assert e4h.E_MARK == 0xE

# digit-rule injectivity over every context (the CSP's runtime shadow)
for oid in range(8):
    for c2 in range(3):
        for half in range(2):
            e4h._digit_to_class(oid, c2, half)

# piece inversion: rot-mod-3 residue -> piece index (reference _ROT3_PIECE)
P3 = [0, 0, 0]
for r, k in e4h._ROT3_PIECE.items():
    P3[r] = k

# ── exact integer forms for the descent ──────────────────────────────
# M_k = 4*conj(a_k): components must be integers in Z[sqrt3].
M_INT = []       # (mrA, mrB, miA, miB): real = mrA + mrB*sqrt3, imag = ...
B_INT = []       # (brA, brB, biA, biB) at fixed scale /4
for a, b in MAPS:
    arA, arB, ark = exact(a.real)
    aiA, aiB, aik = exact(a.imag)
    mr = (4 * arA // (1 << ark), 4 * arB // (1 << ark))
    mi = (4 * aiA // (1 << aik), 4 * aiB // (1 << aik))
    assert (mr[0] * (1 << ark) == 4 * arA and mr[1] * (1 << ark) == 4 * arB
            and mi[0] * (1 << aik) == 4 * aiA and mi[1] * (1 << aik) == 4 * aiB)
    # conj: negate imaginary part
    M_INT.append((mr[0], mr[1], -mi[0], -mi[1]))
    brA, brB, brk = exact(b.real)
    biA, biB, bik = exact(b.imag)
    B_INT.append((brA << (2 - brk), brB << (2 - brk),
                  biA << (2 - bik), biB << (2 - bik)))
    assert brk <= 2 and bik <= 2

# verify M_k(w - b_k) == (w - b_k)/a_k: |a|^2 must be exactly 1/4
for (a, _b), m in zip(MAPS, M_INT):
    v = complex(m[0] + m[1] * S3, m[2] + m[3] * S3)
    assert abs(v * a - 1.0) < 1e-12, '4*conj(a)*a != 1'

# HALVES are +1 / -1
H_SIGN = []
for a, b in HALVES:
    assert b == 0 and a.imag == 0 and abs(abs(a.real) - 1.0) < 1e-15
    H_SIGN.append(int(a.real))
assert H_SIGN == [1, -1]

# EDGES: (corner p, unit direction e) — components at denominator 2.
# Emit 2p and 2e as integer Z[sqrt3] pairs; the score compares
# Im(conj(e)(w - p)) and both scalings are common across candidates.
E_INT = []
for p, e in EDGES:
    prA, prB, prk = exact(p.real); piA, piB, pik = exact(p.imag)
    erA, erB, erk = exact(e.real); eiA, eiB, eik = exact(e.imag)
    assert max(prk, pik, erk, eik) <= 1
    E_INT.append(((prA << (1 - prk), prB << (1 - prk),
                   piA << (1 - pik), piB << (1 - pik)),
                  (erA << (1 - erk), erB << (1 - erk),
                   eiA << (1 - eik), eiB << (1 - eik))))

# ── worst-case component growth bound ────────────────────────────────
# Descent state (xr, xs, yr, ys) at fixed scale 2^-SIGMA. Per level:
# subtract B (|B*2^(SIGMA-2)| bounded by 4*2^SIGMA), then multiply by
# M_k, whose action on the 4-vector has infinity-norm:
def m_norm(m):
    mrA, mrB, miA, miB = m
    #  xr' = mrA*xr + 3*mrB*xs - miA*yr - 3*miB*ys   (and cyclic)
    rows = [(abs(mrA), 3 * abs(mrB), abs(miA), 3 * abs(miB)),
            (abs(mrB), abs(mrA), abs(miB), abs(miA)),
            (abs(miA), 3 * abs(miB), abs(mrA), 3 * abs(mrB)),
            (abs(miB), abs(miA), abs(mrB), abs(mrA))]
    return max(sum(r) for r in rows)


G = max(m_norm(m) for m in M_INT)
MAX_DEPTH = 28                       # deepest tail (layer 0 host)
# Bound over the RUNTIME GUARD BOX (h9e4h_descend admits |w| <= 64, far
# beyond any real seed): seed magnitude <= 64*2^SIGMA (SIGMA+6 bits), then
# per-level growth G, then the score stage (2w, minus P2*2^SIGMA, ring-
# multiplied by E2 with row norm 8) adds the +5 slack.
W0_BITS = SIGMA + 6
g_bits = math.log2(G)
BITS = math.ceil(W0_BITS + MAX_DEPTH * g_bits) + 5
assert BITS <= 120, f'descent bound {BITS} bits exceeds __int128 comfort'

# ── seam tables: per (mode, corner k) unfold + fold affines ──────────
# fold_to_octant's isometry depends only on (mode, k); assert that the
# reference's per-(o,g) unfolds collapse accordingly, then freeze [2][3].
UNF = np.zeros((2, 3, 3, 2))
FLD = np.zeros((2, 3, 3, 2))
AXK = np.zeros((2, 3), dtype=int)
SV = np.array([H9P.sv[0], H9P.sv[1]])
for mode in range(2):
    tri = H9P.sv[mode]
    ax = SV_CORNER_AXIS[mode]
    tri_b = H9P.sv[1 - mode]
    axis_b = SV_CORNER_AXIS[1 - mode]
    for k in range(3):
        AXK[mode, k] = ax[k]
        a, b = tri[(k + 1) % 3], tri[(k + 2) % 3]
        e = b - a
        e = e / np.hypot(e[0], e[1])
        d = tri[k] - a
        refl = a + 2.0 * float(d @ e) * e - d
        src = np.array([a, b, refl])
        axes = (ax[(k + 1) % 3], ax[(k + 2) % 3], ax[k])
        dst = np.array([tri_b[axis_b.index(x)] for x in axes])
        UNF[mode, k] = np.linalg.solve(np.column_stack([dst, np.ones(3)]), src)
        FLD[mode, k] = np.linalg.solve(np.column_stack([src, np.ones(3)]), dst)

# collapse check against the reference's (o, g) table
for (o, g), m in UNFOLDS.items():
    mode = int(H9O.oid_mo[o])
    axis = int(math.log2(o ^ g))
    k = SV_CORNER_AXIS[mode].index(axis)
    assert np.allclose(m, UNF[mode, k], atol=1e-12), f'unfold ({o},{g}) mismatch'

# unfold(fold(p)) == p spot check
rng = np.random.default_rng(0)
for mode in range(2):
    for k in range(3):
        p = rng.uniform(-0.3, 0.3, (5, 2))
        q = np.column_stack([p, np.ones(5)]) @ FLD[mode, k]
        r = np.column_stack([q, np.ones(5)]) @ UNF[mode, k]
        assert np.allclose(r, p, atol=1e-12)

# behavioural cross-check of the frozen fold vs the reference fold
for mode, oid in ((0, 0), (1, 1)):
    pts = rng.uniform(-0.9, 0.9, (200, 2))
    ref_xy, ref_oid = fold_to_octant(pts, oid)
    for i, v in enumerate(pts):
        o, w = oid, v.copy()
        for _ in range(3):
            m = int(H9O.oid_mo[o]); tri = SV[m]
            hit = None
            for kk in range(3):
                aa, bb = tri[(kk + 1) % 3], tri[(kk + 2) % 3]
                cross = ((bb[0] - aa[0]) * (w[1] - aa[1])
                         - (bb[1] - aa[1]) * (w[0] - aa[0]))
                if cross > 1e-12:
                    hit = kk
                    break
            if hit is None:
                break
            w = np.array([w[0], w[1], 1.0]) @ FLD[m, hit]
            o ^= 1 << AXK[m, hit]
        assert o == int(ref_oid[i]) and np.allclose(w, ref_xy[i], atol=1e-12)

# ── emit ─────────────────────────────────────────────────────────────
try:
    commit = subprocess.check_output(
        ['git', '-C', HHG9, 'rev-parse', '--short', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(
        ['git', '-C', HHG9, 'status', '--short', '--', 'hhg9/h9/e4h.py'],
        text=True).strip()
    if dirty:
        commit += '+dirty'
except Exception:
    commit = 'unknown'

L = []
L.append('/* Generated by tools/gen_e4h_tables.py from hhg9 %s — do not edit.'
         % commit)
L.append(' *')
L.append(' * E4H (aperture-4 structural tail) frozen tables. The integer')
L.append(' * tables define the EXACT classifier in Z[1/2, sqrt3]; the FP')
L.append(' * tables are the det-math seed. Regenerating with different')
L.append(' * values is a regime change for E4H addresses (universality.md).')
L.append(' */')
L.append('#pragma once')
L.append('#include <stdint.h>')
L.append('')
L.append('#define E4H_MARK   0xE')
L.append('#define E4H_SIGMA  %d   /* snap: W = llround(w * 2^SIGMA) */' % SIGMA)
L.append('#define E4H_BITS_BOUND %d /* worst-case component bits, rigorous */'
         % BITS)
L.append('#if E4H_BITS_BOUND > 120')
L.append('#error "E4H descent component bound exceeds __int128 comfort"')
L.append('#endif')
L.append('')
L.append('/* rep-4 child rotations (units of 60 deg) and the piece inversion')
L.append(' * P3[(class - s) mod 3] -> piece index (reference _ROT3_PIECE). */')
L.append('static const int E4H_ROT[4] = {%d, %d, %d, %d};' % tuple(ROT))
L.append('static const int E4H_P3[3] = {%d, %d, %d};' % tuple(P3))
L.append('')
L.append('/* five-symbol bases: BASES[oid][slot] (axis->digit gauge baked in) */')
L.append('static const int E4H_BASES[8][3] = {')
for o in range(8):
    L.append('    {%d, %d, %d},' % BASES[o])
L.append('};')
L.append('')
L.append('/* child maps as u = M_k (w - B_k): M integer Z[sqrt3] components')
L.append(' * (rA, rB, iA, iB) with value (rA + rB*s3) + i(iA + iB*s3);')
L.append(' * B at fixed scale /4. Exactly (w - b_k)/a_k. */')
L.append('static const int64_t E4H_M[4][4] = {')
for m in M_INT:
    L.append('    {%d, %d, %d, %d},' % m)
L.append('};')
L.append('static const int64_t E4H_B4[4][4] = {   /* value = B4/4 */')
for b in B_INT:
    L.append('    {%d, %d, %d, %d},' % b)
L.append('};')
L.append('static const int E4H_HALF_SIGN[2] = {1, -1};')
L.append('')
L.append('/* trapezoid edges for the signed-distance score: 2*p and 2*e as')
L.append(' * integer Z[sqrt3] quadruples (score sign/order is scale-free). */')
L.append('static const int64_t E4H_EDGE_P2[4][4] = {')
for p2, _e2 in E_INT:
    L.append('    {%d, %d, %d, %d},' % p2)
L.append('};')
L.append('static const int64_t E4H_EDGE_E2[4][4] = {')
for _p2, e2 in E_INT:
    L.append('    {%d, %d, %d, %d},' % e2)
L.append('};')
L.append('')
L.append('/* state frames A[p_mo][c2] (complex): w = A*3^layer*conj(z - c) */')
L.append('static const double E4H_A_RE[2][3] = {')
for mo in range(2):
    L.append('    {%s, %s, %s},' % tuple(cd(A_FRAMES[mo, c].real)
                                         for c in range(3)))
L.append('};')
L.append('static const double E4H_A_IM[2][3] = {')
for mo in range(2):
    L.append('    {%s, %s, %s},' % tuple(cd(A_FRAMES[mo, c].imag)
                                         for c in range(3)))
L.append('};')
L.append('')
L.append('/* decode probes: canonical trapezoid centroid and its partner */')
L.append('static const double E4H_CEN_RE = %s;   /* -2/(3*sqrt3) */'
         % cd(CEN.real))
L.append('static const double E4H_CEN_IM = %s;' % cd(CEN.imag))
L.append('')
L.append('/* octant seam algebra per (mode, corner k): sv triangles (CW),')
L.append(' * the crossed axis bit, and the unfold/fold affines ([x y 1] @ m).')
L.append(' * Unfold: neighbour-frame point -> this frame; fold: inverse. */')
L.append('static const double E4H_SV[2][3][2] = {')
for mode in range(2):
    L.append('    {{%s, %s}, {%s, %s}, {%s, %s}},'
             % tuple(cd(SV[mode, i, j]) for i in range(3) for j in range(2)))
L.append('};')
L.append('static const int E4H_AXK[2][3] = {{%d, %d, %d}, {%d, %d, %d}};'
         % (*AXK[0], *AXK[1]))
for name, tab in (('E4H_UNF', UNF), ('E4H_FLD', FLD)):
    L.append('static const double %s[2][3][3][2] = {' % name)
    for mode in range(2):
        L.append('  {')
        for k in range(3):
            L.append('    {{%s, %s}, {%s, %s}, {%s, %s}},'
                     % tuple(cd(tab[mode, k, i, j])
                             for i in range(3) for j in range(2)))
        L.append('  },')
    L.append('};')

with open(OUT, 'w') as f:
    f.write('\n'.join(L) + '\n')
print(f'wrote {OUT}  (hhg9 {commit}, sigma={SIGMA}, bits<={BITS})')
