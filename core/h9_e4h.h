/*
 * h9_e4h.h — E4H aperture-4 structural-tail kernel (exact classifier).
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2026, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * Depends on h9_addressing.h (h9a_pack/h9a_unpack, H9_OID_MO) and
 * h9_layout.h; constants from h9_e4h_tables.h (generated from the hhg9
 * reference by tools/gen_e4h_tables.py — hhg9/h9/e4h.py is normative).
 *
 * An E4H address is a hex9 host bin truncated at its attach layer L,
 * an 0xE break marker, one HALF nibble (0/1), and up to (BODYTOP-L-2)
 * TAIL nibbles from {0, 1..5} descending the aperture-4 half-hex
 * rep-4 carrier:
 *
 *     nibbles: [h9 body 0..L] [0xE] [half] [d1..dB] [0xF pad] [key_tail]
 *     label  : <h9-label>E<half><d1..dB>
 *
 * THE EXACT CLASSIFIER. The seed is frozen det-math FP: project the
 * point to the host's canonical state frame (w = A[p_mo][c2] * 3^L *
 * conj(z - centre), seam points unfolded first), then SNAP once:
 * W = llround(w * 2^E4H_SIGMA). Everything after the snap is integer
 * arithmetic in Z[1/2, sqrt3] — each scalar is (A + B*sqrt3)/2^SIGMA
 * with A, B (signed) 128-bit — so half/child classification decisions
 * are exact functions of the snapped seed: they cannot flip across
 * platforms, optimisation levels, or depths. Child steps use
 * u = M_k (w - B_k) with M_k = 4*conj(a_k), which IS (w - b_k)/a_k
 * exactly (|a_k|^2 = 1/4), so the scale 2^-SIGMA never moves; the
 * component bit-growth is bounded by E4H_BITS_BOUND (generator-proven
 * <= 120). Sign of (A + B*sqrt3) with mixed signs is decided by the
 * exact A^2 vs 3B^2 comparison in 256-bit arithmetic.
 *
 * Decode composes the child similarities in doubles (geometry out is
 * advisory; addresses only ever come from the encode side).
 */
#pragma once

#include <stdint.h>
#include <math.h>
#include "h9_layout.h"
#include "h9_e4h_tables.h"

/* sqrt(3) as the correctly-rounded double (frozen literal) */
#define E4H_S3 0x1.bb67ae8584caap+0

/* error codes shared by the ABI wrappers */
#define E4H_EBAD   1   /* bad arguments / malformed input */
#define E4H_EGRAM  2   /* not an E4H address / grammar violation */
#define E4H_ESEAM  3   /* point two seams from host (census: unreachable
                        * for canonically binned points — reject) */

/* ── marker test ────────────────────────────────────────────────────
 * 0xE cannot occur in any h9 uuid (body digits <= 11, key_tail <= 0xD,
 * pad 0xF) nor curve uuid (0xC, slot <= 11, ranks <= 8, pad 0xF), so
 * any-nibble-0xE is a decisive, position-free E4H test. */
static int h9e4h_is_marked(const uint8_t uuid[16]) {
    for (int i = 0; i < 16; ++i) {
        const uint8_t b = uuid[i];
        if ((b >> 4) == E4H_MARK || (b & 0x0Fu) == E4H_MARK) return 1;
    }
    return 0;
}

/* ── grammar: split nibbles ─────────────────────────────────────────
 * nib_in: 32 unpacked nibbles. On success (0) writes the host bin
 * nibbles (marker onwards sentinelled, key_tail kept), the host layer,
 * half, tail digits and count. E4H_EGRAM on any violation. */
static int h9e4h_split_nibs(const uint8_t nib_in[32], uint8_t host_nib[32],
                            int *layer, int *half,
                            uint8_t *digits, int *ndig) {
    int e = -1;
    for (int i = 0; i <= H9_NIB_TAIL - 1; ++i) {
        if (nib_in[i] == E4H_MARK) {
            if (e >= 0) return E4H_EGRAM;      /* exactly one marker */
            e = i;
        }
    }
    if (e < 1 || e + 1 > H9_NIB_TAIL - 1) return E4H_EGRAM;
    const int hf = nib_in[e + 1];
    if (hf > 1) return E4H_EGRAM;
    int nd = 0, stop = 0;
    for (int i = e + 2; i <= H9_NIB_TAIL - 1; ++i) {
        const uint8_t d = nib_in[i];
        if (d == 0x0Fu) { stop = 1; continue; }
        if (stop || d > 5) return E4H_EGRAM;   /* digits then solid pad */
        digits[nd++] = d;
    }
    for (int i = 0; i < e; ++i)                /* host body sanity */
        if (nib_in[i] > 11) return E4H_EGRAM;
    for (int i = 0; i < 32; ++i) host_nib[i] = nib_in[i];
    for (int i = e; i <= H9_NIB_TAIL - 1; ++i) host_nib[i] = 0x0Fu;
    *layer = e - 1;
    *half = hf;
    *ndig = nd;
    return 0;
}

/* ── exact ring arithmetic: (A + B*sqrt3), signed 128-bit ─────────── */
typedef __int128 e4i;
typedef struct { e4i A, B; } e4r;
typedef struct { e4r re, im; } e4c;

static e4r e4r_sub(e4r a, e4r b) { e4r r = { a.A - b.A, a.B - b.B }; return r; }

/* ring multiply where one factor has SMALL (table) components */
static e4r e4r_muls(e4r a, int64_t cA, int64_t cB) {
    e4r r = { a.A * cA + 3 * (a.B * cB), a.A * cB + a.B * cA };
    return r;
}

static void h9e4h_acc256(uint64_t o[4], int idx, uint64_t v) {
    unsigned __int128 t = (unsigned __int128)o[idx] + v;
    o[idx] = (uint64_t)t;
    uint64_t carry = (uint64_t)(t >> 64);
    for (int i = idx + 1; carry && i < 4; ++i) {
        t = (unsigned __int128)o[i] + carry;
        o[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
}

static void h9e4h_sq256(unsigned __int128 v, uint64_t o[4]) {
    const uint64_t lo = (uint64_t)v, hi = (uint64_t)(v >> 64);
    const unsigned __int128 ll = (unsigned __int128)lo * lo;
    const unsigned __int128 lh = (unsigned __int128)lo * hi;
    const unsigned __int128 hh = (unsigned __int128)hi * hi;
    o[0] = (uint64_t)ll; o[1] = (uint64_t)(ll >> 64); o[2] = 0; o[3] = 0;
    h9e4h_acc256(o, 1, (uint64_t)lh); h9e4h_acc256(o, 2, (uint64_t)(lh >> 64));
    h9e4h_acc256(o, 1, (uint64_t)lh); h9e4h_acc256(o, 2, (uint64_t)(lh >> 64));
    h9e4h_acc256(o, 2, (uint64_t)hh); h9e4h_acc256(o, 3, (uint64_t)(hh >> 64));
}

static void h9e4h_mul3_256(uint64_t o[4]) {
    uint64_t carry = 0;
    for (int i = 0; i < 4; ++i) {
        const unsigned __int128 t = (unsigned __int128)o[i] * 3u + carry;
        o[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
}

static int h9e4h_cmp256(const uint64_t a[4], const uint64_t b[4]) {
    for (int i = 3; i >= 0; --i)
        if (a[i] != b[i]) return a[i] > b[i] ? 1 : -1;
    return 0;
}

/* exact sign of A + B*sqrt3 */
static int e4r_sgn(e4r v) {
    const int sA = (v.A > 0) - (v.A < 0);
    const int sB = (v.B > 0) - (v.B < 0);
    if (sA == 0) return sB;
    if (sB == 0 || sA == sB) return sA;
    /* opposite signs: decided by A^2 vs 3B^2, exactly */
    const unsigned __int128 ua =
        v.A < 0 ? (unsigned __int128)(-v.A) : (unsigned __int128)v.A;
    const unsigned __int128 ub =
        v.B < 0 ? (unsigned __int128)(-v.B) : (unsigned __int128)v.B;
    uint64_t a2[4], b2[4];
    h9e4h_sq256(ua, a2);
    h9e4h_sq256(ub, b2);
    h9e4h_mul3_256(b2);
    const int c = h9e4h_cmp256(a2, b2);
    return c > 0 ? sA : (c < 0 ? sB : 0);
}

static int e4r_cmp(e4r a, e4r b) { return e4r_sgn(e4r_sub(a, b)); }

/* ── the signed-distance score, exact ───────────────────────────────
 * min over the 4 trapezoid edges of Im(conj(e)(w - p)), computed on
 * the common positive scale 4*2^SIGMA (2w - 2p, 2e): sign and order
 * free. w is at scale 2^-SIGMA. */
static e4r h9e4h_score_min(const e4c *w) {
    e4r best = { 0, 0 };
    for (int i = 0; i < 4; ++i) {
        const e4i sh = (e4i)1 << E4H_SIGMA;
        e4r tr = { 2 * w->re.A - E4H_EDGE_P2[i][0] * sh,
                   2 * w->re.B - E4H_EDGE_P2[i][1] * sh };
        e4r ti = { 2 * w->im.A - E4H_EDGE_P2[i][2] * sh,
                   2 * w->im.B - E4H_EDGE_P2[i][3] * sh };
        /* Im(conj(E)T) = Er*Ti - Ei*Tr, ring-multiplied */
        const e4r s = e4r_sub(e4r_muls(ti, E4H_EDGE_E2[i][0], E4H_EDGE_E2[i][1]),
                              e4r_muls(tr, E4H_EDGE_E2[i][2], E4H_EDGE_E2[i][3]));
        if (i == 0 || e4r_cmp(s, best) < 0) best = s;
    }
    return best;
}

/* first-max argmax over candidate residuals (numpy argmax tie rule) */
static int h9e4h_argmax(const e4c *cands, int n, e4c *out) {
    int best = 0;
    e4r bs = h9e4h_score_min(&cands[0]);
    for (int k = 1; k < n; ++k) {
        const e4r s = h9e4h_score_min(&cands[k]);
        if (e4r_cmp(s, bs) > 0) { bs = s; best = k; }
    }
    *out = cands[best];
    return best;
}

/* child residual u = M_k (w - B_k); B_k enters at scale /4 */
static e4c h9e4h_child(const e4c *w, int k) {
    const e4i sh2 = (e4i)1 << (E4H_SIGMA - 2);
    e4r tr = { w->re.A - E4H_B4[k][0] * sh2, w->re.B - E4H_B4[k][1] * sh2 };
    e4r ti = { w->im.A - E4H_B4[k][2] * sh2, w->im.B - E4H_B4[k][3] * sh2 };
    /* (Mr + i Mi)(tr + i ti) = (Mr tr - Mi ti) + i(Mr ti + Mi tr) */
    e4c u;
    u.re = e4r_sub(e4r_muls(tr, E4H_M[k][0], E4H_M[k][1]),
                   e4r_muls(ti, E4H_M[k][2], E4H_M[k][3]));
    e4r a = e4r_muls(ti, E4H_M[k][0], E4H_M[k][1]);
    e4r b = e4r_muls(tr, E4H_M[k][2], E4H_M[k][3]);
    u.im.A = a.A + b.A;
    u.im.B = a.B + b.B;
    return u;
}

/* ── the five-symbol digit algebra (pure small integers) ──────────── */
static int h9e4h_cls(int s) { return ((s + 2) % 3) + 1; }

static int h9e4h_digit5(int oid, int c2, int half, int cls) {
    if (cls == 1 + half) return 1;
    const int *b = E4H_BASES[oid];
    if (cls == 3) return b[((1 - c2) % 3 + 3) % 3];
    return half == 0 ? b[((0 - c2) % 3 + 3) % 3]
                     : b[((2 - c2) % 3 + 3) % 3];
}

/* inverse for a context; 0 = digit invalid in this context */
static int h9e4h_digit_class(int oid, int c2, int half, int digit) {
    for (int cl = 1; cl <= 3; ++cl)
        if (h9e4h_digit5(oid, c2, half, cl) == digit) return cl;
    return 0;
}

/* ── the exact descent ──────────────────────────────────────────────
 * (wx, wy): the point in the host's canonical state frame (FP seed).
 * Snaps once, then classifies half and `depth` children exactly.
 * Writes half and the tail digits. */
static int h9e4h_descend(double wx, double wy, int hoid, int c2, int depth,
                         int *half_out, uint8_t *digits) {
    if (!isfinite(wx) || !isfinite(wy) ||
        fabs(wx) > 64.0 || fabs(wy) > 64.0) return E4H_EBAD;
    const double scale = (double)((int64_t)1 << E4H_SIGMA);
    e4c w;
    w.re.A = (e4i)llround(wx * scale); w.re.B = 0;
    w.im.A = (e4i)llround(wy * scale); w.im.B = 0;

    e4c cand[4], u;
    cand[0] = w;
    cand[1].re.A = -w.re.A; cand[1].re.B = -w.re.B;
    cand[1].im.A = -w.im.A; cand[1].im.B = -w.im.B;
    const int half = h9e4h_argmax(cand, 2, &u);
    int s = 3 * half;
    for (int d = 0; d < depth; ++d) {
        for (int k = 0; k < 4; ++k) cand[k] = h9e4h_child(&u, k);
        const int k = h9e4h_argmax(cand, 4, &u);
        s += E4H_ROT[k];
        digits[d] = (uint8_t)(k == 0 ? 0
                              : h9e4h_digit5(hoid, c2, half, h9e4h_cls(s)));
    }
    *half_out = half;
    return 0;
}

/* ── seed helpers (frozen det-math FP) ────────────────────────────── */

/* point (px,py) in octant g -> host octant ho's extended frame.
 * 0 ok; E4H_ESEAM when g is not ho or one of its 3 edge neighbours. */
static int h9e4h_unfold_point(double *px, double *py, int g, int ho) {
    const int x = g ^ ho;
    if (!x) return 0;
    if (x & (x - 1)) return E4H_ESEAM;
    const int axis = (x == 1) ? 0 : (x == 2) ? 1 : 2;
    const int mode = (int)H9_OID_MO[ho];
    int k = -1;
    for (int i = 0; i < 3; ++i) if (E4H_AXK[mode][i] == axis) k = i;
    if (k < 0) return E4H_ESEAM;
    const double x0 = *px, y0 = *py;
    *px = x0 * E4H_UNF[mode][k][0][0] + y0 * E4H_UNF[mode][k][1][0]
        + E4H_UNF[mode][k][2][0];
    *py = x0 * E4H_UNF[mode][k][0][1] + y0 * E4H_UNF[mode][k][1][1]
        + E4H_UNF[mode][k][2][1];
    return 0;
}

/* w = A[p_mo][c2] * 3^layer * conj(z - c): the host state frame */
static void h9e4h_to_frame(double zx, double zy, double cx, double cy,
                           int p_mo, int c2, int layer,
                           double *wx, double *wy) {
    double p3 = 1.0;
    for (int i = 0; i < layer; ++i) p3 *= 3.0;
    const double ar = E4H_A_RE[p_mo][c2] * p3;
    const double ai = E4H_A_IM[p_mo][c2] * p3;
    const double dx = zx - cx, dy = -(zy - cy);   /* conj */
    *wx = ar * dx - ai * dy;
    *wy = ar * dy + ai * dx;
}

/* inverse frame: canonical w -> host chart z (may overhang the octant) */
static void h9e4h_from_frame(double wx, double wy, double cx, double cy,
                             int p_mo, int c2, int layer,
                             double *zx, double *zy) {
    double p3 = 1.0;
    for (int i = 0; i < layer; ++i) p3 *= 3.0;
    const double ar = E4H_A_RE[p_mo][c2] * p3;
    const double ai = E4H_A_IM[p_mo][c2] * p3;
    const double n = ar * ar + ai * ai;
    const double ux = (wx * ar + wy * ai) / n;
    const double uy = (wy * ar - wx * ai) / n;
    /* z = conj(u) + c  (since b = -a*conj(c)) */
    *zx = ux + cx;
    *zy = -uy + cy;
}

/* fold an overhanging host-chart point into its containing octant
 * (frozen port of hhg9 fold_to_octant: <=3 seam reflections) */
static void h9e4h_fold(double *x, double *y, int *oid) {
    for (int it = 0; it < 3; ++it) {
        const int mode = (int)H9_OID_MO[*oid];
        int hit = -1;
        for (int k = 0; k < 3 && hit < 0; ++k) {
            const double ax = E4H_SV[mode][(k + 1) % 3][0];
            const double ay = E4H_SV[mode][(k + 1) % 3][1];
            const double bx = E4H_SV[mode][(k + 2) % 3][0];
            const double by = E4H_SV[mode][(k + 2) % 3][1];
            const double cross = (bx - ax) * (*y - ay) - (by - ay) * (*x - ax);
            if (cross > 1e-12) hit = k;        /* CW polygon: inside <= 0 */
        }
        if (hit < 0) return;
        const double x0 = *x, y0 = *y;
        *x = x0 * E4H_FLD[mode][hit][0][0] + y0 * E4H_FLD[mode][hit][1][0]
           + E4H_FLD[mode][hit][2][0];
        *y = x0 * E4H_FLD[mode][hit][0][1] + y0 * E4H_FLD[mode][hit][1][1]
           + E4H_FLD[mode][hit][2][1];
        *oid ^= 1 << E4H_AXK[mode][hit];
    }
}

/* ── decode composition (advisory FP; digit inversion is integer) ───
 * Composes half+digits into the leaf similarity and applies the probe
 * (E4H_CEN for decode, its long-side mirror for the partner). Returns
 * the representative in the CANONICAL frame; E4H_EGRAM on a digit
 * invalid for the context. */
static int h9e4h_compose(int half, const uint8_t *digits, int nd,
                         int hoid, int c2, double probe_re, double probe_im,
                         double *rep_re, double *rep_im) {
    double ca_re = (double)E4H_HALF_SIGN[half], ca_im = 0.0;
    double cb_re = 0.0, cb_im = 0.0;
    int s = 3 * half;
    for (int i = 0; i < nd; ++i) {
        int k;
        if (digits[i] == 0) {
            k = 0;
        } else {
            const int cl = h9e4h_digit_class(hoid, c2, half, digits[i]);
            if (!cl) return E4H_EGRAM;
            k = E4H_P3[((cl - s) % 3 + 3) % 3];
        }
        /* a_k = conj(M_k)/4, b_k = B4_k/4 (exact table forms) */
        const double a2_re = (E4H_M[k][0] + E4H_M[k][1] * E4H_S3) * 0.25;
        const double a2_im = -(E4H_M[k][2] + E4H_M[k][3] * E4H_S3) * 0.25;
        const double b2_re = (E4H_B4[k][0] + E4H_B4[k][1] * E4H_S3) * 0.25;
        const double b2_im = (E4H_B4[k][2] + E4H_B4[k][3] * E4H_S3) * 0.25;
        const double na_re = ca_re * a2_re - ca_im * a2_im;
        const double na_im = ca_re * a2_im + ca_im * a2_re;
        const double nb_re = ca_re * b2_re - ca_im * b2_im + cb_re;
        const double nb_im = ca_re * b2_im + ca_im * b2_re + cb_im;
        ca_re = na_re; ca_im = na_im; cb_re = nb_re; cb_im = nb_im;
        s += E4H_ROT[k];
    }
    *rep_re = ca_re * probe_re - ca_im * probe_im + cb_re;
    *rep_im = ca_re * probe_im + ca_im * probe_re + cb_im;
    return 0;
}
