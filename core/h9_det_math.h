/* h9_det_math.h — deterministic transcendentals for the hex9 encode/decode chain.
 *
 * WHY THIS FILE EXISTS. A full-depth hex9 uuid is a promise: this lon/lat IS
 * this address, on every platform. libm cannot keep that promise — Apple,
 * glibc and mingw disagree in the last ulps of sin/cos/tan/atan2, and CI
 * measured the consequence: boundary-adjacent points flip deep nibbles per
 * platform (regime_pin, 2026-07-27). Universality therefore requires that the
 * chain's floating-point program be THE definition of the address: same bits
 * in, same operations, same bits out, everywhere. These kernels + pinned FP
 * contraction (-ffp-contract=off, see CMakeLists) are that program.
 *
 * The kernels are musl 1.2.5 (via fdlibm), taken VERBATIM — coefficients,
 * evaluation order, branch structure — renamed h9_* and made header-only.
 * Only three deliberate departures:
 *   1. No __rem_pio2_large: the chain feeds |x| <= ~2*pi. Arguments beyond
 *      2^20*(pi/2) reduce by fmod(x, 2*pi_dbl) — IEEE-exact, deterministic,
 *      slightly WRONG in absolute phase for astronomic args (2*pi_dbl != 2*pi)
 *      and therefore fine: deterministic-but-degraded beats libm-accurate-but-
 *      platform-dependent, and the chain never goes there.
 *   2. FP-exception side effects (FORCE_EVAL) dropped — flags are not results.
 *   3. Word access via memcpy puns (C++-clean), not macro casts.
 *
 * DO NOT "improve" anything here — reordering one addition changes addresses.
 * That is not a style rule, it is the product's central promise. Any change is
 * a regime change: deliberate, versioned, goldens re-derived from source.
 *
 * sqrt and fmod are used freely: IEEE-754 requires them correctly rounded, so
 * they are deterministic on every conforming platform.
 *
 * ---------------------------------------------------------------------------
 * Derived from musl libc 1.2.5 (MIT), in turn from FreeBSD msun / fdlibm:
 *
 *   Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *   Developed at SunSoft, a Sun Microsystems, Inc. business.
 *   Permission to use, copy, modify, and distribute this software is freely
 *   granted, provided that this notice is preserved.
 *
 *   musl as a whole is licensed under the standard MIT license (see the musl
 *   COPYRIGHT file; authors include Rich Felker and contributors).
 * ---------------------------------------------------------------------------
 */
#ifndef H9_DET_MATH_H
#define H9_DET_MATH_H

#include <math.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

/* Universality holds only when double expressions are evaluated in double.
 * True on x86-64 (SSE2), arm64, wasm; false on 32-bit x87 builds — which are
 * therefore not conforming hex9 platforms. */
#if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD != 0
#error "hex9 determinism requires FLT_EVAL_METHOD == 0 (no excess precision)"
#endif

/* Clang honours the standard pragma; GCC needs -ffp-contract=off (pinned in
 * the build — both belts are worn). */
#ifdef __clang__
#pragma STDC FP_CONTRACT OFF
#endif

#ifdef __cplusplus
extern "C" {
#endif

static inline uint64_t h9_det_bits_(double x)  { uint64_t i; memcpy(&i, &x, 8); return i; }
static inline double   h9_det_dbl_(uint64_t i) { double x;   memcpy(&x, &i, 8); return x; }
static inline uint32_t h9_det_hi_(double x)    { return (uint32_t)(h9_det_bits_(x) >> 32); }

/* ── kernel sin on [-pi/4, pi/4]  (musl __sin.c) ─────────────────────────── */
static const double h9_det_S1_ = -1.66666666666666324348e-01, /* 0xBFC55555, 0x55555549 */
                    h9_det_S2_ =  8.33333333332248946124e-03, /* 0x3F811111, 0x1110F8A6 */
                    h9_det_S3_ = -1.98412698298579493134e-04, /* 0xBF2A01A0, 0x19C161D5 */
                    h9_det_S4_ =  2.75573137070700676789e-06, /* 0x3EC71DE3, 0x57B1FE7D */
                    h9_det_S5_ = -2.50507602534068634195e-08, /* 0xBE5AE5E6, 0x8A2B9CEB */
                    h9_det_S6_ =  1.58969099521155010221e-10; /* 0x3DE5D93A, 0x5ACFD57C */

static inline double h9_det_ksin_(double x, double y, int iy)
{
    double z, r, v, w;

    z = x*x;
    w = z*z;
    r = h9_det_S2_ + z*(h9_det_S3_ + z*h9_det_S4_) + z*w*(h9_det_S5_ + z*h9_det_S6_);
    v = z*x;
    if (iy == 0)
        return x + v*(h9_det_S1_ + z*r);
    else
        return x - ((z*(0.5*y - v*r) - y) - v*h9_det_S1_);
}

/* ── kernel cos on [-pi/4, pi/4]  (musl __cos.c) ─────────────────────────── */
static const double h9_det_C1_ =  4.16666666666666019037e-02, /* 0x3FA55555, 0x5555554C */
                    h9_det_C2_ = -1.38888888888741095749e-03, /* 0xBF56C16C, 0x16C15177 */
                    h9_det_C3_ =  2.48015872894767294178e-05, /* 0x3EFA01A0, 0x19CB1590 */
                    h9_det_C4_ = -2.75573143513906633035e-07, /* 0xBE927E4F, 0x809C52AD */
                    h9_det_C5_ =  2.08757232129817482790e-09, /* 0x3E21EE9E, 0xBDB4B1C4 */
                    h9_det_C6_ = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

static inline double h9_det_kcos_(double x, double y)
{
    double hz, z, r, w;

    z  = x*x;
    w  = z*z;
    r  = z*(h9_det_C1_+z*(h9_det_C2_+z*h9_det_C3_)) + w*w*(h9_det_C4_+z*(h9_det_C5_+z*h9_det_C6_));
    hz = 0.5*z;
    w  = 1.0-hz;
    return w + (((1.0-w)-hz) + (z*r-x*y));
}

/* ── kernel tan on ~[-pi/4, pi/4]  (musl __tan.c) ────────────────────────── */
static const double h9_det_T_[] = {
             3.33333333333334091986e-01, /* 3FD55555, 55555563 */
             1.33333333333201242699e-01, /* 3FC11111, 1110FE7A */
             5.39682539762260521377e-02, /* 3FABA1BA, 1BB341FE */
             2.18694882948595424599e-02, /* 3F9664F4, 8406D637 */
             8.86323982359930005737e-03, /* 3F8226E3, E96E8493 */
             3.59207910759131235356e-03, /* 3F6D6D22, C9560328 */
             1.45620945432529025516e-03, /* 3F57DBC8, FEE08315 */
             5.88041240820264096874e-04, /* 3F4344D8, F2F26501 */
             2.46463134818469906812e-04, /* 3F3026F7, 1A8D1068 */
             7.81794442939557092300e-05, /* 3F147E88, A03792A6 */
             7.14072491382608190305e-05, /* 3F12B80F, 32F0A7E9 */
            -1.85586374855275456654e-05, /* BEF375CB, DB605373 */
             2.59073051863633712884e-05, /* 3EFB2A70, 74BF7AD4 */
};
static const double h9_det_pio4_   = 7.85398163397448278999e-01, /* 3FE921FB, 54442D18 */
                    h9_det_pio4lo_ = 3.06161699786838301793e-17; /* 3C81A626, 33145C07 */

static inline double h9_det_ktan_(double x, double y, int odd)
{
    double z, r, v, w, s, a;
    double w0, a0;
    uint32_t hx;
    int big, sign;

    hx = h9_det_hi_(x);
    big = (hx&0x7fffffff) >= 0x3FE59428; /* |x| >= 0.6744 */
    if (big) {
        sign = hx>>31;
        if (sign) {
            x = -x;
            y = -y;
        }
        x = (h9_det_pio4_ - x) + (h9_det_pio4lo_ - y);
        y = 0.0;
    }
    z = x * x;
    w = z * z;
    r = h9_det_T_[1] + w*(h9_det_T_[3] + w*(h9_det_T_[5] + w*(h9_det_T_[7] + w*(h9_det_T_[9] + w*h9_det_T_[11]))));
    v = z*(h9_det_T_[2] + w*(h9_det_T_[4] + w*(h9_det_T_[6] + w*(h9_det_T_[8] + w*(h9_det_T_[10] + w*h9_det_T_[12])))));
    s = z * x;
    r = y + z*(s*(r + v) + y) + s*h9_det_T_[0];
    w = x + r;
    if (big) {
        s = 1 - 2*odd;
        v = s - 2.0 * (x + (r - w*w/(w + s)));
        return sign ? -v : v;
    }
    if (!odd)
        return w;
    /* -1.0/(x+r) has up to 2ulp error, so compute it accurately */
    w0 = w;
    w0 = h9_det_dbl_(h9_det_bits_(w0) & 0xffffffff00000000ULL); /* SET_LOW_WORD(w0, 0) */
    v = r - (w0 - x);       /* w0+v = r+x */
    a0 = a = -1.0 / w;
    a0 = h9_det_dbl_(h9_det_bits_(a0) & 0xffffffff00000000ULL); /* SET_LOW_WORD(a0, 0) */
    return a0 + a*(1.0 + a0*w0 + a0*v);
}

/* ── argument reduction x rem pi/2  (musl __rem_pio2.c, large path replaced
 *    by the deterministic fmod fallback — departure #1 in the header). ───── */
static const double
h9_det_toint_   = 1.5/DBL_EPSILON,
h9_det_invpio2_ = 6.36619772367581382433e-01, /* 0x3FE45F30, 0x6DC9C883 */
h9_det_pio2_1_  = 1.57079632673412561417e+00, /* 0x3FF921FB, 0x54400000 */
h9_det_pio2_1t_ = 6.07710050650619224932e-11, /* 0x3DD0B461, 0x1A626331 */
h9_det_pio2_2_  = 6.07710050630396597660e-11, /* 0x3DD0B461, 0x1A600000 */
h9_det_pio2_2t_ = 2.02226624879595063154e-21, /* 0x3BA3198A, 0x2E037073 */
h9_det_pio2_3_  = 2.02226624871116645580e-21, /* 0x3BA3198A, 0x2E000000 */
h9_det_pio2_3t_ = 8.47842766036889956997e-32; /* 0x397B839A, 0x252049C1 */

static inline int h9_det_rem_pio2_(double x, double *y)
{
    uint64_t bits = h9_det_bits_(x);
    double z, w, t, r, fn;
    uint32_t ix;
    int sign, n, ex, ey;

    sign = (int)(bits>>63);
    ix = (uint32_t)(bits>>32) & 0x7fffffff;
    if (ix <= 0x400f6a7a) {  /* |x| ~<= 5pi/4 */
        if ((ix & 0xfffff) == 0x921fb)  /* |x| ~= pi/2 or 2pi/2 */
            goto medium;  /* cancellation -- use medium case */
        if (ix <= 0x4002d97c) {  /* |x| ~<= 3pi/4 */
            if (!sign) {
                z = x - h9_det_pio2_1_;  /* one round good to 85 bits */
                y[0] = z - h9_det_pio2_1t_;
                y[1] = (z-y[0]) - h9_det_pio2_1t_;
                return 1;
            } else {
                z = x + h9_det_pio2_1_;
                y[0] = z + h9_det_pio2_1t_;
                y[1] = (z-y[0]) + h9_det_pio2_1t_;
                return -1;
            }
        } else {
            if (!sign) {
                z = x - 2*h9_det_pio2_1_;
                y[0] = z - 2*h9_det_pio2_1t_;
                y[1] = (z-y[0]) - 2*h9_det_pio2_1t_;
                return 2;
            } else {
                z = x + 2*h9_det_pio2_1_;
                y[0] = z + 2*h9_det_pio2_1t_;
                y[1] = (z-y[0]) + 2*h9_det_pio2_1t_;
                return -2;
            }
        }
    }
    if (ix <= 0x401c463b) {  /* |x| ~<= 9pi/4 */
        if (ix <= 0x4015fdbc) {  /* |x| ~<= 7pi/4 */
            if (ix == 0x4012d97c)  /* |x| ~= 3pi/2 */
                goto medium;
            if (!sign) {
                z = x - 3*h9_det_pio2_1_;
                y[0] = z - 3*h9_det_pio2_1t_;
                y[1] = (z-y[0]) - 3*h9_det_pio2_1t_;
                return 3;
            } else {
                z = x + 3*h9_det_pio2_1_;
                y[0] = z + 3*h9_det_pio2_1t_;
                y[1] = (z-y[0]) + 3*h9_det_pio2_1t_;
                return -3;
            }
        } else {
            if (ix == 0x401921fb)  /* |x| ~= 4pi/2 */
                goto medium;
            if (!sign) {
                z = x - 4*h9_det_pio2_1_;
                y[0] = z - 4*h9_det_pio2_1t_;
                y[1] = (z-y[0]) - 4*h9_det_pio2_1t_;
                return 4;
            } else {
                z = x + 4*h9_det_pio2_1_;
                y[0] = z + 4*h9_det_pio2_1t_;
                y[1] = (z-y[0]) + 4*h9_det_pio2_1t_;
                return -4;
            }
        }
    }
    if (ix >= 0x7ff00000) {  /* x is inf or NaN */
        y[0] = y[1] = x - x;
        return 0;
    }
    if (ix >= 0x41392000) {
        /* |x| >= 2^20*(pi/2): the chain never sends this. Deterministic
         * fallback — exact fmod into the medium range, degraded absolute
         * phase, identical on every platform (departure #1). */
        x = fmod(x, 2*3.14159265358979311600e+00);
        bits = h9_det_bits_(x);
        ix = (uint32_t)(bits>>32) & 0x7fffffff;
        if (ix <= 0x3fe921fb) { /* landed within pi/4: no reduction needed */
            y[0] = x;
            y[1] = 0.0;
            return 0;
        }
    }
medium:
    /* rint(x/(pi/2)) */
    fn = x*h9_det_invpio2_ + h9_det_toint_ - h9_det_toint_;
    n = (int32_t)fn;
    r = x - fn*h9_det_pio2_1_;
    w = fn*h9_det_pio2_1t_;  /* 1st round, good to 85 bits */
    /* Matters with directed rounding. */
    if (r - w < -h9_det_pio4_) {
        n--;
        fn--;
        r = x - fn*h9_det_pio2_1_;
        w = fn*h9_det_pio2_1t_;
    } else if (r - w > h9_det_pio4_) {
        n++;
        fn++;
        r = x - fn*h9_det_pio2_1_;
        w = fn*h9_det_pio2_1t_;
    }
    y[0] = r - w;
    ey = (int)(h9_det_bits_(y[0])>>52 & 0x7ff);
    ex = (int)(ix>>20);
    if (ex - ey > 16) { /* 2nd round, good to 118 bits */
        t = r;
        w = fn*h9_det_pio2_2_;
        r = t - w;
        w = fn*h9_det_pio2_2t_ - ((t-r)-w);
        y[0] = r - w;
        ey = (int)(h9_det_bits_(y[0])>>52 & 0x7ff);
        if (ex - ey > 49) {  /* 3rd round, good to 151 bits, covers all cases */
            t = r;
            w = fn*h9_det_pio2_3_;
            r = t - w;
            w = fn*h9_det_pio2_3t_ - ((t-r)-w);
            y[0] = r - w;
        }
    }
    y[1] = (r - y[0]) - w;
    return n;
}

/* ── sin / cos / tan drivers  (musl sin.c / cos.c / tan.c) ───────────────── */
static inline double h9_sin(double x)
{
    double y[2];
    uint32_t ix;
    unsigned n;

    ix = h9_det_hi_(x) & 0x7fffffff;

    /* |x| ~< pi/4 */
    if (ix <= 0x3fe921fb) {
        if (ix < 0x3e500000)  /* |x| < 2**-26 */
            return x;
        return h9_det_ksin_(x, 0.0, 0);
    }

    /* sin(Inf or NaN) is NaN */
    if (ix >= 0x7ff00000)
        return x - x;

    /* argument reduction needed */
    n = (unsigned)h9_det_rem_pio2_(x, y);
    switch (n&3) {
    case 0: return  h9_det_ksin_(y[0], y[1], 1);
    case 1: return  h9_det_kcos_(y[0], y[1]);
    case 2: return -h9_det_ksin_(y[0], y[1], 1);
    default:
        return -h9_det_kcos_(y[0], y[1]);
    }
}

static inline double h9_cos(double x)
{
    double y[2];
    uint32_t ix;
    unsigned n;

    ix = h9_det_hi_(x) & 0x7fffffff;

    /* |x| ~< pi/4 */
    if (ix <= 0x3fe921fb) {
        if (ix < 0x3e46a09e)  /* |x| < 2**-27 * sqrt(2) */
            return 1.0;
        return h9_det_kcos_(x, 0);
    }

    /* cos(Inf or NaN) is NaN */
    if (ix >= 0x7ff00000)
        return x-x;

    /* argument reduction */
    n = (unsigned)h9_det_rem_pio2_(x, y);
    switch (n&3) {
    case 0: return  h9_det_kcos_(y[0], y[1]);
    case 1: return -h9_det_ksin_(y[0], y[1], 1);
    case 2: return -h9_det_kcos_(y[0], y[1]);
    default:
        return  h9_det_ksin_(y[0], y[1], 1);
    }
}

static inline double h9_tan(double x)
{
    double y[2];
    uint32_t ix;
    unsigned n;

    ix = h9_det_hi_(x) & 0x7fffffff;

    /* |x| ~< pi/4 */
    if (ix <= 0x3fe921fb) {
        if (ix < 0x3e400000)  /* |x| < 2**-27 */
            return x;
        return h9_det_ktan_(x, 0.0, 0);
    }

    /* tan(Inf or NaN) is NaN */
    if (ix >= 0x7ff00000)
        return x - x;

    /* argument reduction */
    n = (unsigned)h9_det_rem_pio2_(x, y);
    return h9_det_ktan_(y[0], y[1], (int)(n&1));
}

/* ── atan  (musl atan.c) ─────────────────────────────────────────────────── */
static const double h9_det_atanhi_[] = {
  4.63647609000806093515e-01, /* atan(0.5)hi 0x3FDDAC67, 0x0561BB4F */
  7.85398163397448278999e-01, /* atan(1.0)hi 0x3FE921FB, 0x54442D18 */
  9.82793723247329054082e-01, /* atan(1.5)hi 0x3FEF730B, 0xD281F69B */
  1.57079632679489655800e+00, /* atan(inf)hi 0x3FF921FB, 0x54442D18 */
};
static const double h9_det_atanlo_[] = {
  2.26987774529616870924e-17, /* atan(0.5)lo 0x3C7A2B7F, 0x222F65E2 */
  3.06161699786838301793e-17, /* atan(1.0)lo 0x3C81A626, 0x33145C07 */
  1.39033110312309984516e-17, /* atan(1.5)lo 0x3C700788, 0x7AF0CBBD */
  6.12323399573676603587e-17, /* atan(inf)lo 0x3C91A626, 0x33145C07 */
};
static const double h9_det_aT_[] = {
  3.33333333333329318027e-01, /* 0x3FD55555, 0x5555550D */
 -1.99999999998764832476e-01, /* 0xBFC99999, 0x9998EBC4 */
  1.42857142725034663711e-01, /* 0x3FC24924, 0x920083FF */
 -1.11111104054623557880e-01, /* 0xBFBC71C6, 0xFE231671 */
  9.09088713343650656196e-02, /* 0x3FB745CD, 0xC54C206E */
 -7.69187620504482999495e-02, /* 0xBFB3B0F2, 0xAF749A6D */
  6.66107313738753120669e-02, /* 0x3FB10D66, 0xA0D03D51 */
 -5.83357013379057348645e-02, /* 0xBFADDE2D, 0x52DEFD9A */
  4.97687799461593236017e-02, /* 0x3FA97B4B, 0x24760DEB */
 -3.65315727442169155270e-02, /* 0xBFA2B444, 0x2C6A6C2F */
  1.62858201153657823623e-02, /* 0x3F90AD3A, 0xE322DA11 */
};

static inline double h9_atan(double x)
{
    double w, s1, s2, z;
    uint32_t ix, sign;
    int id;

    ix = h9_det_hi_(x);
    sign = ix >> 31;
    ix &= 0x7fffffff;
    if (ix >= 0x44100000) {   /* if |x| >= 2^66 */
        if (isnan(x))
            return x;
        z = h9_det_atanhi_[3] + 0x1p-120f;
        return sign ? -z : z;
    }
    if (ix < 0x3fdc0000) {    /* |x| < 0.4375 */
        if (ix < 0x3e400000)  /* |x| < 2^-27 */
            return x;
        id = -1;
    } else {
        x = fabs(x);
        if (ix < 0x3ff30000) {  /* |x| < 1.1875 */
            if (ix < 0x3fe60000) {  /*  7/16 <= |x| < 11/16 */
                id = 0;
                x = (2.0*x-1.0)/(2.0+x);
            } else {                /* 11/16 <= |x| < 19/16 */
                id = 1;
                x = (x-1.0)/(x+1.0);
            }
        } else {
            if (ix < 0x40038000) {  /* |x| < 2.4375 */
                id = 2;
                x = (x-1.5)/(1.0+1.5*x);
            } else {                /* 2.4375 <= |x| < 2^66 */
                id = 3;
                x = -1.0/x;
            }
        }
    }
    /* end of argument reduction */
    z = x*x;
    w = z*z;
    /* break sum from i=0 to 10 aT[i]z**(i+1) into odd and even poly */
    s1 = z*(h9_det_aT_[0]+w*(h9_det_aT_[2]+w*(h9_det_aT_[4]+w*(h9_det_aT_[6]+w*(h9_det_aT_[8]+w*h9_det_aT_[10])))));
    s2 = w*(h9_det_aT_[1]+w*(h9_det_aT_[3]+w*(h9_det_aT_[5]+w*(h9_det_aT_[7]+w*h9_det_aT_[9]))));
    if (id < 0)
        return x - x*(s1+s2);
    z = h9_det_atanhi_[id] - (x*(s1+s2) - h9_det_atanlo_[id] - x);
    return sign ? -z : z;
}

/* ── atan2  (musl atan2.c) ───────────────────────────────────────────────── */
static const double h9_det_pi_   = 3.1415926535897931160E+00, /* 0x400921FB, 0x54442D18 */
                    h9_det_pilo_ = 1.2246467991473531772E-16; /* 0x3CA1A626, 0x33145C07 */

static inline double h9_atan2(double y, double x)
{
    double z;
    uint32_t m, lx, ly, ix, iy;
    uint64_t bx, by;

    if (isnan(x) || isnan(y))
        return x+y;
    bx = h9_det_bits_(x);
    by = h9_det_bits_(y);
    ix = (uint32_t)(bx>>32); lx = (uint32_t)bx;
    iy = (uint32_t)(by>>32); ly = (uint32_t)by;
    if (((ix-0x3ff00000) | lx) == 0)  /* x = 1.0 */
        return h9_atan(y);
    m = ((iy>>31)&1) | ((ix>>30)&2);  /* 2*sign(x)+sign(y) */
    ix = ix & 0x7fffffff;
    iy = iy & 0x7fffffff;

    /* when y = 0 */
    if ((iy|ly) == 0) {
        switch (m) {
        case 0:
        case 1: return y;             /* atan(+-0,+anything)=+-0 */
        case 2: return  h9_det_pi_;   /* atan(+0,-anything) = pi */
        case 3: return -h9_det_pi_;   /* atan(-0,-anything) =-pi */
        }
    }
    /* when x = 0 */
    if ((ix|lx) == 0)
        return m&1 ? -h9_det_pi_/2 : h9_det_pi_/2;
    /* when x is INF */
    if (ix == 0x7ff00000) {
        if (iy == 0x7ff00000) {
            switch (m) {
            case 0: return  h9_det_pi_/4;   /* atan(+INF,+INF) */
            case 1: return -h9_det_pi_/4;   /* atan(-INF,+INF) */
            case 2: return  3*h9_det_pi_/4; /* atan(+INF,-INF) */
            case 3: return -3*h9_det_pi_/4; /* atan(-INF,-INF) */
            }
        } else {
            switch (m) {
            case 0: return  0.0;          /* atan(+...,+INF) */
            case 1: return -0.0;          /* atan(-...,+INF) */
            case 2: return  h9_det_pi_;   /* atan(+...,-INF) */
            case 3: return -h9_det_pi_;   /* atan(-...,-INF) */
            }
        }
    }
    /* |y/x| > 0x1p64 */
    if (ix+(64<<20) < iy || iy == 0x7ff00000)
        return m&1 ? -h9_det_pi_/2 : h9_det_pi_/2;

    /* z = atan(|y/x|) without spurious underflow */
    if ((m&2) && iy+(64<<20) < ix)  /* |y/x| < 0x1p-64, x<0 */
        z = 0;
    else
        z = h9_atan(fabs(y/x));
    switch (m) {
    case 0: return z;                            /* atan(+,+) */
    case 1: return -z;                           /* atan(-,+) */
    case 2: return h9_det_pi_ - (z-h9_det_pilo_);/* atan(+,-) */
    default: /* case 3 */
        return (z-h9_det_pilo_) - h9_det_pi_;    /* atan(-,-) */
    }
}

/* ── hypot  (musl hypot.c; FLT_EVAL_METHOD==0 ⇒ SPLIT = 2^27+1) ──────────── */
static inline void h9_det_sq_(double *hi, double *lo, double x)
{
    double xh, xl, xc;

    xc = x*(0x1p27 + 1);
    xh = x - xc + xc;
    xl = x - xh;
    *hi = x*x;
    *lo = xh*xh - *hi + 2*xh*xl + xl*xl;
}

static inline double h9_hypot(double x, double y)
{
    uint64_t uxi = h9_det_bits_(x), uyi = h9_det_bits_(y), ut;
    int ex, ey;
    double hx, lx, hy, ly, z;

    /* arrange |x| >= |y| */
    uxi &= (uint64_t)-1>>1;
    uyi &= (uint64_t)-1>>1;
    if (uxi < uyi) {
        ut = uxi;
        uxi = uyi;
        uyi = ut;
    }

    /* special cases */
    ex = (int)(uxi>>52);
    ey = (int)(uyi>>52);
    x = h9_det_dbl_(uxi);
    y = h9_det_dbl_(uyi);
    /* note: hypot(inf,nan) == inf */
    if (ey == 0x7ff)
        return y;
    if (ex == 0x7ff || uyi == 0)
        return x;
    /* note: hypot(x,y) ~= x + y*y/x/2 with inexact for small y/x */
    if (ex - ey > 64)
        return x + y;

    /* precise sqrt argument in nearest rounding mode without overflow */
    z = 1;
    if (ex > 0x3ff+510) {
        z = 0x1p700;
        x *= 0x1p-700;
        y *= 0x1p-700;
    } else if (ey < 0x3ff-450) {
        z = 0x1p-700;
        x *= 0x1p700;
        y *= 0x1p700;
    }
    h9_det_sq_(&hx, &lx, x);
    h9_det_sq_(&hy, &ly, y);
    return z*sqrt(ly+lx+hy+hx);
}

/* ── quarter-root: replaces every pow(x, 0.25) in the chain. Two IEEE-exact
 *    operations — deterministic AND more accurate than libm pow. ─────────── */
static inline double h9_qroot(double x)
{
    return sqrt(sqrt(x));
}

#ifdef __cplusplus
}
#endif

#endif /* H9_DET_MATH_H */
