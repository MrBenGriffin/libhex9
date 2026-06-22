/*
 * h9_math.h — Self-contained Hex9 math for the postgis_hex9 extension.
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2025, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * Ported from h9_boct.cpp (PROJ plugin) and h9_proj/test_h9_standalone.cpp.
 * No external dependencies — pure C++ math, no PROJ, no PostgreSQL.
 *
 * Requires h9_warp_data.h in the same directory (copy from h9_proj/).
 *
 * Public interface
 * ────────────────
 *   H9BOct                   — signed normalised barycentric b_oct point
 *   h9_lonlat_to_boct()      — (lon,lat) radians → H9BOct  [guarded AJ Newton]
 *   h9_lonlat_to_boct_beam() — same, exact beam search [fallback / reference]
 *   h9_boct_to_lonlat()      — H9BOct → (lon,lat) radians  [AK inverse]
 *   h9_lonlatdeg_to_boct()   — degree-input convenience
 *   h9_boct_to_lonlatdeg()   — degree-output convenience
 *   h9_braw_from_boct()      — extract unoriented (fx,fy) b_raw from H9BOct
 *   h9_orient_fwd()          — CCW orient rotation per octant
 *   h9_orient_inv()          — inverse (CW) orient rotation
 *   h9_warp_fwd()            — authalic forward warp (for addressing layer)
 *   h9_warp_inv()            — authalic inverse warp  (Newton-Raphson)
 *
 * Addressing layer interface (implement in h9_addressing.h)
 * ──────────────────────────────────────────────────────────
 *   void   h9_boct_to_uuid(H9BOct, uint8_t uuid[16]);
 *   H9BOct h9_uuid_to_boct(const uint8_t uuid[16]);
 *   void   h9_bin_uuid(const uint8_t in[16], int layer, uint8_t out[16]);
 *
 * UUID layout (32 nibbles = 128 bits) — see h9_addressing.h for full spec.
 *   nibbles  0..29 : body L0..L29  (hierarchy path; valid 0..8 at L1+)
 *   nibble      30 : h_term — terminal RID (0..11)
 *   nibble      31 : key_tail = (p_mo<<3)|(p_c2<<1)|r_mo
 */

#pragma once
#include "h9_layout.h"   /* H9_LMAX, H9_DESC_MAX, ... (HEX9_USE_L29 gate) */
#include <cstdint>
#include <limits>
#include <math.h>
#include <stdint.h>

/* Set to 1 to enable authalic warp (requires the heavy CT mesh headers).
 * Set to 0 during binning / hexgrid work until warp retraining is complete. */
#ifndef H9_WARP_ENABLE
#define H9_WARP_ENABLE 1
#endif

/* Runtime warp:
 *   - load + Bell-Sibson gradient + Alfeld-CT coefficients are built at
 *     extension load from the .h9warp sidecar (embedded or on disk).
 *   - h9_warp_fwd / h9_warp_inv shims live in h9_warp_runtime.h.
 *   - The three former static-data headers (h9_wgs84_warp_data.h /
 *     h9_wgs84_warp_grads.h / h9_wgs84_ct_mesh.h) are gone — see
 *     [[warp-distribution-plan]] for the 620 MB → 2 MB reduction. */

/* ── WGS84 ellipsoid ─────────────────────────────────────────────────────── */
static const double H9_WGS84_A  = 6378137.0;
static const double H9_WGS84_B  = 6356752.3142451793;
static const double H9_WGS84_F  = 0.0033528106647474805;
static const double H9_WGS84_A2 = H9_WGS84_A * H9_WGS84_A;
static const double H9_WGS84_B2 = H9_WGS84_B * H9_WGS84_B;
static const double H9_WGS84_E2 = (2.0 * H9_WGS84_F) - (H9_WGS84_F * H9_WGS84_F);
static const double H9_NORM_B2  = (H9_WGS84_B / H9_WGS84_A) * (H9_WGS84_B / H9_WGS84_A);

/* ── H9_Constants ────────────────────────────────────────────────────────── */
static const double H9_R3 = sqrt(3.0);
static const double H9_W  = sqrt(2.0);
static const double H9_H  = H9_W * H9_R3 * 0.5;
static const double H9_1H_3 = H9_H / 3.0;     // Ḣ Vertical spacing between stacked child cells
static const double H9_2H_3 = H9_1H_3 * 2.0;  // alias of ΛC; dH+dW = H9_H
static const double H9_TR  = H9_W / 2.0;   // Triangle Limit Right.
static const double H9_TL  = -H9_TR;       // Triangle Limit Left.
static const double H9_C1  = H9_2H_3;      // ΛC Mode 1 Ceiling
static const double H9_F1  = -H9_1H_3;     // ΛF Mode 1 Floor
static const double H9_C0  = H9_1H_3;      // VC Mode 0 Ceiling
static const double H9_F0  = -H9_2H_3;     // VF Mode 0 Floor
static const double H9_LU  = H9_W/6;       // U Lattice Unit
static const double H9_LV  = H9_H/9;       // V Lattice Unit

/* ── AK formula constants ────────────────────────────────────────────────── */

static const double H9_ALPHA    = 3.227806237143884260376580;
static const double H9_EPS      = 1e-200;
static const double H9_VERT_EPS = 1e-15;


/* ── Orient rotation per octant (CCW, units of π/3) ─────────────────────── */
/*   Indexed by C++ oct_i = H9BOct.oct_i ^ 7                                  */
/*   C++ oct_i = ((eX≥0)?4:0)|((eY≥0)?2:0)|((eZ≥0)?1:0)                    */
/*   SWP NWP SEP NEP SWA NWA SEA NEA                                          */

static const int H9_ORIENT_TH[8] = {5, 2, 2, 5, 2, 5, 5, 2};

/* ── Child net_mode tables per parent mode ───────────────────────────────── */

static const int H9_UP_MODE[9] = {1, 1, 0, 1, 1, 0, 1, 0, 1};
static const int H9_DN_MODE[9] = {0, 1, 0, 1, 0, 0, 1, 0, 0};

/* ── H9 lattice: child cell offsets in b_raw metric barycentric space ─────── */
/*                                                                             */
/*   All values derived from sqrt(2)/sqrt(3) to preserve ULP consistency       */
/*   with Python's h9_constants() chain.                                       */
/*                                                                             */
/*   W = sqrt(2),  r3 = sqrt(3)                                                */
/*   H = W * r3 / 2    (supercell height)                                      */
/*   U = W / 6         (H9K.U — horizontal lattice step)                       */
/*   V = H / 9         (H9K.V — vertical lattice step)                         */
/*                                                                             */
/*   Up   (mode=1) URIs: [0x16,0x25,0x26,0x2a,0x34,0x35,0x39,0x3a,0x3e]        */
/*   Down (mode=0) URIs: [0x21,0x25,0x26,0x2a,0x2b,0x35,0x39,0x3a,0x49]        */

struct H9Tables {
    double H, W, inv_H, inv_W, inv_2H;
    double UP_X[9], UP_Y[9];
    double DN_X[9], DN_Y[9];

    H9Tables() {
        W = sqrt(2.0);
        H = W * sqrt(3.0) * 0.5;
        const double U  = W / 6.0;
        const double V  = H / 9.0;
        const double U2 = 2.0*U, V2 = 2.0*V, V4 = 4.0*V;
        inv_H  = 1.0 / H;
        inv_W  = 1.0 / W;
        inv_2H = 0.5 / H;
        /* up children */
        const double ux[] = {  0, -U,   0, +U, -U2, -U,   0, +U,  U2};
        const double uy[] = { V4, +V, +V2, +V, -V2, -V, -V2, -V, -V2};
        /* down children */
        const double dx[] = {-U2, -U,   0, +U,  U2, -U,   0, +U,   0};
        const double dy[] = { V2, +V, +V2, +V, +V2, -V, -V2, -V, -V4};
        for (int i = 0; i < 9; ++i) {
            UP_X[i] = ux[i]; UP_Y[i] = uy[i];
            DN_X[i] = dx[i]; DN_Y[i] = dy[i];
        }
    }
};

static const H9Tables H9;

/* ── b_oct point ─────────────────────────────────────────────────────────── */
/*                                                                             */
/*   (u, v, w) — signed normalised barycentric coordinates.                   */
/*     u maps to ECEF X, v to Y, w to Z (via h9_c_oct_to_c_ell).             */
/*     Absolute values sum to ≈ 1 (each ≈ 1/3 at the cell centre).            */
/*                                                                             */
/*   oct_i    — Python oid convention (matches Python H9O.oid):               */
/*              ((eZ<0)?4:0)|((eY<0)?2:0)|((eX<0)?1:0)                       */
/*              bit0=eX<0, bit1=eY<0, bit2=eZ<0 (1 means negative)           */
/*              Extract signs: su=(oct_i&1)?-1:+1, sv=(oct_i&2)?-1:+1,       */
/*                             sw=(oct_i&4)?-1:+1                             */
/*   oct_mode — 0=down (even # negative coords), 1=up (odd # negative)       */
/*              Parity of oct_i bits: mode(oct_i) = parity(oct_i) & 1        */
/*                                                                             */
/*   NOTE: H9BOct stores pre-warp, pre-orient b_raw coordinates expressed     */
/*   as signed barycentric. Python's 'boct' is the post-orient/post-warp      */
/*   form. Apply h9_braw_from_boct → h9_orient_fwd → h9_warp_fwd to get      */
/*   the warped (authalic) coordinates used by the addressing layer.          */

struct H9BOct {
    double u, v, w;
    int    oct_i;
    int    oct_mode;
};

/* ── OID/MODE ────────────────────────────────────────────────────────────── */
/*                                                                             */
/*   oid(eX, eY, eZ) returns the Python oid convention used throughout H9:    */
/*     bit2 = (eZ<0), bit1 = (eY<0), bit0 = (eX<0)   — 1 means negative      */
/*   This matches Python H9O.oid and is stored directly in H9BOct.oct_i.      */
/*                                                                             */
/*   mode(oid) returns parity of oid bits:                                    */
/*     0 = even number of negative axes = down (su·sv·sw > 0)                 */
/*     1 = odd  number of negative axes = up   (su·sv·sw < 0)                 */

template <typename T> uint8_t oid(T x, T y, T z) {
	const T tiny = std::numeric_limits<T>::min();
	if (x == T(0.0)) x = tiny;
	if (y == T(0.0)) y = tiny;
	if (z == T(0.0)) z = tiny;
	return static_cast<uint8_t>(
		((z < T(0.0)) << 2) |
		((y < T(0.0)) << 1) |
		 (x < T(0.0))
	);
}

static uint8_t mode(uint8_t oid) {
	return (oid ^ (oid >> 1) ^ (oid >> 2)) & 1;
}

/* ── AK core: normalised barycentric → un-normalised ECEF ─────────────────── */

static void h9_ak_core(double u, double v, double w,
                        double *x, double *y, double *z) {
    const double a  = H9_ALPHA;
    const double e  = H9_EPS;
    const double tu = tan((M_PI * u + e) * 0.5);
    const double tv = tan((M_PI * v + e) * 0.5);
    const double tw = tan((M_PI * w + e) * 0.5);
    const double u2 = tu*tu, v2 = tv*tv, w2 = tw*tw;
    *x = tu * pow(v2 + w2 + a * w2 * v2, 0.25);
    *y = tv * pow(u2 + w2 + a * u2 * w2, 0.25);
    *z = tw * pow(u2 + v2 + a * u2 * v2, 0.25);
}

/* Normalise raw XYZ to unit-scale ellipsoid (a=1, b=B/A). */
static void h9_ak_normalise(double *x, double *y, double *z) {
    const double n = sqrt(*x * *x + *y * *y + (*z * *z) / H9_NORM_B2);
    *x /= n; *y /= n; *z /= n;
}

/* ── Coordinate conversions ─────────────────────────────────────────────── */

/* Signed normalised barycentric → ECEF (WGS84, metres). */
static void h9_c_oct_to_c_ell(double u, double v, double w,
                               double *X, double *Y, double *Z) {
    const double su = (u >= 0.0) ? 1.0 : -1.0;
    const double sv = (v >= 0.0) ? 1.0 : -1.0;
    const double sw = (w >= 0.0) ? 1.0 : -1.0;
    const double au = fabs(u), av = fabs(v), aw = fabs(w);
    if (av < H9_VERT_EPS && aw < H9_VERT_EPS) {
        *X = H9_WGS84_A * su; *Y = 0.0; *Z = 0.0; return;
    }
    if (au < H9_VERT_EPS && aw < H9_VERT_EPS) {
        *X = 0.0; *Y = H9_WGS84_A * sv; *Z = 0.0; return;
    }
    if (au < H9_VERT_EPS && av < H9_VERT_EPS) {
        *X = 0.0; *Y = 0.0; *Z = H9_WGS84_A * sw; return;
    }
    h9_ak_core(au, av, aw, X, Y, Z);
    h9_ak_normalise(X, Y, Z);
    *X *= H9_WGS84_A * su;
    *Y *= H9_WGS84_A * sv;
    *Z *= H9_WGS84_A * sw;
}

/* ECEF → geodetic lon/lat (radians). Standard surface-Bowring (5 iters).
 * Mirrors Python hhg9/algorithms/wgs84.py::ecef_to_latlon — bit-identical
 * fixed point at h=0 surface points, which is what every caller passes. */
static void h9_c_ell_to_lonlat(double X, double Y, double Z,
                                double *lon, double *lat) {
    const double e2  = 1.0 - H9_WGS84_B2 / H9_WGS84_A2;
    const double p   = sqrt(X*X + Y*Y);
    *lon = atan2(Y, X);
    *lat = atan2(Z, p * (1.0 - e2));
    for (int k = 0; k < 5; ++k) {
        const double sin_lat = sin(*lat);
        const double N       = H9_WGS84_A / sqrt(1.0 - e2 * sin_lat * sin_lat);
        *lat = atan2(Z + e2 * N * sin_lat, p);
    }
}

/* Signed barycentric with explicit octant signs → lon/lat. */
static void h9_coct_to_lonlat(double u, double v, double w,
                               double su, double sv, double sw,
                               double *lon, double *lat) {
    double X, Y, Z;
    h9_c_oct_to_c_ell(su*u, sv*v, sw*w, &X, &Y, &Z);
    h9_c_ell_to_lonlat(X, Y, Z, lon, lat);
}

/* ── Warp: authalic correction ───────────────────────────────────────────── */
/*                                                                             */
/*   Grid stores the b_oct → b_raw displacement (AuthalicWarp.do convention). */
/*   Forward warp b_raw → b_oct requires Newton-Raphson inversion (= Python's */
/*   AuthalicWarp.undo), not direct bilinear application.                      */
/*                                                                             */
/*   Inputs use the braw_y convention: ry0 is positive for mode-1 UP octants. */
/*   The grid is in mode-0 space, so mode-1 y is negated before lookup and    */
/*   the result y is negated back.                                             */
/*                                                                             */
/*   NR iteration (20 steps): cx = rx - dx(cx,cy),  cy = ry - dy(cx,cy)      */
/*     mode-0:  ry = ry0 (no flip);  wy0 = cy        (no flip back)           */
/*     mode-1:  ry = -ry0 (flip);    wy0 = -cy        (flip back)             */

#if H9_WARP_ENABLE
/* h9_warp_fwd / h9_warp_inv live in h9_warp_runtime.h, which is included
 * at the end of this file to keep the legacy call sites unchanged. The
 * pre-port static-data CT walker and edge-blend stubs that used to live
 * here have been replaced by the runtime path.
 *
 * Note: caller MUST have invoked h9::h9_warp_init_embedded() before any
 * warp_fwd/inv call; otherwise the shims fall through to identity. */
#endif /* H9_WARP_ENABLE */

/* Project (fx,fy) onto the fundamental domain via UVW barycentric clamping.
 * Converts to barycentric coords, zeroes negatives, renormalises, converts back.
 * Handles corner violations correctly (sequential edge clamps fail at corners). */
static void h9_clamp_bary(double *fx, double *fy, int oct_mode) {
    const double iH  = H9.inv_H, iW = H9.inv_W, i2H = H9.inv_2H;
    double u, v, w;
    if (oct_mode == 0) {
        u = (1.0/3.0) - (*fy) * iH;
        v = (1.0/3.0) + (*fy) * i2H + (*fx) * iW;
        w = (1.0/3.0) + (*fy) * i2H - (*fx) * iW;
    } else {
        u = (1.0/3.0) + (*fy) * iH;
        v = (1.0/3.0) - (*fy) * i2H + (*fx) * iW;
        w = (1.0/3.0) - (*fy) * i2H - (*fx) * iW;
    }
    if (u < 0.0) u = 0.0;
    if (v < 0.0) v = 0.0;
    if (w < 0.0) w = 0.0;
    const double s = u + v + w;
    if (s > 1e-30) { u /= s; v /= s; w /= s; }
    if (oct_mode == 0) {
        *fy = -H9.H * (u - (1.0/3.0));
        *fx =  H9.W * 0.5 * (v - w);
    } else {
        *fy =  H9.H * (u - (1.0/3.0));
        *fx =  H9.W * 0.5 * (v - w);
    }
}

/* ── Orient rotation ─────────────────────────────────────────────────────── */
/*                                                                             */
/*   CCW rotation by ORIENT_TH[oct_i] * π/3 maps b_raw (fx,fy) to the        */
/*   canonical oriented frame (rx,ry) in which the warp grid is defined.      */

static void h9_orient_fwd(double fx, double fy, int oct_i,
                           double *rx, double *ry) {
    const double th = H9_ORIENT_TH[oct_i ^ 7] * (M_PI / 3.0);  /* ^7: Python oid → C++ oct_i */
    const double fc = cos(th), fs = sin(th);
    *rx = fx * fc - fy * fs;
    *ry = fx * fs + fy * fc;
}

/* Inverse (CW = transpose) rotation: oriented (rx,ry) → b_raw (fx,fy). */
static void h9_orient_inv(double rx, double ry, int oct_i,
                           double *fx, double *fy) {
    const double th = H9_ORIENT_TH[oct_i ^ 7] * (M_PI / 3.0);  /* ^7: Python oid → C++ oct_i */
    const double fc = cos(th), fs = sin(th);
    *fx =  rx * fc + ry * fs;
    *fy = -rx * fs + ry * fc;
}

/* ── Recover unoriented b_raw (fx, fy) from H9BOct ─────────────────────── */
/*                                                                             */
/*   Inverts the (u,v,w) ↔ (fx,fy) c_oct formulas for the given oct_mode.    */
/*   Result is in the same metric-barycentric space as H9Tables offsets.      */
/*                                                                             */
/*   Up   (mode=1): fx = W/2·(|w|−|v|),  fy = H·(|u|−1/3)                   */
/*   Down (mode=0): fx = W/2·(|v|−|w|),  fy = H·(1/3−|u|)                   */

static void h9_braw_from_boct(H9BOct b, double *fx, double *fy) {
    const double au = fabs(b.u), av = fabs(b.v), aw = fabs(b.w);
    if (b.oct_mode == 1) {
        *fy =  H9.H * (au - 1.0/3.0);
        *fx =  H9.W * 0.5 * (aw - av);
    } else {
        *fy = -H9.H * (au - 1.0/3.0);
        *fx =  H9.W * 0.5 * (av - aw);
    }
}

/* ── Encode: (lon, lat) radians → H9BOct ────────────────────────────────── */
/*                                                                             */
/*   h9_lonlat_to_boct_beam — beam search (width=6, depth≈34) in b_raw metric  */
/*   barycentric space. Matches h9_boct_fwd3d in test_h9_standalone.cpp and    */
/*   h9_boct.cpp (without PROJ type wrappers and without authalic warp — that  */
/*   is applied by the addressing layer in warped b_oct space if required).    */
/*                                                                             */
/*   This is now the FALLBACK ORACLE: the production h9_lonlat_to_boct (below) */
/*   solves the same inversion by a guarded analytic-Jacobian Gauss-Newton     */
/*   (~1300× faster), routing only the thin seam/vertex skin and any           */
/*   non-convergence here. Kept as the exact reference + skin solver.          */

static void h9_rad_lonlat_to_ecef(double lon_rad, double lat_rad, double *x, double *y, double *z) {
	const double sin_lat = sin(lat_rad);
	const double cos_lat = cos(lat_rad);
	const double n = H9_WGS84_A / sqrt(1.0 - H9_WGS84_E2 * (sin_lat * sin_lat));
	*x = n * cos_lat * cos(lon_rad);
	*y = n * cos_lat * sin(lon_rad);
	*z = n * (1.0 - H9_WGS84_E2) * sin_lat;
}

static void h9_ecef_to_rad_lonlat(double x, double y, double z, double *lon_rad, double *lat_rad) {
	const double p = sqrt(x * x + y * y);
	*lon_rad = atan2(y, x);
	*lat_rad = atan2(z, p * H9_NORM_B2);
}


/* Axis-vertex snap threshold (shared by the beam and the production solver).
 * The vertex shortcut snaps a point to the exact axis vertex; the proximity test
 * must use the TANGENTIAL (off-axis) distance², which is LINEAR near the vertex,
 * NOT |e_axis|−A, which is quadratic (`≈ −A·θ²/2`) and below its own ~1.4 nm ULP
 * for the whole region — so the old `|e|−A < 1e-9` captured a ~0.11 m disk and
 * snapped points up to 0.1 m away (visible as the ex0064 annulus / ex0059 seam
 * spikes). Tangential `e_b²+e_c² < TAN2` gives a clean √TAN2 capture radius.
 * Default (100 nm)²; set 0 to snap only the float-exact vertex. */
#ifndef H9_BOCT_VERTEX_TAN2
#define H9_BOCT_VERTEX_TAN2 1e-14
#endif

static H9BOct h9_lonlat_to_boct_beam(double lon_rad, double lat_rad) {
	H9BOct result;
	double eX, eY, eZ;
	h9_rad_lonlat_to_ecef(lon_rad, lat_rad, &eX, &eY, &eZ);

	const double cl = cos(lat_rad);
	result.oct_i = oid(eX, eY, eZ);
	result.oct_mode = mode(result.oct_i);

	/* Axis sign values — used in beam search (h9_coct_to_lonlat) and final
	 * result assignment. Separate from oct_i bit encoding. */
	const double su = (eX >= 0.0) ? 1.0 : -1.0;
	const double sv = (eY >= 0.0) ? 1.0 : -1.0;
	const double sw = (eZ >= 0.0) ? 1.0 : -1.0;

    /* Axis-vertex shortcut: snap to the exact vertex when the TANGENTIAL
     * (off-axis) distance² is within H9_BOCT_VERTEX_TAN2 — linear capture, see
     * the macro note. (Was |e_axis|−A < 1e-9: quadratic, ~0.11 m capture.) */
	if (eX*eX + eY*eY < H9_BOCT_VERTEX_TAN2) {          /* pole (Z axis) */
		result.w = eZ > 0 ? 1.0 : -1.0;
		result.u = 0.0; result.v = 0.0;
		return result;
	}

	if (eY*eY + eZ*eZ < H9_BOCT_VERTEX_TAN2) {          /* X axis */
		result.u = eX > 0 ? 1.0 : -1.0;
		result.v = 0.0; result.w = 0.0;
		return result;
	}

	if (eX*eX + eZ*eZ < H9_BOCT_VERTEX_TAN2) {          /* Y axis */
		result.v = eY > 0 ? 1.0 : -1.0;
		result.u = 0.0; result.w = 0.0;
		return result;
	}

	const int oct_mode = result.oct_mode;
    const double cos_lat  = cl;
    const double inv_H    = H9.inv_H, inv_W = H9.inv_W, inv_2H = H9.inv_2H;

    /* Beam search. BEAM=6 is empirically required (tricky nearest-triangle
     * regions need the extra candidate); DEPTH reduced 40→34 — still well past
     * layer-29 lattice resolution, validated output-preserving vs DEPTH=40. */
    const int BEAM = 6, DEPTH = 34;
    struct Cand { int mode; double x, y, dist; };
    Cand cur[BEAM], nxt[BEAM * 9];
    int bw = 1;
    cur[0] = {oct_mode, 0.0, 0.0, 1e300};

    double scale = 1.0;
    for (int d = 0; d < DEPTH; ++d) {
        int nn = 0;
        for (int k = 0; k < bw; ++k) {
            const double *off_x     = (cur[k].mode == 1) ? H9.UP_X : H9.DN_X;
            const double *off_y     = (cur[k].mode == 1) ? H9.UP_Y : H9.DN_Y;
            const int    *chld_mode = (cur[k].mode == 1) ? H9_UP_MODE : H9_DN_MODE;
            for (int j = 0; j < 9; ++j) {
                const double cx = cur[k].x + off_x[j] * scale;
                const double cy = cur[k].y + off_y[j] * scale;
                double cu, cv, cw;
                if (oct_mode == 1) {
                    cu = cy * inv_H  + (1.0/3.0);
                    cv = (1.0/3.0) - cy * inv_2H - cx * inv_W;
                    cw = (1.0/3.0) - cy * inv_2H + cx * inv_W;
                } else {
                    cu = (1.0/3.0) - cy * inv_H;
                    cv = (1.0/3.0) + cy * inv_2H + cx * inv_W;
                    cw = (1.0/3.0) + cy * inv_2H - cx * inv_W;
                }
                if (cu < 0.0) cu = 0.0;
                if (cv < 0.0) cv = 0.0;
                if (cw < 0.0) cw = 0.0;
                const double s = cu + cv + cw;
                if (s > 1e-30) { cu /= s; cv /= s; cw /= s; }
                double lon0, lat0;
                h9_coct_to_lonlat(cu, cv, cw, su, sv, sw, &lon0, &lat0);
                double dlat = lat0 - lat_rad;
                double dlon = lon0 - lon_rad;
                if (dlon >  M_PI) dlon -= 2.0*M_PI;
                if (dlon < -M_PI) dlon += 2.0*M_PI;
                const double sl = sin(dlat * 0.5), sd = sin(dlon * 0.5);
                nxt[nn++] = {chld_mode[j], cx, cy,
                             sl*sl + cos_lat * cos(lat0) * sd*sd};
            }
        }
        int bw_new = (nn < BEAM) ? nn : BEAM;
        /* Partial selection sort: put the best bw_new candidates first. */
        for (int i = 0; i < bw_new; ++i) {
            int best = i;
            for (int j = i+1; j < nn; ++j)
                if (nxt[j].dist < nxt[best].dist) best = j;
            Cand tmp = nxt[i]; nxt[i] = nxt[best]; nxt[best] = tmp;
        }
        for (int i = 0; i < bw_new; ++i) cur[i] = nxt[i];
        bw = bw_new;
        if (cur[0].dist < 3e-31) break;   /* tight: NR basin is sub-cell; looser eps
                                             corrupts ~19%+ of UUIDs (probed) */
        scale /= 3.0;
    }

    /* NR post-refinement: converge (fx,fy) to machine precision.
     * Skipped when the beam lands near a vertex: ∂lat/∂fy → 0 there, making
     * the Jacobian ill-conditioned and the NR step diverge despite det > 1e-30. */
    double fx = cur[0].x, fy = cur[0].y;
    bool skip_nr = true;
    {
        /* Vertex proximity check: skip NR if any barycentric coord > 0.85. */
        double vt_u, vt_v, vt_w;
        if (oct_mode == 0) {
            vt_u = (1.0/3.0) - fy * inv_H;
            vt_v = (1.0/3.0) + fy * inv_2H + fx * inv_W;
            vt_w = (1.0/3.0) + fy * inv_2H - fx * inv_W;
        } else {
            vt_u = (1.0/3.0) + fy * inv_H;
            vt_v = (1.0/3.0) - fy * inv_2H + fx * inv_W;
            vt_w = (1.0/3.0) - fy * inv_2H - fx * inv_W;
        }
        skip_nr = (vt_u > 0.85 || vt_v > 0.85 || vt_w > 0.85);
    }
    if (!skip_nr)
    {
        /* NR usually exits at iter 0 since beam search already converges
         * below the lon/lat residual threshold; nr_eps only matters in
         * the rare cases where beam plateaus above threshold. */
        const double nr_eps = 1e-7;
        auto to_lonlat = [&](double qx, double qy, double *qlon, double *qlat) {
            double qu, qv, qw;
            if (oct_mode == 1) {
                qu = qy * inv_H  + (1.0/3.0);
                qv = (1.0/3.0) - qy * inv_2H - qx * inv_W;
                qw = (1.0/3.0) - qy * inv_2H + qx * inv_W;
            } else {
                qu = (1.0/3.0) - qy * inv_H;
                qv = (1.0/3.0) + qy * inv_2H + qx * inv_W;
                qw = (1.0/3.0) + qy * inv_2H - qx * inv_W;
            }
            if (qu < 0.0) qu = 0.0;
            if (qv < 0.0) qv = 0.0;
            if (qw < 0.0) qw = 0.0;
            const double qs = qu + qv + qw;
            if (qs > 1e-30) { qu /= qs; qv /= qs; qw /= qs; }
            h9_coct_to_lonlat(qu, qv, qw, su, sv, sw, qlon, qlat);
        };
        for (int nr = 0; nr < 8; ++nr) {
            double lon0, lat0;
            to_lonlat(fx, fy, &lon0, &lat0);
            double dlon = lon_rad - lon0, dlat = lat_rad - lat0;
            if (dlon >  M_PI) dlon -= 2.0 * M_PI;
            if (dlon < -M_PI) dlon += 2.0 * M_PI;
            if (fabs(dlon) < 1e-15 && fabs(dlat) < 1e-15) break;
            double lp, lm, ap, am;
            to_lonlat(fx + nr_eps, fy, &lp, &ap);
            to_lonlat(fx - nr_eps, fy, &lm, &am);
            const double lon_fx = (lp - lm) / (2.0 * nr_eps);
            const double lat_fx = (ap - am) / (2.0 * nr_eps);
            to_lonlat(fx, fy + nr_eps, &lp, &ap);
            to_lonlat(fx, fy - nr_eps, &lm, &am);
            const double lon_fy = (lp - lm) / (2.0 * nr_eps);
            const double lat_fy = (ap - am) / (2.0 * nr_eps);
            const double det = lon_fx * lat_fy - lon_fy * lat_fx;
            if (fabs(det) < 1e-30) break;
            fx += ( lat_fy * dlon - lon_fy * dlat) / det;
            fy += (-lat_fx * dlon + lon_fx * dlat) / det;
        }
    }

    /* Convert final (fx,fy) to normalised barycentric (u,v,w). */
    double u, v, w;
    if (oct_mode == 1) {
        u = fy * inv_H  + (1.0/3.0);
        v = (1.0/3.0) - fy * inv_2H - fx * inv_W;
        w = (1.0/3.0) - fy * inv_2H + fx * inv_W;
    } else {
        u = (1.0/3.0) - fy * inv_H;
        v = (1.0/3.0) + fy * inv_2H + fx * inv_W;
        w = (1.0/3.0) + fy * inv_2H - fx * inv_W;
    }
    if (u < 0.0) u = 0.0;
    if (v < 0.0) v = 0.0;
    if (w < 0.0) w = 0.0;
    const double s = u + v + w;
    if (s > 1e-30) { u /= s; v /= s; w /= s; }

    result.u = su * u;
    result.v = sv * v;
    result.w = sw * w;
    return result;
}

/* ── Encode (fast): guarded analytic-Jacobian Gauss-Newton inversion ──────── */
/*                                                                             */
/*   Replaces the beam's ~1,800 trig projections with a few analytic steps     */
/*   (~1300× faster; the inversion is now sub-µs, on par with decode). It      */
/*   solves forward(fx,fy)=target in ECEF SURFACE space: the AK forward P      */
/*   lands exactly on the WGS84 ellipsoid (X²/A²+Y²/A²+Z²/B²=1) and the target */
/*   E lands on it too, so r3 = P−E has a reachable zero at the SAME root as   */
/*   the geodetic match — no Bowring in the loop, fully analytic Jacobian.     */
/*                                                                             */
/*   A-PRIORI GUARD: the FD/clamp/vertex conditioning degrades in a thin       */
/*   seam/vertex skin (validated via tools/newton_hotzones). There — and on    */
/*   any non-convergence — it falls back to the exact beam, which is clean.    */
/*   Output is bit-identical to the beam at every bin layer; only the          */
/*   non-address L29 leaf tail differs (sub-ULP representative-point noise).   */
/*                                                                             */
/*   NOTE: guard thresholds (1e-6 seam, 1e-4 vertex) are conditioning-derived  */
/*   and ELLIPSOID-DEPENDENT — re-derive (re-run newton_hotzones) if the       */
/*   reference ellipsoid changes.                                              */

/* analytic forward: nonneg barycentric (u,v,w), sum 1 → ECEF P on the         */
/* ellipsoid, plus the 3×3 raw-ak Jacobian dPdc[i][j]=∂P_i/∂coord_j (already   */
/* propagated through the normalise + A·s scaling). Mirrors h9_ak_core +       */
/* h9_ak_normalise + h9_c_oct_to_c_ell, fused with their derivatives.          */
static inline void h9_aj_ak_fwd_jac(double u, double v, double w,
                                    double su, double sv, double sw,
                                    double P[3], double dPdc[3][3]) {
    const double a = H9_ALPHA;
    const double hp = M_PI * 0.5;
    const double tu = tan(hp*u + 0.5*H9_EPS);
    const double tv = tan(hp*v + 0.5*H9_EPS);
    const double tw = tan(hp*w + 0.5*H9_EPS);
    const double u2 = tu*tu, v2 = tv*tv, w2 = tw*tw;
    const double pu = hp*(1.0+u2), pv = hp*(1.0+v2), pw = hp*(1.0+w2);   /* dt/dcoord */

    const double Sx = v2 + w2 + a*v2*w2;
    const double Sy = u2 + w2 + a*u2*w2;
    const double Sz = u2 + v2 + a*u2*v2;
    const double Px = pow(Sx, 0.25), Py = pow(Sy, 0.25), Pz = pow(Sz, 0.25);
    const double xr = tu*Px, yr = tv*Py, zr = tw*Pz;

    const double qx = 0.25*Px/Sx, qy = 0.25*Py/Sy, qz = 0.25*Pz/Sz;     /* ¼ S^{-¾} */
    dPdc[0][0] = pu*Px;
    dPdc[0][1] = tu*qx*(2.0*tv*pv*(1.0 + a*w2));
    dPdc[0][2] = tu*qx*(2.0*tw*pw*(1.0 + a*v2));
    dPdc[1][0] = tv*qy*(2.0*tu*pu*(1.0 + a*w2));
    dPdc[1][1] = pv*Py;
    dPdc[1][2] = tv*qy*(2.0*tw*pw*(1.0 + a*u2));
    dPdc[2][0] = tw*qz*(2.0*tu*pu*(1.0 + a*v2));
    dPdc[2][1] = tw*qz*(2.0*tv*pv*(1.0 + a*u2));
    dPdc[2][2] = pw*Pz;

    const double NB2 = H9_NORM_B2;
    const double n2  = xr*xr + yr*yr + zr*zr/NB2;
    const double n   = sqrt(n2);
    const double A   = H9_WGS84_A;
    P[0] = A*su*xr/n; P[1] = A*sv*yr/n; P[2] = A*sw*zr/n;
    for (int j = 0; j < 3; ++j) {
        const double dxr = dPdc[0][j], dyr = dPdc[1][j], dzr = dPdc[2][j];
        const double dn = (xr*dxr + yr*dyr + zr*dzr/NB2)/n;
        dPdc[0][j] = A*su*(dxr*n - xr*dn)/n2;
        dPdc[1][j] = A*sv*(dyr*n - yr*dn)/n2;
        dPdc[2][j] = A*sw*(dzr*n - zr*dn)/n2;
    }
}

/* analytic J3 = ∂P/∂(fx,fy) [3×2] and P at (fx,fy). */
static inline void h9_aj_jac3(double fx, double fy, int oct_mode,
                              double su, double sv, double sw,
                              double P[3], double J3[3][2]) {
    const double iH = H9.inv_H, iW = H9.inv_W, i2H = H9.inv_2H;
    double u, v, w;
    if (oct_mode == 1) {
        u = fy*iH + (1.0/3.0); v = (1.0/3.0) - fy*i2H - fx*iW; w = (1.0/3.0) - fy*i2H + fx*iW;
    } else {
        u = (1.0/3.0) - fy*iH; v = (1.0/3.0) + fy*i2H + fx*iW; w = (1.0/3.0) + fy*i2H - fx*iW;
    }
    /* affine ∂(u,v,w)/∂(fx,fy) (sum=1 ⇒ no renorm term) */
    double dcdf[3][2];
    if (oct_mode == 1) { dcdf[0][0]=0; dcdf[0][1]=iH;  dcdf[1][0]=-iW; dcdf[1][1]=-i2H; dcdf[2][0]=iW;  dcdf[2][1]=-i2H; }
    else               { dcdf[0][0]=0; dcdf[0][1]=-iH; dcdf[1][0]=iW;  dcdf[1][1]=i2H;  dcdf[2][0]=-iW; dcdf[2][1]=i2H;  }
    double dPdc[3][3]; h9_aj_ak_fwd_jac(u, v, w, su, sv, sw, P, dPdc);
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 2; ++k)
            J3[i][k] = dPdc[i][0]*dcdf[0][k] + dPdc[i][1]*dcdf[1][k] + dPdc[i][2]*dcdf[2][k];
}

static H9BOct h9_lonlat_to_boct(double lon_rad, double lat_rad) {
    H9BOct result;
    double eX, eY, eZ;
    h9_rad_lonlat_to_ecef(lon_rad, lat_rad, &eX, &eY, &eZ);
    result.oct_i    = oid(eX, eY, eZ);
    result.oct_mode = mode(result.oct_i);
    const double su = (eX >= 0.0) ? 1.0 : -1.0;
    const double sv = (eY >= 0.0) ? 1.0 : -1.0;
    const double sw = (eZ >= 0.0) ? 1.0 : -1.0;

    /* pole / axis shortcuts — tangential-distance² capture, identical to the
     * beam (see H9_BOCT_VERTEX_TAN2 note; was the quadratic |e|−A < 1e-9). */
    if (eX*eX + eY*eY < H9_BOCT_VERTEX_TAN2) { result.w = sw; result.u = 0; result.v = 0; return result; }  /* pole */
    if (eY*eY + eZ*eZ < H9_BOCT_VERTEX_TAN2) { result.u = su; result.v = 0; result.w = 0; return result; }  /* X */
    if (eX*eX + eZ*eZ < H9_BOCT_VERTEX_TAN2) { result.v = sv; result.u = 0; result.w = 0; return result; }  /* Y */

    const int    oct_mode = result.oct_mode;
    const double inv_H = H9.inv_H, inv_W = H9.inv_W, inv_2H = H9.inv_2H;

    /* gnomonic seed + a-priori seam/vertex guard (route the thin skin to beam) */
    const double a_u = fabs(eX), a_v = fabs(eY), a_w = fabs(eZ);
    const double s1  = a_u + a_v + a_w;
    const double u0 = a_u/s1, v0 = a_v/s1, w0 = a_w/s1;
    const double mn = fmin(u0, fmin(v0, w0));
    const double mx = fmax(u0, fmax(v0, w0));
    if (mn < 1e-6 || mx > 1.0 - 1e-4)
        return h9_lonlat_to_boct_beam(lon_rad, lat_rad);

    double fx, fy;
    if (oct_mode == 1) { fy = (u0 - 1.0/3.0) * H9.H; fx = (w0 - v0) * H9.W * 0.5; }
    else               { fy = (1.0/3.0 - u0) * H9.H; fx = (v0 - w0) * H9.W * 0.5; }

    /* Gauss-Newton on r3 = P − E in ECEF metres; ‖r3‖² is the veto metric.
     * TOL2 is the early-out: at 1e-16 (‖r3‖<10 nm) the solve stopped a layer
     * short of the float64 floor (the ECEF ULP is ~1.4 nm/coord ⇒ ‖r3‖ ≈ 2–3 nm).
     * Default now drives to that floor (it exits via the line-search stall when no
     * step improves r2); override -DH9_BOCT_TOL2=1e-16 to restore the legacy ~10 nm
     * early-out, or promote to a runtime fineness knob later. FAIL2 (beam route)
     * is unchanged, and the floor (~4e-18) is far below it, so no spurious routing. */
#ifndef H9_BOCT_TOL2
#define H9_BOCT_TOL2 1e-20
#endif
    const double TOL2  = H9_BOCT_TOL2;
    const double FAIL2 = 1e-14;   /* ‖r3‖ > 100 nm ⇒ genuine non-convergence */
    double P[3], J3[3][2];
    h9_aj_jac3(fx, fy, oct_mode, su, sv, sw, P, J3);
    double r2 = (P[0]-eX)*(P[0]-eX) + (P[1]-eY)*(P[1]-eY) + (P[2]-eZ)*(P[2]-eZ);
    for (int it = 0; it < 12; ++it) {
        if (r2 < TOL2) break;
        const double rx = P[0]-eX, ry = P[1]-eY, rz = P[2]-eZ;
        const double a00 = J3[0][0]*J3[0][0]+J3[1][0]*J3[1][0]+J3[2][0]*J3[2][0];
        const double a01 = J3[0][0]*J3[0][1]+J3[1][0]*J3[1][1]+J3[2][0]*J3[2][1];
        const double a11 = J3[0][1]*J3[0][1]+J3[1][1]*J3[1][1]+J3[2][1]*J3[2][1];
        const double b0  = -(J3[0][0]*rx + J3[1][0]*ry + J3[2][0]*rz);
        const double b1  = -(J3[0][1]*rx + J3[1][1]*ry + J3[2][1]*rz);
        const double det = a00*a11 - a01*a01;
        if (fabs(det) < 1e-30) break;
        const double sx = ( a11*b0 - a01*b1)/det;
        const double sy = (-a01*b0 + a00*b1)/det;
        double t = 1.0, rn2 = r2, nx = fx, ny = fy, Pn[3], Jn[3][2]; int damp = 0;
        for (damp = 0; damp < 10; ++damp) {
            nx = fx + t*sx; ny = fy + t*sy;
            h9_aj_jac3(nx, ny, oct_mode, su, sv, sw, Pn, Jn);
            rn2 = (Pn[0]-eX)*(Pn[0]-eX) + (Pn[1]-eY)*(Pn[1]-eY) + (Pn[2]-eZ)*(Pn[2]-eZ);
            if (rn2 < r2) break;
            t *= 0.5;
        }
        if (damp == 10) break;
        fx = nx; fy = ny; r2 = rn2;
        for (int i = 0; i < 3; ++i) { P[i]=Pn[i]; J3[i][0]=Jn[i][0]; J3[i][1]=Jn[i][1]; }
    }
    if (r2 > FAIL2)                                   /* true non-convergence → beam */
        return h9_lonlat_to_boct_beam(lon_rad, lat_rad);

    double u, v, w;
    if (oct_mode == 1) {
        u = fy*inv_H + (1.0/3.0); v = (1.0/3.0) - fy*inv_2H - fx*inv_W; w = (1.0/3.0) - fy*inv_2H + fx*inv_W;
    } else {
        u = (1.0/3.0) - fy*inv_H; v = (1.0/3.0) + fy*inv_2H + fx*inv_W; w = (1.0/3.0) + fy*inv_2H - fx*inv_W;
    }
    if (u < 0.0) u = 0.0;
    if (v < 0.0) v = 0.0;
    if (w < 0.0) w = 0.0;
    const double s = u + v + w;
    if (s > 1e-30) { u /= s; v /= s; w /= s; }
    result.u = su*u; result.v = sv*v; result.w = sw*w;
    return result;
}

/* ── Decode: H9BOct → (lon, lat) radians ────────────────────────────────── */
/*                                                                             */
/*   AK analytical inverse — no beam search. The signed (u,v,w) carry enough  */
/*   information for the full inverse: h9_c_oct_to_c_ell handles the signs    */
/*   internally, so no oct_i or oct_mode is needed here.                      */

static void h9_boct_to_lonlat(H9BOct b, double *lon_rad, double *lat_rad) {
    double X, Y, Z;
    h9_c_oct_to_c_ell(b.u, b.v, b.w, &X, &Y, &Z);
    h9_c_ell_to_lonlat(X, Y, Z, lon_rad, lat_rad);
}

/* ── Degree convenience wrappers ─────────────────────────────────────────── */

static H9BOct h9_lonlatdeg_to_boct(double lon_deg, double lat_deg) {
    return h9_lonlat_to_boct(lon_deg * (M_PI / 180.0), lat_deg * (M_PI / 180.0));
}

static void h9_boct_to_lonlatdeg(H9BOct b, double *lon_deg, double *lat_deg) {
    double lon_r, lat_r;
    h9_boct_to_lonlat(b, &lon_r, &lat_r);
    *lon_deg = lon_r * (180.0 / M_PI);
    *lat_deg = lat_r * (180.0 / M_PI);
}


/* Pull in the runtime warp shims (h9_warp_fwd / h9_warp_inv). The header
 * is guarded; including from h9_math.h means every legacy call site
 * keeps its current signature without source edits. */
#include "h9_warp_runtime.h"
