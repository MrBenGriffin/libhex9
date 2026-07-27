/* h9_authalic.h — geodetic ↔ authalic latitude (Karney 6th-order series).
 *
 * Port of hhg9/algorithms/authalic.py (production series path only). The
 * authalic latitude ξ maps the ellipsoid to the equal-area sphere; it is
 * the ellipsoid-specific part of the via-sphere chain and nothing else.
 *
 *     ξ − φ = Σₖ Fₖ·sin(2kφ),  k = 1..6,  Fₖ = polynomial in n (order 6)
 *
 * with n = f/(2−f) the third flattening, coefficients verbatim from PROJ
 * src/latitudes.cpp (pj_auxlat_coeffs; Karney, "On auxiliary latitudes",
 * Survey Review 56 (2024), Eqs. A19/A20), evaluated by Clenshaw
 * summation. Full float64 accuracy for |n| < 0.01 — every geodetic
 * ellipsoid. The exact closed-form fallback for |n| ≥ 0.01 is NOT ported
 * (hhg9 keeps it as the test oracle; the C library only ever runs the
 * series, and h9_authalic_init reports validity so callers can refuse).
 *
 * Operation order mirrors the Python exactly (same Horner fill, same
 * Clenshaw recurrence), so agreement is limited only by libm-vs-numpy
 * sin/cos differences: ≤ ~2 ULP in practice, poles exact.
 *
 * Header-only, C-compatible subset of C++ (matches h9_math.h style).
 */
#ifndef H9_AUTHALIC_H
#define H9_AUTHALIC_H

#include <math.h>
#include "h9_det_math.h" /* deterministic sin/cos — universality, not libm */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double F_fwd[6];   /* C[ξ,φ]: φ → ξ Fourier coefficients */
    double F_inv[6];   /* C[φ,ξ]: ξ → φ Fourier coefficients */
    int    valid;      /* 1 iff |n| < 0.01 (series in range) and e2 > 0 */
    int    sphere;     /* 1 iff e2 == 0 (identity both directions) */
} H9Authalic;

/* Taylor polynomials in n (ascending powers) whose value, times n^(l+1),
 * gives Fourier coefficient F[l] of sin(2(l+1)ζ). Row lengths 6-l. */
static const double H9_AUX_C_XI_PHI[6][6] = {
    {-4.0 / 3.0, -4.0 / 45.0, 88.0 / 315.0, 538.0 / 4725.0,
     20824.0 / 467775.0, -44732.0 / 2837835.0},
    {34.0 / 45.0, 8.0 / 105.0, -2482.0 / 14175.0, -37192.0 / 467775.0,
     -12467764.0 / 212837625.0, 0.0},
    {-1532.0 / 2835.0, -898.0 / 14175.0, 54968.0 / 467775.0,
     100320856.0 / 1915538625.0, 0.0, 0.0},
    {6007.0 / 14175.0, 24496.0 / 467775.0, -5884124.0 / 70945875.0,
     0.0, 0.0, 0.0},
    {-23356.0 / 66825.0, -839792.0 / 19348875.0, 0.0, 0.0, 0.0, 0.0},
    {570284222.0 / 1915538625.0, 0.0, 0.0, 0.0, 0.0, 0.0},
};

static const double H9_AUX_C_PHI_XI[6][6] = {
    {4.0 / 3.0, 4.0 / 45.0, -16.0 / 35.0, -2582.0 / 14175.0,
     60136.0 / 467775.0, 28112932.0 / 212837625.0},
    {46.0 / 45.0, 152.0 / 945.0, -11966.0 / 14175.0, -21016.0 / 51975.0,
     251310128.0 / 638512875.0, 0.0},
    {3044.0 / 2835.0, 3802.0 / 14175.0, -94388.0 / 66825.0,
     -8797648.0 / 10945935.0, 0.0, 0.0},
    {6059.0 / 4725.0, 41072.0 / 93555.0, -1472637812.0 / 638512875.0,
     0.0, 0.0, 0.0},
    {768272.0 / 467775.0, 455935736.0 / 638512875.0, 0.0, 0.0, 0.0, 0.0},
    {4210684958.0 / 1915538625.0, 0.0, 0.0, 0.0, 0.0, 0.0},
};

/* Fill the six Fourier coefficients from the row polynomials: Horner in
 * ascending-power storage, F[l] = n^(l+1) · poly_l(n) — mirrors
 * authalic.py::_series_coeffs::fill (row length 6-l). */
static void h9_authalic_fill(const double rows[6][6], double n, double F[6])
{
    double d = n;
    for (int l = 0; l < 6; ++l) {
        const int len = 6 - l;
        double acc = rows[l][len - 1];
        for (int c = len - 2; c >= 0; --c)
            acc = acc * n + rows[l][c];
        F[l] = d * acc;
        d *= n;
    }
}

/* Initialise for squared eccentricity e2. Returns aux.valid. */
static int h9_authalic_init(H9Authalic *aux, double e2)
{
    aux->sphere = (e2 == 0.0);
    aux->valid  = 0;
    if (aux->sphere) return 1;
    const double b_a = sqrt(1.0 - e2);
    const double n   = (1.0 - b_a) / (1.0 + b_a);   /* third flattening */
    if (fabs(n) >= 0.01) return 0;                   /* series out of range */
    h9_authalic_fill(H9_AUX_C_XI_PHI, n, aux->F_fwd);
    h9_authalic_fill(H9_AUX_C_PHI_XI, n, aux->F_inv);
    aux->valid = 1;
    return 1;
}

/* Σ F[k]·sin((2k+2)ζ) by Clenshaw summation (mirrors PROJ pj_clenshaw
 * and authalic.py::_clenshaw — identical recurrence order). */
static double h9_authalic_clenshaw(double szeta, double czeta,
                                   const double F[6])
{
    double u0 = 0.0, u1 = 0.0;
    const double X = 2.0 * (czeta - szeta) * (czeta + szeta);  /* 2cos(2ζ) */
    for (int k = 5; k >= 0; --k) {
        const double t = X * u0 - u1 + F[k];
        u1 = u0;
        u0 = t;
    }
    return 2.0 * szeta * czeta * u0;                            /* sin(2ζ)·u0 */
}

/* Geodetic latitude φ (radians) → authalic ξ (radians). */
static double h9_geodetic_to_authalic(const H9Authalic *aux, double phi)
{
    if (aux->sphere) return phi;
    return phi + h9_authalic_clenshaw(h9_sin(phi), h9_cos(phi), aux->F_fwd);
}

/* Authalic latitude ξ (radians) → geodetic φ (radians). */
static double h9_authalic_to_geodetic(const H9Authalic *aux, double xi)
{
    if (aux->sphere) return xi;
    return xi + h9_authalic_clenshaw(h9_sin(xi), h9_cos(xi), aux->F_inv);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* H9_AUTHALIC_H */
