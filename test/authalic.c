/* authalic.c — h9_authalic.h vs the hhg9 Python oracle.
 *
 * Reference values generated from hhg9.algorithms.authalic (the Karney
 * series, itself machine-verified against the exact closed form by
 * hex9 tests/test_authalic.py). Agreement bound 3e-16 rad (~2 ULP —
 * libm vs numpy sin/cos); poles must be exact; series round-trip must
 * stay at the same bound (Python's own RT worst on this table is 0.0).
 */
#include "../core/h9_authalic.h"

#include <math.h>
#include <stdio.h>

/* Generated from hhg9.algorithms.authalic (WGS84 e2 below). */
static const double E2 = 0.0066943799901413199;
static const double REF_PHI[19] = {
    -1.57079632679489656e+00,
    -1.57079458146564455e+00,
    -1.04719755119659763e+00,
    -7.85398163397448279e-01,
    -5.23598775598298816e-01,
    -1.74532925199432954e-09,
    0.00000000000000000e+00,
    1.74532925199432954e-09,
    2.61799387799149408e-01,
    5.23598775598298816e-01,
    7.85398145944155801e-01,
    7.85398163397448279e-01,
    7.85398180850740757e-01,
    1.04719755119659763e+00,
    1.30899693899574721e+00,
    1.55334303427495324e+00,
    1.57077887350237666e+00,
    1.57079630934160419e+00,
    1.57079632679489656e+00,
};
static const double REF_XI[19] = {
    -1.57079632679489656e+00,
    -1.57079457363442465e+00,
    -1.04525649321536473e+00,
    -7.83158956118022664e-01,
    -5.21661408370140411e-01,
    -1.73752778458892934e-09,
    0.00000000000000000e+00,
    1.73752778458892934e-09,
    2.60681625702924746e-01,
    5.21661408370140411e-01,
    7.83158938664878956e-01,
    7.83158956118022664e-01,
    7.83158973571166372e-01,
    1.04525649321536473e+00,
    1.30787548613477322e+00,
    1.55326473806926191e+00,
    1.57077879519017749e+00,
    1.57079630926329195e+00,
    1.57079632679489656e+00,
};
static const int N_REF = 19;

int main(void)
{
    int fails = 0;
    H9Authalic aux;
    if (!h9_authalic_init(&aux, E2) || !aux.valid) {
        printf("FAIL: init rejected WGS84 e2\n");
        return 1;
    }

    /* 1. forward vs oracle */
    double worst_f = 0.0;
    for (int i = 0; i < N_REF; ++i) {
        const double xi = h9_geodetic_to_authalic(&aux, REF_PHI[i]);
        const double d = fabs(xi - REF_XI[i]);
        if (d > worst_f) worst_f = d;
    }
    printf("fwd vs python oracle: worst %.3e rad %s\n", worst_f,
           worst_f <= 3e-16 ? "PASS" : "FAIL");
    if (worst_f > 3e-16) ++fails;

    /* 2. inverse vs oracle */
    double worst_i = 0.0;
    for (int i = 0; i < N_REF; ++i) {
        const double phi = h9_authalic_to_geodetic(&aux, REF_XI[i]);
        const double d = fabs(phi - REF_PHI[i]);
        if (d > worst_i) worst_i = d;
    }
    printf("inv vs python oracle: worst %.3e rad %s\n", worst_i,
           worst_i <= 3e-16 ? "PASS" : "FAIL");
    if (worst_i > 3e-16) ++fails;

    /* 3. poles exact both ways */
    const double hp = 1.57079632679489656e+00;   /* pi/2 as stored */
    if (h9_geodetic_to_authalic(&aux, hp) != hp ||
        h9_authalic_to_geodetic(&aux, hp) != hp ||
        h9_geodetic_to_authalic(&aux, -hp) != -hp ||
        h9_authalic_to_geodetic(&aux, -hp) != -hp) {
        printf("FAIL: pole not exact\n");
        ++fails;
    } else {
        printf("poles exact: PASS\n");
    }

    /* 4. dense series round-trip */
    double worst_rt = 0.0;
    for (int i = 0; i <= 100000; ++i) {
        const double phi = -hp + (2.0 * hp) * (double)i / 100000.0;
        const double back =
            h9_authalic_to_geodetic(&aux, h9_geodetic_to_authalic(&aux, phi));
        const double d = fabs(back - phi);
        if (d > worst_rt) worst_rt = d;
    }
    printf("series RT dense 100k: worst %.3e rad %s\n", worst_rt,
           worst_rt <= 3e-16 ? "PASS" : "FAIL");
    if (worst_rt > 3e-16) ++fails;

    /* 5. sphere short-circuits to identity */
    H9Authalic sph;
    h9_authalic_init(&sph, 0.0);
    if (!sph.sphere || h9_geodetic_to_authalic(&sph, 0.7) != 0.7 ||
        h9_authalic_to_geodetic(&sph, -0.3) != -0.3) {
        printf("FAIL: sphere identity\n");
        ++fails;
    } else {
        printf("sphere identity: PASS\n");
    }

    /* 6. out-of-range n refused (e.g. e2=0.5 → n≈0.17) */
    H9Authalic big;
    if (h9_authalic_init(&big, 0.5) || big.valid) {
        printf("FAIL: |n|>=0.01 not refused\n");
        ++fails;
    } else {
        printf("series range gate: PASS\n");
    }

    printf(fails ? "=== authalic: FAIL (%d) ===\n" : "=== authalic: PASS ===\n",
           fails);
    return fails ? 1 : 0;
}
