/* newton_invert.h — PROTOTYPE guarded 2D-Newton inversion, shared by the
 * newton_probe (global accuracy/timing) and newton_hotzones (stress zones)
 * probes. Drop-in alternative to h9_lonlat_to_boct (core/h9_math.h): same
 * pole/axis shortcuts, same forward-map kernel, same signed-barycentric tail —
 * only the BEAM=6×DEPTH=34 search is replaced by an analytic gnomonic seed +
 * guarded FD-Newton. NOT wired into the encoder; measurement only.
 *
 * Include AFTER h9_grid.h / h9_math.h (uses H9BOct, H9, h9_coct_to_lonlat,
 * h9_rad_lonlat_to_ecef, oid/mode, H9_WGS84_*).
 */
#ifndef NEWTON_INVERT_H
#define NEWTON_INVERT_H
#include <cmath>

struct NStat { int iters; int fellback; double final_resid; };

static H9BOct lonlat_to_boct_newton(double lon_rad, double lat_rad, NStat *st) {
    H9BOct result;
    double eX, eY, eZ;
    h9_rad_lonlat_to_ecef(lon_rad, lat_rad, &eX, &eY, &eZ);
    const double cl = cos(lat_rad);
    result.oct_i    = oid(eX, eY, eZ);
    result.oct_mode = mode(result.oct_i);

    const double su = (eX >= 0.0) ? 1.0 : -1.0;
    const double sv = (eY >= 0.0) ? 1.0 : -1.0;
    const double sw = (eZ >= 0.0) ? 1.0 : -1.0;

    /* pole / axis shortcuts — identical to the beam */
    const double VT = 1e-9;
    if (fabs(fabs(eZ) - H9_WGS84_B) < VT) { result.w = sw; result.u = 0; result.v = 0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }
    if (fabs(fabs(eX) - H9_WGS84_A) < VT) { result.u = su; result.v = 0; result.w = 0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }
    if (fabs(fabs(eY) - H9_WGS84_A) < VT) { result.v = sv; result.u = 0; result.w = 0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }

    const int    oct_mode = result.oct_mode;
    const double inv_H = H9.inv_H, inv_W = H9.inv_W, inv_2H = H9.inv_2H;

    /* forward map (fx,fy) → lon/lat — the ORACLE (same kernel as the beam) */
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
    auto resid = [&](double qx, double qy) {
        double l0, a0; to_lonlat(qx, qy, &l0, &a0);
        double dlat = a0 - lat_rad, dlon = l0 - lon_rad;
        if (dlon >  M_PI) dlon -= 2.0*M_PI;
        if (dlon < -M_PI) dlon += 2.0*M_PI;
        const double sl = sin(dlat*0.5), sd = sin(dlon*0.5);
        return sl*sl + cl*cos(a0)*sd*sd;          /* haversine² (same metric as beam) */
    };

    /* ── analytic gnomonic seed: octant face vertices ARE the axis dirs, so the
     * linear barycentric of the target is just |eX|,|eY|,|eZ| renormalised ── */
    const double a_u = fabs(eX), a_v = fabs(eY), a_w = fabs(eZ);
    const double s1  = a_u + a_v + a_w;
    const double u0 = a_u/s1, v0 = a_v/s1, w0 = a_w/s1;

    /* ── A-PRIORI seam/vertex guard (the design memo's "beam fallback on
     * seam-cross", made predictive) ─────────────────────────────────────────
     * The FD Jacobian is corrupted in a thin skin the hot-zone probe pins down:
     *   - SEAM (one bary coord → 0): the forward map's clamp (qu,qv,qw ≥ 0) is a
     *     kink; the ±eps stencil straddles it → garbage derivative. Failures at
     *     min coord ≲ 1.6e-8 (≈10 cm); clean by ≈1.6e-7 (1 m).
     *   - EQUATORIAL VERTEX (∂lat/∂fy → 0, det → 0): failures at max coord
     *     ≳ 1−1.6e-5 (≈100 m); clean by 1 km. (Poles are handled by the axis
     *     shortcut above and never reach here.)
     * The skin is a negligible fraction of the globe (~0.003%), so we route it
     * to the exact beam with generous margin. Both tests are needed: near a
     * vertex two coords are small, so min alone misses the 6 m–100 m band.
     *
     * NOTE (fold-in / generalisation): these thresholds are tied to THIS map's
     * conditioning (WGS84 AK forward + clamp). They are ELLIPSOID-DEPENDENT and
     * Jacobian/warp-dependent — a different reference ellipsoid (or wiring the
     * warp's own Jacobian in) shifts where the FD derivative degrades, so the
     * skin must be re-derived (re-run newton_hotzones) per ellipsoid, or the
     * guard generalised to a conditioning estimate rather than fixed constants. */
    const double mn = fmin(u0, fmin(v0, w0));
    const double mx = fmax(u0, fmax(v0, w0));
    if (mn < 1e-6 || mx > 1.0 - 1e-4) {
        H9BOct b = h9_lonlat_to_boct_beam(lon_rad, lat_rad);   /* exact, ~45 nm clean here */
        st->iters = 0; st->fellback = 1; st->final_resid = 0.0;
        return b;
    }

    double fx, fy;
    if (oct_mode == 1) { fy = (u0 - 1.0/3.0) * H9.H; fx = (w0 - v0) * H9.W * 0.5; }
    else               { fy = (1.0/3.0 - u0) * H9.H; fx = (v0 - w0) * H9.W * 0.5; }

    /* ── guarded Newton (FD Jacobian, residual-vetoed steps) ─────────────── */
    const double eps = 1e-7;
    double r = resid(fx, fy);
    int it = 0, fell = 0;
    for (it = 0; it < 12; ++it) {
        if (r < 3e-31) break;                     /* beam's tightness target */
        double lon0, lat0; to_lonlat(fx, fy, &lon0, &lat0);
        double dlon = lon_rad - lon0, dlat = lat_rad - lat0;
        if (dlon >  M_PI) dlon -= 2.0*M_PI;
        if (dlon < -M_PI) dlon += 2.0*M_PI;
        if (fabs(dlon) < 1e-15 && fabs(dlat) < 1e-15) break;

        double lp, lm, ap, am;
        to_lonlat(fx + eps, fy, &lp, &ap); to_lonlat(fx - eps, fy, &lm, &am);
        const double lon_fx = (lp - lm) / (2.0*eps), lat_fx = (ap - am) / (2.0*eps);
        to_lonlat(fx, fy + eps, &lp, &ap); to_lonlat(fx, fy - eps, &lm, &am);
        const double lon_fy = (lp - lm) / (2.0*eps), lat_fy = (ap - am) / (2.0*eps);
        const double det = lon_fx * lat_fy - lon_fy * lat_fx;
        if (fabs(det) < 1e-30) { fell = 1; break; }

        const double sx = ( lat_fy * dlon - lon_fy * dlat) / det;
        const double sy = (-lat_fx * dlon + lon_fx * dlat) / det;

        /* guard: damp until the oracle residual actually shrinks */
        double t = 1.0, rn = r, nx = fx, ny = fy; int damp = 0;
        for (damp = 0; damp < 10; ++damp) {
            nx = fx + t*sx; ny = fy + t*sy;
            rn = resid(nx, ny);
            if (rn < r) break;
            t *= 0.5;
        }
        if (damp == 10) { fell = 1; break; }      /* cannot make progress → punt */
        fx = nx; fy = ny; r = rn;
    }

    st->iters = it;
    st->final_resid = r;
    /* guard verdict: Newton "owns" this point only if it reached lattice
     * precision; otherwise the production path would fall back to the beam. */
    st->fellback = fell || (r > 1e-26);

    /* final (fx,fy) → signed barycentric (same tail as the beam) */
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
    result.u = su * u; result.v = sv * v; result.w = sw * w;
    return result;
}

#endif /* NEWTON_INVERT_H */
