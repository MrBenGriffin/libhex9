/* newton_invert_aj.h — PROTOTYPE 2D-Newton inversion, ANALYTIC-JACOBIAN variant.
 *
 * The FD variant (newton_invert.h) spends 4 forward evals/iter, each running the
 * 5-iter Bowring ecef→geodetic. This variant instead does GAUSS-NEWTON in ECEF
 * SURFACE SPACE: the forward P(fx,fy) lands exactly on the WGS84 ellipsoid (the
 * ak normalise makes X²/A²+Y²/A²+Z²/B² = 1), and the target E = rad_lonlat_to_ecef
 * is on it too, so the residual r3 = P − E has a reachable zero at the SAME root
 * as the geodetic match — but with NO Bowring in the loop and a fully ANALYTIC
 * Jacobian (ak_core derivatives + the normalise quotient + the affine map).
 *
 * Cost: 1 analytic forward+Jacobian per iter instead of 4 FD evals × Bowring.
 * Output is root-identical to the FD variant / beam in the smooth interior; the
 * same a-priori seam/vertex guard routes the thin skin to the exact beam.
 *
 * `aj_jacobian_max_relerr()` cross-checks the analytic J3 against central
 * differences (verification for "test and verify"); call it from a probe.
 *
 * Include AFTER h9_grid.h / h9_math.h and newton_invert.h (reuses NStat).
 */
#ifndef NEWTON_INVERT_AJ_H
#define NEWTON_INVERT_AJ_H
#include <cmath>

/* analytic forward: (u,v,w) nonneg, sum 1 → ECEF P on the ellipsoid, plus the
 * 3×3 raw-ak Jacobian dR[i][j] = ∂R_i/∂coord_j (coords: 0=u,1=v,2=w). */
static inline void aj_ak_fwd_jac(double u, double v, double w,
                                 double su, double sv, double sw,
                                 double P[3], double dPdc[3][3]) {
    const double a = H9_ALPHA;
    const double hp = M_PI * 0.5;
    const double tu = std::tan(hp*u + 0.5*H9_EPS);
    const double tv = std::tan(hp*v + 0.5*H9_EPS);
    const double tw = std::tan(hp*w + 0.5*H9_EPS);
    const double u2 = tu*tu, v2 = tv*tv, w2 = tw*tw;
    const double pu = hp*(1.0+u2), pv = hp*(1.0+v2), pw = hp*(1.0+w2);   /* dt/dcoord */

    const double Sx = v2 + w2 + a*v2*w2;          /* Px = Sx^¼ */
    const double Sy = u2 + w2 + a*u2*w2;
    const double Sz = u2 + v2 + a*u2*v2;
    const double Px = std::pow(Sx, 0.25), Py = std::pow(Sy, 0.25), Pz = std::pow(Sz, 0.25);
    const double xr = tu*Px, yr = tv*Py, zr = tw*Pz;

    /* ∂(2-coord part)/∂coord via Px = ¼ S^{-¾} ∂S/∂coord, ∂t²/∂coord = 2 t p */
    const double qx = 0.25*Px/Sx, qy = 0.25*Py/Sy, qz = 0.25*Pz/Sz;   /* ¼ S^{-¾} */
    /* xr = tu·Px(v,w) */
    dPdc[0][0] = pu*Px;                               /* ∂xr/∂u */
    dPdc[0][1] = tu*qx*(2.0*tv*pv*(1.0 + a*w2));      /* ∂xr/∂v */
    dPdc[0][2] = tu*qx*(2.0*tw*pw*(1.0 + a*v2));      /* ∂xr/∂w */
    /* yr = tv·Py(u,w) */
    dPdc[1][0] = tv*qy*(2.0*tu*pu*(1.0 + a*w2));      /* ∂yr/∂u */
    dPdc[1][1] = pv*Py;                               /* ∂yr/∂v */
    dPdc[1][2] = tv*qy*(2.0*tw*pw*(1.0 + a*u2));      /* ∂yr/∂w */
    /* zr = tw·Pz(u,v) */
    dPdc[2][0] = tw*qz*(2.0*tu*pu*(1.0 + a*v2));      /* ∂zr/∂u */
    dPdc[2][1] = tw*qz*(2.0*tv*pv*(1.0 + a*u2));      /* ∂zr/∂v */
    dPdc[2][2] = pw*Pz;                               /* ∂zr/∂w */

    /* normalise R → R̂ (ellipsoid), then scale A·s. Propagate the quotient. */
    const double NB2 = H9_NORM_B2;
    const double n2  = xr*xr + yr*yr + zr*zr/NB2;
    const double n   = std::sqrt(n2);
    const double A   = H9_WGS84_A;
    P[0] = A*su*xr/n; P[1] = A*sv*yr/n; P[2] = A*sw*zr/n;

    /* ∂(R̂)/∂coord = (∂R/∂coord · n − R · ∂n/∂coord)/n², ∂n/∂coord = (R·wR ∂R/∂coord)/n
     * with weight wR = (1,1,1/NB2). Then × A·s. */
    for (int j = 0; j < 3; ++j) {
        const double dxr = dPdc[0][j], dyr = dPdc[1][j], dzr = dPdc[2][j];
        const double dn = (xr*dxr + yr*dyr + zr*dzr/NB2)/n;
        dPdc[0][j] = A*su*(dxr*n - xr*dn)/n2;
        dPdc[1][j] = A*sv*(dyr*n - yr*dn)/n2;
        dPdc[2][j] = A*sw*(dzr*n - zr*dn)/n2;
    }
}

/* affine ∂(u,v,w)/∂(fx,fy) for the given oct_mode (sum=1 ⇒ no renorm term). */
static inline void aj_affine_dcdf(int oct_mode, double dcdf[3][2]) {
    const double iH = H9.inv_H, iW = H9.inv_W, i2H = H9.inv_2H;
    if (oct_mode == 1) {
        dcdf[0][0]=0.0;   dcdf[0][1]= iH;
        dcdf[1][0]=-iW;   dcdf[1][1]=-i2H;
        dcdf[2][0]= iW;   dcdf[2][1]=-i2H;
    } else {
        dcdf[0][0]=0.0;   dcdf[0][1]=-iH;
        dcdf[1][0]= iW;   dcdf[1][1]= i2H;
        dcdf[2][0]=-iW;   dcdf[2][1]= i2H;
    }
}

/* analytic J3 = ∂P/∂(fx,fy) [3×2] at (fx,fy). */
static inline void aj_jac3(double fx, double fy, int oct_mode,
                           double su, double sv, double sw,
                           double P[3], double J3[3][2]) {
    const double iH = H9.inv_H, iW = H9.inv_W, i2H = H9.inv_2H;
    double u, v, w;
    if (oct_mode == 1) {
        u = fy*iH + (1.0/3.0); v = (1.0/3.0) - fy*i2H - fx*iW; w = (1.0/3.0) - fy*i2H + fx*iW;
    } else {
        u = (1.0/3.0) - fy*iH; v = (1.0/3.0) + fy*i2H + fx*iW; w = (1.0/3.0) + fy*i2H - fx*iW;
    }
    double dPdc[3][3]; aj_ak_fwd_jac(u, v, w, su, sv, sw, P, dPdc);
    double dcdf[3][2]; aj_affine_dcdf(oct_mode, dcdf);
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 2; ++k)
            J3[i][k] = dPdc[i][0]*dcdf[0][k] + dPdc[i][1]*dcdf[1][k] + dPdc[i][2]*dcdf[2][k];
}

/* verification: max relative error of analytic J3 vs central differences over
 * `n` random INTERIOR points. Each Jacobian COLUMN is normalised by its own 2-norm
 * (~1e6 m), so the FD roundoff on near-zero components isn't spuriously inflated
 * (dividing a ~7e-3 m FD roundoff by max(1,|fd|)=1 fakes a 1e-3 "error"). `eps` is
 * the near-optimal central-difference step for this O(1e6 m) scale. */
static inline double aj_jacobian_max_relerr(long n, uint64_t seed) {
    uint64_t s = seed; auto rnd = [&]{ s += 0x9e3779b97f4a7c15ULL; uint64_t z=s;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; z^=z>>31;
        return (double)(z>>11)/9007199254740992.0; };
    double worst = 0.0;
    const double eps = 3e-6;
    for (long i = 0; i < n; ++i) {
        const int oct_mode = (rnd() < 0.5) ? 0 : 1;
        const double su = (rnd()<0.5?1:-1), sv=(rnd()<0.5?1:-1), sw=(rnd()<0.5?1:-1);
        const double fx = (rnd()-0.5)*0.5*H9.W, fy = (rnd()-0.4)*0.5*H9.H;   /* interior */
        double P0[3], Ja[3][2]; aj_jac3(fx, fy, oct_mode, su, sv, sw, P0, Ja);
        for (int k = 0; k < 2; ++k) {
            double Pp[3], Pm[3], dum[3][2];
            aj_jac3(fx + (k==0?eps:0), fy + (k==1?eps:0), oct_mode, su, sv, sw, Pp, dum);
            aj_jac3(fx - (k==0?eps:0), fy - (k==1?eps:0), oct_mode, su, sv, sw, Pm, dum);
            double fd[3], cn = 0.0;
            for (int r = 0; r < 3; ++r) { fd[r] = (Pp[r]-Pm[r])/(2*eps); cn += fd[r]*fd[r]; }
            cn = std::sqrt(cn);                       /* column 2-norm (~1e6 m) */
            if (cn < 1e-12) continue;
            for (int r = 0; r < 3; ++r) {
                const double re = std::fabs(Ja[r][k]-fd[r])/cn;
                if (re > worst) worst = re;
            }
        }
    }
    return worst;
}

/* ── analytic-Jacobian guarded Gauss-Newton inversion ────────────────────── */
static H9BOct lonlat_to_boct_newton_aj(double lon_rad, double lat_rad, NStat *st) {
    H9BOct result;
    double eX, eY, eZ;
    h9_rad_lonlat_to_ecef(lon_rad, lat_rad, &eX, &eY, &eZ);
    result.oct_i    = oid(eX, eY, eZ);
    result.oct_mode = mode(result.oct_i);
    const double su = (eX >= 0.0) ? 1.0 : -1.0;
    const double sv = (eY >= 0.0) ? 1.0 : -1.0;
    const double sw = (eZ >= 0.0) ? 1.0 : -1.0;

    const double VT = 1e-9;
    if (fabs(fabs(eZ) - H9_WGS84_B) < VT) { result.w = sw; result.u=0; result.v=0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }
    if (fabs(fabs(eX) - H9_WGS84_A) < VT) { result.u = su; result.v=0; result.w=0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }
    if (fabs(fabs(eY) - H9_WGS84_A) < VT) { result.v = sv; result.u=0; result.w=0; st->iters=0; st->fellback=0; st->final_resid=0; return result; }

    const int oct_mode = result.oct_mode;
    const double inv_H = H9.inv_H, inv_W = H9.inv_W, inv_2H = H9.inv_2H;

    /* gnomonic seed + a-priori guard (identical policy to the FD variant) */
    const double a_u = fabs(eX), a_v = fabs(eY), a_w = fabs(eZ);
    const double s1  = a_u + a_v + a_w;
    const double u0 = a_u/s1, v0 = a_v/s1, w0 = a_w/s1;
    const double mn = fmin(u0, fmin(v0, w0));
    const double mx = fmax(u0, fmax(v0, w0));
    if (mn < 1e-6 || mx > 1.0 - 1e-4) {
        H9BOct b = h9_lonlat_to_boct_beam(lon_rad, lat_rad);
        st->iters=0; st->fellback=1; st->final_resid=0.0; return b;
    }
    double fx, fy;
    if (oct_mode == 1) { fy = (u0 - 1.0/3.0) * H9.H; fx = (w0 - v0) * H9.W * 0.5; }
    else               { fy = (1.0/3.0 - u0) * H9.H; fx = (v0 - w0) * H9.W * 0.5; }

    /* residual in ECEF metres; ‖r3‖² is the veto metric (= ground error²) */
    auto res2 = [&](const double P[3]) {
        const double dx=P[0]-eX, dy=P[1]-eY, dz=P[2]-eZ; return dx*dx+dy*dy+dz*dz;
    };
    const double TOL2  = 1e-16;   /* converged: ‖r3‖ < 10 nm (below the ~16 nm
                                     L29 roundtrip floor; quadratic GN reaches it) */
    const double FAIL2 = 1e-14;   /* genuine non-convergence: ‖r3‖ > 100 nm */
    double P[3], J3[3][2];
    aj_jac3(fx, fy, oct_mode, su, sv, sw, P, J3);
    double r2 = res2(P);
    int it = 0;
    for (it = 0; it < 12; ++it) {
        if (r2 < TOL2) break;
        const double rx = P[0]-eX, ry = P[1]-eY, rz = P[2]-eZ;
        /* normal equations: (JᵀJ) δ = −Jᵀ r   (2×2) */
        const double a00 = J3[0][0]*J3[0][0]+J3[1][0]*J3[1][0]+J3[2][0]*J3[2][0];
        const double a01 = J3[0][0]*J3[0][1]+J3[1][0]*J3[1][1]+J3[2][0]*J3[2][1];
        const double a11 = J3[0][1]*J3[0][1]+J3[1][1]*J3[1][1]+J3[2][1]*J3[2][1];
        const double b0  = -(J3[0][0]*rx + J3[1][0]*ry + J3[2][0]*rz);
        const double b1  = -(J3[0][1]*rx + J3[1][1]*ry + J3[2][1]*rz);
        const double det = a00*a11 - a01*a01;
        if (fabs(det) < 1e-30) break;
        const double sx = ( a11*b0 - a01*b1)/det;
        const double sy = (-a01*b0 + a00*b1)/det;
        /* guard: damp until the ECEF residual shrinks */
        double t = 1.0, rn2 = r2, nx = fx, ny = fy, Pn[3], Jn[3][2]; int damp = 0;
        for (damp = 0; damp < 10; ++damp) {
            nx = fx + t*sx; ny = fy + t*sy;
            aj_jac3(nx, ny, oct_mode, su, sv, sw, Pn, Jn);
            rn2 = res2(Pn);
            if (rn2 < r2) break;
            t *= 0.5;
        }
        if (damp == 10) break;        /* can't reduce further — at the floor or stuck;
                                         the FAIL2 test below decides if it's a failure */
        fx = nx; fy = ny; r2 = rn2;
        for (int i=0;i<3;++i){ P[i]=Pn[i]; J3[i][0]=Jn[i][0]; J3[i][1]=Jn[i][1]; }
    }
    st->iters = it; st->final_resid = r2;
    st->fellback = (r2 > FAIL2);      /* only a TRUE non-convergence punts to beam */

    double u, v, w;
    if (oct_mode == 1) {
        u = fy*inv_H + (1.0/3.0); v = (1.0/3.0) - fy*inv_2H - fx*inv_W; w = (1.0/3.0) - fy*inv_2H + fx*inv_W;
    } else {
        u = (1.0/3.0) - fy*inv_H; v = (1.0/3.0) + fy*inv_2H + fx*inv_W; w = (1.0/3.0) + fy*inv_2H - fx*inv_W;
    }
    if (u<0) u=0; if (v<0) v=0; if (w<0) w=0;
    const double s = u+v+w; if (s>1e-30){u/=s;v/=s;w/=s;}
    result.u = su*u; result.v = sv*v; result.w = sw*w;
    return result;
}

#endif /* NEWTON_INVERT_AJ_H */
