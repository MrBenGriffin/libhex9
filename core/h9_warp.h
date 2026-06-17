/* h9_warp.h — public WarpState API for forward and inverse authalic warp.
 *
 * Port of hhg9/domains/octahedral_barycentric.py :: AuthalicWarp.do /
 * AuthalicWarp.undo (lines 316-475). Composed of:
 *
 *   load .h9warp        : h9_warp_io.h
 *   ghost padding       : h9_warp_mesh.h
 *   vertex gradients    : h9_ct.h :: estimate_gradients_2d_global
 *   Alfeld-CT coeffs    : h9_ct.h :: ct_coeffs_for_tri
 *   729² grid index     : h9_ct.h :: build_ct_state
 *
 * `WarpState` holds all of the above plus the runtime knobs (Newton
 * iteration count, edge-band tolerance). Build once at extension load.
 */
#ifndef H9_WARP_H
#define H9_WARP_H

#include "h9_ct.h"
#include "h9_warp_io.h"
#include "h9_warp_mesh.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace h9 {

struct WarpState {
    H9WarpData      data;
    WarpMesh        mesh;
    VertexNeighbors vn;
    CTState         ct;
    int             newton_iter   = 25;
    double          edge_tol      = 1e-7;     /* lateral-edge band */
    bool            edge_bypass   = true;     /* LATERAL_EDGE_BYPASS (legacy) */
    double          newton_h      = 1e-7;     /* finite-diff step for J */
    double          newton_safe   = 1e-12;    /* min |det| for Cramer */

    /* Lateral-edge FEATHER (F6 fix, 2026-06-11). The legacy hard bypass
     * makes the warp discontinuous at the band edge by the local warp
     * delta (~3 m near the poles): cells straddling it render as slivers
     * and, near the poles, polygons stop containing the points that
     * encode into them. The feather instead scales the delta smoothly
     * from 0 on the lateral edge to full at line-value F: on-edge points
     * remain EXACT identity (the carve-out — both octant frames agree by
     * construction, which the trained field itself does not guarantee:
     * its on-edge delta has a ~3.2e-7 cross-edge component near the
     * apex), and the warp is C1 everywhere else. Width is measured on
     * the line value y ∓ √3·x + W (2× the perpendicular distance):
     * F = 1e-4 ≈ 350 m half-zone, worst shear ~1.4 % at the poles.
     * edge_feather=false restores the legacy bypass for hhg9 parity. */
    bool   edge_feather = true;
    double feather_w    = 1e-4;

    /* F6 edge-tangent field (v3 blob, gradients shipped): the field is
     * boundary-preserving by construction and MUST run raw — feather or
     * bypass would mangle its legitimate km-scale tangential seam slide.
     * Set by finish_warp_state(); the runtime binds the edge mode to it. */
    bool field_has_grads = false;

    /* Cached constants. Mode-0 lateral edges: y = ±√3·x − √6/3. */
    double R3   = std::sqrt(3.0);
    double W_eq = std::sqrt(6.0) / 3.0;
};

/* Build mesh, gradients and CT state from an already-loaded `out.data`.
 * v3 fields ship hhg9's final per-vertex gradients (ghost orbits, x-mirror
 * symmetrisation and the edge-tangent projection already applied): use
 * them directly — interior triangles then reproduce hhg9's field
 * bit-exactly. Equator ghost vertices get the y-mirror transform of their
 * band twin's gradients (dx is even across the equator, dy odd — only
 * out-of-face queries ever see ghost triangles). v1/v2 fields keep the
 * legacy global gradient estimation. */
inline bool finish_warp_state(WarpState& out,
                              int    grad_maxiter = 2000,
                              double grad_tol     = 1e-12)
{
    out.mesh = build_warp_mesh(out.data.deltas, out.data.header.level,
                               out.data.header.mode);
    const std::size_t n     = out.mesh.padded_src.size();
    const std::size_t n_int = out.mesh.n_interior_verts;

    std::vector<double> fdx(n), fdy(n);
    for (std::size_t i = 0; i < n; ++i) {
        fdx[i] = out.mesh.padded_diff[i][0];
        fdy[i] = out.mesh.padded_diff[i][1];
    }

    std::vector<std::array<double, 2>> gdx, gdy;
    if (out.data.has_gradients()) {
        if (out.data.grad_dx.size() != n_int) return false;
        gdx.resize(n);
        gdy.resize(n);
        for (std::size_t i = 0; i < n_int; ++i) {
            gdx[i] = out.data.grad_dx[i];
            gdy[i] = out.data.grad_dy[i];
        }
        for (std::size_t k = 0; k < out.mesh.band_idx.size(); ++k) {
            const std::uint32_t t = out.mesh.band_idx[k];
            gdx[n_int + k] = { out.data.grad_dx[t][0], -out.data.grad_dx[t][1] };
            gdy[n_int + k] = { -out.data.grad_dy[t][0], out.data.grad_dy[t][1] };
        }
        out.field_has_grads = true;
    } else {
        out.vn = build_vertex_neighbors(out.mesh.padded_src.size(), out.mesh.tris);
        estimate_gradients_2d_global(out.mesh.padded_src, out.vn, fdx.data(), gdx,
                                     grad_maxiter, grad_tol);
        estimate_gradients_2d_global(out.mesh.padded_src, out.vn, fdy.data(), gdy,
                                     grad_maxiter, grad_tol);
        out.field_has_grads = false;
    }
    out.ct = build_ct_state(out.mesh, fdx.data(), fdy.data(), gdx, gdy);
    return true;
}

/* Smoothstep ramp of the lateral-edge feather at mode-0 point (cx, cy):
 * 0 exactly on either lateral edge → 1 beyond line-value feather_w. */
inline double warp_edge_scale(const WarpState& ws, double cx, double cy)
{
    if (!ws.edge_feather) return 1.0;
    double s = 1.0;
    const double f[2] = { std::fabs(cy - ws.R3 * cx + ws.W_eq),
                          std::fabs(cy + ws.R3 * cx + ws.W_eq) };
    for (int e = 0; e < 2; e++) {
        if (f[e] >= ws.feather_w) continue;
        const double t = f[e] / ws.feather_w;
        s *= t * t * (3.0 - 2.0 * t);
    }
    return s;
}

/* Build a fully-loaded WarpState from a .h9warp file. */
inline bool build_warp_state(const std::string& path,
                             WarpState&         out,
                             std::string*       err = nullptr,
                             int                grad_maxiter = 2000,
                             double             grad_tol     = 1e-12)
{
    if (!load_h9warp(path, out.data, err)) return false;
    if (!finish_warp_state(out, grad_maxiter, grad_tol)) {
        if (err) *err = "h9warp: gradient count mismatch vs mesh";
        return false;
    }
    return true;
}

/* Forward warp: b_raw → b_oct in mode-`mo`.
 *
 * mo=0 keeps y as-is; mo=1 mirrors y on entry and exit so the warp
 * mesh (built for mode-0) can be reused for the antipodal octant. */
inline void warp_do(const WarpState& ws,
                    double rx, double ry, int mo,
                    double* wx, double* wy)
{
    const double sign = (mo == 0) ? 1.0 : -1.0;
    const double cx = rx;
    const double cy = ry * sign;

    /* Lateral edges in mode-0 frame: y − R3·x + Ẇ = 0 (right) and
     * y + R3·x + Ẇ = 0 (left). Feathered (default): delta scales
     * smoothly to exact identity on the edge. Legacy hard bypass
     * (edge_feather=false; matches AuthalicWarp.do() lines 322-329)
     * kept for hhg9 parity runs. */
    if (!ws.edge_feather) {
        const bool on_right = std::fabs(cy - ws.R3 * cx + ws.W_eq) < ws.edge_tol;
        const bool on_left  = std::fabs(cy + ws.R3 * cx + ws.W_eq) < ws.edge_tol;
        if (ws.edge_bypass && (on_right || on_left)) {
            *wx = cx;
            *wy = cy * sign;
            return;
        }
    }

    const double scale = warp_edge_scale(ws, cx, cy);
    double b[3];
    const std::int32_t ti = ct_find_tri(ws.ct, cx, cy, b);
    double dx = 0.0, dy = 0.0;
    if (scale > 0.0 && ti >= 0) {
        dx = scale * ct_eval_with(ws.ct.coeffs_dx[ti], b);
        dy = scale * ct_eval_with(ws.ct.coeffs_dy[ti], b);
    }
    *wx = cx + dx;
    *wy = (cy + dy) * sign;
}

/* Inverse warp: b_oct (in mode-`mo`) → b_raw.
 *
 * Newton in 2D with finite-difference Jacobian (2×2 Cramer solve) —
 * matches AuthalicWarp.undo() lines 408-445. Seeded with identity
 * (curr = target); the warp's max magnitude (~1e-2 b_oct) is well
 * inside Newton's quadratic basin. Edge points are identity-bypassed. */
inline void warp_undo(const WarpState& ws,
                      double wx_in, double wy_in, int mo,
                      double* rx, double* ry)
{
    const double sign = (mo == 0) ? 1.0 : -1.0;
    const double tx = wx_in;
    const double ty = wy_in * sign;

    if (!ws.edge_feather) {
        const bool on_right = std::fabs(ty - ws.R3 * tx + ws.W_eq) < ws.edge_tol;
        const bool on_left  = std::fabs(ty + ws.R3 * tx + ws.W_eq) < ws.edge_tol;
        if (ws.edge_bypass && (on_right || on_left)) {
            *rx = tx;
            *ry = ty * sign;
            return;
        }
    }

    /* Feathered forward delta at (x, y) — the function Newton inverts.
     * The feather scale moves with the iterate, so it lives inside the
     * evaluation (and therefore inside the FD Jacobian). */
    auto eval_d = [&ws](double x, double y, double* dx, double* dy) {
        const double scale = warp_edge_scale(ws, x, y);
        double b[3];
        const std::int32_t ti = (scale > 0.0) ? ct_find_tri(ws.ct, x, y, b) : -1;
        *dx = (ti < 0) ? 0.0 : scale * ct_eval_with(ws.ct.coeffs_dx[ti], b);
        *dy = (ti < 0) ? 0.0 : scale * ct_eval_with(ws.ct.coeffs_dy[ti], b);
    };

    /* Identity seed. */
    double cx = tx, cy = ty;
    const double h = ws.newton_h;
    for (int it = 0; it < ws.newton_iter; ++it) {
        double dx, dy;
        eval_d(cx, cy, &dx, &dy);
        const double ex = cx + dx - tx;
        const double ey = cy + dy - ty;

        /* Converged — typically after 2-4 iterations (identity seed,
         * |delta| ≲ 1e-2, quadratic basin). 1e-15 b_raw ≈ 7 nm, below
         * any representable address distinction; checking BEFORE the
         * Jacobian saves its two FD evaluations on the final pass.
         * (encode was ~680 µs/pt with the historic fixed 25×3 CT
         * evaluations — this is the dominant encode cost.) */
        if (std::fabs(ex) + std::fabs(ey) < 1e-15) break;

        /* Finite-diff Jacobian of (curr + warp(curr)) wrt curr. */
        double dx_px, dy_px, dx_py, dy_py;
        eval_d(cx + h, cy, &dx_px, &dy_px);
        eval_d(cx, cy + h, &dx_py, &dy_py);

        const double a =  1.0 + (dx_px - dx) / h;   /* 1 + ∂dx/∂x */
        const double bj =        (dx_py - dx) / h;  /* ∂dx/∂y     */
        const double c  =        (dy_px - dy) / h;  /* ∂dy/∂x     */
        const double d  = 1.0 + (dy_py - dy) / h;   /* 1 + ∂dy/∂y */
        const double det = a * d - bj * c;

        double delta_x, delta_y;
        if (std::fabs(det) > ws.newton_safe) {
            delta_x = -(d * ex - bj * ey) / det;
            delta_y = -(a * ey - c  * ex) / det;
        } else {
            /* Singular Jacobian — fall back to simple Picard step. */
            delta_x = -ex;
            delta_y = -ey;
        }
        cx += delta_x;
        cy += delta_y;
    }

    *rx = cx;
    *ry = cy * sign;
}

} /* namespace h9 */

#endif /* H9_WARP_H */
