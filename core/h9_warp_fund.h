/* h9_warp_fund.h — wedge (fundamental-domain) warp mesh for the v4 blob.
 *
 * Port of hhg9/domains/octahedral_barycentric.py :: AuthalicWarp._build_fold
 * (the wedge-CT build, h9_proj L6_MOBIUS_WARP.md §3.7b). The v4 sidecar
 * ships deltas + FINAL gradients for the (2 3 6) fundamental wedge only
 * (399,675 vertices at L6); this header rebuilds everything else:
 *
 *   - the wedge vertex list: tri_mesh(level, 0) filtered by the closed
 *     wedge test (same order as the exporter — no coordinates shipped);
 *   - the mirror halos: every halo image is q = F + A·(p − F) with
 *     δ → A·δ and ∇ → A·(A_row·∇) — the same float ops as Python, so
 *     halo triangle geometry matches hhg9 to the ULP. Internal mirrors
 *     (C–G at x=0, G–M median) use the D3 symmetry; the C–M seam edge
 *     reflection equals the F6 ghost construction; corner orbits close
 *     the two incident reflections into their dihedral group (order 12
 *     at C, 6 at G, 4 at M);
 *   - the triangulation: all wedge+halo points are lattice points (each
 *     halo op is a lattice symmetry) and the equilateral lattice has a
 *     UNIQUE Delaunay, so the mesh is enumerated deterministically in
 *     integer lattice coordinates — no Qhull, no cocircular ties.
 *
 * Because the shipped gradients are hhg9's final (halo-padded Bell-Sibson
 * + fold-C1 / seam-tangency projection) values, no gradient estimation
 * runs here at all: interior wedge triangles reproduce hhg9's field
 * bit-near-exactly by construction (same values + same gradients + same
 * per-triangle CT formula), exactly as the v3 face blob does.
 *
 * Evaluation folds queries into the wedge by the D3 group and unfolds
 * the displacement — see the fold branch in h9_warp.h.
 *
 * Header-only; depends on h9_warp_mesh.h (WarpMesh, build_neighbors) and
 * h9_warp_io.h (H9WarpData).
 */
#ifndef H9_WARP_FUND_H
#define H9_WARP_FUND_H

#include "h9_warp_io.h"
#include "h9_warp_mesh.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace h9 {

/* Fold geometry: wedge corners, side normals, halo reflections and the
 * D3 fold group (same op order as Python — first wedge test that passes
 * wins, so on-line ties resolve identically). All 2×2 matrices row-major
 * {m00, m01, m10, m11}; all products/normalisations mirror the numpy
 * expressions so the float values agree to the ULP. */
struct FundGeom {
    double C[2], B[2], M[2];          /* wedge corners (G = origin) */
    double u2[2], n2[2];              /* G–M median direction / normal */
    double u3[2], n3[2];              /* seam-edge direction / normal */
    double s1[4], s2[4], s3[4];       /* side reflections */
    double T[6][4];                   /* fold group [I, R120, R120ᵀ, S1,
                                         S1·R120, S1·R120ᵀ] */
};

inline FundGeom fund_geom()
{
    FundGeom g;
    /* Frame constants via the same W/H chain as h9_math.h / Python H9K. */
    const double W  = std::sqrt(2.0);
    const double H  = W * std::sqrt(3.0) * 0.5;
    const double TR = W / 2.0;
    const double VC = H / 3.0;
    const double VF = -(H / 3.0) * 2.0;
    const double r3 = std::sqrt(3.0);

    g.C[0] = 0.0; g.C[1] = VF;
    g.B[0] = TR;  g.B[1] = VC;
    g.M[0] = 0.5 * (g.C[0] + g.B[0]);
    g.M[1] = 0.5 * (g.C[1] + g.B[1]);

    const double mn = std::sqrt(g.M[0] * g.M[0] + g.M[1] * g.M[1]);
    g.u2[0] = g.M[0] / mn; g.u2[1] = g.M[1] / mn;
    g.n2[0] = -g.u2[1];    g.n2[1] = g.u2[0];
    if (g.n2[0] * g.C[0] + g.n2[1] * g.C[1] < 0.0) {
        g.n2[0] = -g.n2[0]; g.n2[1] = -g.n2[1];   /* wedge side: p·n2 ≥ 0 */
    }
    const double e3x = g.B[0] - g.C[0], e3y = g.B[1] - g.C[1];
    const double en = std::sqrt(e3x * e3x + e3y * e3y);
    g.u3[0] = e3x / en; g.u3[1] = e3y / en;
    g.n3[0] = -g.u3[1]; g.n3[1] = g.u3[0];
    if (g.n3[0] * (0.0 - g.C[0]) + g.n3[1] * (0.0 - g.C[1]) < 0.0) {
        g.n3[0] = -g.n3[0]; g.n3[1] = -g.n3[1];   /* interior: (p−C)·n3 ≥ 0 */
    }

    /* Side reflections: s1 = diag(−1, 1); s_k = 2·u_k u_kᵀ − I. */
    g.s1[0] = -1.0; g.s1[1] = 0.0; g.s1[2] = 0.0; g.s1[3] = 1.0;
    g.s2[0] = 2.0 * g.u2[0] * g.u2[0] - 1.0;
    g.s2[1] = 2.0 * g.u2[0] * g.u2[1];
    g.s2[2] = g.s2[1];
    g.s2[3] = 2.0 * g.u2[1] * g.u2[1] - 1.0;
    g.s3[0] = 2.0 * g.u3[0] * g.u3[0] - 1.0;
    g.s3[1] = 2.0 * g.u3[0] * g.u3[1];
    g.s3[2] = g.s3[1];
    g.s3[3] = 2.0 * g.u3[1] * g.u3[1] - 1.0;

    /* Fold group, same order as Python GROUP. */
    const double R120[4]  = {-0.5, -r3 / 2.0, r3 / 2.0, -0.5};
    const double R120T[4] = {-0.5, r3 / 2.0, -r3 / 2.0, -0.5};
    const double* base[3] = {nullptr, R120, R120T};
    g.T[0][0] = 1.0; g.T[0][1] = 0.0; g.T[0][2] = 0.0; g.T[0][3] = 1.0;
    for (int i = 1; i < 3; ++i)
        for (int k = 0; k < 4; ++k) g.T[i][k] = base[i][k];
    /* S1·X negates the top row of X (S1 = diag(−1, 1)). */
    for (int i = 0; i < 3; ++i) {
        g.T[3 + i][0] = -g.T[i][0];
        g.T[3 + i][1] = -g.T[i][1];
        g.T[3 + i][2] = g.T[i][2];
        g.T[3 + i][3] = g.T[i][3];
    }
    return g;
}

/* Fold (x, y) into the closed wedge cone: first group op whose image
 * passes the wedge test wins (matches Python's op order and 1e-9 tol).
 * Returns the op index; non-finite input keeps identity. */
inline int fund_fold(const FundGeom& g, double x, double y,
                     double* fx, double* fy)
{
    for (int i = 0; i < 6; ++i) {
        const double* T = g.T[i];
        const double cx = x * T[0] + y * T[1];
        const double cy = x * T[2] + y * T[3];
        if (cx >= -1e-9 && cx * g.n2[0] + cy * g.n2[1] >= -1e-9) {
            *fx = cx; *fy = cy;
            return i;
        }
    }
    *fx = x; *fy = y;
    return 0;
}

namespace detail {

/* Dihedral closure of two reflections about a common fixed point,
 * mirroring Python's BFS exactly (products q = r·e, dedup at 1e-9).
 * Emits the elements EXCLUDING identity and the two generators — the
 * side bands already cover those sectors. Max order 12 (corner C). */
inline int fund_corner_extras(const double sa[4], const double sb[4],
                              double out[12][4])
{
    double elems[16][4];
    int n_elems = 1;
    elems[0][0] = 1.0; elems[0][1] = 0.0; elems[0][2] = 0.0; elems[0][3] = 1.0;
    int frontier[16], n_frontier = 1;
    frontier[0] = 0;
    auto same = [](const double a[4], const double b[4]) {
        double w = 0.0;
        for (int k = 0; k < 4; ++k) {
            const double d = std::fabs(a[k] - b[k]);
            if (d > w) w = d;
        }
        return w < 1e-9;
    };
    while (n_frontier > 0) {
        int next[16], n_next = 0;
        for (int fi = 0; fi < n_frontier; ++fi) {
            const double* e = elems[frontier[fi]];
            for (int ri = 0; ri < 2; ++ri) {
                const double* r = (ri == 0) ? sa : sb;
                double q[4];
                q[0] = r[0] * e[0] + r[1] * e[2];
                q[1] = r[0] * e[1] + r[1] * e[3];
                q[2] = r[2] * e[0] + r[3] * e[2];
                q[3] = r[2] * e[1] + r[3] * e[3];
                bool seen = false;
                for (int k = 0; k < n_elems; ++k)
                    if (same(q, elems[k])) { seen = true; break; }
                if (!seen && n_elems < 16) {
                    for (int k = 0; k < 4; ++k) elems[n_elems][k] = q[k];
                    next[n_next++] = n_elems;
                    ++n_elems;
                }
            }
        }
        n_frontier = n_next;
        for (int k = 0; k < n_next; ++k) frontier[k] = next[k];
    }
    int n_out = 0;
    for (int k = 1; k < n_elems; ++k) {
        if (same(elems[k], sa) || same(elems[k], sb)) continue;
        for (int j = 0; j < 4; ++j) out[n_out][j] = elems[k][j];
        ++n_out;
    }
    return n_out;
}

} /* namespace detail */

/* Build the wedge+halo mesh from a v4 fund H9WarpData. On success fills
 * `m` (padded verts/deltas/tris/neighbors; n_interior_verts = wedge
 * count) and the padded per-vertex gradients, and sets `geom`. */
inline bool build_fund_mesh(const H9WarpData& data,
                            WarpMesh& m,
                            std::vector<std::array<double, 2>>& gdx,
                            std::vector<std::array<double, 2>>& gdy,
                            FundGeom& geom,
                            std::string* err = nullptr)
{
    const int level = data.header.level;
    geom = fund_geom();
    const double GHOST_BAND = 0.05;

    /* 1. wedge vertex list — tri_mesh order filtered by the wedge test
     * (the exporter's canonical order; no coordinates in the blob). */
    const LatticeMesh face = tri_mesh(level, 0);
    std::vector<std::uint32_t> wsel;
    wsel.reserve(face.verts.size() / 6 + 4096);
    for (std::uint32_t i = 0; i < face.verts.size(); ++i) {
        const double x = face.verts[i].x, y = face.verts[i].y;
        if (x >= -1e-9 && x * geom.n2[0] + y * geom.n2[1] >= -1e-9)
            wsel.push_back(i);
    }
    const std::size_t n_w = wsel.size();
    if (n_w != data.deltas.size() || n_w != data.grad_dx.size()) {
        if (err) *err = "h9warp v4: wedge vertex count mismatch ("
                        + std::to_string(n_w) + " derived vs "
                        + std::to_string(data.deltas.size()) + " shipped)";
        return false;
    }

    m.padded_src.clear();  m.padded_src.reserve(n_w * 8 / 5);
    m.padded_diff.clear(); m.padded_diff.reserve(n_w * 8 / 5);
    gdx.clear(); gdx.reserve(n_w * 8 / 5);
    gdy.clear(); gdy.reserve(n_w * 8 / 5);
    for (std::size_t k = 0; k < n_w; ++k) {
        const auto& v = face.verts[wsel[k]];
        m.padded_src.push_back({v.x, v.y});
        m.padded_diff.push_back(data.deltas[k]);
        gdx.push_back(data.grad_dx[k]);
        gdy.push_back(data.grad_dy[k]);
    }
    m.n_interior_verts = n_w;
    m.band_idx.clear();

    /* Integer lattice coords (x-pitch a/2, row-pitch H/3^(L+1)): every
     * wedge/halo point is a lattice point (halo ops are lattice
     * symmetries); the incremental map both DEDUPES halo generation (the
     * on-line orbit completion below — first op wins, matching Python's
     * position-keyed dedupe) and drives the triangulation. */
    const double W  = std::sqrt(2.0);
    const double H  = W * std::sqrt(3.0) * 0.5;
    const double VF = -(H / 3.0) * 2.0;
    double nseg = 1.0;
    for (int l = 0; l < level + 1; ++l) nseg *= 3.0;
    const double ux = (W / nseg) / 2.0;
    const double uy = H / nseg;
    std::unordered_map<std::int64_t, std::uint32_t> at;
    at.reserve(n_w * 4);
    std::vector<std::int32_t> vix, viy;
    vix.reserve(n_w * 8 / 5); viy.reserve(n_w * 8 / 5);
    auto int_key = [&](double x, double y, std::int32_t* ix, std::int32_t* iy) {
        *ix = static_cast<std::int32_t>(std::nearbyint(x / ux));
        *iy = static_cast<std::int32_t>(std::nearbyint((y - VF) / uy));
        return (static_cast<std::int64_t>(*ix) << 32) |
               static_cast<std::uint32_t>(*iy);
    };
    bool dup_seen = false;
    auto reg_last = [&]() {   /* register the most recently pushed vertex */
        std::int32_t ix, iy;
        const std::size_t k = m.padded_src.size() - 1;
        const std::int64_t key = int_key(m.padded_src[k][0],
                                         m.padded_src[k][1], &ix, &iy);
        vix.push_back(ix); viy.push_back(iy);
        if (!at.emplace(key, static_cast<std::uint32_t>(k)).second)
            dup_seen = true;
    };
    for (std::size_t k = 0; k < n_w; ++k) {
        std::int32_t ix, iy;
        const std::int64_t key = int_key(m.padded_src[k][0],
                                         m.padded_src[k][1], &ix, &iy);
        vix.push_back(ix); viy.push_back(iy);
        if (!at.emplace(key, static_cast<std::uint32_t>(k)).second)
            dup_seen = true;
    }

    /* Halo append helper: q = F + A·(p − F), δ → A·δ, ∇δ_a → A·(A_a0·∇δx
     * + A_a1·∇δy) — the uniform rule (sides A symmetric; corner ops A
     * orthogonal). Source index is a WEDGE index. */
    auto add_halo = [&](std::size_t src_i, const double A[4],
                        double Fx, double Fy) {
        const double px = m.padded_src[src_i][0] - Fx;
        const double py = m.padded_src[src_i][1] - Fy;
        m.padded_src.push_back({Fx + (px * A[0] + py * A[1]),
                                Fy + (px * A[2] + py * A[3])});
        const double dx = m.padded_diff[src_i][0];
        const double dy = m.padded_diff[src_i][1];
        m.padded_diff.push_back({dx * A[0] + dy * A[1],
                                 dx * A[2] + dy * A[3]});
        const double gxx = gdx[src_i][0], gxy = gdx[src_i][1];
        const double gyx = gdy[src_i][0], gyy = gdy[src_i][1];
        const double uxx = A[0] * gxx + A[1] * gyx;   /* A00·∇δx + A01·∇δy */
        const double uxy = A[0] * gxy + A[1] * gyy;
        const double uyx = A[2] * gxx + A[3] * gyx;   /* A10·∇δx + A11·∇δy */
        const double uyy = A[2] * gxy + A[3] * gyy;
        gdx.push_back({A[0] * uxx + A[1] * uxy, A[2] * uxx + A[3] * uxy});
        gdy.push_back({A[0] * uyx + A[1] * uyy, A[2] * uyx + A[3] * uyy});
        reg_last();
    };

    /* 2. side halos (strict bands — on-line points excluded). */
    std::vector<double> d1(n_w), d2(n_w), d3(n_w);
    for (std::size_t k = 0; k < n_w; ++k) {
        const double x = m.padded_src[k][0], y = m.padded_src[k][1];
        d1[k] = x;
        d2[k] = x * geom.n2[0] + y * geom.n2[1];
        d3[k] = (x - geom.C[0]) * geom.n3[0] + (y - geom.C[1]) * geom.n3[1];
    }
    struct Side { const double* d; const double* A; double Fx, Fy; };
    const Side sides[3] = {
        {d1.data(), geom.s1, 0.0, 0.0},
        {d2.data(), geom.s2, 0.0, 0.0},
        {d3.data(), geom.s3, geom.C[0], geom.C[1]},
    };
    for (const auto& s : sides)
        for (std::size_t k = 0; k < n_w; ++k)
            if (s.d[k] > 1e-12 && s.d[k] < GHOST_BAND)
                add_halo(k, s.A, s.Fx, s.Fy);

    /* 3. corner orbits. Strict-interior points cannot collide; the
     * ON-LINE points' orbits complete the back-sector boundary rays and
     * DO collide (stabiliser + first-generation images already in the
     * side bands), so they are position-deduped — first op wins, same
     * order and dedupe decisions as the Python build. Without the ray
     * points the in-wedge triangles at the corner would face gaps (in
     * hhg9: Qhull filler slivers) instead of proper lattice neighbours. */
    struct Corner { double Fx, Fy; const double *sa, *sb; };
    const Corner corners[3] = {
        {geom.C[0], geom.C[1], geom.s1, geom.s3},
        {0.0, 0.0, geom.s1, geom.s2},
        {geom.M[0], geom.M[1], geom.s2, geom.s3},
    };
    for (const auto& c : corners) {
        double extras[12][4];
        const int n_ex = detail::fund_corner_extras(c.sa, c.sb, extras);
        for (int e = 0; e < n_ex; ++e) {
            for (std::size_t k = 0; k < n_w; ++k) {
                if (!(d1[k] > 1e-12 && d2[k] > 1e-12 && d3[k] > 1e-12))
                    continue;
                const double rx = m.padded_src[k][0] - c.Fx;
                const double ry = m.padded_src[k][1] - c.Fy;
                if (std::sqrt(rx * rx + ry * ry) >= GHOST_BAND) continue;
                add_halo(k, extras[e], c.Fx, c.Fy);
            }
        }
        for (int e = 0; e < n_ex; ++e) {
            const double* A = extras[e];
            for (std::size_t k = 0; k < n_w; ++k) {
                if (d1[k] > 1e-12 && d2[k] > 1e-12 && d3[k] > 1e-12)
                    continue;                        /* interior handled above */
                const double rx = m.padded_src[k][0] - c.Fx;
                const double ry = m.padded_src[k][1] - c.Fy;
                const double r = std::sqrt(rx * rx + ry * ry);
                if (r <= 1e-12 || r >= GHOST_BAND) continue;
                const double qx = c.Fx + (rx * A[0] + ry * A[1]);
                const double qy = c.Fy + (rx * A[2] + ry * A[3]);
                std::int32_t ix, iy;
                const std::int64_t key = int_key(qx, qy, &ix, &iy);
                if (at.find(key) != at.end()) continue;   /* dup: first wins */
                add_halo(k, A, c.Fx, c.Fy);
            }
        }
    }
    if (dup_seen) {
        if (err) *err = "h9warp v4: duplicate lattice point in wedge+halo";
        return false;
    }

    /* 4. deterministic lattice triangulation off the incremental map.
     * Unit triangles: up {(ix,iy),(ix+2,iy),(ix+1,iy+1)} from the
     * bottom-left vertex and down {(ix,iy),(ix+2,iy),(ix+1,iy−1)} from
     * the top-left vertex — each triangle generated exactly once. */
    const std::size_t n_pad = m.padded_src.size();
    auto lookup = [&](std::int32_t ix, std::int32_t iy) -> std::int64_t {
        const std::int64_t key =
            (static_cast<std::int64_t>(ix) << 32) |
            static_cast<std::uint32_t>(iy);
        const auto it = at.find(key);
        return it == at.end() ? -1 : static_cast<std::int64_t>(it->second);
    };
    m.tris.clear();
    m.tris.reserve(n_pad * 2);
    for (std::size_t k = 0; k < n_pad; ++k) {
        const std::int32_t ix = vix[k], iy = viy[k];
        const std::int64_t r = lookup(ix + 2, iy);
        if (r < 0) continue;
        const std::int64_t up = lookup(ix + 1, iy + 1);
        if (up >= 0)
            m.tris.push_back({static_cast<std::uint32_t>(k),
                              static_cast<std::uint32_t>(r),
                              static_cast<std::uint32_t>(up)});
        const std::int64_t dn = lookup(ix + 1, iy - 1);
        if (dn >= 0)
            m.tris.push_back({static_cast<std::uint32_t>(k),
                              static_cast<std::uint32_t>(r),
                              static_cast<std::uint32_t>(dn)});
    }
    m.n_interior_tris = m.tris.size();
    detail::build_neighbors(m.tris, m.neighbors);
    return true;
}

} /* namespace h9 */

#endif /* H9_WARP_FUND_H */
