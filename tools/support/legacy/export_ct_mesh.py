#!/usr/bin/env python3
"""
export_ct_mesh.py — export exact scipy CloughTocher2D mesh for C++ evaluation.

Extracts the Delaunay triangulation and per-triangle Bernstein coefficients
that reproduce the scipy CT cubic exactly.  Precomputes the 19 Alfeld-CT
control points for each triangle, using affine-invariant g-values derived
from neighbour centroids (g = -½ at convex-hull boundary).

g-value formula (affine-invariant, from scipy source):
    For vertex k, neighbour centroid V4' and current centroid V4:
        d = V4' - V4
        a = V4  - V2        (V2 = vertex (k+1)%3)
        b = V3  - V2        (V3 = vertex (k+2)%3)
        g[k] = (d.y*a.x - d.x*a.y) / (d.x*b.y - d.y*b.x)
    No neighbour → g[k] = -½

Inner control-point formula (from scipy _interpnd.pyx):
    c0111 = (g[0]*(-c0300 + 3*c0210 - 3*c0120 + c0030)
             + (-c0300 + 2*c0210 - c0120 + c0021 + c0201)) / 2
    c1011 = (g[1]*(-c0030 + 3*c1020 - 3*c2010 + c3000)
             + (-c0030 + 2*c1020 - c2010 + c2001 + c0021)) / 2
    c1101 = (g[2]*(-c3000 + 3*c2100 - 3*c1200 + c0300)
             + (-c3000 + 2*c2100 - c1200 + c2001 + c0201)) / 2

Output: h9_wgs84_ct_mesh.h  (same directory)
Run from hex9_cli source directory:
    python export_ct_mesh.py
"""

import os, math, datetime
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator, NearestNDInterpolator
from scipy.spatial import cKDTree

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH   = os.path.join(SCRIPT_DIR, 'WGS84_l5_warp_data.npz')
OUT_PATH   = os.path.join(SCRIPT_DIR, 'h9_wgs84_ct_mesh.h')

N     = 729
X_MIN = -math.sqrt(2.0) / 2.0
X_MAX =  math.sqrt(2.0) / 2.0
Y_MIN = -math.sqrt(6.0) / 3.0
Y_MAX =  math.sqrt(6.0) / 6.0

# ── Load NPZ ──────────────────────────────────────────────────────────────────
print(f"Loading {os.path.normpath(NPZ_PATH)} ...")
repo = np.load(NPZ_PATH, allow_pickle=True)
src  = repo['source_pts']
dst  = repo['target_pts']
diff_arr = dst - src
print(f"  {len(src)} source/target pairs")

# ── Ghost-row padding (strict band — must match octahedral_barycentric.py) ────
# Excludes y == Y_EQ: their mirror is the point itself with identical (dx, dy),
# so including them only feeds Qhull redundant duplicates that perturb its
# tie-break on cocircular equator-strip rhombi.
Y_EQ = math.sqrt(6.0) / 6.0
delta_y = Y_EQ - src[:, 1]
eq_band_mask = (delta_y > 0.0) & (delta_y < 0.05)
if np.any(eq_band_mask):
    ghost_src        = src[eq_band_mask].copy()
    ghost_src[:, 1]  = 2.0 * Y_EQ - ghost_src[:, 1]
    ghost_diff       = diff_arr[eq_band_mask].copy()
    ghost_diff[:, 1] *= -1.0
    padded_src  = np.vstack([src,  ghost_src])
    padded_diff = np.vstack([diff_arr, ghost_diff])
    print(f"  Ghost padding: {np.sum(eq_band_mask)} → {len(padded_src)} total")
else:
    padded_src  = src
    padded_diff = diff_arr

# ── Build CT interpolators ────────────────────────────────────────────────────
print("Building CloughTocher2D interpolators ...")
ct_dx = CloughTocher2DInterpolator(padded_src, padded_diff[:, 0])
ct_dy = CloughTocher2DInterpolator(padded_src, padded_diff[:, 1])

# ── Extract Delaunay mesh ──────────────────────────────────────────────────────
tri       = ct_dx.tri
pts       = tri.points                         # (npts, 2)  float64
simplices = tri.simplices.astype(np.int32)     # (ntri, 3)  int32
neighbors = tri.neighbors.astype(np.int32)     # (ntri, 3)  int32

npts = len(pts)
ntri = len(simplices)
print(f"  Mesh: {npts} vertices, {ntri} triangles")

# Values: (npts, 1) → (npts,)
val_dx = ct_dx.values.ravel().astype(np.float64)
val_dy = ct_dy.values.ravel().astype(np.float64)

# Gradients: (npts, 1, 2) → (npts, 2)
grad_dx = ct_dx.grad.reshape(npts, 2).astype(np.float64)
grad_dy = ct_dy.grad.reshape(npts, 2).astype(np.float64)

# ── Per-triangle coefficient computation ──────────────────────────────────────
def g_values(i_tri):
    """Affine-invariant g[k] for each vertex k of triangle i_tri."""
    s  = simplices[i_tri]
    P  = pts[s]                        # (3, 2)
    V4 = P.mean(axis=0)                # centroid of current triangle
    g  = np.full(3, -0.5)
    for k in range(3):
        n = neighbors[i_tri, k]
        if n < 0:
            continue                    # boundary → g[k] = -½
        V4p = pts[simplices[n]].mean(axis=0)   # centroid of neighbour
        d   = V4p - V4
        V2  = P[(k + 1) % 3]
        V3  = P[(k + 2) % 3]
        a   = V4 - V2
        b   = V3 - V2
        denom = d[0]*b[1] - d[1]*b[0]
        if abs(denom) > 1e-30:
            g[k] = (d[1]*a[0] - d[0]*a[1]) / denom
    return g


def ct_coeffs(i_tri, val, grad):
    """19 Alfeld-CT Bernstein control points for triangle i_tri.

    Coefficient order (matches C++ evaluation loop):
        c3000, c2100, c2010, c2001, c1200, c1101, c1020, c1011, c1002,
        c0300, c0210, c0201, c0120, c0111, c0102,
        c0030, c0021, c0012, c0003
    """
    s  = simplices[i_tri]
    P  = pts[s]
    v  = val[s]
    gv = grad[s]                       # (3, 2) gradient at each vertex

    # Edge vectors from each vertex
    e01 = P[1]-P[0]; e02 = P[2]-P[0]
    e10 = P[0]-P[1]; e12 = P[2]-P[1]
    e20 = P[0]-P[2]; e21 = P[1]-P[2]

    # Directional derivatives at each corner
    d01 = np.dot(gv[0], e01); d02 = np.dot(gv[0], e02)
    d10 = np.dot(gv[1], e10); d12 = np.dot(gv[1], e12)
    d20 = np.dot(gv[2], e20); d21 = np.dot(gv[2], e21)

    # Corner control points
    c3000 = v[0]; c0300 = v[1]; c0030 = v[2]

    # Edge control points
    c2100 = (d01 + 3*c3000) / 3
    c1200 = (d10 + 3*c0300) / 3
    c2010 = (d02 + 3*c3000) / 3
    c0210 = (d12 + 3*c0300) / 3
    c1020 = (d20 + 3*c0030) / 3
    c0120 = (d21 + 3*c0030) / 3

    # Trisector control points — centroid formula (no g dependency)
    c2001 = (c2100 + c2010 + c3000) / 3
    c0201 = (c1200 + c0300 + c0210) / 3
    c0021 = (c1020 + c0120 + c0030) / 3

    # Inner control points — depend on affine-invariant g-values
    g = g_values(i_tri)
    c0111 = (g[0]*(-c0300 + 3*c0210 - 3*c0120 + c0030)
             + (-c0300 + 2*c0210 - c0120 + c0021 + c0201)) / 2
    c1011 = (g[1]*(-c0030 + 3*c1020 - 3*c2010 + c3000)
             + (-c0030 + 2*c1020 - c2010 + c2001 + c0021)) / 2
    c1101 = (g[2]*(-c3000 + 3*c2100 - 3*c1200 + c0300)
             + (-c3000 + 2*c2100 - c1200 + c2001 + c0201)) / 2

    # Spine control points
    c1002 = (c1101 + c1011 + c2001) / 3
    c0102 = (c1101 + c0111 + c0201) / 3
    c0012 = (c1011 + c0111 + c0021) / 3

    # Centroid control point
    c0003 = (c1002 + c0102 + c0012) / 3

    return np.array([c3000, c2100, c2010, c2001, c1200, c1101, c1020, c1011, c1002,
                     c0300, c0210, c0201, c0120, c0111, c0102,
                     c0030, c0021, c0012, c0003])


print(f"Computing {ntri} × 19 CT coefficients ...")
coeffs_dx = np.zeros((ntri, 19), dtype=np.float64)
coeffs_dy = np.zeros((ntri, 19), dtype=np.float64)
for i in range(ntri):
    if i % 2000 == 0:
        print(f"  {i}/{ntri}")
    coeffs_dx[i] = ct_coeffs(i, val_dx, grad_dx)
    coeffs_dy[i] = ct_coeffs(i, val_dy, grad_dy)
print(f"  done.")

# ── Verify coefficients against scipy CT evaluator ────────────────────────────
def eval_poly(c, P, x, y):
    """Evaluate precomputed 19-coeff CT polynomial at (x,y) inside triangle P."""
    T  = np.array([[P[0,0]-P[2,0], P[1,0]-P[2,0]],
                   [P[0,1]-P[2,1], P[1,1]-P[2,1]]])
    b01 = np.linalg.solve(T, np.array([x - P[2,0], y - P[2,1]]))
    b0, b1, b2 = b01[0], b01[1], 1.0 - b01[0] - b01[1]
    mn  = min(b0, b1, b2)
    b1s, b2s, b3s, b4s = b0-mn, b1-mn, b2-mn, 3*mn
    c3000,c2100,c2010,c2001,c1200,c1101,c1020,c1011,c1002, \
    c0300,c0210,c0201,c0120,c0111,c0102, \
    c0030,c0021,c0012,c0003 = c
    return (b1s**3*c3000 + 3*b1s**2*b2s*c2100 + 3*b1s**2*b3s*c2010 +
            3*b1s**2*b4s*c2001 + 3*b1s*b2s**2*c1200 +
            6*b1s*b2s*b4s*c1101 + 3*b1s*b3s**2*c1020 + 6*b1s*b3s*b4s*c1011 +
            3*b1s*b4s**2*c1002 + b2s**3*c0300 + 3*b2s**2*b3s*c0210 +
            3*b2s**2*b4s*c0201 + 3*b2s*b3s**2*c0120 + 6*b2s*b3s*b4s*c0111 +
            3*b2s*b4s**2*c0102 + b3s**3*c0030 + 3*b3s**2*b4s*c0021 +
            3*b3s*b4s**2*c0012 + b4s**3*c0003)

print("Verifying against scipy CT evaluator (200 random interior points) ...")
rng = np.random.default_rng(42)
cands = rng.uniform([X_MIN, Y_MIN], [X_MAX, Y_MAX], (4000, 2))
cand_tri = tri.find_simplex(cands)
interior = cands[cand_tri >= 0][:200]
int_tri  = tri.find_simplex(interior)

scipy_dx_v = ct_dx(interior)
scipy_dy_v = ct_dy(interior)

errs_dx = []
errs_dy = []
for i, pt in enumerate(interior):
    ti = int_tri[i]
    P  = pts[simplices[ti]]
    my_dx = eval_poly(coeffs_dx[ti], P, pt[0], pt[1])
    my_dy = eval_poly(coeffs_dy[ti], P, pt[0], pt[1])
    errs_dx.append(abs(my_dx - scipy_dx_v[i]))
    errs_dy.append(abs(my_dy - scipy_dy_v[i]))

max_err_dx = max(errs_dx)
max_err_dy = max(errs_dy)
mean_err_dx = sum(errs_dx) / len(errs_dx)
mean_err_dy = sum(errs_dy) / len(errs_dy)
print(f"  dx: max={max_err_dx:.2e}  mean={mean_err_dx:.2e}")
print(f"  dy: max={max_err_dy:.2e}  mean={mean_err_dy:.2e}")
if max(max_err_dx, max_err_dy) > 1e-10:
    print("  WARNING: coefficient error larger than expected!")
else:
    print("  OK: machine-precision match with scipy CT")

# ── Build N×N grid-to-triangle index ──────────────────────────────────────────
print(f"\nBuilding {N}×{N} grid-to-triangle acceleration index ...")
xs = np.linspace(X_MIN, X_MAX, N)
ys = np.linspace(Y_MIN, Y_MAX, N)
gx_grid, gy_grid = np.meshgrid(xs, ys)
grid_pts = np.stack([gx_grid.ravel(), gy_grid.ravel()], axis=1)  # (N*N, 2)
grid_tri_idx = tri.find_simplex(grid_pts).astype(np.int32)

n_outside = (grid_tri_idx < 0).sum()
print(f"  Outside hull: {n_outside}")
if n_outside > 0:
    centroids = pts[simplices].mean(axis=1)
    tree = cKDTree(centroids)
    outside_mask = grid_tri_idx < 0
    _, nearest = tree.query(grid_pts[outside_mask])
    grid_tri_idx[outside_mask] = nearest.astype(np.int32)
    print(f"  Filled with nearest-centroid triangle")

# ── Write header ───────────────────────────────────────────────────────────────
def fmt_dbl_arr(arr, name, per_row=8):
    flat = arr.ravel()
    rows = ["    " + ", ".join(f"{v:.17g}" for v in flat[i:i+per_row])
            for i in range(0, len(flat), per_row)]
    return f"static const double {name}[{len(flat)}] = {{\n" + ",\n".join(rows) + "\n};"

def fmt_int_arr(arr, name, per_row=16):
    flat = arr.ravel().astype(np.int32)
    rows = ["    " + ", ".join(str(v) for v in flat[i:i+per_row])
            for i in range(0, len(flat), per_row)]
    return f"static const int {name}[{len(flat)}] = {{\n" + ",\n".join(rows) + "\n};"

print(f"\nWriting {os.path.normpath(OUT_PATH)} ...")
hdr = f"""\
/* h9_wgs84_ct_mesh.h — AUTO-GENERATED by export_ct_mesh.py — DO NOT EDIT
 *
 * Generated : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}
 * Source    : WGS84_l5_warp_data.npz
 * Mesh      : {npts} vertices, {ntri} triangles (Delaunay of padded source pts)
 * Grid      : H9_CT_GRID_N = {N}  (O(1) triangle lookup)
 *
 * Exact scipy CloughTocher2DInterpolator mesh for h9_ct_eval() in h9_math.h.
 * 19 Alfeld-CT Bernstein coefficients per triangle, with affine-invariant
 * g-values from neighbour centroids (g = -½ at convex-hull boundary).
 *
 * Verification vs scipy: dx_max={max_err_dx:.2e}  dy_max={max_err_dy:.2e}
 */
#ifndef H9_CT_MESH_H
#define H9_CT_MESH_H

#define H9_CT_NPTS     {npts}
#define H9_CT_NTRI     {ntri}
#define H9_CT_GRID_N   {N}
#define H9_CT_XMIN     ({X_MIN:.17g})
#define H9_CT_XMAX     ({X_MAX:.17g})
#define H9_CT_YMIN     ({Y_MIN:.17g})
#define H9_CT_YMAX     ({Y_MAX:.17g})

/* Delaunay vertex positions — x0,y0, x1,y1, ...  [H9_CT_NPTS*2] */
{fmt_dbl_arr(pts, 'h9_ct_pts', 6)}

/* Delaunay simplices — vertex indices per triangle  [H9_CT_NTRI*3] */
{fmt_int_arr(simplices, 'h9_ct_simp', 12)}

/* Delaunay neighbours — triangle opposite vertex k (-1=boundary)  [H9_CT_NTRI*3] */
{fmt_int_arr(neighbors, 'h9_ct_neigh', 12)}

/* CT Bernstein coefficients for warp_dx  [H9_CT_NTRI*19]
 * Per triangle: c3000,c2100,c2010,c2001,c1200,c1101,c1020,c1011,c1002,
 *               c0300,c0210,c0201,c0120,c0111,c0102,c0030,c0021,c0012,c0003 */
{fmt_dbl_arr(coeffs_dx, 'h9_ct_dx_c', 4)}

/* CT Bernstein coefficients for warp_dy  [H9_CT_NTRI*19] */
{fmt_dbl_arr(coeffs_dy, 'h9_ct_dy_c', 4)}

/* Grid-to-triangle acceleration index  [H9_CT_GRID_N * H9_CT_GRID_N]
 * Index: iy*H9_CT_GRID_N + ix → starting Delaunay triangle for walk.
 * Lay-out: row-major, iy=0 → y=Y_MIN, ix=0 → x=X_MIN. */
{fmt_int_arr(grid_tri_idx, 'h9_ct_grid_tri', 20)}

#endif /* H9_CT_MESH_H */
"""
with open(OUT_PATH, 'w') as f:
    f.write(hdr)

size_mb = os.path.getsize(OUT_PATH) / (1024 * 1024)
print(f"Done.  {size_mb:.1f} MB")
