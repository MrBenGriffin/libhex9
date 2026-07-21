# F6 C-side port brief — edge-tangent warp field + tangent interpolant
## What happened on the Python side (hhg9)

Three changes, all deployed in `hex9`:

1. **Retrained field** — `hhg9/data/WGS84_l5_warp_data.npz` (266,815
   vertices, L5, float64, max |delta| 9.905e-3 units, training area MAE
   0.00066). Unlike the old field, lateral-edge vertices are NOT pinned:
   they **slide tangentially along the seam** (up to 8.35 km), with edge-
   normal component exactly zero and exact x-mirror symmetry
   (dx antisymmetric, dy symmetric). Equator vertices keep dy = 0.
2. **Interpolant construction** (reference implementation:
   `domains/octahedral_barycentric.py`, `AuthalicWarp.__init__`,
   blocks labelled "LATERAL GHOST PADDING", "CORNER ORBIT PADDING",
   "EDGE-TANGENT GRADIENT PROJECTION") — see §Construction below.
3. **`raw` is now the hhg9 default edge mode.** With this field the
   feather is actively WRONG: it would compress km-scale legitimate
   tangential slide into the 350 m ramp (massive shear). Bypass likewise.
   Both are retained in hhg9 only for old-field comparison.

Python validation (`experimental/sinkhorn/validate_f6.py`, raw mode):
on-edge CT normal component **2e-16 units** (brief asked ~1e-9);
cross-frame seam agreement 2.3e-16; equator dy exactly 0; round-trip
0.4 nm; L22 Greenwich and L15 lat ±89.1 hex meshes fully regular
(side-length max/min 1.36 incl. seam-crossing cells — no slivers).

## Why the data alone is not enough (the key insight)
Tangent vertex deltas do NOT give a tangent interpolant: the CT/ACT
gradient estimation couples the steep near-corner tangential slide
(0 → 2.7 km in the first cell off the apex) into an edge-normal bulge
BETWEEN on-edge vertices (~135 m worst, near the apex). Ghost padding
helps but plateaus ~1e-5 units (the global gradient solve always sees
the asymmetric boundary of any finite pad).

The exact fix is structural: **CT restricted to a mesh edge is the cubic
Hermite of the endpoint values and the ALONG-EDGE gradient components
only.** On-edge values are exactly tangent, so forcing the along-edge
derivative of the delta's edge-normal component to zero at every
boundary vertex makes the on-edge normal identically zero, for any
interior gradients. C1 continuity is preserved (CT is C1 for any choice
of vertex gradients).

## Selected approach — ship the gradients (when possible)
Skip Bell & Sibson (Bell's boundary constraints with Sibson's gradient 
estimations) for this field entirely. Extend `.h9warp` (or add a
sidecar) to carry the **final per-vertex gradients** computed by hhg9:
per vertex, delta (2×f64) + ∇dx (2×f64) + ∇dy (2×f64) = 48 B/vertex
(~12.8 MB total). ACT evaluation already consumes (value, gradient) per
vertex — just load instead of estimate.

Why: bit-exact field parity with hhg9 **by construction** (same values,
same gradients, same per-triangle CT formula ⇒ same field), no need to
port ghost orbits or the projection, faster load. The Delaunay interior
of the padded scattered set coincides with the fixed L5 lattice mesh, so
in-triangle queries never see Python's ghost triangles.

Exporter (run from the hhg9 repo root; the loader has already done the
padding + symmetrisation + projection by the time `grad` is read, and now
exposes the interpolators as `self.ct_dx` / `self.ct_dy`):

```python
import numpy as np
from hhg9.domains.octahedral_barycentric import AuthalicWarp
w = AuthalicWarp('hhg9/data/WGS84_l5_warp_data.npz')  # H9_WARP_EDGE irrelevant here
n = len(w.src)
# scipy stores grad as (n_padded, 1, 2); the first n entries are the
# real vertices, in source_pts order (ghosts follow).
np.savez('WGS84_l5_warp_grads.npz',
         source_pts=w.src, target_pts=w.dst,
         grad_dx=np.asarray(w.ct_dx.grad)[:n, 0, :],
         grad_dy=np.asarray(w.ct_dy.grad)[:n, 0, :])
```

Then regenerate the mirrored 12-octant blob from delta+gradients
(gradients transform with the same mirror matrices as the deltas).

## Fallback — port the construction (please avoid if possible).

If you must keep Bell & Sibson, reproduce all three steps **in order**:

1. **Ghost padding** (improves near-edge gradient quality; needed for
   value parity with hhg9 near seams, not for exactness):
   - Equator band (already exists): reflect points with
     0 < (√6/6 − y) < 0.05 across y = √6/6, negate the dy of the ghost.
   - Lateral bands: for each edge line `sgn·√3·x − y − Ẇ = 0`
     (sgn = ±1, Ẇ = √6/3), unit normal n̂ = (sgn·√3, −1)/2, signed
     distance s = p·n̂ − Ẇ/2 (negative inside). Reflect points with
     −0.05 < s < −1e-12: p' = p − 2s·n̂, d' = d − 2(d·n̂)n̂.
     **Strict band** — never include on-edge points (self-duplicates
     break the triangulation).
   - Corner orbits: at each corner the two adjacent edge mirrors are 60°
     apart and generate a dihedral group (3 mirror lines + ±120°
     rotations). For strictly-interior points within 0.05 of a corner,
     add rotation(±120° about the corner) images and the third-mirror
     image. **Pitfall:** the third mirror line is `u_a` rotated +120°
     ONLY when the other edge sits at +60° from `u_a`. Correct choices
     (line directions): apex (0,−√6/3): the 60° line (towards the
     equator-R corner); equator-R (√2/2, √6/6): the 0° line (equator
     direction); equator-L: the 120° line (towards apex). Get this wrong
     and you duplicate the equator ghosts (degenerate points) and leave
     a wedge unfilled — I hit exactly this.
2. **X-mirror gradient symmetrisation** (real vertices only):
   with mirror partner m(i) = vertex at (−x, y):
   ∂dx/∂x ← avg even, ∂dx/∂y ← avg odd, ∂dy/∂x ← avg odd,
   ∂dy/∂y ← avg even. (odd: ½(g_i − g_m); even: ½(g_i + g_m).)
3. **Edge-tangent projection** at boundary vertices (after 2):
   per boundary line with unit normal n̂, unit tangent t̂
   (right: n̂=(√3,−1)/2, t̂=(1,√3)/2; left: n̂=(−√3,−1)/2, t̂=(1,−√3)/2;
   equator: n̂=(0,1), t̂=(1,0)):
   `s = n̂ₓ·(t̂·∇dx) + n̂_y·(t̂·∇dy);  ∇dx −= s·n̂ₓ·t̂;  ∇dy −= s·n̂_y·t̂`
   This subtracts s·t̂ from ∇f_n and leaves ∇f_t untouched.
   The three corners lie on TWO lines: apply a joint null-space
   projection of both constraints on (∇dx, ∇dy) ∈ R⁴
   (rows a_e = [n̂ₓt̂ₓ, n̂ₓt̂_y, n̂_yt̂ₓ, n̂_yt̂_y]; G −= Aᵀ(AAᵀ)⁻¹AG).

## Inference-side changes

- Default `H9_WARP_EDGE=raw`. Feather/bypass must never be applied to
  the new field (document, or hard-bind mode to field version).
- Newton inverse on/near edges is fine in raw with this field (it is
  boundary-preserving); hhg9 measures 0.4 nm round-trip including
  exact-on-edge vertices.

## Blob conversion spot-checks (verify before anything else)

Right lateral edge, apex end, tangential slide (t̂=(1,√3)/2 direction,
1 unit = 7.076e6 m): y=−0.81482 → +2746.9 m; y=−0.81314 → +4808.1 m;
y=−0.81146 → +5498.0 m; y=−0.80978 → +4474.9 m. Equator end:
y=+0.40321 → −5347.7 m. Normal components ≤ 1e-16 units everywhere.
730 on-edge vertices per lateral edge.

## Acceptance (all in raw mode)

1. `seam_transect 16 89.1` → zero slivers; same at the L22 Greenwich
   band edges (x = ±0.434 m was the old failure).
2. `seam_profile` → on-edge delta tangential ≤ 1e-9 (hhg9: 2e-16).
3. `ctest` 9/9, `gc_kring 5000` green.
4. Cross-check field values against hhg9 at shared probe points
   (hhg9 oracle: `experimental/sinkhorn/validate_f6.py`).

## Blast radius
The field changed **globally**, not just near seams: ~4.9 km shift 1 km
off the seam meridians (intended slide), but also Westminster 169 m,
Edinburgh 47 m, mid-octant 81 m. Fine-level addresses change everywhere;
anchor fixtures (Westminster, Edinburgh, 1255-cell) need re-baselining.

## Pointers
- Reference implementation + inline rationale:
  `hhg9/domains/octahedral_barycentric.py` (AuthalicWarp.__init__).
- Python oracle: `experimental/sinkhorn/validate_f6.py`.
- Retrain provenance: `experimental/sinkhorn/retrain_f6.sh`,
  `output/retrain_f6.log`, stage-3 best
  `output/stage3/l05_WGS84_best.npz`; old field archived at
  `experimental/sinkhorn/archive/pre_f6_retrain/`.
- Now-obsolete in hhg9 (do not copy semantics from):
  `experimental/sinkhorn/warp.py`, `test_warp.py`, `config.edge_blend`
  — these encode the abandoned pinned-edge design.
