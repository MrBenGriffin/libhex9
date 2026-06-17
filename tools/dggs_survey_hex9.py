#!/usr/bin/env python3
"""DGGS aspect-ratio survey, extended with hex9.

Mirrors ajfriend/skar_py notebooks/dggs_survey.ipynb: for N random cells at the
finest resolution of each available DGGS, compute the tightest enclosing-cone
aspect ratio with `skar`, then report the per-system distribution and the
best/worst cell.

hex9 is the point of this script and is always surveyed (finest = layer 29).
H3 / S2 / A5 are surveyed too *if* importable, so the numbers sit alongside the
notebook's three reference systems. `skar` is required (it is the solver).

Run from the repo root after building the python ext:
    PYTHONPATH=build python3 tools/dggs_survey_hex9.py
"""
import sys
import time

import numpy as np

import skar  # required: the enclosing-cone aspect-ratio solver

try:
    import hex9_ext as hex9
except ImportError:
    sys.exit("hex9_ext not importable — build it and run with PYTHONPATH=build")

# CRITICAL: enable the authalic warp. It is OFF by default; without it, cell()
# returns the RAW (un-warped) lattice — the very thing the warp exists to correct.
# Measuring raw cells gives a bogus area CV of ~7.4% (max/min ~1.4); with the warp
# on, hex9 is equal-area to MAE ~0.001% (area CV ~0.12%). Aspect ratio is
# warp-independent (~1.37 either way), but area is meaningless without it.
hex9.warp_init()
hex9.set_use_warp(True)

# Optional reference systems (mirror the notebook when present).
try:
    import h3
except ImportError:
    h3 = None
try:
    import s2sphere
except ImportError:
    s2sphere = None
try:
    import a5_fast as a5
except ImportError:
    a5 = None

# ── config ────────────────────────────────────────────────────────────────
N = 5000
SEED = 0xC0FFEE
GAP_TOL = 1e-3        # sub-metre cells floor skar's strict 1e-6 duality gap

EARTH_R_M = 6_371_008.8  # mean Earth radius — scale unit-sphere areas to metres²

# Per-system finest *numerically clean* survey resolution. Both metrics are
# resolution-invariant, so these are chosen to sit clear of float-precision
# floors: H3/S2 at their maxima; hex9 at L15 (its L29 leaf is ~60 nm, vertices
# collapse in float64); A5 at r20 (at r30 its boundary coords quantize and floor
# the area to a bogus 72% CV — aspect ratio survives, but area does not).
# A whole-survey coarsening offset (argv[1], default 0) drops every system by N
# levels — used to confirm the rankings are transitive across resolution.
H3_RES = 15          # h3 supports 0..15
S2_LEVEL = 30        # s2sphere supports 0..30
A5_RES = 20          # a5 supports 0..30 (r30 floors area; see above)
HEX9_LAYER = 15      # hex9 supports 0..29 (L29 unsurveyable; see above)

SYS_LABEL = {}       # filled from the resolutions in main()
PCTL = (50, 90, 95)  # authalicity percentiles, reported as ratio to mean area


def sample_uniform_lonlat(n, rng):
    """Uniform-on-sphere samples as (lon, lat) degrees, shape (n, 2)."""
    lon = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lon, lat])


def cell_area_m2(verts):
    """Area (m²) of a cell from its vec3 boundary (unit vectors, shape (n,3)).

    Uses the planar 3D-polygon area 0.5·‖Σ wᵢ×wᵢ₊₁‖ rather than spherical excess:
    at DGGS finest resolutions a cell spans micro-radians, where Girard's
    (Σangles − (n−2)π) loses all precision to cancellation. For such cells the
    chord polygon area matches the spherical area to ~1e-15 relative.

    The vertices wᵢ are taken relative to the cell centroid before the cross
    products. This is mathematically identical for a closed polygon but vital
    numerically: the absolute vectors are O(1) unit vectors, so for a sub-metre
    cell their cross products cancel in the ~15th digit and the area floors to
    noise (S2 L30, ~mm², lost ~9% of cells this way). Centred, the edge vectors
    are O(1e-9) and the area is computed at full relative precision. Scaled by R².
    """
    v = np.asarray(verts, float)
    w = v - v.mean(axis=0)
    area_vec = np.cross(w, np.roll(w, -1, axis=0)).sum(axis=0)
    return 0.5 * np.linalg.norm(area_vec) * EARTH_R_M ** 2


try:
    from pyproj import Geod
    _GEOD = Geod(ellps='WGS84')          # geographiclib-backed geodesic areas
except ImportError:
    _GEOD = None


def cell_area_wgs84(verts):
    """True WGS84 *ellipsoidal* area (m²) of a cell from its vec3 boundary.

    skar's vec3 and the DGGS lon/lat are spherical, so cell_area_m2 measures on a
    sphere. But the grids' I/O is WGS84 geodetic lon/lat, and some are equal-area
    on the *ellipsoid* (A5: CV 0.40% on the sphere collapses to 0.01% here),
    while others target the sphere (hex9: unchanged). We recover lon/lat from the
    vec3 — exact, since to_vec3 treated the input lat/lng as spherical — and take
    geographiclib's geodesic polygon area, the physically meaningful real-Earth
    equal-area test and the common datum across all four systems.
    """
    if _GEOD is None:
        return cell_area_m2(verts)       # graceful fallback: sphere
    v = np.asarray(verts, float)
    lon = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(v[:, 2], -1.0, 1.0)))
    area, _ = _GEOD.polygon_area_perimeter(lon, lat)
    return abs(area)


# ── per-system cell streams: yield (cell_id, boundary-as-vec3) ─────────────
def iter_h3(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = h3.latlng_to_cell(float(lat), float(lon), H3_RES)
        if cid in seen:
            continue
        seen.add(cid)
        yield cid, skar.to_vec3(h3.cell_to_boundary(cid), geo='latlng')


def iter_s2(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = s2sphere.CellId.from_lat_lng(
            s2sphere.LatLng.from_degrees(float(lat), float(lon))).parent(S2_LEVEL)
        if cid.id() in seen:
            continue
        seen.add(cid.id())
        cell = s2sphere.Cell(cid)
        verts = []
        for k in range(4):
            ll = s2sphere.LatLng.from_point(cell.get_vertex(k))
            verts.append((ll.lat().degrees, ll.lng().degrees))
        yield cid.id(), skar.to_vec3(verts, geo='latlng_deg')


def iter_a5(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = a5.lonlat_to_cell(float(lon), float(lat), A5_RES)
        if cid in seen:
            continue
        seen.add(cid)
        ring = a5.cell_to_boundary(cid)            # closed ring of (lon, lat)
        if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
            ring = ring[:-1]
        latlng = [(lat_, lon_) for lon_, lat_ in ring]
        yield a5.u64_to_hex(cid), skar.to_vec3(latlng, geo='latlng_deg')


def iter_hex9(n, seed):
    """hex9 leaf cells: encode points to L29, dedupe, yield the hex boundary.

    `hex9.cell(uuid, layer, 0)` returns a closed 7-point (lon, lat) ring; drop
    the duplicate closing vertex and re-order to (lat, lon) like the a5 stream.
    """
    rng = np.random.default_rng(seed)
    pts = sample_uniform_lonlat(n, rng)            # (n, 2) lon, lat
    uuids = hex9.encode(pts[:, 0].copy(), pts[:, 1].copy())  # (n, 16) L29 cells
    seen = set()
    for u in uuids:
        cid = u.tobytes()
        if cid in seen:
            continue
        seen.add(cid)
        ring = np.asarray(hex9.cell(u, HEX9_LAYER, 0))       # (7, 2) lon, lat
        if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
            ring = ring[:-1]
        latlng = [(float(la), float(lo)) for lo, la in ring]
        yield cid.hex(), skar.to_vec3(latlng, geo='latlng_deg')


ITERATORS = {'h3': iter_h3, 's2': iter_s2, 'a5': iter_a5, 'hex9': iter_hex9}
AVAILABLE = (['h3'] if h3 else []) + (['s2'] if s2sphere else []) + \
            (['a5'] if a5 else []) + ['hex9']


def run_system(name):
    """Stream every cell, record its area, solve aspect ratio + the extremes.

    Area is geometry-only (no skar), so it is kept for every cell even when the
    aspect-ratio solve does not converge.
    """
    ars, ell, sph, dnc, best, worst = [], [], [], 0, None, None
    for cid, verts in ITERATORS[name](N, SEED):
        ell.append(cell_area_wgs84(verts))   # primary: real-Earth ellipsoidal
        sph.append(cell_area_m2(verts))       # secondary: spherical, for compare
        r = skar.solve(verts, geo='vec3', gap_tol=GAP_TOL)
        if not isinstance(r, skar.Converged):
            dnc += 1
            continue
        ar = r.aspect_ratio
        ars.append(ar)
        if best is None or ar < best[1]:
            best = (cid, ar)
        if worst is None or ar > worst[1]:
            worst = (cid, ar)
    return np.array(ars), np.array(ell), np.array(sph), dnc, best, worst


def main():
    global H3_RES, S2_LEVEL, A5_RES, HEX9_LAYER
    delta = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # coarsen every system
    H3_RES -= delta
    S2_LEVEL -= delta
    A5_RES -= delta
    HEX9_LAYER -= delta
    SYS_LABEL.update({'h3': f'H3 r{H3_RES}', 's2': f'S2 L{S2_LEVEL}',
                      'a5': f'A5 r{A5_RES}', 'hex9': f'hex9 L{HEX9_LAYER}'})
    print(f"N={N}  seed={hex(SEED)}  gap_tol={GAP_TOL}  coarsen={delta}")
    print(f"systems: {', '.join(SYS_LABEL[s] for s in AVAILABLE)}\n")

    res = {}
    for name in AVAILABLE:
        t0 = time.time()
        res[name] = run_system(name) + (time.time() - t0,)

    # ── aspect ratio: cell circularity (1.0 = perfect circle; lower better) ──
    print("ASPECT RATIO  (enclosing-cone circularity; 1.0 = circle, lower better)")
    print(f"{'system':<10}{'cells':>7}{'dnc':>6}{'mean':>9}{'median':>9}"
          f"{'p99':>9}{'max':>9}{'sec':>8}")
    for name in AVAILABLE:
        ars, ell, sph, dnc, best, worst, dt = res[name]
        if len(ars) == 0:
            print(f"{SYS_LABEL[name]:<10}{0:>7}{dnc:>6}  (no converged cells)")
            continue
        print(f"{SYS_LABEL[name]:<10}{len(ars):>7}{dnc:>6}"
              f"{ars.mean():>9.4f}{np.median(ars):>9.4f}"
              f"{np.percentile(ars, 99):>9.4f}{ars.max():>9.4f}{dt:>8.2f}")
        print(f"           best  {best[1]:.4f}  {best[0]}")
        print(f"           worst {worst[1]:.4f}  {worst[0]}")

    # ── authalicity: cell-area uniformity. CV% & max/min plus the area
    #    distribution at p50/p90/p95, all as a ratio to the mean so the columns
    #    are comparable across systems with very different cell sizes (ideal=1).
    #    Primary area is WGS84 ellipsoidal (real Earth); the spherical CV is shown
    #    alongside — they agree except for A5, which is ellipsoid-native.
    datum = "WGS84 ellipsoidal" if _GEOD else "SPHERICAL (pyproj absent)"
    print(f"\nAUTHALICITY  ({datum} cell-area uniformity; CV%→0, ratios→1.0 = ideal)")
    pc = ''.join(f"{f'p{p}/μ':>9}" for p in PCTL)
    print(f"{'system':<10}{'cells':>7}{'mean m²':>12}{'cv%':>8}{pc}"
          f"{'max/min':>9}{'sphCV%':>9}")
    for name in AVAILABLE:
        _, ell, sph, _, _, _, _ = res[name]
        m = np.isfinite(ell) & (ell > 0)
        a, s = ell[m], sph[np.isfinite(sph) & (sph > 0)]
        if len(a) == 0:
            print(f"{SYS_LABEL[name]:<10}  (no measurable cells)")
            continue
        mu = a.mean()
        cv = 100.0 * a.std() / mu
        pcols = ''.join(f"{np.percentile(a, p) / mu:>9.4f}" for p in PCTL)
        sphcv = 100.0 * s.std() / s.mean()
        print(f"{SYS_LABEL[name]:<10}{len(a):>7}{mu:>12.4g}{cv:>8.2f}"
              f"{pcols}{a.max() / a.min():>9.4f}{sphcv:>9.2f}")


if __name__ == '__main__':
    main()
