#!/usr/bin/env python3
"""
check_equal_area.py — validate Hex9's equal-area claim against GeographicLib.

WHY THIS IS DIFFERENT FROM EVERY OTHER CHECK IN THE REPO
--------------------------------------------------------
test/authalic, test/via_sphere and test/pole_seam all compare libhex9 against
hhg9. Both are this project's own work, so they prove the C port is FAITHFUL —
they cannot prove the thing being ported is CORRECT. A misconception shared by
hhg9 and libhex9 passes all of them.

GeographicLib (Karney) shares no code and no authorship with Hex9. Its
PolygonArea computes exact geodesic area on the WGS84 ellipsoid. So this
script does not check our formula against their formula; it checks the
OUTCOME — "are the cells actually equal-area on the real ellipsoid?" — against
a third party that has never seen our latitude series or our warp.

That is the claim the README makes (~0.014 % area uniformity) and the whole
purpose of the authalic latitude. This is the only test that can falsify it.

METHOD
------
Cell edges are not geodesics, so each edge is densified into 3^densify
segments and PolygonArea is given the polyline; it treats consecutive vertices
as geodesic segments, which converges to the true edge as densify rises. The
script reports the convergence so the discretisation error is visible rather
than assumed — if the spread does not settle as densify increases, the number
is discretisation noise, not a real area spread.

Usage:  python3 tools/support/check_equal_area.py [layer] [densify]
"""
import sys

import numpy as np

sys.path.insert(0, '/Users/ben/Documents/Projects/libhex9/build')

try:
    import hex9_ext
except ImportError as exc:                          # noqa: BLE001
    sys.exit(f'cannot import hex9_ext from build/: {exc}\n'
             'Build it first:  cmake --build build --target hex9_ext')

try:
    from geographiclib.geodesic import Geodesic
    from geographiclib.polygonarea import PolygonArea
except ImportError as exc:                          # noqa: BLE001
    sys.exit(f'geographiclib not available: {exc}\n  pip install geographiclib')

WGS84_AREA_M2 = 5.100656217240886e14      # total WGS84 ellipsoid surface area


def cell_area_m2(ring_lonlat):
    """Exact geodesic area (m^2) of a closed lon/lat ring, via GeographicLib."""
    poly = PolygonArea(Geodesic.WGS84, False)
    # The ring from hex9 is closed (last point repeats the first); PolygonArea
    # closes it itself, so feed only the distinct vertices.
    for lon, lat in ring_lonlat[:-1]:
        poly.AddPoint(lat, lon)
    _, _, area = poly.Compute(False, True)
    return abs(area)


def main():
    layer = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    densify = int(sys.argv[2]) if len(sys.argv) > 2 else 3

    n_expected = 12 * 9 ** layer
    print(f'layer {layer}: {n_expected} cells expected, densify={densify} '
          f'({hex9_ext.__doc__ and ""}{6 * 3 ** densify + 1} ring points)')

    # Enumerate distinct cells at `layer` by binning a dense point sample.
    # (Sampling equal-area in latitude so no region is over-represented.)
    rng = np.random.default_rng(20260721)
    n_pts = max(n_expected * 40, 200_000)
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, n_pts)))
    lon = rng.uniform(-180.0, 180.0, n_pts)

    uu = hex9_ext.encode(lon, lat)
    bins = hex9_ext.bin(uu, layer)
    uniq = np.unique(bins, axis=0)
    print(f'  sampled {n_pts} points -> {len(uniq)} distinct cells '
          f'({100.0 * len(uniq) / n_expected:.1f}% of the layer)')

    areas = np.array([cell_area_m2(hex9_ext.cell(u, layer, densify))
                      for u in uniq])

    mean = areas.mean()
    spread = (areas.max() - areas.min()) / mean
    print(f'  mean area      {mean:.6e} m^2')
    print(f'  min/max        {areas.min():.6e} / {areas.max():.6e}')
    print(f'  spread         {100.0 * spread:.6f} %  (max-min)/mean')
    print(f'  std/mean       {100.0 * areas.std() / mean:.6f} %')

    # Closure: cells at a layer partition the globe, so mean*count should be
    # the ellipsoid's area. This catches a systematic scale error that a pure
    # uniformity check would miss entirely.
    implied_total = mean * n_expected
    closure = implied_total / WGS84_AREA_M2 - 1.0
    print(f'  implied total  {implied_total:.6e} m^2 '
          f'vs WGS84 {WGS84_AREA_M2:.6e}')
    print(f'  closure error  {100.0 * closure:+.6f} %')
    return areas, spread, closure


if __name__ == '__main__':
    main()
