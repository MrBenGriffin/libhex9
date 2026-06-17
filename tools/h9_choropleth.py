#!/usr/bin/env python3
"""h9_choropleth — point CSV -> adaptive Hex9 choropleth, no PostGIS/QGIS needed.

Reads a CSV of points (lon/lat, optional weight), runs the Hex9 population
digest (h9_adaptive: dense areas resolve to fine cells, sparse areas aggregate
upward, the cell set partitions the SAMPLE and values sum exactly to the input
total), and writes a GeoJSON FeatureCollection of hexagon polygons with
per-cell properties — drop it straight into geojson.io, kepler.gl, a Leaflet
map, or geopandas. Optional --png renders a quick matplotlib choropleth.

Properties per feature: layer, value, npoints, density (persons/km^2), grade
(log9 graduation; +1 == 9x denser). Geometry is the canonical cell polygon
(== h9_grid / h9_cell), so it lines up with anything you later do in PostGIS.

Examples
--------
  h9_choropleth.py pts.csv > cells.geojson
  h9_choropleth.py --weight pop --max-layer 12 --ceiling 100 pts.csv out.geojson
  h9_choropleth.py --png map.png --color density pts.csv out.geojson
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, os.pardir, "build"), os.path.join(_HERE, os.pardir)):
    if os.path.isdir(_cand):
        sys.path.insert(0, os.path.abspath(_cand))
try:
    import numpy as np
    import hex9_ext as h9
except ImportError as e:  # pragma: no cover
    sys.exit(f"h9_choropleth: cannot import hex9_ext / numpy ({e}). "
             f"Build the extension (cmake --build build) or set PYTHONPATH.")

EARTH_KM2 = 510065622.0
LON_NAMES = ("longitude", "lon", "lng", "long", "x")
LAT_NAMES = ("latitude", "lat", "y")


def _find_col(header, wanted, names, role, required=True):
    lower = [c.strip().lower() for c in header]
    if wanted is not None:
        if wanted in header:
            return header.index(wanted)
        if wanted.lower() in lower:
            return lower.index(wanted.lower())
        sys.exit(f"h9_choropleth: {role} column {wanted!r} not in header {header}")
    for n in names:
        if n in lower:
            return lower.index(n)
    if required:
        sys.exit(f"h9_choropleth: could not autodetect the {role} column "
                 f"(looked for {', '.join(names)}); pass --{role}. Header: {header}")
    return None


def _read_points(path, lon_c, lat_c, wt_c, delim):
    """Load lon/lat (+ optional weight) from a CSV; returns float arrays."""
    fin = sys.stdin if path == "-" else open(path, newline="")
    try:
        reader = csv.reader(fin, delimiter=delim)
        header = next(reader, None)
        if header is None:
            sys.exit("h9_choropleth: empty input")
        li = _find_col(header, lon_c, LON_NAMES, "lon")
        ai = _find_col(header, lat_c, LAT_NAMES, "lat")
        wi = _find_col(header, wt_c, (), "weight", required=False) if wt_c else None
        lon, lat, wt = [], [], ([] if wi is not None else None)
        skipped = 0
        for row in reader:
            try:
                lo, la = float(row[li]), float(row[ai])
                w = float(row[wi]) if wi is not None else 1.0
            except (ValueError, IndexError):
                skipped += 1
                continue
            lon.append(lo); lat.append(la)
            if wt is not None:
                wt.append(w)
    finally:
        if fin is not sys.stdin:
            fin.close()
    return (np.asarray(lon), np.asarray(lat),
            (np.asarray(wt) if wt is not None else None), skipped)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="h9_choropleth",
        description="Point CSV -> adaptive Hex9 choropleth (GeoJSON, optional PNG).")
    p.add_argument("input", nargs="?", default="-", help="input CSV (- = stdin)")
    p.add_argument("output", nargs="?", default="-", help="output GeoJSON (- = stdout)")
    p.add_argument("--lon"); p.add_argument("--lat")
    p.add_argument("--weight", help="weight column (default: each point counts 1)")
    p.add_argument("--min-layer", type=int, default=6)
    p.add_argument("--max-layer", type=int, default=12)
    p.add_argument("--ceiling", type=float, default=100.0,
                   help="per-cell value ceiling (digest splits finer above it)")
    p.add_argument("--floor", type=float, default=20.0,
                   help="per-cell value floor (digest aggregates up below it)")
    p.add_argument("--densify", type=int, default=0,
                   help="polygon edge subdivision 0..9 (default 0 = 6 corners)")
    p.add_argument("--delimiter", default=",")
    p.add_argument("--png", metavar="PATH", help="also render a matplotlib PNG")
    p.add_argument("--dpi", type=int, default=150,
                   help="PNG resolution in pixels/inch (default: 150)")
    p.add_argument("--size", type=float, default=8.0,
                   help="map height in inches; width follows the bbox (default: 8)")
    p.add_argument("--color", default="grade",
                   choices=("density", "value", "grade", "layer", "npoints"),
                   help="PNG fill field (default: grade — the log9 graduation "
                        "designed for multi-layer display; density is log-scaled)")
    args = p.parse_args(argv)

    if not (0 <= args.min_layer <= args.max_layer <= 29):
        sys.exit("h9_choropleth: need 0 <= min-layer <= max-layer <= 29")

    h9.warp_init()
    lon, lat, wt, skipped = _read_points(args.input, args.lon, args.lat,
                                         args.weight, args.delimiter)
    if lon.size == 0:
        sys.exit("h9_choropleth: no usable points")

    full = h9.encode(lon, lat)
    uu, layers, values, npoints, _assign = h9.adaptive(
        full, args.min_layer, args.max_layer, args.ceiling, args.floor, wt)

    features = []
    for i in range(uu.shape[0]):
        L = int(layers[i]); v = float(values[i])
        ring = h9.cell(uu[i], L, args.densify)              # (P,2) lon/lat, closed
        density = v * 12.0 * (9.0 ** L) / EARTH_KM2
        grade = (L + math.log(v) / math.log(9.0)) if v > 0 else None
        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon",
                         "coordinates": [[[round(x, 7), round(y, 7)] for x, y in ring]]},
            "properties": {"layer": L, "value": v, "npoints": int(npoints[i]),
                           "density": density, "grade": grade},
        })

    fc = {"type": "FeatureCollection", "features": features}
    fout = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        json.dump(fc, fout)
        fout.write("\n")
    finally:
        if fout is not sys.stdout:
            fout.close()

    if args.png:
        _render_png(features, args.color, args.png, args.size, args.dpi)

    total_v = float(values.sum())
    sys.stderr.write(
        f"h9_choropleth: {lon.size} points -> {len(features)} cells "
        f"(L{int(layers.min())}..L{int(layers.max())}), value sum {total_v:g}"
        + (f", {skipped} skipped" if skipped else "") + "\n")
    return 0


def _render_png(features, field, path, size_in, dpi):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        from matplotlib.colors import Normalize, LogNorm
    except ImportError:
        sys.stderr.write("h9_choropleth: --png needs matplotlib (pip install matplotlib); "
                         "GeoJSON was still written.\n")
        return
    # Cells of different layers OVERLAP (they partition the sample, not the
    # surface), so draw coarse layers first and fine cells last — otherwise a
    # big coarse hex painted on top hides the fine detail underneath.
    rows = []
    for f in features:
        c = f["properties"].get(field)
        if c is None:
            continue
        rows.append((f["properties"]["layer"], f["geometry"]["coordinates"][0], c))
    if not rows:
        sys.stderr.write("h9_choropleth: nothing to plot for --color field.\n")
        return
    rows.sort(key=lambda r: r[0])              # ascending layer => fine cells on top
    polys = [r[1] for r in rows]
    vals = [r[2] for r in rows]
    vals = np.asarray(vals, dtype=float)
    if field == "density":
        # density spans orders of magnitude across layers -> log scale
        pos = vals[vals > 0]
        norm = LogNorm(vmin=(pos.min() if pos.size else 1.0), vmax=vals.max() or 1.0)
    else:
        vmax = np.quantile(vals, 0.98) or vals.max() or 1.0   # clip outliers
        norm = Normalize(vmin=vals.min(), vmax=vmax)
    pc = PolyCollection(polys, array=vals, cmap="viridis", norm=norm,
                        edgecolors="none")

    # Figure sized to the data bbox so there's no wasted whitespace, with a
    # latitude correction (plate-carrée) so the hexes aren't stretched, plus a
    # fixed strip for the colorbar.
    allxy = np.concatenate([np.asarray(pl) for pl in polys])
    minx, miny = allxy.min(axis=0); maxx, maxy = allxy.max(axis=0)
    lat_corr = 1.0 / max(math.cos(math.radians((miny + maxy) / 2.0)), 0.1)
    disp_w = max(maxx - minx, 1e-9)
    disp_h = max((maxy - miny) * lat_corr, 1e-9)
    legend_in = 1.4
    fig_w = size_in * (disp_w / disp_h) + legend_in
    fig, ax = plt.subplots(figsize=(fig_w, size_in), layout="constrained")
    ax.add_collection(pc)
    ax.set_xlim(minx, maxx); ax.set_ylim(miny, maxy)
    ax.set_aspect(lat_corr); ax.set_axis_off()
    fig.colorbar(pc, ax=ax, label=field, shrink=0.85)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    sys.stderr.write(f"h9_choropleth: wrote {path} "
                     f"({fig_w * dpi:.0f}x{size_in * dpi:.0f}px @ {dpi}dpi)\n")


if __name__ == "__main__":
    raise SystemExit(main())
