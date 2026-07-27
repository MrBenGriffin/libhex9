#!/usr/bin/env python3
"""h9_net.py — render Hex9 digests on octahedral-net layouts, ownership-direct.

Ready-rolled hexagons on a flat net (rhombus / windmill / ...): every cell is
a per-layer template hexagon placed by ADDRESS, not by geometry —

    cell → L0 owner (h9 cell_ancestor)
         → which of the owner's two d-cells (HX_OC2, by the cell's octant)
         → that (octant, c2) piece's static affine from the layout table.

No face classification, no projection engine, and — since cell_uv (2.1.x) —
no geometry at all until the final affine: each cell's polygon is its EXACT
integer lattice vertices (canonical shared keys, hex9_cell_uv), each vertex
mapped by its own side's piece affine, so polygons bend correctly at octant
folds and adjacent cells share vertices to the last bit. No chevrons, no
templates, no lon/lat anywhere: the whole pipeline is address-side, which
also makes it DATUM-FREE — a WGS84-minted and a sphere-minted digest render
identically, so there is no datum flag.

Layout tables and the d-cell incidence are vendored from the Hex9 project's
hhg9 (Ben Griffin, Apache-2.0: hhg9/domains/nets.py `net_layouts`, grid units
U=W/6 · V=H/9 with placements in 3-unit steps; H9O.l0hex_by_id). hhg9 remains
authoritative; re-derive with the session's derive_hex_table.py if they move.

Usage:
  h9_net.py digest.csv [--layout rhombus] [--color grade] [--clip 2,98]
            [--output out.png] [--title "..."]
  h9_net.py --placeholder 2 --output l2.png   # 12 L0 owners, colour-coded

The astronomy recipe (flux digests are bimodal — a few star cells sit many
grades above the field): --clip 2,90 with the default --color grade.
"""
from __future__ import annotations

import argparse
import csv
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
    sys.exit(f"h9_net: cannot import hex9_ext / numpy ({e}). "
             f"Build the extension (cmake --build build) or set PYTHONPATH.")

GW = math.sqrt(2.0) / 2.0     # net grid x-unit (W/2 = 3·U)
GH = math.sqrt(6.0) / 6.0     # net grid y-unit (H/3 = 3·V)

# face sign triple -> (gx, gy, rot60); split faces: [primary, c2_0, c2_1, c2_2]
NET_LAYOUTS = {
    'rhombus': {
        (+1, +1, +1): (3., 2., 3),
        (-1, +1, +1): (4., 3., 5),
        (+1, -1, +1): (2., 3., 1),
        (-1, -1, +1): [(2., 5., 5), (4., 0., 4), (2., 0, -4), (0, 0, 0)],
        (+1, +1, -1): [(3., 0., 3), (0., 0., 0), (2., 0, -2), (-2., 0, -4)],
        (-1, +1, -1): (5., 2., 5),
        (+1, -1, -1): (1., 2., 1),
        (-1, -1, -1): [(6., 3., 3), (0., 0., 0), (-6., 0, 0), (0, 0, 0)],
    },
    'windmill': {
        (-1, -1, +1): (3., 7., 5),
        (-1, -1, -1): (1., 5., 3),
        (+1, -1, +1): (3., 5., 1),
        (+1, -1, -1): (2., 4., 1),
        (+1, +1, +1): (4., 4., 3),
        (+1, +1, -1): (4., 2., 3),
        (-1, +1, +1): (5., 5., 5),
        (-1, +1, -1): (6., 4., 5),
    },
    'windmill_pacific': {
        (-1, -1, +1): (4., 5., 2),
        (-1, -1, -1): (5., 4., 2),
        (+1, -1, +1): (4., 7., 4),
        (+1, -1, -1): (6., 5., 0),
        (+1, +1, +1): (2., 5., 4),
        (+1, +1, -1): (1., 4., 4),
        (-1, +1, +1): (3., 4., 0),
        (-1, +1, -1): (3., 2., 0),
    },
}

# L0 hex -> ((oid, c2) of the mode-0 d-cell, (oid, c2) of the mode-1 d-cell).
# Derived from hhg9 H9O.l0hex_by_id / H9O.oid_mo.
HX_OC2 = {
    0x0: ((0, 0), (4, 0)), 0x1: ((5, 0), (1, 0)), 0x2: ((6, 0), (2, 0)),
    0x3: ((3, 0), (7, 0)), 0x4: ((0, 1), (2, 2)), 0x5: ((0, 2), (1, 1)),
    0x6: ((3, 2), (2, 1)), 0x7: ((3, 1), (1, 2)), 0x8: ((5, 2), (4, 1)),
    0x9: ((5, 1), (7, 2)), 0xa: ((6, 1), (4, 2)), 0xb: ((6, 2), (7, 1)),
}


def rot60(th):
    t = (th % 6) * math.pi / 3.0
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s], [s, c]])


def piece_affines(layout_name):
    """(oid, c2) -> (M, t): the 24 rigid placements of one layout."""
    layout = NET_LAYOUTS[layout_name]
    out = {}
    for sign, val in layout.items():
        sx, sy, sz = sign
        oid = ((sz < 0) << 2) | ((sy < 0) << 1) | (sx < 0)
        if isinstance(val, list):
            gx, gy, th = val[0]
            Mf, tf = rot60(th), np.array([gx * GW, gy * GH])
            for c2, (x, y, t) in enumerate(val[1:]):
                out[(oid, c2)] = (Mf @ rot60(t),
                                  np.array([x * GW, y * GH]) + tf)
        else:
            gx, gy, th = val
            M, t = rot60(th), np.array([gx * GW, gy * GH])
            for c2 in (0, 1, 2):
                out[(oid, c2)] = (M, t)
    return out


# The universal integer hexagon: every H9 cell's ring is its centre key plus
# a rotation of these offsets (verified: one distinct ring across all
# interior cells at any layer). Used to reconstruct, in the cell's governing
# frame, the rare vertex whose canonical key lives in a third octant (cone
# corners) — everything else renders from its canonical key directly.
OFF = np.array([(-2, 0), (-1, 1), (1, 1), (2, 0), (1, -1), (-1, -1)],
               dtype=np.int64)


def render(bins, layers, cvals, layout_name, out_png, title, by_owner=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPoly
    from matplotlib.collections import PatchCollection

    aff = piece_affines(layout_name)
    u1, v3 = h9.uv_units()
    layers = np.asarray(layers, dtype=np.int32)
    # Adaptive digests report cells max_layer first (digestion order), and
    # cells of different layers OVERLAP (a parent holds what its descendants
    # did not digest) — paint coarse first so fine cells stay on top.
    order = np.argsort(layers, kind="stable")
    bins = bins[order]
    layers = layers[order]
    cvals = np.asarray(cvals)[order]
    cia, cib, coid, via, vib, void_, _ = h9.cell_uv(bins, layers)
    owners = h9.cell_ancestor(bins, 0)
    hexid = np.array([int(h9.label(u, 0, False), 16) for u in owners])

    cmap = plt.get_cmap('tab20' if by_owner else 'viridis')
    patches, cols = [], []
    for i in range(len(bins)):
        s0, s1 = HX_OC2[int(hexid[i])]
        div = 3.0 ** int(layers[i])
        # phase-lock the offset ring from any governing-frame vertex, so
        # third-octant vertices can be reconstructed in the governing frame
        k = None
        for v in range(6):
            if int(void_[i, v]) == int(coid[i]):
                du, dv = int(via[i, v] - cia[i]), int(vib[i, v] - cib[i])
                hits = np.where((OFF[:, 0] == du) & (OFF[:, 1] == dv))[0]
                if len(hits) == 1:
                    k = (int(hits[0]) - v) % 6
                    break
        poly = np.empty((6, 2))
        for v in range(6):
            o = int(void_[i, v])
            if o == s0[0]:
                piece, u, w = s0, via[i, v], vib[i, v]
            elif o == s1[0]:
                piece, u, w = s1, via[i, v], vib[i, v]
            else:
                # cone corner: reconstruct in the governing frame
                piece = s0 if int(coid[i]) == s0[0] else s1
                u = cia[i] + OFF[(v + (k or 0)) % 6, 0]
                w = cib[i] + OFF[(v + (k or 0)) % 6, 1]
            M, t = aff[piece]
            poly[v] = np.array([u * u1 / div, w * v3 / div]) @ M + t
        patches.append(MplPoly(poly, closed=True))
        cols.append(cmap(int(hexid[i]) % 20 if by_owner else float(cvals[i])))

    allv = np.vstack([p.get_xy() for p in patches])
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.add_collection(PatchCollection(patches, facecolors=cols,
                                      edgecolors='none'))
    ax.set_xlim(allv[:, 0].min() - .1, allv[:, 0].max() + .1)
    ax.set_ylim(allv[:, 1].min() - .1, allv[:, 1].max() + .1)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, color='0.25')
    fig.savefig(out_png, dpi=400, bbox_inches='tight', facecolor='white')
    print(f"h9_net: wrote {out_png} ({len(patches)} cells, "
          f"layout {layout_name})")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Render an h9_choropleth --csv digest on an octahedral "
                    "net, placed by address ownership.")
    p.add_argument("input", nargs="?",
                   help="digest CSV (h9_bin, layer, ..., density)",
                   default=f"{_HERE}/astro/out/sky_flux.csv")
    p.add_argument("--output", help="output PNG", default="h9_net.png")
    p.add_argument("--title", nargs="?", default="")
    p.add_argument("--layout", default="rhombus",
                   choices=sorted(NET_LAYOUTS))
    p.add_argument("--color", default="grade",
                   choices=("grade", "density", "value", "npoints", "layer"),
                   help="digest column to colour by (default grade — the "
                        "digest's log9 graduation, built for multi-layer "
                        "display; density/value are log-scaled)")
    p.add_argument("--clip", default="2,98", metavar="LO,HI",
                   help="percentile range mapped to the colour ramp (default "
                        "2,98). Cells outside saturate — with a few extreme "
                        "cells (bright stars), tighten HI (e.g. 2,90) to give "
                        "the ramp back to the background field")
    p.add_argument("--placeholder", metavar="L",
                   help="ignore input; render all layer-L cells coloured by "
                        "their L0 owner (L = 0..3)")
    args = p.parse_args(argv)
    # No datum flag: rendering is address-side (cell_uv), identical for
    # WGS84- and sphere-minted digests.

    if args.placeholder is not None:
        L = int(args.placeholder)
        uu, _, _ = h9.grid(-180.0, -90.0, 180.0, 90.0, max(L, 1),
                           max_cells=200000)
        if L == 0:
            uu = np.unique(h9.bin(uu, 0), axis=0)
        render(uu, np.full(len(uu), L), np.zeros(len(uu)), args.layout,
               args.output,
               args.title or f"L{L} cells on '{args.layout}' by L0 owner",
               by_owner=True)
        return 0

    if not args.input:
        p.error("digest CSV required (or use --placeholder)")
    rows = list(csv.DictReader(open(args.input)))
    bins = np.array([np.frombuffer(
        bytes.fromhex(r["h9_bin"].replace("-", "")), dtype=np.uint8)
        for r in rows])
    layers = np.array([int(r["layer"]) for r in rows])
    # grade is already the digest's log9 graduation (linear scale); density
    # and value span orders of magnitude (log scale). Empty grade (value<=0)
    # sinks to the bottom of the ramp.
    vals = np.array([float(r[args.color]) if r[args.color] != "" else np.nan
                     for r in rows])
    if args.color in ("density", "value"):
        vals = np.log10(np.maximum(vals, 1e-300))
    try:
        lo_p, hi_p = (float(x) for x in args.clip.split(","))
    except ValueError:
        p.error("--clip wants LO,HI percentiles, e.g. 2,90")
    if not (0 <= lo_p < hi_p <= 100):
        p.error("--clip wants 0 <= LO < HI <= 100")
    lo, hi = np.nanpercentile(vals, [lo_p, hi_p])
    if hi <= lo:
        hi = lo + 1e-12
    cvals = np.clip((np.nan_to_num(vals, nan=lo) - lo) / (hi - lo), 0, 1)
    render(bins, layers, cvals, args.layout, args.output, args.title)
    return 0


if __name__ == "__main__":
    sys.exit(main())
