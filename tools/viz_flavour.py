# viz_flavour.py — plot the dump from viz_flavour.c.
#
# Part of the Hex9 (H9) Project
# Copyright ©2026, Ben Griffin
# Licensed under the Apache License, Version 2.0
#
# Two panels over the same grid of encoded sample points:
#   left  — coloured by WALK BIN (h9_bin; byte-identical to Python = truth)
#   right — coloured by IDENTITY resolution (ring rendered from the full uuid)
# Black outlines: the disk(2) cell polygons (canonical kring emission).
# Where the right panel's colour regions cut THROUGH outlines, identity is
# resolving those points to the wrong cell — the mode-1-half defect.
#
# Usage: viz_flavour <dump.csv> <out.png> [title]

import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dump, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else ""

pts = []            # (lon, lat, walkbin, ringhash)
rings = {}          # cellidx -> [(lon, lat), ...]
with open(dump) as f:
    for row in csv.reader(f):
        if row[0] == "P":
            pts.append((float(row[1]), float(row[2]), row[3], row[4]))
        elif row[0] == "R":
            rings.setdefault(int(row[1]), []).append((float(row[3]), float(row[4])))

def colour_map(keys):
    # one visually-distinct colour per key — collisions in a small cycled
    # palette previously made adjacent same-coloured cells read as split
    # hexagons (the phantom "striping")
    uniq = sorted(set(keys))
    cmap = plt.get_cmap("gist_ncar")
    n = max(len(uniq), 2)
    return {k: cmap(0.05 + 0.9 * i / (n - 1)) for i, k in enumerate(uniq)}

lons = [p[0] for p in pts]
lats = [p[1] for p in pts]
wcol = colour_map([p[2] for p in pts])
icol = colour_map([p[3] for p in pts])

fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharex=True, sharey=True)
for ax, key, name in ((axes[0], 2, "walk bins (h9_bin == Python: truth)"),
                      (axes[1], 3, "identity resolution (ring of full uuid)")):
    ax.scatter(lons, lats, c=[(wcol if key == 2 else icol)[p[key]] for p in pts],
               s=4, marker="s", linewidths=0)
    for poly in rings.values():
        xs = [v[0] for v in poly]
        ys = [v[1] for v in poly]
        ax.plot(xs, ys, "k-", linewidth=1.2)
    # annotate each colour region with the trailing digits of its key at the
    # region's mean position (walk panel: last 3 body digits of the bin label)
    regions = {}
    for p in pts:
        regions.setdefault(p[key], []).append((p[0], p[1]))
    for k, ps in regions.items():
        if len(ps) < 12:
            continue
        mx = sum(v[0] for v in ps) / len(ps)
        my = sum(v[1] for v in ps) / len(ps)
        if key == 2:
            # bin hex: body nibbles up to sentinel + tail nibble
            body = k.split("f" * 6)[0]
            tag = body[-3:] + "." + k[-1]
        else:
            tag = k[-3:]
        ax.text(mx, my, tag, fontsize=6, ha="center", va="center")
    ax.set_title(name)
    ax.set_aspect("equal", adjustable="box")

fig.suptitle(title)
fig.tight_layout()
fig.savefig(out, dpi=140)
print(f"wrote {out}  ({len(pts)} samples, {len(rings)} outlines)")
