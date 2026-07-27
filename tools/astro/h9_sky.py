#!/usr/bin/env python3
"""h9_sky.py — adaptive Hex9 digests of the celestial sphere (HYG catalogue).

The sphere-datum demonstration: star positions are RA/dec on the celestial
sphere — already-spherical coordinates, no ellipsoid anywhere — so the whole
pipeline runs through the `sphere=True` entry points with zero special-casing.
Density is per STERADIAN (the unit sphere's cell area is intrinsic; the
sphere datum carries no radius — see docs/warp-regimes.md, "Two datums, one
regime").

Two digests over the ~120k-star HYG catalogue:

  counts — every star weighs 1. The density map is star counts / steradian:
           the Milky Way's galactic plane draws itself.
  flux   — weight = 10^(-0.4·mag), apparent brightness in Vega units. The
           digest isolates dominant point sources into deep cells: the top
           cells BY VALUE are the brightest stars of the night sky, in order
           (Sirius 3.77, Canopus 1.77, Arcturus 1.05, ...), while the diffuse
           background aggregates at the floor layer.

Conventions: HYG's `ra` is in HOURS — converted to degrees and normalised to
[-180, 180) as the spherical lon; `dec` is the spherical lat. The Sun is
excluded. Addresses minted here are SPHERE-datum: decode them with
`decode(..., sphere=True)` / `h9_decode_sphere` only.

The catalogue (~35 MB) is not in the repo: this script offers to download it
into data/ (gitignored) on first run — or fetch it yourself, see README.md.

Everything downstream of the prep is ../h9_choropleth.py — this script just
prepares (lon, lat, flux) and drives it.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, os.pardir))   # for h9_choropleth
import h9_choropleth  # noqa: E402

HYG_URL = ("https://raw.githubusercontent.com/astronexus/HYG-Database/"
           "main/hyg/CURRENT/hygdata_v41.csv")
DATA_DIR = os.path.join(_HERE, "data")
OUT_DIR = os.path.join(_HERE, "out")


def fetch_catalogue(path: str, assume_yes: bool) -> None:
    print(f"h9_sky: catalogue not found at {path}", file=sys.stderr)
    print(f"h9_sky: source: {HYG_URL} (~35 MB, CC BY-SA 4.0)", file=sys.stderr)
    if not assume_yes:
        if not sys.stdin.isatty():
            sys.exit("h9_sky: refusing to download without a terminal — "
                     "re-run with --download, or fetch it yourself (README.md)")
        if input("h9_sky: download now? [y/N] ").strip().lower() != "y":
            sys.exit("h9_sky: aborted")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    urllib.request.urlretrieve(HYG_URL, tmp)
    os.replace(tmp, path)
    print(f"h9_sky: downloaded {path}", file=sys.stderr)


def prepare(cat_path: str, stars_path: str) -> tuple[int, float]:
    """HYG → (lon, lat, flux) CSV. Returns (n_stars, total_flux)."""
    n = 0
    tot = 0.0
    with open(cat_path, newline="") as f, \
         open(stars_path, "w", newline="") as g:
        r = csv.DictReader(f)
        w = csv.writer(g)
        w.writerow(["lon", "lat", "flux"])
        for row in r:
            try:
                ra_h = float(row["ra"])
                dec = float(row["dec"])
                mag = float(row["mag"])
            except (ValueError, KeyError):
                continue
            if row.get("proper") == "Sol" or mag < -20:   # no Sun
                continue
            lon = (ra_h * 15.0 + 180.0) % 360.0 - 180.0   # hours → [-180, 180)°
            flux = 10.0 ** (-0.4 * mag)                   # Vega units
            tot += flux
            w.writerow([f"{lon:.6f}", f"{dec:.6f}", f"{flux:.8g}"])
            n += 1
    return n, tot


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Adaptive Hex9 digests of the HYG star catalogue "
                    "(sphere datum; density per steradian).")
    p.add_argument("--data", default=os.path.join(DATA_DIR, "hyg.csv"),
                   help="HYG catalogue CSV (default: data/hyg.csv; offered "
                        "for download when missing)")
    p.add_argument("--download", action="store_true",
                   help="fetch the catalogue without asking")
    p.add_argument("--out", default=OUT_DIR,
                   help="output directory (default: out/)")
    p.add_argument("--which", choices=("counts", "flux", "both"),
                   default="both")
    p.add_argument("--png", action="store_true",
                   help="also render density PNGs (matplotlib)")
    p.add_argument("--counts-layers", default="1:9", metavar="MIN:MAX")
    p.add_argument("--counts-ceiling", type=float, default=60.0,
                   help="stars per cell before refining (default 60)")
    p.add_argument("--flux-layers", default="1:9", metavar="MIN:MAX")
    p.add_argument("--flux-ceiling", type=float, default=None,
                   help="Vega-flux per cell before refining "
                        "(default: total/200, printed)")
    args = p.parse_args(argv)

    if not os.path.exists(args.data):
        fetch_catalogue(args.data, args.download)

    os.makedirs(args.out, exist_ok=True)
    stars = os.path.join(args.out, "stars.csv")
    n, tot = prepare(args.data, stars)
    print(f"h9_sky: {n} stars, total flux {tot:.2f} Vega units",
          file=sys.stderr)

    def run(tag, layers, ceiling, floor, weight):
        lo, hi = (int(x) for x in layers.split(":"))
        argv = ["--sphere", "--csv", "--color", "density",
                "--min-layer", str(lo), "--max-layer", str(hi),
                "--ceiling", str(ceiling), "--floor", str(floor)]
        if weight:
            argv += ["--weight", weight]
        if args.png:
            argv += ["--png", os.path.join(args.out, f"sky_{tag}.png")]
        argv += [stars, os.path.join(args.out, f"sky_{tag}.csv")]
        print(f"h9_sky: {tag}: ceiling {ceiling:g}, floor {floor:g}, "
              f"layers {lo}..{hi}", file=sys.stderr)
        rc = h9_choropleth.main(argv)
        if rc:
            sys.exit(f"h9_sky: {tag} digest failed ({rc})")

    if args.which in ("counts", "both"):
        run("counts", args.counts_layers,
            args.counts_ceiling, args.counts_ceiling / 4.0, None)
    if args.which in ("flux", "both"):
        ceiling = args.flux_ceiling if args.flux_ceiling is not None \
            else tot / 200.0
        run("flux", args.flux_layers, ceiling, ceiling / 5.0, "flux")
    return 0


if __name__ == "__main__":
    sys.exit(main())
