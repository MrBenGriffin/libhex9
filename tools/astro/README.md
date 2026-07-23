# tools/astro — Hex9 on the celestial sphere

The sphere-datum pipeline (2.1.0) applied to real sky: adaptive Hex9 digests
of the [HYG star catalogue](https://github.com/astronexus/HYG-Database)
(~120k stars, CC BY-SA 4.0).

Star positions are RA/dec on the celestial sphere — already-spherical
coordinates with no ellipsoid anywhere — so everything runs through the
`sphere=True` entry points untouched, and density comes out per **steradian**
(the unit sphere's cell area is intrinsic; the sphere datum carries no
radius). Doctrine: `docs/warp-regimes.md`, "Two datums, one regime".

## Run

```sh
python3 h9_sky.py --png            # offers to download the catalogue (~35 MB)
python3 h9_sky.py --download       # non-interactive fetch
```

Outputs in `out/` (gitignored): `sky_counts.csv/png` and `sky_flux.csv/png` —
the `h9_adaptive` row shape (`h9_bin, layer, value, npoints, density, grade`).

The catalogue lands in `data/hyg.csv` (gitignored — not redistributed here).
To fetch it yourself:

```sh
mkdir -p data
curl -o data/hyg.csv https://raw.githubusercontent.com/astronexus/HYG-Database/main/hyg/CURRENT/hygdata_v41.csv
```

## The two digests

- **counts** — every star weighs 1; the density map is star counts per
  steradian. The galactic plane draws itself as the classic sinusoidal band.
- **flux** — weight = 10^(−0.4·mag) (apparent brightness, Vega units). The
  digest isolates dominant point sources into deep cells; verified 2026-07-22:
  the top eight cells by value are the eight brightest stars of the night
  sky, **in order** — Sirius (3.767), Canopus (1.770), Arcturus (1.047),
  Rigil Kentaurus (1.009, at L7 because α Cen A+B share a cell), Vega
  (0.973), Capella, Rigel, Procyon — with decoded positions matching
  published RA/dec to ~3 decimals.

## Conventions (important)

- HYG's `ra` is in **hours**; the prep converts to degrees and normalises to
  [−180, 180) as spherical lon. `dec` is spherical lat. The Sun is excluded.
- Addresses minted here are **sphere-datum**. Decode them only with
  `decode(..., sphere=True)` / `h9_decode_sphere`; the 16 bytes do not record
  the datum, so keeping that straight is this dataset's metadata contract.
- Per-km² on an actual body would be `density × 4π / body_area` — not
  meaningful for the sky; per-steradian is the natural unit here (sources/sr
  and flux/sr, i.e. surface brightness).
