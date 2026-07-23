# Warp regimes: why 2.0.0 changed every address, and what to do about it

## The one-sentence version

A Hex9 address is only meaningful together with the projection that produced
it, nothing in the 16 bytes records which projection that was, and 2.0.0
changed the projection — so addresses written by 1.x and addresses written by
2.x are different things that look identical.

## What changed

Up to 1.5.0 the core ran on the WGS84 ellipsoid against a WGS84-trained warp
field. Since 2.0.0 it reduces geodetic latitude to the **authalic** latitude
(Karney series), runs on the **unit sphere**, and applies the **Sphere-L6
wedge-fold** field. One trained field now serves every ellipsoid; the
ellipsoid-specific part of the chain is the latitude series and nothing else.

The old regime is **removed**, not deprecated — no GUC, no build flag, no
function argument reaches it. That is deliberate. While both existed, one point
had two valid addresses with no way to tell them apart, which is a trap that
only ever springs long after the mistake.

## How much addresses move

Point-dependent, and shallower than is comfortable. Measured at the commit
where both regimes still coexisted:

| Point | First differing nibble |
|---|---|
| Westminster, Tokyo, New York, Cape Town, octant corner | 7 |
| Sydney | 8 |
| Edinburgh | 9 |
| Near the poles | 17 |
| Null island | identical |

So coarse bins are frequently unchanged and **that is exactly what makes this
dangerous**: a spot-check on a handful of cells can easily show no difference
at all. The project's own two anchors demonstrate the trap — Edinburgh agrees
to nibble 9, so the entire regression suite passed through the regime change
untouched, while Westminster's layer-8 bin moved:

```
westminster  435878503  →  435878530
edinburgh    432177478  →  432177468
```

## The doctrine: re-derive, never remap

**Re-derive addresses from the source geometry you encoded from.**

```sql
UPDATE my_table SET h9 = h9_encode(geom);
```

**Do not decode-and-re-encode existing addresses.** It is technically easy and
it is wrong, for two independent reasons.

### It cannot detect itself

Remapping is not idempotent, and there is no way to tell a converted address
from an unconverted one — every 16-byte pattern is a valid address, so a second
pass "succeeds" and silently displaces every cell again. An interrupted or
retried migration leaves data in a state nothing can diagnose afterwards.

### It invalidates everything downstream

This is the deeper reason, and it holds even if the migration runs exactly
once. A population density, an aggregate, a published figure, a rendered map —
anything computed against 1.x addresses — is no longer reproducible from
remapped data. The addresses change underneath the derived work, and nothing
announces it. The cells that moved are precisely the ones whose contents
changed, so the derived numbers are not merely stale, they are unverifiable.

Treat 1.x addresses as historical values belonging to a prior regime — the way
you would treat measurements from a since-recalibrated instrument. Keep them
and know what they are, or re-derive from source. Do not rewrite them in place.

## If the source geometry is gone

Then the addresses are the only record of position, and **only the old library
can interpret them**. Recover coordinates *before* upgrading:

```sql
-- while still on postgis_hex9 1.5.0 / libhex9 1.x
ALTER TABLE my_table ADD COLUMN geom_recovered geometry(Point, 4326);
UPDATE my_table SET geom_recovered = h9_centroid(h9);
```

Then upgrade and re-encode from `geom_recovered`.

Mind the precision ceiling: a centroid is the centre of its cell, not the
original point. The round trip is lossless *at that cell's layer* and discards
any finer precision the original carried. Prefer real source geometry wherever
it still exists.

If you have already upgraded, `git checkout v1.5.0` still builds the old
library — that is part of why the tag exists.

## What guards this now

The 2.0.0 change slipped through the entire test suite without a single
failure. Three guards exist so that cannot recur:

1. **`test/regime_pin`** — 230 frozen `lon/lat → uuid → curve` goldens plus the
   layer-0..29 ownership ladder, in `test_data/regime_pin*.tsv`. Any projection
   change fails here loudly. Regenerating is a deliberate act
   (`tools/gen_regime_pin.c`), and the goldens say so in their header.
2. **`test/via_sphere` and `test/pole_seam`** — the chain against the hhg9
   Python reference, the second sampling only the worst-conditioned points
   (polar approach, straddled seams into the polar band, octant triple points).
3. **`_PG_init` version check** — `postgis_hex9` refuses to load against a
   libhex9 whose version differs from the header it was built against. Because
   the shared library carries an SOVERSION, a major upgrade installs under a
   *new filename* and leaves the old one in place, so an un-rebuilt module will
   happily keep loading 1.x and emitting the previous regime's addresses while
   every version string in sight says 2.0.0.

For independent confirmation that the equal-area property actually holds —
rather than merely agreeing with our own Python — see
`tools/support/check_equal_area.py`, which measures cell areas with
GeographicLib's exact geodesic `PolygonArea`.

## Two datums, one regime (2.1.0)

2.1.0 added `*_sphere` twins of the lon/lat entry points (`hex9_encode_sphere`,
`hex9_decode_sphere`, `hex9_project_sphere`, `hex9_unproject_sphere`, their
`_many` forms, `hex9_grid_create_sphere`, `hex9_cell_ring_sphere`; in Python a
`sphere=True` keyword). **This is not a second regime**, and it must never
become one. Since 2.0.0 every address is minted on the unit sphere; WGS84
enters the chain at exactly one place — the authalic latitude reduction at the
lon/lat boundary. The sphere twins run the *identical* chain minus that
reduction: their lon/lat are already-spherical degrees. They exist for callers
that own their datum — another body's authalic frame, celestial RA/dec — and
they make libhex9 usable as a pure spherical addressing engine.

Everything the rest of this document says about regimes applies to datums,
scaled down one level:

- The 16 bytes still record nothing. A WGS84-minted and a sphere-minted
  address for the same numeric lon/lat differ below roughly layer 5 (the
  authalic shift peaks at ~21 km near 45°; coarse bins agree, deep ones do
  not). `test/sphere_mode.c` pins both facts — the parity (sphere twin ==
  WGS84 chain given pre-reduced input) and the divergence.
- Therefore the datum is **dataset metadata**, owned by the caller, exactly
  like the body the data belongs to. Never mix datums within one dataset.
- The choice is carried by *which function you call* — part of the function's
  identity, immutable, index-safe — never by ambient state. No GUC, no
  setting, no push/pop. That is why they are twins and not a mode: the
  pre-2.0.0 `hex9.use_warp` GUC is the cautionary tale (see
  `extension/postgis_hex9/lwgeom_hex9.cpp`).
- The grid handle records its datum at create and emits all subsequent
  lon/lat in it; it is fixed for the handle's lifetime.
- One `hex9_init()` serves both datums (the sphere path needs no series, both
  need the warp field).

The deliberately small twin surface: anything composable from symbolic ops
plus one datum-crossing primitive got **no** twin (w_oct is a pure rotation
from b_oct; label centroids = `hex9_parse_label` + `hex9_decode_sphere`; grid
centroids = `hex9_grid_cell_id` + `hex9_decode_sphere`). The adaptive digest
takes and returns addresses and never crosses the boundary at all — its SQL
twin `h9_adaptive_sphere` differs only in rendering (geom SRID 0) and in the
density unit: value per **steradian**, because a layer-L cell's area on the
unit sphere (4π/(12·9^L) sr) is intrinsic, and the sphere datum deliberately
carries no radius. Per-km² on a specific body is density × 4π / body_area,
caller-side — the body, like the datum, is dataset metadata.

## If the projection changes again

It might. Treat it as a major version, always, and:

- regenerate `test_data/regime_pin*.tsv` **as a separate, reviewable commit**,
  so the address movement is visible in the diff rather than buried;
- bump `HEX9_VERSION` in `hex9_c.h` (CMake fails the build if it disagrees with
  `project()`), so the `_PG_init` guard trips for stale consumers;
- say so in `CHANGELOG.md` in the terms above — re-derive, never remap.

The warp *mechanism* was kept general when the second regime was removed:
`WarpState`, the blob loader and `h9_warp_init_from_path` are all still there.
Swapping the field is a blob change, not an architecture change.
