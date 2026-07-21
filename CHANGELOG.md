# Changelog

All notable changes to libhex9, the `postgis_hex9` extension, the PDAL filters
and the `hex9-sys` crate. These share one version number.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning is [semantic](https://semver.org/), with the addition that **any
change to the projection is a major version**, because it changes every stored
address.

---

## [2.0.0] — 2026-07-21

### ⚠️ Addresses change. Read this first.

The addressing chain was replaced. **The same point generally addresses to a
different cell from about layer 7 downward**, and nothing in a 16-byte address
records which regime produced it — 1.x and 2.x addresses are different things
that look identical.

Re-derive addresses from source geometry (`UPDATE t SET h9 = h9_encode(geom)`).
**Do not decode-and-re-encode stored addresses**: it is not idempotent, cannot
detect a second run, and silently invalidates anything already derived from
those addresses. Full reasoning, and what to do when the source geometry is
gone, in [`docs/warp-regimes.md`](docs/warp-regimes.md).

The movement is point-dependent and often invisible in a spot check — coarse
bins frequently do not move at all. Both project anchors illustrate it:

```
westminster  435878503  →  435878530
edinburgh    432177478  →  432177468     (unchanged through layer 8's prefix)
```

### Added

- **Regime pin tests** — `test/regime_pin` freezes 230 `lon/lat → uuid → curve`
  goldens and the layer-0..29 ownership ladder (`test_data/regime_pin*.tsv`), so
  a projection change can never again pass silently. Produced deliberately by
  `tools/gen_regime_pin.c`.
- **`test/pole_seam`** — 3,918 adversarial points against the hhg9 Python
  reference: polar approach from 1e-12° out, both poles from every seam
  meridian, seams straddled at ±1e-13/±1e-9 into the polar band, and the octant
  triple points. Covers the 88°–90° pole×seam band the Python suite reaches from
  neither side.
- **`tools/support/check_equal_area.py`** — independent validation of the
  equal-area claim using GeographicLib's exact geodesic `PolygonArea`. Measures
  ≤0.0013 % cell-area spread (measurement-limited) and exact closure against the
  WGS84 ellipsoid's surface area. Every other test compares this project against
  itself; this one does not.
- **`HEX9_VERSION` / `_MAJOR` / `_MINOR` / `_PATCH`** in `hex9_c.h`, so consumers
  can compare compile-time against runtime `hex9_version()`.
- **`postgis_hex9` refuses to load** against a libhex9 it was not built for
  (`_PG_init`). A major SOVERSION bump installs under a new filename and leaves
  its predecessor in place, so an un-rebuilt module would otherwise keep loading
  1.x and emitting the old regime's addresses indefinitely.
- **Build-time environment guards** in the extension Makefile: errors if
  `PG_CONFIG` points at a different PostgreSQL than the running server (a
  `make install` that succeeds into a prefix nothing reads), and if
  `POSTGIS_SRC` is a different PostGIS *minor* release than the server loaded
  (liblwgeom structs are resolved by `dlsym`, so a mismatch corrupts memory
  rather than failing to load).
- **`Makefile.local`** (gitignored, optional) for pinning machine-specific
  `PG_CONFIG` / `POSTGIS_SRC`.
- CMake asserts `hex9_c.h`'s `HEX9_VERSION` matches `project()`, so they cannot
  drift.

### Changed

- **Addressing chain** — geodetic latitude is now reduced to the **authalic**
  latitude (Karney series), the core runs on the **unit sphere**, and the warp
  is the **Sphere-L6 wedge-fold** field. One trained field serves every
  ellipsoid; the ellipsoid-specific part of the chain is the latitude series and
  nothing else.
- Library is **~40 % smaller** (32.5 MB → 19.6 MB) and warp init drops from
  ~13 s to ~1 s — the new blob ships its own gradients, so there is no
  gradient-estimation pass at startup.
- `hex9_version()` reports the real version. It returned a hardcoded `0.1.0`.
- `h9_version()` (SQL) derives from `postgis_hex9.control`. It returned a
  hardcoded `1.2.0` across three releases.
- Shared library `VERSION`/`SOVERSION` derive from `project()` instead of a
  stale hardcoded `0.1.0`/`0`.
- Warp-init failure in `postgis_hex9` is now an **ERROR**, not a warning. The
  identity fallback produces plausible values that are not Hex9 addresses;
  refusing to start is better than emitting them.
- macOS builds default to `CMAKE_OSX_DEPLOYMENT_TARGET=11.0`, matching common
  PostgreSQL distributions, so the extension no longer declares support for a
  macOS version its dependency cannot run on.
- `POSTGIS_SRC` has no default path. It is discovered from the running server's
  PostGIS version where possible, and otherwise reports what to set.

### Removed

- **The WGS84-trained warp regime**, its blob (`WGS84_l5_warp_f6.full.f64g`,
  12.8 MB), and the v2 mirrored ancestor. A second regime meant one point had
  two valid addresses with nothing to distinguish them.
- **`hex9_set_via_sphere()`** — there is one regime; there is nothing to select.
- **The `hex9.use_warp` GUC.** It was `PGC_USERSET` and changed the output of
  functions declared `IMMUTABLE`, which PostgreSQL may cache in functional
  indexes and generated columns — so two sessions could disagree about what a
  stored index entry meant. The capability remains in libhex9's C ABI
  (`hex9_set_use_warp`, for the hhg9 A/B parity harness) where it cannot reach
  SQL, and in no build is it reachable from a query.
- `set_use_warp` from the Python module.

### Fixed

- The `hex9.sql` regression test hard-coded Westminster's layer-8 label as an
  *input*, so under the new regime two assertions silently evaluated false.
  Literals updated; the assertions are live again. (A wholesale re-baseline
  would have frozen them false forever.)
- Edinburgh's literal in the same test was stale but never compared against an
  encoded point, so it failed nothing while being wrong.

### Migration

1. Rebuild and reinstall libhex9, **then rebuild and reinstall the extension** —
   this is not optional (see the `_PG_init` guard above).
2. Remove superseded `libhex9.0*` shared libraries. While they exist, a rebuild
   may silently link 1.x again, because `-L/usr/local/lib` can win the search
   order.
3. `ALTER EXTENSION postgis_hex9 UPDATE TO '2.0.0';`
4. Re-derive addresses from source geometry. Do **not** remap.

---

## [1.0.0] — 2026-07-08

First tagged release. Equal-area hexagonal DGGS with a geometry-free C ABI,
self-contained 16-byte UUID addresses, the L30 reclaimed layout, PostGIS
extension, PDAL filters, Python module and Rust FFI crate.

Releases between `v1.0.0` and 2.0.0 were extension-only version bumps
(`postgis_hex9` 1.1.0–1.5.0) and were not tagged in this repository; they added
the canonical parent/ancestor pair, Hamiltonian curve addressing, the `w_oct`
storage CRS and the PDAL filters.

[2.0.0]: https://github.com/MrBenGriffin/libhex9/releases/tag/v2.0.0
[1.0.0]: https://github.com/MrBenGriffin/libhex9/releases/tag/v1.0.0
