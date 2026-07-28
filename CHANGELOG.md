# Changelog

All notable changes to libhex9, the `postgis_hex9` extension, the PDAL filters
and the `hex9-sys` crate. These share one version number.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning is [semantic](https://semver.org/), with the addition that **any
change to the projection is a major version**, because it changes every stored
address.

---

## [2.2.0] — 2026-07-28

### Added — face-coordinate addressing (bring your own ⟨polyhedron⟩)

The ⟨face⟩⟨cell⟩ layers exposed directly at the warped chart — the
(cx, cy, oid) coordinates `hex9_project` / `hex9_woct_to_boct` emit — for
callers who supply their **own** sphere→octahedron projection (research
projections, AKW-comparison studies within the hex9 frame — and the frame
measuring them back). Descent is projection-blind; bins, lineage, k-rings
and the curve all work unchanged on the resulting addresses.

- `hex9_encode_boct` (+`_many`): mint addresses from face coordinates.
  Always the grid-canonical containment descent; entered at the warped
  chart, so `project → encode_boct` is **bit-identical** to `encode`
  (internally: `uuid_from_boct_full`, the post-warp half of the descent —
  the BRAW-in path now bridges to it via `h9_warp_inv`, byte-identical
  refactor, regime pins unchanged).
- `hex9_decode_boct` (+`_many`): the cell's lattice centroid in face
  coordinates — the same centroid `hex9_decode` unprojects, emitted before
  unprojection. Projection-free.
- `hex9_cell_ring_boct`: cell boundary in face coordinates, per-vertex
  frames (each vertex in the chart of its own oid) — render through your
  own inverse.
- Python: `encode_boct` / `decode_boct` / `cell_ring_boct` (vectorized).
- Doctrine (see `docs/projection.md`): addresses minted from a non-hex9
  projection are **non-canonical** — same 128-bit space, different meaning,
  no marker bits (deliberately). Projection identity is dataset metadata,
  exactly as the sphere/WGS84 datum is: never mix within one dataset.
  Guarantees split at the seam: ⟨face⟩⟨cell⟩ keeps the universality
  regime; everything upstream is the caller's.
- `test/boct_io.c` pins the seam: byte-equality with the canonical chain
  (2D route), working-depth agreement (3D w_oct route — the chart bridge is
  a rotation out and back, ULP-exact rather than bit-exact), centroid and
  ring coherence through `unproject`. Note for embedders: the boct surface
  speaks the *warped* chart, so `hex9_init` must have run — an
  uninitialised warp collapses BRAW == b_oct and hides frame errors.

## [2.1.0] — 2026-07-22

### Changed — deterministic encode chain (universality) — 2026-07-27

The floating-point program of the encode/decode chain is now **the definition
of the address**: the same lon/lat yields the bit-identical uuid (and curve)
on every conforming platform. The first multi-platform CI run proved that
libm/FMA variance flips deep nibbles for boundary-adjacent points (Linux and
Windows vs the macOS-generated goldens; primary source: compiler FMA
contraction on arm64, libm differences the residual). Universality cannot be
bought with accuracy — distance-to-boundary has no floor — only with
bit-identical arithmetic.

- `core/h9_det_math.h` (new): `sin`/`cos`/`tan`/`atan2`/`hypot` vendored
  verbatim from musl 1.2.5 (fdlibm lineage — see COPYRIGHT), header-only,
  `h9_`-prefixed. `sqrt`/`fmod` remain libm: IEEE requires them correctly
  rounded, so they are already deterministic. Guards: `FLT_EVAL_METHOD == 0`
  enforced at compile time; clang `FP_CONTRACT OFF` pragma in the header.
- `pow(x, 0.25)` eliminated from the chain → `h9_qroot(x)` =
  `sqrt(sqrt(x))` — two IEEE-exact operations, deterministic and faster.
- `-ffp-contract=off` pinned project-wide, and separately in the DuckDB
  extension tree (which compiles the core itself).
- `regime_pin` goldens re-derived once, deliberately, from a
  via_sphere-green pinned build. Every moved row diverges at nibble ≥ 26
  (leaf tail); the ownership ladder is byte-identical. Determinism evidence:
  `gen_regime_pin` output is byte-identical across `-O0` / `-O2` /
  `-O3 -march=native`; kernels are ≤ 1 ULP of Apple libm over the chain's
  domains (tan ≤ 3, near π/2); the via_sphere oracle and full suite pass;
  encode wall-time is unchanged. Cross-platform proof: the CI matrix runs
  `regime_pin` strict — its green **is** the universality claim.
- **Spec freeze.** Any change to the chain's arithmetic — kernels, operation
  order, contraction policy, warp blob — is a regime change: deliberate,
  versioned, goldens re-derived from source. Never "improved" in place.

*Versioning note.* Full-depth addresses of boundary-adjacent points move at
nibble ≥ 26 relative to pre-pin macOS builds. This stays within 2.1.0 rather
than forcing a major bump: no 2.x address has been published (v1.0.0 is the
only external tag), working-depth bins and both anchor addresses are
unchanged, and pre-pin full-depth tails were platform-dependent — there was
no single address to preserve. The universality invariant starts here.

### Added — sphere-datum entry points (additive; no existing address moves)

Ten `*_sphere` twins of the lon/lat entry points, for callers that own their
own datum (planetary authalic frames, celestial RA/dec): `hex9_encode_sphere`,
`hex9_decode_sphere`, `hex9_encode_many_sphere`, `hex9_decode_many_sphere`,
`hex9_project_sphere`, `hex9_unproject_sphere`, `hex9_project_many_sphere`,
`hex9_unproject_many_sphere`, `hex9_grid_create_sphere`,
`hex9_cell_ring_sphere`. Identical chain minus the WGS84 authalic reduction —
**not a second regime**: every address has been minted on the unit sphere
since 2.0.0, and the WGS84 functions are unchanged bit-for-bit (regime pins
untouched). The datum is part of the function's identity, never process
state; it is dataset metadata — never mix datums within one dataset. Doctrine
in `docs/warp-regimes.md` ("Two datums, one regime"); pinned by
`test/sphere_mode.c`.

The grid handle records its datum at create; `hex9_grid_cell_centroid` /
`hex9_grid_cell_ring` emit in the handle's datum with no new accessors.

Python: `sphere=False` keyword on `encode`, `decode`, `cell`, `grid` — no new
function names at that layer.

PostGIS: six SQL twins — `h9_encode_sphere`, `h9_encode_many_sphere`,
`h9_decode_sphere`, `h9_cell_sphere`, `h9_grid_sphere`, `h9_adaptive_sphere`
— as distinct IMMUTABLE functions (datum in the function identity; nothing
address-affecting is ever a GUC). Sphere geometries are emitted with
**SRID 0** (spherical degrees have no EPSG identity; the datum/body is
dataset metadata — `ST_SetSRID` yours if it has one; never tag them 4326).
`h9_adaptive_sphere` runs the identical digest (addresses never cross the
datum boundary); only the rendering differs — geom SRID 0, and **density in
value per steradian**: a layer-L cell's unit-sphere area (4π/(12·9^L) sr) is
intrinsic, so no body radius is needed; per-km² on a body = density × 4π /
body_area, caller-side (`h9_choropleth.py --sphere` follows the same rule).
No twins where SQL composition serves: label centroids =
`h9_decode_sphere(h9_parse_label(t))`. Extension packaged as
2.1.0 with a `2.0.0--2.1.0` upgrade script (additive, apply in place);
regress suite extended and green, upgrade path tested.

Deliberately **no** twins for composable surfaces: w_oct (pure rotation from
b_oct — compose `hex9_project_sphere` + `hex9_boct_to_woct`), label centroids
(`hex9_parse_label` + `hex9_decode_sphere`), grid centroids
(`hex9_grid_cell_id` + `hex9_decode_sphere`).

### Added — integer lattice identity (`hex9_cell_uv`)

`hex9_cell_uv` / `hex9_cell_uv_many` / `hex9_uv_units`: a cell's EXACT
integer lattice identity — centre key (the lattice hex centre, the mesh
anchor; deliberately distinct from `h9_decode`'s representative point) and
six canonical vertex pool keys (`(ia, ib, oid)`; on-boundary vertices get a
deterministic representative across the seam maps, so the same physical
vertex keys identically from every touching cell — shared-vertex meshes,
cf. hhg9 HexMesh). Datum-free by construction: upstream of b_oct, the warp,
and lon/lat — no `*_sphere` twins exist or ever will. `ext` flags the
12·3^layer seam-chain cells. Removes the decode→re-project round-trip
(through the warp and a Newton solve, twice) that address-side consumers
previously paid. Python: `cell_uv(uuid, layers)`, `uv_units()`. Pinned by
`test/cell_uv.c` (arithmetic parity with the geometric chain, canonical
sharing incl. the cone-point doubled neighbour, ext census).
`tools/h9_net.py` now renders from exact lattice vertices — seamless,
template-free, and datum-flag-free.

### Changed

- `hex9_init()` is the canonical init name (it always built the whole chain,
  not just the warp); `hex9_warp_init()` remains forever as an alias. All
  in-tree consumers now call `hex9_init`.

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
