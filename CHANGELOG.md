# Changelog

All notable changes to libhex9, the `postgis_hex9` extension, the PDAL filters
and the `hex9-sys` crate. These share one version number.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning is [semantic](https://semver.org/), with the addition that **any
change to the projection is a major version**, because it changes every stored
address.

---

## [2.3.0] — 2026-08-05

### Added — universality made citable, and curve labels reach the wheel

No address moves: everything here documents, checks, or exposes what the
chain already does.

- `docs/universality.md`: the normative statement of the bit-exact address
  contract — previously enforced by, but scattered across, `h9_det_math.h`,
  the build pins, `regime_pin`, `wheel_pin` and the CI comments. Defines the
  conforming platform, blesses `test_data/regime_pin.tsv` as the conformance
  interchange corpus (third parties conform by reproducing it bit-for-bit),
  and states the regime-change procedure (major version = regime name).
- Python `hex9.selftest()`: the installed wheel carries nine regime pins
  (landmarks, poles, antimeridian, an octant seam, s005/s030 — the historic
  FMA-flippers) and re-mints them on demand, so an embedder can prove *their*
  environment conforms — a nonconforming rebuild fails at first check, not as
  silent address drift in production. Gated per platform by `wheel_pin.py`.
- Python `curve_label` / `curve_parse_label`: bindings over the existing
  (and `test/curve.c`-pinned) `hex9_curve_label` / `hex9_curve_parse_label` —
  the curve side of the label surface was C-only; hhg9 and the wheel now
  agree. A full-depth curve label is the curve-uuid's hex spelling.
- Python `label` and `curve_label` are overloaded for batch input: a (16,)
  uuid returns one string, an (n,16) array returns a list of n — matching
  the batch shape of `encode`/`curve` they are fed from.
- Python `label()` now rejects curve-uuids (nibble 0 = 0xC) instead of
  silently minting a plausible-looking but meaningless cell label from curve
  ranks — the addresses-are-not-labels failure mode, made loud.
- `wheel_pin.py` additionally pins curve uuids, curve labels, the parse
  round-trip, and the marker guard on every wheel platform.
- Python `hex9.verbs`: cell-first grid verbs over the validated neighbour
  algebra — `aim` (bearing → neighbour, batchable), `walk_to` (A* shortest
  hex path, obstacles impassable), `vision_cone` (k-ring cone of visibility
  with great-circle occlusion), plus the bind helpers `encode_keys` /
  `encode_keys_multi` / `centroids` and spherical helpers `gc_bearing` /
  `gc_angle`. Ported from the hhg9 fox-and-rabbits PoC (`ex_0510`); bearings
  are advisory FP that steer but never mint — every returned key comes from
  the canonical encode/bin/neighbour chain, so the universality contract is
  untouched. Python-only in this take; SQL/DuckDB verb twins are a later
  take. Behavioural test `test/verbs.py` (ctest `verbs_py`) runs the field
  across the equator octant seam.

### Added — E4H: the aperture-4 structural tail (a new address surface)

E4H extends a hex9 host bin at attach layer L with an `0xE` break marker,
one HALF nibble (which state-cut trapezoid of the host hexagon), and up to
`28−L` tail nibbles from `{0, 1..5}` descending the aperture-4 half-hex
rep-4 carrier — exact nesting, straight edges, suffix-local truncation
(truncation IS binning, unlike the a9 body), globally matched pairs (the
two halves of every fine hexagon share their final digit). Digit semantics
are the five-symbol enumeration proven minimal-and-necessary in hhg9
(`docs/dggs-transport-tilings.md` §4b–§4d, CSP-verified). E4H uuids are
ADDRESSES (ruling 2026-08-04), and no existing h9/curve address moves.

- **The exact classifier.** The canonical E4H program is a frozen det-math
  FP seed (project → host lattice-centre frame → seam unfold), ONE snap at
  `2^-46`, then pure integer descent in ℤ[½, √3] (`core/h9_e4h.h`;
  128-bit components, generator-proven ≤109 bits; sign decisions by exact
  `A² vs 3B²` in 256-bit). Every half/child decision after the snap is
  exact, so E4H addresses cannot drift with depth, platform, or optimiser —
  the reference's open "all-depths stability" question is closed by
  construction. Constants are frozen from the hhg9 reference by
  `tools/gen_e4h_tables.py` (provenance commit recorded in the header).
- C ABI: `hex9_e4h_encode` / `_decode` / `_partner` (+ `_sphere` and
  `_many` twins), `_split`, `_bin` (tail truncation), `_depth`,
  `hex9_is_e4h`, `_label` / `_parse_label`
  (`<h9-label>E<half><digits>`). Full depth budget admitted from day one
  (Ben's ruling 2026-08-05): layer 0 carries a 28-digit tail.
- Marker guards: `0xE` cannot occur in any h9 or curve uuid, so
  `hex9_is_e4h` is decisive and EVERY h9/curve uuid-consuming entry point
  (decode, bin, parent/ancestor, curve family, k-rings, labels, cell_uv,
  rings, adaptive…) now rejects E4H input rather than mis-reading a tail.
- Conformance: `test_data/e4h_pin.tsv` — 504 hhg9-minted rows over all
  octants, seams, cone-point rings and poles, depths probing far past the
  CSP-verified range up to the full budget, knife-edge points margin-
  filtered — reproduced byte-exactly by `test/e4h_parity.c`. Structural
  laws (count law `2·4^d`, decode→encode identity, truncation=binning,
  budget boundary, grammar rejects, matched-pair involution, guards) in
  `test/e4h.c`, plus a GLOBAL census at host layer 1: all 108 hosts × all
  tail combinations at depths 0/1 — exactly 216/864 unique addresses
  (2·4^d per host), every one decode→encode stable, and the partner
  involution partitions them into exactly 108/432 fine hexagons (108×4 at
  depth 1 — addresses are half-trapezoids, hexagons are matched pairs;
  both counts pinned so they cannot be conflated). Universality doctrine
  addendum in `docs/universality.md`.
- Python: `hex9.e4h_encode/e4h_decode/e4h_partner/e4h_split/e4h_bin/`
  `e4h_depth/e4h_label/e4h_parse_label/is_e4h` (batch numpy shapes, sphere
  kwarg, label batch overload). `wheel_pin.py` gains four E4H pins (to
  depth 28) plus guard checks on every wheel platform; `hex9.selftest()`
  now re-mints 10 pins (9 h9/curve + 1 full-budget E4H).

### Added — hexagon binning: `e4h_hex`, the GIS surface (ruling 2026-08-05)

E4H addresses are HALF-hex trapezoids; the fine HEXAGONS are the matched
pairs, and most straddle parent/host boundaries. `hex9_e4h_hex` (+`_many`,
python `e4h_hex`, SQL/DuckDB `h9e_hex`) returns the CANONICAL HEXAGON KEY:
the pair's **mode-0** member, where mode is the parity of the descent's
rotation accumulator — pure integer arithmetic on the address, exposed as
`hex9_e4h_mode`/`e4h_mode`/`h9e_mode`. Both halves of every pair map to
one key (idempotent), and every host owns exactly 4^depth hexagons — the
aperture count, so host-level roll-ups of hexagon bins partition evenly.
Datum-free (the pairing and the mode are structural; no `_sphere` twins).
The mode formula was derived from hhg9's verified transport-mode ownership
doctrine and machine-checked against the geometric triangle-parity oracle
on 20,608 census addresses (11,840 pairs, 820 cross-seam): mode ≡ s mod 2
exactly, with no base term. Pinned by the layer-1 census in `test/e4h.c`
(108/432 canonical keys, 4-per-host balance), byte-pinned hexagon keys per
platform in `wheel_pin.py`, and law rows in pg_regress + sqllogictest.

### Fixed — 2.3.0 pre-release review (three-agent adversarial pass)

All three reviews returned positive verdicts (the exact classifier "sound
to freeze"; C verbs bit-equivalent to the python module under randomized
fuzz; SQL surfaces solid). Remediations:

- Curve-uuids are now rejected by `hex9_decode(_sphere/_many)`,
  `hex9_bin(_many)`, `hex9_label(_key)` — previously they fell through the
  resolution machinery and silently minted wrong-hemisphere points and
  plausible-looking fake bins (the marker-guard doctrine now covers BOTH
  foreign kinds; the verbs inherit the guard through their centroid path).
- `hex9_disk_ncells` caps k (int64 overflow at k ≳ 1.75e9 could bypass
  callers' 60M-cell guards — UB, though provably no memory corruption).
- `hex9_walk_to` rejects wrong-form src/dest (a full uuid or wrong-layer
  bin previously burned max_expand and read as "no path"), and caps
  max_expand at 1e6 (heap-allocation exposure at absurd values).
- E4H batch forms report 0/1 (per-item codes blended under OR: 1|2 read
  as the seam code 3).
- The generator's descent bit bound now covers the runtime guard box
  (109 → 113 bits; still 14 bits of __int128 margin — comment-level).
- `hex9.selftest()` round-trip now checks longitude (cos-scaled, wrapped),
  not just latitude; one exact-FP bearing pin in each SQL test relaxed to
  1e-9 tolerance (libm-portability).

### Added — take 2: the verbs and E4H reach SQL (PostGIS + DuckDB)

- **C verbs ABI** — `hex9_aim`, `hex9_walk_to`, `hex9_vision_cone`
  (hex9_c.h): the C twins of the python wheel's `hex9.verbs`, cell-first
  (bin keys in and out, full-uuid re-derivation internal), bearings advisory
  FP. A* matches the python module's heap order exactly (f, then g, then key
  bytes); occlusion samples sight lines at a third of the local hex pitch.
  Behavioural laws in `test/verbs_abi.c` (aim-inverse, cone ⊆ disk,
  full-circle == disk, seam crossing, occlusion front-face, walk adjacency /
  detour / sealed-destination, E4H guards).
- **PostGIS `postgis_hex9` 2.3.0** (upgrade path 2.1.0 → 2.3.0; 2.2.x had
  no SQL surface): `h9e_encode(/_sphere)`, `h9e_decode(/_sphere)`,
  `h9e_partner(/_sphere)`, `h9e_bin`, `h9e_depth`, `h9_is_e4h`, `h9e_host`,
  `h9e_label`, `h9e_parse_label`; verbs `h9_aim`, `h9_walk_to` (uuid[] path,
  NULL when no path), `h9_vision_cone` (SETOF uuid). pg_regress rows pin the
  Edinburgh E4H mint byte-exactly and assert the structural/verb laws; the
  0xE guard errors are pinned as errors. `postgis_hex9--2.1.0.sql` is now a
  frozen artifact; `hex9.sql.in` generates `--2.3.0.sql`.
- **True-bearing correction** for physical/material headings, on every
  surface: python `hex9.verbs.from_true_bearing`/`to_true_bearing`, C
  `hex9_bearing_from_true`/`_to_true`, SQL/DuckDB
  `h9_bearing_from_true`/`_to_true`. A measured true-north heading is a
  WGS84 geodesic azimuth; the verbs' bearings are spherical trig over
  geodetic coordinates — the difference is the flattening term
  `tan(b_verb) = (M/N)·tan(b_true)` at the heading's latitude (identity at
  poles and cardinal bearings, ~0.19° max at equatorial diagonals). Apply
  before `aim`/`vision_cone` when cone tolerance is that tight. Magnetic
  headings must already be reduced to true north — declination is
  epoch/model-dependent, i.e. dataset metadata, deliberately out of scope.
- **DuckDB `hex9` extension**: `h9e_encode/h9e_decode(_wkb)/h9e_partner/`
  `h9e_bin/h9e_depth/h9_is_e4h/h9e_host/h9e_label/h9e_parse_label` +
  `h9_aim`, `h9_walk_to(src, dest, layer[, obstacles])` → LIST(UUID) (NULL
  when no path), `h9_vision_cone(key, layer, bearing, half, k[, obstacles])`
  → LIST(UUID). sqllogictest `test/sql/h9_e4h.test` pins the same Edinburgh
  mint and laws.

## [2.2.1] — 2026-07-28

### Fixed — packaging and metadata only; no code, no address moves

- PyPI storefront: the `hex9` package now carries its own README
  (`python/README.pypi.md` — the Python package's contract, including what
  it is *not*: it contains no hhg9 code) instead of libhex9's C-library
  README, which 2.2.0 shipped by accident.
- `hex9_c.h`: corrected a stale header comment claiming GPL-2.0-or-later —
  the work is Apache-2.0 (LICENSE/COPYRIGHT are canonical); Apache-2.0 is
  one-way GPL-compatible, so the PostGIS extension never needed the core to
  be GPL.
- TestPyPI publishes tolerate re-dispatch at an existing version
  (`skip-existing`); real PyPI stays strict.

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
