# postgis_hex9

A PostGIS extension for the **Hex9 (H9) Discrete Global Grid System** — an
equal-area hexagonal global grid with self-contained UUID cell identifiers and a
9-fold (3×3) hierarchy from layer 0 (12 base cells) down to layer 30
(~32 nm cell diameter). The legacy `USE_L29` build layout stops at layer 29
(~95 nm); the default reclaims the old terminal nibble as a real body digit, so
a full UUID *is* its own deepest-layer bin. The deepest layer is reported by
`h9_version()`/`hex9_lmax()`.

Each cell ID is a standard `uuid`, so it indexes, joins, and `GROUP BY`s with
zero custom types. An authalic (equal-area) warp keeps cell areas uniform to
~0.014 % across the globe.

---

## Install

```sql
CREATE EXTENSION postgis;          -- required dependency
CREATE EXTENSION postgis_hex9;
```

Build from source (standalone / development):

```sh
cd extension/postgis_hex9
make && sudo make install
make installcheck        # needs a PostgreSQL with PostGIS available
```

> The module resolves liblwgeom symbols at load time and works in any backend
> regardless of whether a PostGIS function was called first — no `shared_preload_libraries`
> or `LOAD` needed.

---

## Concepts

| Term | Meaning |
|------|---------|
| **layer** | Resolution level `0`–`30` (`0`–`29` on the legacy `USE_L29` layout). Each step is a 3×3 (9×) subdivision. Layer 30 ≈ 32 nm. |
| **full UUID** | From `h9_encode()`. Carries the exact point context — losslessly reversible via `h9_decode()`. |
| **bin UUID** | From `h9_bin(uuid, layer)`. All points in the same layer-`L` cell collapse to one UUID — the cell key for aggregation. (Lossy below `L`.) |
| **cell** | The hexagon for one UUID at a layer — `h9_cell()` (one) or `h9_grid()` (many over an area). |

---

## Quick start

```sql
-- 1. Encode a point to its layer-29 cell UUID, and round-trip back.
SELECT h9_encode('SRID=4326;POINT(-3.19 55.95)'::geometry);   -- → uuid
SELECT ST_AsText(h9_decode(h9_encode('SRID=4326;POINT(-3.19 55.95)'::geometry)));

-- 2. Bin points into layer-8 cells and count (spatial GROUP BY).
SELECT h9_bin(h9_encode(geom), 8) AS cell, count(*)
FROM   my_points
GROUP  BY cell;

-- 3. Draw the hexagon for a cell.
SELECT h9_cell(h9_bin(h9_encode(geom), 8), 8) FROM my_points;

-- 4. Tile an area with a layer-8 grid: full identity, bin key, polygon, centroid.
SELECT h9_id, h9_bin, ST_AsText(geom), ST_AsText(centroid)
FROM   h9_grid(ST_MakeEnvelope(-0.2, 51.4, 0.0, 51.6, 4326), 8);
```

A bin UUID makes a clean **functional index** or **generated column** (all
functions are `IMMUTABLE STRICT`):

```sql
ALTER TABLE my_points
  ADD COLUMN cell8 uuid GENERATED ALWAYS AS (h9_bin(h9_encode(geom), 8)) STORED;
CREATE INDEX ON my_points (cell8);
```

---

## Function reference

| Function | Returns | Purpose |
|----------|---------|---------|
| `h9_encode(geometry)` | `uuid` | Encode a POINT (any SRID transformable to 4326) to its deepest-layer cell UUID (layer 30; 29 on the legacy layout). |
| `h9_encode_many(geometry[])` | `uuid[]` | Batch `h9_encode`, one OpenMP-parallel pass, input order preserved; NULL elements → NULL. *(1.3.0)* |
| `h9_decode(uuid)` | `geometry` | Decode a UUID back to a POINT (SRID 4326). Round-trips `h9_encode` to ~32 nm (~95 nm on the legacy layout). |
| `h9_bin(uuid, layer)` | `uuid` | Cell-key UUID at `layer` 0–30. Equal for all points sharing that cell — use for `GROUP BY` / indexes. |
| `h9_cell(uuid, layer, densify DEFAULT 0)` | `geometry` | Hexagon polygon (SRID 4326) for a UUID at `layer` 1–30. |
| `h9_grid(bounds, layer, densify DEFAULT 0)` | `TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)` | Every cell at `layer` 1–30 whose centre lies in `bounds`. `h9_id` = full reversible identity; `h9_bin` = grouping key. *(h9_id + rename, 1.4.0)* |
| `h9_label(uuid, layer)` | `text` | Compact body-nibble label, e.g. `478232778`. |
| `h9_label_key(uuid, layer)` | `text` | `h9_label` plus the `.key_tail` nibble, e.g. `478232778.9`. |
| `h9_parse_label(text)` | `uuid` | Canonical bin UUID from a bare (`h9_label`) or keyed (`h9_label_key`) label. *(1.2.0)* |
| `h9_label_centroid(text)` | `geometry` | Exact geographic centroid of the labelled cell. *(1.2.0)* |
| `h9_common_ancestor(uuid[], layer)` | `(label, h9_bin, layer)` | Deepest common address ancestor of the cells; NULL when they span L0 hexes. *(1.2.0; `hex9`→`h9_bin` 1.4.0)* |
| `h9_neighbors(uuid, layer)` | `SETOF uuid` | The (up to 6) edge-adjacent cells. *(1.2.0)* |
| `h9_kring(uuid, layer, k)` | `SETOF uuid` | Cells at graph distance exactly `k` (6k for an interior cell). *(1.2.0)* |
| `h9_kdisk(uuid, layer, k)` | `SETOF uuid` | Cells at distance ≤ `k`, centre included (nominally `1+3k(k+1)`). *(1.2.0)* |
| `h9_adaptive(uuid[], float8[], min_layer, max_layer, ceiling, floor DEFAULT 0)` | `TABLE(h9_bin, layer, value, npoints, density, grade, geom)` | Population-digest multi-layer grid: bottom-up digestion under a per-cell value ceiling; `sum(value)` = total input weight. `weights` may be NULL (all 1). `density`/`grade` are cartography columns (see below). *(1.2.0; uuid[]+geom 1.2.x; density/grade + `hex9`→`h9_bin` 1.4.0)* |
| `h9_version()` | `text` | Extension version + build stamp. |

Adjacency notes: the k-family is symbolic (exact integer lattice arithmetic,
no floating point). **Input must be a full UUID from `h9_encode`** — bin
input is rejected (a bin is a layer-scoped key, not an address; resolving
one is approximate, and the approximation risks failed conclusions from an
unknown underlying error). Outputs are sorted canonical bin UUIDs at the
requested layer — keys for joining (e.g. to `h9_grid`), not addresses for
further traversal; to walk onward, encode a point and pass the full UUID.
The 12 half-hex cells per layer at the octahedron vertices have 5
neighbours, and rings/disks covering them are correspondingly smaller than
the nominal counts.

### `h9_adaptive` density & grade

The digest equalises *value*, and its cells are equal-area per layer, so two
convenience columns are derived from `value` and `layer`:

- **`density`** = `value · 12 · 9^layer / 510065622` — persons/km² (510 065 622 km²
  is the Earth's surface; `12·9^layer` is the cell count at that layer). It is
  **exact for the digest**.
- **`grade`** = `layer + ln(value)/ln(9)` — a log₉ graduation for symbolising
  (`+1` ⇒ 9× denser); `NULL` when `value ≤ 0`.

The numbers are exact arithmetic, but their **real-world reading carries the
source's error bars**, not the formula's: `density` is only as meaningful as the
sample behind `value`. In particular `npoints = 1` is a **point-mass reading**
(one source point fell in the cell, areally unsupported) — treat it as an
honesty flag and symbolise/outline those cells. Quantised inputs (e.g. census
allocated over footprints) and cells where the `floor`/`ceiling` binds likewise
propagate into `density`. Render with `ORDER BY layer ASC` so coarser parents
draw under their finer children.

### `densify` (rendering large cells)

`h9_cell` and `h9_grid` take an optional `densify` (default `0`) that subdivides
each of the 6 hex edges into `3^densify` segments, so edges follow the H9 lattice
instead of straight `(lon, lat)` chords. Ring size is `6·3^densify + 1`
(`0`→7, `1`→19, `2`→55, `3`→163 points). Useful at coarse layers (≤ 3); at
layer ≥ 5 the corner-only form is normally fine. Constraints: `densify ≥ 0`,
`layer + densify ≤ 30` (≤ 29 on the legacy layout), `densify ≤ 9`.

```sql
-- A layer-2 grid over Europe with smooth (densified) hex edges:
SELECT h9_bin, geom
FROM   h9_grid(ST_MakeEnvelope(-10, 35, 30, 60, 4326), 2, 3);
```

---

## Caveats

- **Bin UUIDs are lossy below their layer.** `h9_cell` on a bin UUID may collapse
  distinct cells to the same polygon — for per-cell rendering over an area prefer
  `h9_grid`, which keeps non-lossy state internally.
- **`h9_grid` is capped** by the `hex9.grid_max_cells` GUC (default 708 588). For a
  large area, use a coarser layer/smaller bounds or raise the cap:
  `SET hex9.grid_max_cells = 5000000;`
- **Antimeridian / poles:** cell polygons spanning the antimeridian are returned
  unsplit; split with `ST_Split`/your own logic if your renderer needs it.
- **Ancestry is address ancestry, not geometric containment.** Hexagon children
  straddle their parent's boundary, so a point's layer-`L` bin cell (and
  `h9_common_ancestor` results) need not geometrically contain the point —
  the cell that *contains* it at `L` may be an address cousin. Use `h9_grid`
  when you need the geometric cover.
- **Bins and labels are not addresses.** Only the full UUID from `h9_encode`
  is a guaranteed-reversible address. A bin UUID guarantees uniqueness per
  (hexagon, layer) — a grouping/join key, nothing more — and a label is a
  human-readable name. Operations that treat bins or bare labels as
  addresses (`h9_decode` on a bin, re-binning a bin coarser, `h9_parse_label`
  on a bare body, `h9_cell`/k-family resolution of meta-bearing bins) *often*
  work but are **not guaranteed**: the split hexagons (digits 6/7/8) need
  meta the bin/label no longer carries. These paths are kept as documented
  fossils; the failure modes, worked examples, and the guaranteed patterns
  (re-bin from the full UUID; join bins to `h9_grid` for geometry) are
  catalogued in `docs/addressing-doctrine.md`.

## Settings (GUCs)

| GUC | Default | Effect |
|-----|---------|--------|
| `hex9.grid_max_cells` | `708588` | Max cells `h9_grid()` will emit before erroring. |
| `hex9.use_warp` | `on` | Apply the authalic (equal-area) warp. `off` = raw octahedral projection. |

`h9_encode` uses the **containment** (grid-canonical) encoder — the one that
agrees with `h9_grid`/`h9_id` at every layer. It is deliberately **not** a GUC:
`h9_encode` is `IMMUTABLE` (used in functional indexes and generated columns),
so its result must not depend on a session setting. (The legacy NN encoder
differs only in the sub-metre leaf tail, ~3.7% of full UUIDs and 0% of bins
through L20; it remains in libhex9 for research but is unreachable from SQL.)
