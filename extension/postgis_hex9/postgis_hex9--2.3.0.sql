-- hex9.sql.in — SQL declarations for the postgis_hex9 extension.
--
-- Part of the Hex9 (H9) Project
-- Copyright ©2025, Ben Griffin
-- Licensed under the Apache License, Version 2.0
--
-- Processed by the build system: MODULE_PATHNAME is replaced with the
-- actual shared-library path before installation.

-- NOTE: no `LOAD '$libdir/postgis-3'` is needed. The module resolves liblwgeom
-- entry points at runtime via h9_lwgeom_resolve() in _PG_init (see
-- h9_lwgeom_shim.h), which load_file()s PostGIS itself — so the extension works
-- in any backend regardless of whether a PostGIS function was called first.

-- ── Version ──────────────────────────────────────────────────────────────────

-- h9_version() → text
--   Returns the extension version and build timestamp.
--   Use this to confirm which binary is loaded after an install.
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_version()
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9_version'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

-- h9_lmax() → integer
--   Deepest addressable layer of the loaded libhex9: 30 for the default
--   (reclaimed) layout, 29 on a legacy HEX9_USE_L29 build. The value is
--   compiled into the linked library, so this reports the layout the extension
--   is actually running — use it to confirm an L30 install:
--       SELECT h9_lmax();   -- 30
--
-- Availability: Hex9 1.4.0
CREATE OR REPLACE FUNCTION h9_lmax()
    RETURNS integer
    AS 'MODULE_PATHNAME', 'h9_lmax'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

-- ── Encode ───────────────────────────────────────────────────────────────────

-- h9_encode(geometry) → uuid
--   Encodes a POINT to a self-contained UUID at the deepest layer (30; 29 on
--   the legacy USE_L29 layout). See h9_lmax().
--   The UUID carries all decode context internally; no companion byte needed.
--   Suitable for exact round-trip via h9_decode(), indexing, and h9_bin().
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_encode(geometry)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_encode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_encode_sphere(geometry) → uuid
--   Sphere-datum twin of h9_encode: the POINT's coordinates are
--   ALREADY-SPHERICAL lon/lat degrees (another body's authalic frame,
--   celestial RA/dec) — the WGS84 authalic reduction is skipped; everything
--   else is the identical chain. The datum is part of the function's
--   identity, never a setting: a WGS84-minted and a sphere-minted address
--   for the same numeric lon/lat differ below ~layer 5, and the 16 bytes do
--   not record which. The datum (and body) is DATASET metadata — never mix
--   datums within one dataset. See docs/warp-regimes.md,
--   "Two datums, one regime".
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_encode_sphere(geometry)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_encode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_encode_many(geometry[]) → uuid[]
--   Batch form of h9_encode: encodes an array of POINTs to their layer-29
--   UUIDs in input order, in a single OpenMP-parallel pass. Crosses the
--   SQL/C boundary once and parallelises the (independent) point work — the
--   fast path for encoding a whole column. NULL elements yield NULL UUIDs
--   (position preserved); a non-POINT element raises an error.
--
--   Composes with the digest: h9_adaptive(h9_encode_many(...), ...).
--
--   Example (encode a column, then digest it):
--     SELECT a.h9_bin, a.layer, a.value, a.density, a.geom
--     FROM h9_adaptive(
--              h9_encode_many((SELECT array_agg(geom) FROM dwellings)),
--              (SELECT array_agg(occupants::float8) FROM dwellings),
--              6, 12, 500, 100) AS a;
--
-- Availability: Hex9 1.3.0
CREATE OR REPLACE FUNCTION h9_encode_many(geometry[])
    RETURNS uuid[]
    AS 'MODULE_PATHNAME', 'h9_encode_many'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_encode_many_sphere(geometry[]) → uuid[]
--   Batch form of h9_encode_sphere (see its datum doctrine).
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_encode_many_sphere(geometry[])
    RETURNS uuid[]
    AS 'MODULE_PATHNAME', 'h9_encode_many_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- ── Decode ───────────────────────────────────────────────────────────────────

-- h9_decode(uuid) → geometry
--   Decode a UUID to a POINT geometry (SRID 4326).
--   Exact round-trip: h9_decode(h9_encode(pt)) recovers the input point
--   to within the resolution of the deepest layer (layer 30 ≈ 32 nm; layer 29
--   ≈ 95 nm on the legacy layout).
--
--   GUARANTEED for full UUIDs only. Bins are layer-scoped grouping keys,
--   not addresses; decoding one is a fossil convenience that mis-locates
--   meta-bearing (split-hex 6/7/8 / seam-flavour) cells — see
--   docs/addressing-doctrine.md (F2). For bin-keyed geometry join to
--   h9_grid instead.
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_decode(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_decode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

-- h9_decode_sphere(uuid) → geometry
--   Sphere-datum twin of h9_decode: emits the representative point in
--   spherical degrees, SRID 0 (unknown) — spherical coordinates have no
--   EPSG identity here; the datum/body is dataset metadata. ST_SetSRID it
--   yourself if your dataset has one. Only meaningful for addresses MINTED
--   with h9_encode_sphere: the address does not record its datum, so
--   decoding a WGS84-minted uuid here returns coordinates in the wrong
--   frame with no error — keeping the datums straight is the caller's
--   contract. (h9_label_centroid has no twin: compose
--   h9_decode_sphere(h9_parse_label(t)).)
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_decode_sphere(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_decode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

-- ── Hierarchy ────────────────────────────────────────────────────────────────

-- h9_bin(uuid, integer) → uuid
--   Returns the bin-key UUID at the given layer (0..lmax; lmax = 30, or 29 on
--   the legacy USE_L29 layout).
--   All points in the same H9 cell at layer L produce the same output UUID,
--   making this safe for GROUP BY aggregation and spatial binning.
--   Body nibbles above layer are replaced with the 0xF sentinel; nibble[31]
--   stores the backward-pass context needed by h9_cell().
--   IMMUTABLE and STRICT: usable in GENERATED STORED columns and functional
--   indexes without any performance penalty.
--
--   GUARANTEED for full-UUID input (byte-identical to the Python
--   reference). Re-binning a bin at its own layer is the identity;
--   re-binning a bin COARSER is a fossil — unguaranteed at split-hex
--   ancestry, and can emit an invalid all-sentinel UUID. Always re-bin
--   from the full UUID. See docs/addressing-doctrine.md (F3).
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_bin(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_bin'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- ── Labels ───────────────────────────────────────────────────────────────────

-- h9_label(uuid, integer) → text
--   Returns the body nibbles 0..layer as a compact text label.
--   Digits 0..9 map to '0'..'9'; 10 → 'a', 11 → 'b'.
--   Example: h9_label(h9_bin(h9_encode(pt), 8), 8) → '478232778'
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_label(uuid, integer)
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- h9_label_key(uuid, integer) → text
--   Like h9_label but appends '.' and the key_tail nibble (nibble[31]).
--   The key_tail encodes (p_mo, p_c2, r_mo) — enough to uniquely identify
--   a point within its layer-L cell for exact decode.
--   Example: h9_label_key(h9_encode(pt), 8) → '478232778.9'
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_label_key(uuid, integer)
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9_label_key'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- ── Cell geometry ────────────────────────────────────────────────────────────

-- h9_cell(uuid, integer, integer DEFAULT 0) → geometry
--   Returns the hexagonal cell polygon (SRID 4326) for the H9 UUID at the
--   given layer (1..29). Accepts either:
--     * full UUIDs from h9_encode() — fully reversible to (fx, fy);
--     * bin UUIDs from h9_bin()    — LOSSY: different bin UUIDs in the same
--       equivalence class decode to the same cell. For per-cell rendering
--       use h9_grid (which keeps non-lossy state internally) — h9_cell on
--       bin UUIDs may collapse distinct grid rows to the same polygon.
--
--   Third arg `densify` (default 0) subdivides each of the 6 hex edges into
--   3^densify segments using the recursive 1/3-step rule from
--   HexMesh._get_verts. Intermediate vertices land on the H9 lattice at
--   (layer + densify); output ring has 6·3^densify + 1 points.
--
--     densify = 0 → 7 points (corners only — fast path, matches the 1.0.0 form).
--     densify = 1 → 19; densify = 2 → 55; densify = 3 → 163; ...
--
--   densify > 0 is intended for rendering large hexes (layer ≤ 3) where
--   straight-line edges in (lon, lat) visibly diverge from the lattice
--   geometry. For layer ≥ 5 the corner-only output is usually sufficient.
--
--   Errors: layer outside 1..29; densify < 0; layer + densify > 29; densify > 9
--   (last is a soft cap on ring size, 6·3^9+1 = 118 099 points).
--
--   Example:
--     SELECT h9_cell(h9_bin(h9_encode(geom), 8), 8) FROM points;      -- corners only
--     SELECT h9_cell(h9_bin(h9_encode(geom), 2), 2, 3) FROM points;   -- 163-pt L2 hex
--
-- Availability: Hex9 1.1.0 (densify arg added; supersedes the 1.0.0 two-arg form)
CREATE OR REPLACE FUNCTION h9_cell(uuid, integer, integer DEFAULT 0)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_cell'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 200;

-- h9_cell_sphere(uuid, integer, integer DEFAULT 0) → geometry
--   Sphere-datum twin of h9_cell: the polygon's vertices are spherical
--   degrees, SRID 0 (see h9_decode_sphere for the datum contract).
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_cell_sphere(uuid, integer, integer DEFAULT 0)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_cell_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 200;

-- h9_grid(geometry, integer, integer DEFAULT 0)
--     → TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)
--   Set-returning function: yields one row per H9 cell at the given layer
--   (1..29) whose geographic centre falls within `bounds`.
--
--   Two UUID columns:
--     * h9_id  — the cell's FULL reversible identity (layer-29 UUID). This is
--                the address: h9_bin(h9_id, L') recovers the correct bin at
--                every L' ≤ layer, so it is the safe key to persist and to
--                re-bin to any coarser layer (see docs/addressing-doctrine.md
--                F4). decode(h9_id) lands ~0.1 cell-radius off centre — use
--                the `centroid` column for the geographic centre.
--     * h9_bin — the layer-scoped grouping key (matches
--                h9_bin(h9_encode(pt), layer) for any point in the cell).
--                Equal/GROUP BY/join key, NOT an address: re-binning a bin to
--                a coarser layer is unguaranteed (the F3 fossil). Re-bin from
--                h9_id instead.
--
--   Argument order changed in 1.1.0: (bounds, layer, [densify]) — was
--   (layer, bounds) in 1.0.0. Aligns with PostGIS convention (geometry
--   first). The 1.0.0 two-arg form is dropped in this version.
--
--   Third arg `densify` (default 0) applies the same structural densification
--   as h9_cell: each cell's 6 edges are subdivided into 3^densify segments.
--   Polygon ring size: 6·3^densify + 1.
--
--     densify = 0 → 7-point polygons (fast path; pre-computed BFS corners).
--     densify > 0 → re-decodes each cell's UUID through h9_bary_to_lonlat for
--                   every intermediate vertex. Per-row cost rises with 3^densify.
--
--   Implementation: integer-UV supercell BFS with strict octant-region
--   pruning + v4/v5 seam reflection — same algorithm Python
--   HexMesh.create_clipped uses. UUIDs are derived via the containment-based
--   xy_regions encoder. Centroid-in-polygon test clips to the actual bounds
--   geometry (not just its bbox).
--
--   Safety cap: hex9.grid_max_cells GUC (default 708 588 = 12 × 9⁵). Raises
--   an error if the estimated cell count exceeds this; use a coarser layer
--   or smaller bounds for large areas.
--
--   The `centroid` column is the cell's geographic centre as a POINT
--   (SRID 4326), computed by the BFS in the cell's own frame — useful for
--   labelling, clustering, or as a join key without re-decoding the UUID.
--   (Unlike decode(h9_id), it is the exact geometric centre.)
--
--   Example:
--     SELECT h9_id, h9_bin, ST_AsText(geom), ST_AsText(centroid)
--     FROM h9_grid(ST_MakeEnvelope(-0.2, 51.4, 0.0, 51.6, 4326), 8);
--
--     -- Render an L2 grid over Europe with densified polygons:
--     SELECT h9_bin, ST_AsText(geom)
--     FROM h9_grid(ST_MakeEnvelope(-10, 35, 30, 60, 4326), 2, 3);
--
-- Availability: Hex9 1.4.0 (h9_id identity column added; bin column renamed
-- hex9 → h9_bin). Was 1.1.0 (hex9 uuid, geom, centroid).
CREATE OR REPLACE FUNCTION h9_grid(geometry, integer, integer DEFAULT 0)
    RETURNS TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)
    AS 'MODULE_PATHNAME', 'h9_grid'
    LANGUAGE 'c' STABLE STRICT PARALLEL SAFE
    ROWS 1000 COST 5000;

-- h9_grid_sphere(geometry, integer, integer DEFAULT 0)
--   Sphere-datum twin of h9_grid: the bounds geometry is read as spherical
--   lon/lat degrees, and every emitted geom/centroid is spherical, SRID 0.
--   The enumeration is the identical chain minus the WGS84 authalic
--   reduction — same cells h9_encode_sphere addresses into (see
--   h9_encode_sphere for the datum doctrine).
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_grid_sphere(geometry, integer, integer DEFAULT 0)
    RETURNS TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)
    AS 'MODULE_PATHNAME', 'h9_grid_sphere'
    LANGUAGE 'c' STABLE STRICT PARALLEL SAFE
    ROWS 1000 COST 5000;

-- ── Adjacency: neighbours / k-ring / k-disk ──────────────────────────────────
-- Symbolic adjacency on the H9 mesh — exact integer lattice arithmetic, no
-- floating point. INPUT MUST BE A FULL UUID from h9_encode — bin input is
-- rejected with an error: a bin is a layer-scoped key, not an address (its
-- key tail cannot carry the meta the resolution needs, and the approximation
-- risks failed conclusions from an unknown underlying error). Output cells
-- are canonical bin UUIDs at the requested layer, sorted — keys for joining
-- (e.g. to h9_grid), not addresses for further traversal; to walk onward,
-- encode a point and pass the full UUID. Every cell has 6 neighbours except
-- the 12 half-hex cells per layer that make up the 6 octahedron-vertex
-- hexagons; those have 5 (4 sides + the partner half-hex, which shares the
-- same hexagonal outline).

-- h9_neighbors(uuid, integer) → SETOF uuid
--   The (up to 6) edge-adjacent cells of the given cell at layer (1..29).
--
--   Example:
--     SELECT h9_cell(n, 8) FROM h9_neighbors(h9_encode(pt), 8) AS n;
--
-- Availability: Hex9 1.2.0 (full-UUID input enforced since 1.2.x)
CREATE OR REPLACE FUNCTION h9_neighbors(uuid, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_neighbors'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 6 COST 50;

-- h9_kring(uuid, integer, integer) → SETOF uuid
--   Cells at graph distance exactly k from the given cell at layer (1..29).
--   k = 0 yields the cell itself; a full ring has 6k cells (fewer when the
--   ring covers an octahedron-vertex hexagon).
--
-- Availability: Hex9 1.2.0
CREATE OR REPLACE FUNCTION h9_kring(uuid, integer, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_kring'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 12 COST 100;

-- h9_kdisk(uuid, integer, integer) → SETOF uuid
--   Cells at graph distance at most k, including the centre cell — the
--   filled disk; nominally 1 + 3k(k+1) cells.
--
--   Example (all cells within 2 steps of a point's L8 cell):
--     SELECT d FROM h9_kdisk(h9_encode(pt), 8, 2) AS d;
--
-- Availability: Hex9 1.2.0
CREATE OR REPLACE FUNCTION h9_kdisk(uuid, integer, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_kdisk'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 19 COST 100;

-- ── Label parsing / mesh prefix ops ──────────────────────────────────────────

-- h9_parse_label(text) → uuid
--   Recover a bin UUID from a label. Accepts both forms:
--     * keyed labels from h9_label_key() ('.k' tail) — parsed SYNTACTICALLY
--       as the exact inverse: h9_parse_label(h9_label_key(u, L)) =
--       h9_bin(u, L) is an identity (the keyed label IS the bin in text
--       form; no resolution is performed);
--     * bare labels from h9_label() — resolved to the canonical bin by
--       tail search + canonical re-encode verification.
--   The layer is implicit in the label length (body chars = layer + 1).
--   Errors on strings that are not valid H9 labels.
--
--   FOSSIL CAVEAT: bare labels are NOT unique at split-hex (6/7/8) bodies —
--   the same body names two cells and the parse silently returns one of
--   them. Labels are names, not addresses; only the full UUID is
--   guaranteed recoverable. See docs/addressing-doctrine.md (F1).
--
--   Example:
--     SELECT h9_cell(h9_parse_label('435878133'), 8);
--
-- Availability: Hex9 1.2.0
CREATE OR REPLACE FUNCTION h9_parse_label(text)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_parse_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_label_centroid(text) → geometry
--   Geographic centroid (POINT, SRID 4326) of the labelled cell — the same
--   convention as h9_grid's centroid column and h9_decode on a bin UUID
--   (including the half-hex 4-vertex mean at octahedron vertices).
--
-- Availability: Hex9 1.2.0
CREATE OR REPLACE FUNCTION h9_label_centroid(text)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_label_centroid'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_common_ancestor(uuid[], integer) → (label text, h9_bin uuid, layer integer)
--   Deepest common ancestor (in the ADDRESS hierarchy) of the given cells,
--   all treated at the given layer (0..29). Returns the ancestor's label,
--   bin UUID, and layer; NULL when the cells span L0 hexes. The common
--   prefix lets a mesh be stored as one ancestor label plus per-cell
--   suffixes (label chars ancestor_layer+1 .. layer).
--
--   NOTE: ancestry is descent containment — hexagon children straddle their
--   parent's geometric boundary, so this is not exact geometric containment.
--   The returned UUID is recovered by bare-prefix parse and is only
--   trustworthy when the ancestor body avoids split-hex (6/7/8) ambiguity;
--   the label and layer are always right (docs/addressing-doctrine.md, F1).
--
--   Example:
--     SELECT (h9_common_ancestor(array_agg(h9_bin), 8)).*
--     FROM h9_grid(ST_MakeEnvelope(-0.14, 51.49, -0.11, 51.52, 4326), 8);
--
-- Availability: Hex9 1.2.0
CREATE OR REPLACE FUNCTION h9_common_ancestor(
        uuid[], integer,
        OUT label text, OUT h9_bin uuid, OUT layer integer)
    RETURNS record
    AS 'MODULE_PATHNAME', 'h9_common_ancestor'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 100;

-- ── Adaptive multi-layer grid (population digest) ────────────────────────────

-- h9_adaptive(uuids uuid[], weights float8[],
--             min_layer integer, max_layer integer,
--             ceiling float8, floor float8 DEFAULT 0)
--     → TABLE(h9_bin uuid, layer integer, value float8, npoints bigint,
--             density float8, grade float8, geom geometry)
--   Aggregate weighted addresses into a mixed-layer cell set whose cell
--   values respect a population ceiling, by bottom-up digestion: addresses
--   are binned at max_layer; a cell whose accumulated value reaches `floor`
--   emits itself, consuming whole points until the next would push past
--   `ceiling` (always at least one, so a single overweight point still
--   lands); unconsumed excess re-bins to the parent layer and the pass
--   repeats down to min_layer, which emits everything that remains.
--
--   INPUT IS FULL UUIDs (from h9_encode) — bin input is REJECTED, like the
--   k-family: the digest re-bins across layers, which is guaranteed only
--   from the full uuid. Addresses, not coordinates, are the digest's
--   natural input — they pair directly with full-uuid sample tables, and
--   starting from geometry is just array_agg(h9_encode(geom)).
--
--   Dense areas resolve into fine cells, sparse areas aggregate upward, and
--   the sample is captured exactly: sum(value) = sum(weights). `weights`
--   may be NULL (every point weighs 1). Cells of different layers overlap
--   geometrically — a parent holds only what its descendants did not
--   digest; the cell set partitions the SAMPLE, not the surface. Cells are
--   reported max_layer first, UUID-sorted within each layer.
--
--   `h9_bin` is the cell's layer-scoped bin key (renamed from `hex9` in
--   1.4.0). No h9_id identity column: the digest is a terminal aggregate, not
--   a re-binnable source.
--
--   `density` and `grade` are the cartography columns (added 1.4.0), derived
--   from value and layer so callers need not re-derive (and re-invert) them:
--     density = value · 12 · 9^layer / 510065622   -- persons/km², EXACT for
--               the digest (cells are equal-area per layer; 510065622 km² =
--               Earth area, 12·9^layer = cells at the layer).
--     grade   = layer + ln(value)/ln(9)            -- log₉ graduation; +1 ⇒ 9×
--               denser. NULL when value ≤ 0.
--   The numbers are exact arithmetic; their real-world READING carries the
--   source's caveats — see the API README. In particular npoints = 1 is a
--   point-mass reading (areally unsupported), so symbolise/outline those.
--
--   Each row carries its hexagon (`geom`, SRID 4326, corners-only — the
--   same identity render as h9_cell), so the digest is directly mappable:
--   no companion h9_grid is needed, and the displayable max_layer is not
--   bounded by grid-enumeration size.
--
--   Example (≤ 500 people per displayed cell, no cell finer than L12;
--   colour by density, order so parents draw under children):
--     SELECT a.h9_bin, a.layer, a.value, a.density, a.grade, a.geom
--     FROM h9_adaptive(
--              (SELECT array_agg(h9_encode(geom)) FROM dwellings),
--              (SELECT array_agg(occupants::float8) FROM dwellings),
--              6, 12, 500, 100) AS a
--     ORDER BY a.layer ASC;
--
-- Availability: Hex9 1.4.0 (density + grade columns; bin column renamed
-- hex9 → h9_bin). Was 1.2.0 (uuid[] input + geom column since 2026-06-12).
CREATE OR REPLACE FUNCTION h9_adaptive(
        uuid[], float8[],
        integer, integer,
        float8, float8 DEFAULT 0)
    RETURNS TABLE(h9_bin uuid, layer integer, value float8, npoints bigint,
                  density float8, grade float8, geom geometry)
    AS 'MODULE_PATHNAME', 'h9_adaptive'
    LANGUAGE 'c' IMMUTABLE CALLED ON NULL INPUT PARALLEL SAFE
    ROWS 100 COST 1000;

-- h9_adaptive_sphere(uuid[], float8[], min_layer, max_layer, ceiling,
--                    floor DEFAULT 0)
--   Sphere-datum twin of h9_adaptive. The DIGEST is identical — its input
--   and output are addresses, which never cross the datum boundary. What
--   changes is the rendering and the density unit:
--     * geom is the cell polygon in spherical degrees, SRID 0 (see
--       h9_decode_sphere for the datum contract; input uuids should be
--       sphere-minted — the digest cannot check).
--     * density is value per STERADIAN: a layer-L cell's area on the unit
--       sphere is intrinsic (4π/(12·9^L) sr), so no body radius is needed —
--       the sphere datum deliberately carries none. Per-km² on a specific
--       body = density × 4π / body_area_km2, caller-side.
--     * grade is unchanged (datum-free: derived from value and layer only).
--
-- Availability: Hex9 2.1.0
CREATE OR REPLACE FUNCTION h9_adaptive_sphere(
        uuid[], float8[],
        integer, integer,
        float8, float8 DEFAULT 0)
    RETURNS TABLE(h9_bin uuid, layer integer, value float8, npoints bigint,
                  density float8, grade float8, geom geometry)
    AS 'MODULE_PATHNAME', 'h9_adaptive_sphere'
    LANGUAGE 'c' IMMUTABLE CALLED ON NULL INPUT PARALLEL SAFE
    ROWS 100 COST 1000;

-- ── Hamiltonian curve addressing (space-filling curve) ───────────────────────
--
--   The H9 space-filling curve visits every layer-L cell exactly once, in an
--   order that is edge-adjacent between consecutive indices and REFINES: a
--   cell's 9 lineage children occupy curve indices index*9 .. index*9+8.
--   Its two jobs are SORTING (curve order is locality order) and GENERATION
--   (enumerate a region's cells in locality order).
--
--   The packed CURVE-UUID is the sortable form (nibble 0 = 'c' marks it; the
--   rest is the axiom slot + one base-9 rank nibble per layer, f-padded).
--   Native uuid btree order over curve-uuids AT ONE LAYER is exactly curve
--   order, so locality clustering needs no new types or operators:
--
--       ALTER TABLE pings ADD COLUMN curve uuid
--           GENERATED ALWAYS AS (h9_curve(h9_bin(h9_id, 12))) STORED;
--       CREATE INDEX ON pings (curve);
--       CLUSTER pings USING pings_curve_idx;      -- locality-ordered heap
--
--   NOTE in MIXED-layer collections the f-padding makes an ancestor sort
--   AFTER its descendants (post-order); group or filter by h9_curve_layer
--   when mixing layers. Unlike an h9-uuid body, a curve-uuid truncates
--   EXACTLY: h9_curve_bin is pure prefix truncation and yields the LINEAGE
--   ancestor's curve address (the curve tree is iterated one-generation
--   canonical parents — cf. h9_cell_parent, not the deep-fold
--   h9_cell_ancestor).

-- h9_cell_parent(uuid) → uuid
--   Canonical CELL parent of a layer-L bin (mode-0 convention): the single
--   layer-(L-1) cell containing the cell's mode-0 d_cell. This is the
--   curve's generation step — the curve/lineage tree is ITERATED
--   one-generation parents. Distinct from h9_bin, which answers the POINT
--   question (and is only guaranteed from a full uuid). Every parent has
--   exactly 9 canonical children. Errors on an L0 input.
--
-- Availability: Hex9 1.5.0 (ABI since 1.4.0)
CREATE OR REPLACE FUNCTION h9_cell_parent(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_cell_parent'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

-- h9_cell_ancestor(uuid, integer) → uuid
--   Canonical layer-L ancestor of a bin — the DIRECT deep fold (exact at
--   any depth), NOT the iterated one-generation parent: the two relations
--   differ on hexagon-band cells beyond one generation (composing parents
--   re-adjudicates splits every layer). Pass-through when already at L.
--
-- Availability: Hex9 1.5.0 (ABI since 1.4.0)
CREATE OR REPLACE FUNCTION h9_cell_ancestor(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_cell_ancestor'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

-- h9_curve(uuid) → uuid
--   Curve-uuid of a cell. Accepts a canonical bin uuid at any layer (from
--   h9_bin / h9_grid / h9_adaptive), a full uuid from h9_encode (re-binned
--   canonically at its own depth first), or a curve-uuid (returned as-is).
--   Encode is pure address-table arithmetic — exact at any depth.
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_curve_decode(uuid) → uuid
--   Constructive inverse: the canonical bin uuid of a curve address, at the
--   curve address's own layer. An h9-uuid input passes through unchanged.
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_decode(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_decode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 500;

-- h9_is_curve(uuid) → boolean
--   True when the uuid is a packed curve address (nibble 0 = 'c').
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_is_curve(uuid)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'h9_is_curve'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

-- h9_curve_layer(uuid) → integer
--   Layer of a curve-uuid (0..h9_lmax()).
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_layer(uuid)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'h9_curve_layer'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

-- h9_curve_bin(uuid, integer) → uuid
--   Layer-L curve ancestor — PURE PREFIX TRUNCATION of the rank nibbles,
--   exact on the curve/lineage tree (no fold, no geometry). Input must be a
--   curve-uuid at a layer >= L. This is the safe cross-layer coarsening the
--   h9-uuid body famously lacks (docs/addressing-doctrine.md F3): group
--   deep curve keys by h9_curve_bin(curve, L) freely.
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_bin(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_bin'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- h9_curve_index(uuid) → numeric
--   Position on the curve as a number, 0 .. 12*9^L - 1 (numeric: it exceeds
--   bigint above layer 18). Accepts an h9 uuid (encoded first) or a
--   curve-uuid (pure arithmetic). The index is the same numeral the
--   curve-uuid packs, so ORDER BY h9_curve_index == ORDER BY h9_curve at
--   one layer — use the uuid for storage/indexes, the numeric for
--   arithmetic (ranges, striding, partitioning).
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_index(uuid)
    RETURNS numeric
    AS 'MODULE_PATHNAME', 'h9_curve_index'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_curve_pack(numeric, integer) → uuid
--   Curve-uuid from a plain curve index at a known layer (the inverse of
--   h9_curve_index on the representation; a bare index does not
--   self-describe its layer).
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_pack(numeric, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_pack'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- h9_curve_label(uuid) → text
--   Human-readable curve address: 'c' + slot hex char + base-9 rank chars,
--   e.g. 'c112504' (length carries the layer; empty ranks at L0). Accepts
--   an h9 uuid (encoded first) or a curve-uuid.
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_label(uuid)
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9_curve_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_curve_from_label(text) → uuid
--   Parse a curve label back to a curve-uuid.
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_from_label(text)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_from_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

-- h9_curve_cells(uuid, integer) → TABLE(h9_bin uuid, h9_curve uuid)
--   The GENERATION primitive: every layer-L LINEAGE descendant of the
--   cell (iterated one-generation canonical children — the curve's own
--   tree), emitted IN CURVE ORDER (ascending curve index), as (bin key,
--   curve key) pairs — 9^(L - cell_layer) rows, capped by
--   hex9.grid_max_cells. Input is an h9 bin/full uuid or a curve-uuid.
--   Lineage is transitive and prefix-exact on curve-uuids, but its
--   descendant sets displace 1/6 of the ancestor's area — for the
--   geometrically-bounded aggregation partition use h9_owned_cells.
--
--   Example — synthesize an L6 dataset over Edinburgh's L2 cell, written
--   in locality order (neighbouring rows are neighbouring cells):
--     SELECT c.h9_bin, c.h9_curve, h9_cell(c.h9_bin, 6) AS geom
--     FROM h9_curve_cells(h9_bin(h9_encode(
--              ST_SetSRID(ST_MakePoint(-3.19, 55.95), 4326)), 2), 6) AS c;
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_curve_cells(uuid, integer)
    RETURNS TABLE(h9_bin uuid, h9_curve uuid)
    AS 'MODULE_PATHNAME', 'h9_curve_cells'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 729 COST 5000;

-- h9_cell_children(uuid) → SETOF uuid
--   The 9 canonical children of a cell, in CURVE-RANK order (one
--   generation — the lineage and ownership readings coincide here).
--   Input is an h9 bin/full uuid or a curve-uuid; errors at h9_lmax().
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_cell_children(uuid)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_cell_children'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 9 COST 500;

-- h9_owned_cells(uuid, integer) → TABLE(h9_bin uuid, h9_curve uuid)
--   The AGGREGATION primitive: every layer-L OWNED sub-zone of the cell —
--   the cells whose h9_cell_ancestor at the cell's layer IS the cell.
--   Exactly 9^(L - cell_layer) rows, and the owned sets over all zones at
--   one layer PARTITION the globe: aggregation by owned sub-zone is
--   loss-free and double-count-free (the "owned sub-zone" relation of
--   OGC API-DGGS issue #108 — sub-zones give coverage, identifiers give
--   lineage, aggregation needs ownership). Ownership is geometrically
--   bounded (only rim splits protrude, by their far half) but NOT
--   transitive; the lineage set (h9_curve_cells) is the transitive one.
--   Rows are CURVE-SORTED; the h9 bin is the primary key form (h9→curve
--   is the cheap direction). Capped by hex9.grid_max_cells.
--
--   Example — loss-free roll-up of an L6 population column to L4 zones:
--     SELECT o.h9_bin AS zone, sum(p.pop)
--     FROM zones z, h9_owned_cells(z.h9_bin, 6) o
--     JOIN pop_l6 p ON p.h9_bin = o.h9_bin
--     GROUP BY 1;
--
-- Availability: Hex9 1.5.0
CREATE OR REPLACE FUNCTION h9_owned_cells(uuid, integer)
    RETURNS TABLE(h9_bin uuid, h9_curve uuid)
    AS 'MODULE_PATHNAME', 'h9_owned_cells'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 729 COST 5000;

-- ═══ E4H: the aperture-4 structural tail ═════════════════════════════════════
--
-- An E4H address extends a hex9 HOST bin at attach layer L with an 0xE break
-- marker, one HALF digit (which state-cut trapezoid of the host hexagon) and
-- up to (h9_lmax()-2-L) tail digits from {0,1..5} descending the aperture-4
-- half-hex rep-4 carrier: exact nesting, straight edges, and suffix-local
-- truncation (h9e_bin IS binning — unlike the a9 body). E4H uuids are
-- ADDRESSES under the universality contract (libhex9 docs/universality.md):
-- the classifier is exact — the same point mints the bit-identical uuid on
-- every platform at every depth. The 0xE marker cannot occur in h9 or curve
-- uuids, so h9_is_e4h is decisive, and the h9_*/curve machinery REJECTS E4H
-- input (tails have their own operations — the functions below).
--
--   label form: <h9-label>E<half><digits>       e.g. 0031586E1213

-- h9e_encode(geometry, layer, depth) → uuid
--   Mint the E4H address of a POINT: host bin at `layer`, then `depth`
--   aperture-4 digits below the half cut. layer + 2 + depth <= h9_lmax().
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_encode(geometry, integer, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9e_encode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

-- Sphere-datum twin (already-spherical lon/lat degrees; datum doctrine as
-- h9_encode_sphere). Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_encode_sphere(geometry, integer, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9e_encode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

-- h9e_decode(uuid) → geometry
--   Representative POINT of the leaf trapezoid (its centroid), SRID 4326.
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_decode(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9e_decode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

CREATE OR REPLACE FUNCTION h9e_decode_sphere(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9e_decode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

-- h9e_partner(uuid) → geometry
--   Representative POINT of the leaf's PARTNER half — the mirror trapezoid
--   across the leaf's long side, i.e. the other half of the same fine
--   hexagon. Matched pairs are global: h9e_encode of this point yields an
--   address sharing the final digit (for centre children the partner is the
--   host's other half). Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_partner(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9e_partner'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

CREATE OR REPLACE FUNCTION h9e_partner_sphere(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9e_partner_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

-- h9e_bin(uuid, depth) → uuid
--   Truncate the tail to `depth` digits — suffix-local, exact (this is the
--   one place truncation IS binning; deepening is rejected).
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_bin(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9e_bin'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 10;

-- h9e_depth(uuid) → integer — tail digit count; -1 if not E4H.
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_depth(uuid)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'h9e_depth'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 1;

-- h9_is_e4h(uuid) → boolean — decisive marker test (any nibble = 0xE).
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_is_e4h(uuid)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'h9_is_e4h'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 1;

-- h9e_host(uuid) → uuid — the host bin (the h9 part above the 0xE marker).
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_host(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9e_host'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 10;

-- h9e_label(uuid) → text / h9e_parse_label(text) → uuid
--   '<h9-label>E<half><digits>'. Bare host labels carry the F1 split-hex
--   ambiguity caveat of h9_parse_label. Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9e_label(uuid)
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9e_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 10;

CREATE OR REPLACE FUNCTION h9e_parse_label(text)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9e_parse_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 100;

-- ═══ Grid verbs: aim / walk_to / vision_cone ═════════════════════════════════
--
-- Cell-first field primitives over the neighbour algebra — the SQL twins of
-- the python wheel's hex9.verbs (hhg9 fox-and-rabbits PoC). Every verb takes
-- and returns layer-scoped BIN KEYS (h9_bin / h9_grid keys) — deliberately
-- unlike the k-family, whose inputs must be full uuids: a verb's output is
-- the next verb's input, and the guaranteed full-uuid re-derivation happens
-- inside. Bearings are degrees clockwise from north, advisory FP over cell
-- centroids: they steer, they never mint. The octant seam changes nothing.

-- h9_aim(key, layer, bearing) → uuid
--   The neighbour key whose centroid bearing best matches.
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_aim(uuid, integer, double precision)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_aim'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 500;

-- h9_walk_to(src, dest, layer, obstacles, max_expand) → uuid[]
--   Shortest hex path src → dest (A* over the neighbour algebra; obstacle
--   keys impassable, dest passable even if listed). Returns the key path
--   (first = src, last = dest), or NULL when no path exists within
--   max_expand node expansions. Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_walk_to(uuid, uuid, integer,
                                      uuid[] DEFAULT '{}',
                                      integer DEFAULT 5000)
    RETURNS uuid[]
    AS 'MODULE_PATHNAME', 'h9_walk_to'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 5000;

-- h9_vision_cone(key, layer, bearing, half_angle, k, obstacles) → SETOF uuid
--   Keys within k rings, within half_angle of bearing, not occluded by an
--   obstacle key (great-circle sight lines sampled at a third of the local
--   hex pitch; an obstacle never occludes itself — you can see the hedge,
--   just not through it). The source key is always included; rows are
--   UUID-sorted. Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_vision_cone(uuid, integer, double precision,
                                          double precision, integer,
                                          uuid[] DEFAULT '{}')
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_vision_cone'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 37 COST 5000;

-- h9_bearing_from_true(bearing, lat) → double precision
--   TRUE (geodesic, WGS84) azimuth → the verbs' bearing convention, at the
--   geodetic latitude where the heading holds (the source cell — compose
--   ST_Y(h9_decode(key))). The verbs compute bearings by spherical trig
--   over geodetic coordinates; a measured true-north heading differs by
--   the flattening term tan(b_verb) = (M/N)·tan(b_true): identity at the
--   poles and cardinal bearings, largest ~0.19° at equatorial diagonals —
--   apply it when that is inside your cone tolerance. Magnetic headings
--   must already be reduced to true north (declination is epoch/model-
--   dependent — dataset metadata). Meaningless on the sphere datum.
--
--   Example — a drone heading feeding a vision cone:
--     SELECT v FROM h9_vision_cone(key, 8,
--         h9_bearing_from_true(hdg, ST_Y(h9_decode(key))), 10.0, 6) AS v;
--
-- Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_bearing_from_true(double precision, double precision)
    RETURNS double precision
    AS 'MODULE_PATHNAME', 'h9_bearing_from_true'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 1;

-- Inverse: verb bearing → true azimuth. Availability: Hex9 2.3.0
CREATE OR REPLACE FUNCTION h9_bearing_to_true(double precision, double precision)
    RETURNS double precision
    AS 'MODULE_PATHNAME', 'h9_bearing_to_true'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE COST 1;
