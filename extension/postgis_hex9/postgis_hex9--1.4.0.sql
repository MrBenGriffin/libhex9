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

-- ── Encode ───────────────────────────────────────────────────────────────────

-- h9_encode(geometry) → uuid
--   Encodes a POINT to a self-contained UUID at layer 29 resolution.
--   The UUID carries all decode context internally; no companion byte needed.
--   Suitable for exact round-trip via h9_decode(), indexing, and h9_bin().
--
-- Availability: Hex9 1.0.0
CREATE OR REPLACE FUNCTION h9_encode(geometry)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_encode'
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

-- ── Decode ───────────────────────────────────────────────────────────────────

-- h9_decode(uuid) → geometry
--   Decode a UUID to a POINT geometry (SRID 4326).
--   Exact round-trip: h9_decode(h9_encode(pt)) recovers the input point
--   to within the resolution of layer 29 (~95 nm cell diameter).
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

-- ── Hierarchy ────────────────────────────────────────────────────────────────

-- h9_bin(uuid, integer) → uuid
--   Returns the bin-key UUID at the given layer (0..29).
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
