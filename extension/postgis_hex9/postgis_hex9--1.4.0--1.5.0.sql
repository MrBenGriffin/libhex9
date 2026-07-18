-- postgis_hex9 1.4.0 → 1.5.0 upgrade
--
-- Hamiltonian curve addressing (the H9 space-filling curve): sortable
-- curve-uuids, numeric indices, labels, exact prefix coarsening, and the
-- curve-ordered descendant generator. All functions are NEW — no existing
-- signature changes, no persisted values change.
--
--   * h9_cell_parent(uuid) → uuid          canonical cell parent (SQL surface
--   * h9_cell_ancestor(uuid, int) → uuid     for the 1.4.0 ABI pair; the
--                                            curve's lineage step)
--   * h9_curve(uuid) → uuid                curve-uuid (the sortable key)
--   * h9_curve_decode(uuid) → uuid         constructive inverse (bin uuid)
--   * h9_is_curve(uuid) → boolean          type test (nibble 0 = 'c')
--   * h9_curve_layer(uuid) → integer       layer of a curve-uuid
--   * h9_curve_bin(uuid, int) → uuid       EXACT prefix coarsening
--   * h9_curve_index(uuid) → numeric       position 0..12*9^L-1
--   * h9_curve_pack(numeric, int) → uuid   index + layer → curve-uuid
--   * h9_curve_label(uuid) → text          'c<slot><ranks>'
--   * h9_curve_from_label(text) → uuid
--   * h9_curve_cells(uuid, int)
--         → TABLE(h9_bin uuid, h9_curve uuid)
--     LINEAGE descendants in curve order (the generation primitive)
--   * h9_cell_children(uuid) → SETOF uuid   9 canonical children,
--                                            curve-rank order
--   * h9_owned_cells(uuid, int)
--         → TABLE(h9_bin uuid, h9_curve uuid)
--     OWNED sub-zones, curve-sorted (the aggregation partition —
--     exactly 9^d per zone; OGC API-DGGS issue #108)
--
-- Doctrine notes (hex9_c.h, core/h9_curve.h): the curve tree is the
-- LINEAGE tree (iterated one-generation canonical parents), curve order at
-- one layer == native uuid btree order over curve-uuids, and in mixed-layer
-- collections an ancestor sorts AFTER its descendants (post-order).
--
-- Run via:
--   ALTER EXTENSION postgis_hex9 UPDATE TO '1.5.0';

\echo Use "ALTER EXTENSION postgis_hex9 UPDATE TO '1.5.0'" to load this file. \quit

CREATE OR REPLACE FUNCTION h9_cell_parent(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_cell_parent'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

CREATE OR REPLACE FUNCTION h9_cell_ancestor(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_cell_ancestor'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

CREATE OR REPLACE FUNCTION h9_curve(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

CREATE OR REPLACE FUNCTION h9_curve_decode(uuid)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_decode'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 500;

CREATE OR REPLACE FUNCTION h9_is_curve(uuid)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'h9_is_curve'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

CREATE OR REPLACE FUNCTION h9_curve_layer(uuid)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'h9_curve_layer'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 1;

CREATE OR REPLACE FUNCTION h9_curve_bin(uuid, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_bin'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

CREATE OR REPLACE FUNCTION h9_curve_index(uuid)
    RETURNS numeric
    AS 'MODULE_PATHNAME', 'h9_curve_index'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

CREATE OR REPLACE FUNCTION h9_curve_pack(numeric, integer)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_pack'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

CREATE OR REPLACE FUNCTION h9_curve_label(uuid)
    RETURNS text
    AS 'MODULE_PATHNAME', 'h9_curve_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

CREATE OR REPLACE FUNCTION h9_curve_from_label(text)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_curve_from_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 5;

CREATE OR REPLACE FUNCTION h9_curve_cells(uuid, integer)
    RETURNS TABLE(h9_bin uuid, h9_curve uuid)
    AS 'MODULE_PATHNAME', 'h9_curve_cells'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 729 COST 5000;

CREATE OR REPLACE FUNCTION h9_cell_children(uuid)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_cell_children'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 9 COST 500;

CREATE OR REPLACE FUNCTION h9_owned_cells(uuid, integer)
    RETURNS TABLE(h9_bin uuid, h9_curve uuid)
    AS 'MODULE_PATHNAME', 'h9_owned_cells'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 729 COST 5000;
