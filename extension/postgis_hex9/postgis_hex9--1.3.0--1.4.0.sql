-- postgis_hex9 1.3.0 → 1.4.0 upgrade
--
-- Identity columns + cartography columns + UUID-column naming standardisation.
-- Both set-returning functions change their result rowtype, so each must be
-- DROPped and recreated (CREATE OR REPLACE cannot change a function's output
-- columns).
--
--   * h9_grid → (h9_id, h9_bin, geom, centroid)
--       - h9_id  ADDED: the cell's full reversible identity (layer-29 UUID).
--                h9_bin(h9_id, L') is exact at every L' ≤ layer — the safe key
--                to persist and re-bin (a bin's own coarsening is the F3
--                fossil). See docs/addressing-doctrine.md F4.
--       - hex9 RENAMED → h9_bin (it always was the layer-scoped bin key).
--
--   * h9_adaptive → (h9_bin, layer, value, npoints, density, grade, geom)
--       - hex9 RENAMED → h9_bin (no h9_id: the digest is a terminal aggregate).
--       - density, grade ADDED (cartography; derived from value and layer).
--
-- Run via:
--   ALTER EXTENSION postgis_hex9 UPDATE TO '1.4.0';

\echo Use "ALTER EXTENSION postgis_hex9 UPDATE TO '1.4.0'" to load this file. \quit

DROP FUNCTION IF EXISTS h9_grid(geometry, integer, integer);
CREATE OR REPLACE FUNCTION h9_grid(geometry, integer, integer DEFAULT 0)
    RETURNS TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)
    AS 'MODULE_PATHNAME', 'h9_grid'
    LANGUAGE 'c' STABLE STRICT PARALLEL SAFE
    ROWS 1000 COST 5000;

DROP FUNCTION IF EXISTS h9_adaptive(uuid[], float8[], integer, integer, float8, float8);
CREATE OR REPLACE FUNCTION h9_adaptive(
        uuid[], float8[],
        integer, integer,
        float8, float8 DEFAULT 0)
    RETURNS TABLE(h9_bin uuid, layer integer, value float8, npoints bigint,
                  density float8, grade float8, geom geometry)
    AS 'MODULE_PATHNAME', 'h9_adaptive'
    LANGUAGE 'c' IMMUTABLE CALLED ON NULL INPUT PARALLEL SAFE
    ROWS 100 COST 1000;

-- h9_common_ancestor: rename OUT column hex9 → h9_bin (it returns a bin).
DROP FUNCTION IF EXISTS h9_common_ancestor(uuid[], integer);
CREATE OR REPLACE FUNCTION h9_common_ancestor(
        uuid[], integer,
        OUT label text, OUT h9_bin uuid, OUT layer integer)
    RETURNS record
    AS 'MODULE_PATHNAME', 'h9_common_ancestor'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 100;
