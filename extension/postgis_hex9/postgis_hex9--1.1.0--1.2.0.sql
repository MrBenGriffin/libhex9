-- postgis_hex9 1.1.0 → 1.2.0 upgrade
--
-- 1.2.0 is additive — no existing signature changes. New functions:
--   * h9_neighbors(uuid, integer)            → SETOF uuid
--   * h9_kring(uuid, integer, integer)       → SETOF uuid
--   * h9_kdisk(uuid, integer, integer)       → SETOF uuid
--   * h9_parse_label(text)                   → uuid
--   * h9_label_centroid(text)                → geometry
--   * h9_common_ancestor(uuid[], integer)    → (label text, hex9 uuid, layer integer)
--   * h9_adaptive(geometry[], float8[], …)   → TABLE(hex9 uuid, layer int, value float8, npoints bigint)
--
-- Behavioural fix carried by the module (no SQL change): h9_bin on a bin
-- UUID is now exact — re-binning at the same layer is the identity and
-- coarsening a bin agrees with coarsening the full UUID (previously
-- undefined behaviour).
--
-- Run via:
--   ALTER EXTENSION postgis_hex9 UPDATE TO '1.2.0';

\echo Use "ALTER EXTENSION postgis_hex9 UPDATE TO '1.2.0'" to load this file. \quit

-- Adjacency (see hex9.sql.in for full documentation).

CREATE OR REPLACE FUNCTION h9_neighbors(uuid, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_neighbors'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 6 COST 50;

CREATE OR REPLACE FUNCTION h9_kring(uuid, integer, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_kring'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 12 COST 100;

CREATE OR REPLACE FUNCTION h9_kdisk(uuid, integer, integer)
    RETURNS SETOF uuid
    AS 'MODULE_PATHNAME', 'h9_kdisk'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    ROWS 19 COST 100;

-- Label parsing / mesh prefix ops.

CREATE OR REPLACE FUNCTION h9_parse_label(text)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_parse_label'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

CREATE OR REPLACE FUNCTION h9_label_centroid(text)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_label_centroid'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

CREATE OR REPLACE FUNCTION h9_common_ancestor(
        uuid[], integer,
        OUT label text, OUT hex9 uuid, OUT layer integer)
    RETURNS record
    AS 'MODULE_PATHNAME', 'h9_common_ancestor'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 100;

-- Adaptive multi-layer grid (population digest).

CREATE OR REPLACE FUNCTION h9_adaptive(
        geometry[], float8[],
        integer, integer,
        float8, float8 DEFAULT 0)
    RETURNS TABLE(hex9 uuid, layer integer, value float8, npoints bigint)
    AS 'MODULE_PATHNAME', 'h9_adaptive'
    LANGUAGE 'c' IMMUTABLE CALLED ON NULL INPUT PARALLEL SAFE
    ROWS 100 COST 1000;
