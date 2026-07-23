-- postgis_hex9 2.0.0 → 2.1.0 upgrade
--
-- ADDITIVE ONLY. No existing function changes behaviour; no stored address
-- moves. Safe to apply in place.
--
-- WHAT IS NEW: the sphere-datum twins. Since 2.0.0 every address is minted
-- on the unit sphere; WGS84 enters the chain at exactly one place — the
-- authalic latitude reduction at the lon/lat boundary. The *_sphere twins
-- run the IDENTICAL chain minus that reduction: their lon/lat are
-- already-spherical degrees, for datasets that own their own datum
-- (another body's authalic frame, celestial RA/dec).
--
-- THE DATUM CONTRACT (docs/warp-regimes.md, "Two datums, one regime"):
--   * The datum is part of the FUNCTION's identity — that is why these are
--     distinct functions and not a setting. Nothing address-affecting is a
--     GUC, ever (the pre-2.0.0 hex9.use_warp GUC is the cautionary tale).
--   * The 16 bytes do not record the datum. A WGS84-minted and a
--     sphere-minted address for the same numeric lon/lat differ below
--     roughly layer 5, with no way to tell them apart afterwards. The datum
--     is DATASET metadata, owned by you. Never mix datums in one dataset.
--   * Sphere geometries are emitted with SRID 0 (unknown): spherical
--     degrees have no EPSG identity here. ST_SetSRID if your dataset has
--     one. Never tag them 4326.
--
-- No twin exists where composition already serves: label centroids are
-- h9_decode_sphere(h9_parse_label(t)).

-- h9_encode_sphere(geometry) → uuid  (see h9_encode; input POINT is
-- spherical lon/lat degrees)
CREATE OR REPLACE FUNCTION h9_encode_sphere(geometry)
    RETURNS uuid
    AS 'MODULE_PATHNAME', 'h9_encode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_encode_many_sphere(geometry[]) → uuid[]  (batch form)
CREATE OR REPLACE FUNCTION h9_encode_many_sphere(geometry[])
    RETURNS uuid[]
    AS 'MODULE_PATHNAME', 'h9_encode_many_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;

-- h9_decode_sphere(uuid) → geometry  (POINT, spherical degrees, SRID 0)
CREATE OR REPLACE FUNCTION h9_decode_sphere(uuid)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_decode_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 10;

-- h9_cell_sphere(uuid, layer, densify DEFAULT 0) → geometry  (POLYGON, SRID 0)
CREATE OR REPLACE FUNCTION h9_cell_sphere(uuid, integer, integer DEFAULT 0)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'h9_cell_sphere'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 200;

-- h9_grid_sphere(geometry, layer, densify DEFAULT 0)  (bounds and all
-- emitted geometry are spherical degrees, SRID 0)
CREATE OR REPLACE FUNCTION h9_grid_sphere(geometry, integer, integer DEFAULT 0)
    RETURNS TABLE(h9_id uuid, h9_bin uuid, geom geometry, centroid geometry)
    AS 'MODULE_PATHNAME', 'h9_grid_sphere'
    LANGUAGE 'c' STABLE STRICT PARALLEL SAFE
    ROWS 1000 COST 5000;

-- h9_adaptive_sphere(uuid[], float8[], min_layer, max_layer, ceiling,
-- floor DEFAULT 0) — same digest as h9_adaptive (addresses never cross the
-- datum boundary); geom in spherical degrees SRID 0; density is value per
-- STERADIAN (a layer-L cell is intrinsically 4π/(12·9^L) sr on the unit
-- sphere, so no body radius is needed — per-km² on a body is
-- density × 4π / body_area, caller-side); grade unchanged (datum-free).
CREATE OR REPLACE FUNCTION h9_adaptive_sphere(
        uuid[], float8[],
        integer, integer,
        float8, float8 DEFAULT 0)
    RETURNS TABLE(h9_bin uuid, layer integer, value float8, npoints bigint,
                  density float8, grade float8, geom geometry)
    AS 'MODULE_PATHNAME', 'h9_adaptive_sphere'
    LANGUAGE 'c' IMMUTABLE CALLED ON NULL INPUT PARALLEL SAFE
    ROWS 100 COST 1000;
