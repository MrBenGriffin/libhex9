-- postgis_hex9 upgrade: 2.1.0 → 2.3.0 (2.2.x had no SQL surface changes).
-- Adds the E4H aperture-4 tail functions and the grid verbs — additive
-- only; no existing function, address, or bin value changes.

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
