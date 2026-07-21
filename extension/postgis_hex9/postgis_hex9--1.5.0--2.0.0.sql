-- postgis_hex9 1.5.0 → 2.0.0 upgrade
--
-- ============================================================================
--  THIS UPGRADE CHANGES WHAT EVERY ADDRESSING FUNCTION RETURNS.
--  Addresses already stored in your tables are NOT updated and are NOT
--  comparable with addresses produced after it.
-- ============================================================================
--
-- No function signature changes. No new objects. Nothing in this script
-- touches your data — and that is deliberate, not an oversight. Read on.
--
-- WHAT CHANGED
--
-- libhex9 2.0.0 replaced the addressing chain. Up to 1.5.0 the core ran on
-- the WGS84 ellipsoid against a WGS84-trained warp field. It now reduces
-- geodetic latitude to the AUTHALIC latitude (Karney series), runs on the
-- unit sphere, and applies the Sphere-L6 wedge-fold field. One trained field
-- serves every ellipsoid; the ellipsoid-specific part of the chain is the
-- latitude series and nothing else.
--
-- The old regime has been removed outright, not deprecated. It is not
-- reachable by a GUC, a build flag, or a function argument.
--
-- WHAT THAT MEANS FOR STORED ADDRESSES
--
-- The same point now generally addresses to a DIFFERENT cell from around
-- layer 7 downward. It is point-dependent: coarse bins are often unchanged,
-- some points are identical at every layer, and some diverge as shallow as
-- layer 7. For example, Westminster's layer-8 bin moves:
--
--     1.5.0:  435878503
--     2.0.0:  435878530
--
-- Nothing in a 16-byte address records which regime produced it. The two
-- address sets are indistinguishable by inspection.
--
-- WHAT TO DO — AND WHAT NOT TO
--
-- DO re-derive addresses from the SOURCE GEOMETRY you encoded from:
--
--     UPDATE my_table SET h9 = h9_encode(geom);
--
-- DO NOT decode-and-re-encode existing addresses to "migrate" them. It is
-- technically possible and it is a trap, for two reasons:
--
--   1. It is not idempotent and it cannot detect itself. Running it twice
--      silently displaces every cell again — the second pass decodes an
--      already-converted address as though it were a 1.5.0 one, which
--      "succeeds" because every 16-byte pattern is a valid address. There is
--      no way afterwards to tell converted data from doubly-converted data.
--
--   2. It quietly invalidates everything derived from those addresses.
--      A population density, an aggregate, a published figure or a map
--      computed against 1.5.0 addresses is no longer reproducible from
--      migrated data, and nothing announces that.
--
-- Treat 1.5.0 addresses as historical values belonging to a prior regime, the
-- way you would treat measurements taken with a since-recalibrated
-- instrument. Keep them, or re-derive from source. Do not rewrite them.
--
-- IF THE SOURCE GEOMETRY IS GONE
--
-- The addresses themselves are then the only record of position, and only the
-- OLD library can interpret them. Recover coordinates BEFORE upgrading:
--
--     -- while still on 1.5.0:
--     ALTER TABLE my_table ADD COLUMN geom_recovered geometry(Point, 4326);
--     UPDATE my_table SET geom_recovered = h9_centroid(h9);
--
-- then upgrade and re-encode from geom_recovered. Note the precision ceiling:
-- a centroid is the centre of its cell, not the original point, so this is
-- lossless only at that cell's layer and discards any finer precision the
-- original point carried. Prefer real source geometry whenever it exists.
--
-- ALSO REMOVED IN 2.0.0
--
--   * The `hex9.use_warp` GUC. It could change the output of functions
--     declared IMMUTABLE, which PostgreSQL may cache in functional indexes
--     and generated columns — so two sessions could disagree about what a
--     stored index entry meant. It was only ever intended for development
--     A/B testing; that capability still exists in libhex9's C ABI, where it
--     cannot reach SQL.
--
-- Run via:
--   ALTER EXTENSION postgis_hex9 UPDATE TO '2.0.0';

\echo Use "ALTER EXTENSION postgis_hex9 UPDATE TO '2.0.0'" to load this file. \quit

DO $$
BEGIN
    RAISE WARNING 'postgis_hex9 2.0.0: the addressing chain changed. Hex9 '
                  'addresses stored under 1.5.0 are NOT comparable with '
                  'addresses produced from now on (they diverge from about '
                  'layer 7).';
    RAISE WARNING 'postgis_hex9 2.0.0: re-derive affected columns from their '
                  'SOURCE geometry, e.g. UPDATE t SET h9 = h9_encode(geom). '
                  'Do NOT decode-and-re-encode stored addresses: that is not '
                  'idempotent, cannot detect a second run, and invalidates '
                  'anything already derived from them.';
END
$$;
