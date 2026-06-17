-- postgis_hex9 1.2.0 → 1.3.0 upgrade
--
-- Two surface changes:
--   * h9_diag(geometry) → text is REMOVED. It was an internal frame-parity
--     diagnostic ("temporary aid"); encoder parity has long been confirmed,
--     so it no longer ships. (The core C ABI keeps hex9_diag for tooling;
--     only the SQL entry point is dropped.)
--   * h9_encode_many(geometry[]) → uuid[] is ADDED — the batch form of
--     h9_encode, one OpenMP-parallel pass over an array of POINTs, returning
--     the UUIDs in input order. Composes with h9_adaptive(uuid[], …).
--
-- Run via:
--   ALTER EXTENSION postgis_hex9 UPDATE TO '1.3.0';

\echo Use "ALTER EXTENSION postgis_hex9 UPDATE TO '1.3.0'" to load this file. \quit

-- Drop the retired diagnostic (IF EXISTS: a DB may already lack it).
DROP FUNCTION IF EXISTS h9_diag(geometry);

-- Batch encode (see hex9.sql.in for full documentation).
CREATE OR REPLACE FUNCTION h9_encode_many(geometry[])
    RETURNS uuid[]
    AS 'MODULE_PATHNAME', 'h9_encode_many'
    LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
    COST 50;
