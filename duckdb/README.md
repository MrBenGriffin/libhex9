# hex9 DuckDB extension (`h9`)

A DuckDB extension exposing the Hex9 DGGS over the libhex9 C ABI (`hex9_c.h`):
vectorized scalar functions for encode/decode/bin, cell geometry (WKB),
labels, lineage, adjacency, and the Hamiltonian curve. Built with the
standard DuckDB C++ extension-template flow so it lifts unchanged into a
community-extensions wrapper repo.

## Layout

This directory is a self-contained extension-template tree that compiles the
three libhex9 core sources directly (see `CMakeLists.txt`; `HEX9_ROOT` is the
only lift-out knob). The `duckdb/` and `extension-ci-tools/` subdirectories
are shallow clones (not submodules, gitignored):

```sh
git clone --depth 1 --branch v1.5.4 https://github.com/duckdb/duckdb.git duckdb
git clone --depth 1 --branch v1.5.4 https://github.com/duckdb/extension-ci-tools.git extension-ci-tools
```

## Build & test

```sh
GEN=ninja make release   # builds duckdb + the h9 extension
make test                # runs test/sql/*.test (sqllogictest)
./build/release/duckdb   # shell with h9 statically loaded
```

The loadable artifact is `build/release/extension/h9/h9.duckdb_extension`.

## Types

| hex9 concept | DuckDB type | Note |
|---|---|---|
| cell / bin / curve id | `UUID` | sorts bytewise, same as Postgres `uuid` |
| id as integer | `UHUGEINT` via `h9_id_int` / `h9_id_from_int` | unsigned so ordering matches UUID order |
| curve index | `HUGEINT` | exceeds `BIGINT` above L18; fits int128 |
| geometry | `BLOB` (WKB, lon/lat WGS84) | `ST_GeomFromWKB(...)` with the `spatial` extension |
| decoded point | `STRUCT(lon DOUBLE, lat DOUBLE)` | `h9_decode(u).lon` |

Load policy: `LOAD h9` fails hard if the warp field cannot initialise
(a degraded field would not produce Hex9 addresses), and if the compiled
`HEX9_VERSION` differs from the linked core's `hex9_version()`.

Doctrine: bins are layer-scoped keys, not addresses — adjacency functions
(`h9_neighbors`, `h9_kring`, `h9_kdisk`) reject bin input; re-derive from the
full UUID (`docs/addressing-doctrine.md`).

**Leading zeros are load-bearing** in hex9 digit strings: `'00101'` and
`'101'` name different cells. Never cast labels (or uuid text) to a numeric
type. The integer forms this extension provides are the only value-safe ones:
`h9_id_int` is a full-width 128-bit reinterpretation (reconstruction always
writes all 16 bytes back), and `h9_curve_index` is an ordinal whose digit
width is re-derived from the separate `layer` argument in `h9_curve_pack`.

## Lift-out (community submission)

1. New repo from the extension-template skeleton; add libhex9 as a submodule.
2. Copy `src/`, `test/`, `CMakeLists.txt`, `extension_config.cmake`, `Makefile`.
3. Point `HEX9_ROOT` at the submodule path.
4. Submit a descriptor to duckdb/community-extensions; its CI builds all
   platforms including wasm.

wasm note: the warp blob is embedded via `__asm__ .incbin`
(`core/h9_warp_embedded.cpp`). If Emscripten rejects that, substitute a
generated C byte-array TU exposing the same `h9_warp_blob`/`h9_warp_blob_end`
symbols (not yet needed — verify with `make wasm_mvp` first).
