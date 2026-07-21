# Installing libhex9

libhex9 is a geometry-free DGGS core with a stable C ABI (`hex9_c.h`). It is
consumed five ways — pick the section that matches your purpose. The C library
is the foundation; the PostGIS extension, the PDAL plugins, and (optionally)
the Python module build *on top of* it.

> **Addressing layout.** All builds default to the **L30** layout (a full UUID
> is a max-depth bin). The legacy 29-layer format is opt-in with
> `-DHEX9_USE_L29=ON` — only needed to reproduce `test/golden_l29.py`.

---

## 1. The C library (`libhex9`) — the foundation

Builds the shared + static libraries and installs the public header.

```sh
cmake -S . -B build
cmake --build build -j
cmake --install build            # add: sudo, if installing to a system prefix
```

Installs into `CMAKE_INSTALL_PREFIX` (default `/usr/local`):

| Artifact | Destination | Used by |
|---|---|---|
| `libhex9.{so,dylib}` | `lib/` | the PostGIS extension (dynamic) |
| `libhex9.a` | `lib/` | static linkers (e.g. the Rust crate) |
| `hex9_c.h` | `include/` | every consumer |

Choose a non-system prefix with `-DCMAKE_INSTALL_PREFIX=/path` (no sudo needed).
The shared library links OpenMP when present (parallel batch loops); the static
archive deliberately omits it (self-contained, serial-correct).

Run the test suite with `ctest --test-dir build`.

---

## 2. The PostGIS extension (`postgis_hex9`)

A PGXS-built PostgreSQL/PostGIS extension that **links the installed
`libhex9`** — so install the C library (section 1) first.

**Prerequisites**

- `libhex9` installed (section 1) — note its prefix.
- The **`pg_config` of the server you actually connect to** (see the pitfall
  below) — a PostgreSQL dev install.
- A PostGIS **source** tree matching your loaded PostGIS (liblwgeom headers are
  not shipped by packagers) — point `POSTGIS_SRC` at it. Find your loaded
  version inside psql with `SELECT postgis_full_version();`, then unpack that
  exact version and run `./configure` in it (this *generates* `liblwgeom.h`
  from `liblwgeom.h.in` — an unconfigured tarball has only the `.in`). Match the
  full version including patch level where you can; across a patch release
  (e.g. 3.6.2 vs 3.6.3) the headers are usually ABI-compatible, but exact is
  safest.
  - **`POSTGIS_SRC` must be an absolute path — no `~`.** It is passed straight
    into a compiler `-I` flag, and `-I` does not perform tilde expansion, so a
    `~` yields `fatal error: 'liblwgeom.h' file not found`. Use `$HOME/...`.

```sh
cd extension/postgis_hex9
make PG_CONFIG=/path/to/pg_config HEX9_PREFIX=/usr/local POSTGIS_SRC=/path/to/postgis-3.x.y
sudo make install
```

> **Pitfall — which `pg_config`?** The Makefile defaults to `pg_config` on
> `PATH` (often Homebrew's `/usr/local/bin/pg_config`). If your running server
> is a *different* install — e.g. **Postgres.app** — `make install` will drop
> the module and SQL files into the wrong tree and the server won't see them.
> Two installs can even share a major version (Postgres.app 18.4 vs Homebrew
> 18.3) yet be entirely separate.
>
> Confirm they match before building:
>
> ```sql
> -- in psql, connected to your target server:
> SELECT version();        -- e.g. PostgreSQL 18.4 (Postgres.app)
> SHOW config_file;        -- ~/Library/Application Support/Postgres/var-18/... => Postgres.app
> ```
> ```sh
> pg_config --version      # must match SELECT version() above
> ```
>
> If they differ, pass the right one explicitly, e.g. for Postgres.app:
> `PG_CONFIG=/Applications/Postgres.app/Contents/Versions/18/bin/pg_config`.
>
> Since 2.0.0 you do not have to remember: `make installcheck` compares the
> `pg_config` you built with against the running server's own `PKGLIBDIR` and
> refuses to proceed on a mismatch, printing the correct `PG_CONFIG=` line.
> `POSTGIS_SRC` is likewise checked against the server's loaded PostGIS, and
> discovered automatically where possible.

**Tip — pin it once.** Both variables are machine-specific and easy to get
silently wrong. Put them in `extension/postgis_hex9/Makefile.local`
(gitignored, included automatically) and then plain `make` does the right
thing:

```make
# Makefile.local
PG_CONFIG   = /Applications/Postgres.app/Contents/Versions/18/bin/pg_config
POSTGIS_SRC = /Users/you/src/postgis-3.6.3
```

Then, in the database:

```sql
CREATE EXTENSION postgis;          -- required first
CREATE EXTENSION postgis_hex9;
ALTER EXTENSION postgis_hex9 UPDATE TO '2.0.0'; -- latest version
```

- `HEX9_PREFIX` must match the prefix you installed `libhex9` to (default
  `/usr/local`). The module records an rpath to `$(HEX9_PREFIX)/lib`, so
  `libhex9.{so,dylib}` must be present on the **database host** at runtime.
- OpenMP acceleration comes for free — it lives inside the linked `libhex9`.
- Verify against a running server: `make installcheck` (needs postgis loaded).

### Upgrading from 1.x — rebuild the extension, and delete the old library

Upgrading libhex9 alone is **not enough**, and the failure is silent.

The shared library carries an `SOVERSION`, so 2.0.0 installs as
`libhex9.2.dylib` / `libhex9.so.2` — a *different filename*. It does not
replace `libhex9.0.*`, which stays on disk. An extension module built against
1.x goes on resolving and loading 1.x quite happily, so the database keeps
producing **1.x addresses** while the control file, `h9_version()` and every
release note say 2.0.0. Since 2.0.0 changed the projection, that is wrong data,
not a cosmetic version skew.

`postgis_hex9` 2.0.0 detects this and refuses to load, telling you what to do.
The fix:

```sh
# 1. install the new library
cmake --install build

# 2. remove superseded libhex9 shared libraries — do NOT skip this. While they
#    remain, a rebuild may link 1.x again, because -L/usr/local/lib can win the
#    linker's search order over your build tree.
rm /usr/local/lib/libhex9.0.*        # adjust for your prefix / .so on Linux

# 3. rebuild and reinstall the extension against the new library
cd extension/postgis_hex9 && make clean && make && sudo make install
```

Then in the database:

```sql
ALTER EXTENSION postgis_hex9 UPDATE TO '2.0.0';
SELECT h9_version();   -- should report 2.0.0 and libhex9 2.0.0
```

> **Your stored addresses are now stale.** They are not updated by the upgrade
> and must be re-derived from source geometry — never remapped by
> decode-and-re-encode. Read [`docs/warp-regimes.md`](docs/warp-regimes.md)
> before touching any table that stores Hex9 addresses.

---

## 3. The Python module (`hex9_ext`, nanobind)

Built by the same CMake build when nanobind is available — it links the
in-tree `hex9` target, so no separate install of the C library is required.

```sh
pip install nanobind            # enables the module
cmake -S . -B build             # configures with HEX9_PYTHON=ON (default)
cmake --build build -j
```

The compiled module lands in the build tree (it is **not** installed into a
system prefix). Use it from there:

```sh
PYTHONPATH=build python -c "import hex9_ext; print(hex9_ext)"
```

Disable with `-DHEX9_PYTHON=OFF`. Smoke test: `ctest --test-dir build -R smoke_py`.

---

## 4. The PDAL plugins (`filters.hex9`, `filters.hex9bin`)

Two PDAL filter plugins built by the same CMake build, **opt-in** with
`-DHEX9_PDAL=ON` (requires PDAL — dev headers and the `pdal` CLI):

- **`filters.hex9`** — per-point `WGS84 lon/lat → w_oct` (the 3D storage CRS).
- **`filters.hex9bin`** — aggregate points into Hex9 cells (population-ceiling
  adaptive digest) with `mean/min/max/mode` rollups.

```sh
cmake -S . -B build -DHEX9_PDAL=ON
cmake --build build -j
```

PDAL discovers a plugin by filename (`libpdal_plugin_filter_<name>.{so,dylib}`,
`pdal_plugin_filter_<name>.dll` on Windows — CMake gets the naming right) in
one of two places:

1. **A directory on `PDAL_DRIVER_PATH`** — simplest for development, no
   install needed:

   ```sh
   export PDAL_DRIVER_PATH=$PWD/build/pdal
   ```

2. **PDAL's own plugin dir** — zero-config discovery. Install the plugins
   there:

   ```sh
   cmake -S . -B build -DHEX9_PDAL=ON \
         -DHEX9_PDAL_INSTALL_DIR=$(pdal-config --prefix)/lib
   cmake --build build -j
   cmake --install build            # installs the plugins + libhex9
   ```

   The plugins link `libhex9`; `cmake --install` installs it too and sets the
   plugins' rpath so they resolve it. (If you install `libhex9` to a
   non-standard prefix, that dir must be on the runtime library path.)

Verify either way:

```sh
pdal --drivers | grep hex9          # -> filters.hex9, filters.hex9bin
pdal --options filters.hex9bin
```

When the `pdal` CLI is present, the build adds guarded ctests (`pdal_filter`,
`pdal_hex9bin`) that run real pipelines through the plugins and verify against
the C ABI: `ctest --test-dir build -R pdal`.

See `pdal/examples/README.md` for an end-to-end LiDAR → Hex9 → maps walkthrough.

---

## 5. The Rust crate (`geoplegma/`, static FFI)

The `hex9-sys` crate at `geoplegma/` builds `libhex9` **itself** (via the
`cmake` crate, `hex9Static` target) and generates FFI bindings from `hex9_c.h` —
nothing needs to be pre-installed.

```sh
cd geoplegma
cargo build
```

It links the static archive (no OpenMP, no external runtime dependency; the
warp blob is embedded at compile time), so the resulting binary is
self-contained.
