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

Then, in the database:

```sql
CREATE EXTENSION postgis;          -- required first
CREATE EXTENSION postgis_hex9;
ALTER EXTENSION postgis_hex9 UPDATE TO '1.5.0'; -- latest version
```

- `HEX9_PREFIX` must match the prefix you installed `libhex9` to (default
  `/usr/local`). The module records an rpath to `$(HEX9_PREFIX)/lib`, so
  `libhex9.{so,dylib}` must be present on the **database host** at runtime.
- OpenMP acceleration comes for free — it lives inside the linked `libhex9`.
- Verify against a running server: `make installcheck` (needs postgis loaded).

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
