# Hex9 PDAL examples

End-to-end demo: real survey LiDAR → reproject → **Hex9** cells → maps.

The plugins (`filters.hex9`, `filters.hex9bin`) are built with the libhex9 tree
when configured with `-DHEX9_PDAL=ON`:

```sh
cmake -S . -B build -DHEX9_PDAL=ON
cmake --build build
export PDAL_DRIVER_PATH=$PWD/build/pdal     # so `pdal` finds the plugins
```

- **`filters.hex9`** — per-point `WGS84 lon/lat → w_oct` (the 3D storage CRS).
- **`filters.hex9bin`** — aggregate points into Hex9 cells (population-ceiling
  adaptive digest) with `mean/min/max/mode` rollups. See `bin_to_las.json`.

## Making PDAL find the plugins

PDAL discovers a plugin from a shared library named `libpdal_plugin_filter_<name>`
in one of:

1. a directory on **`PDAL_DRIVER_PATH`** (simplest for dev — no install):
   ```sh
   export PDAL_DRIVER_PATH=$PWD/build/pdal
   ```
2. **PDAL's own plugin dir** (zero-config — where its bundled plugins live, e.g.
   `$(pdal-config --prefix)/lib`). Install there:
   ```sh
   cmake -S . -B build -DHEX9_PDAL=ON \
         -DHEX9_PDAL_INSTALL_DIR=$(pdal-config --prefix)/lib
   cmake --build build
   cmake --install build            # installs the plugins + libhex9
   ```
   The plugins link `libhex9`; `cmake --install` also installs it and sets the
   plugins' rpath so they resolve it. (If you install libhex9 to a non-standard
   prefix, that dir must be on the runtime library path.)

Verify either way:
```sh
pdal --drivers | grep hex9          # -> filters.hex9, filters.hex9bin
pdal --options filters.hex9bin
```

## Platforms

The plugin builds the same way on all three; only the shared-library naming and
the runtime resolution of `libhex9` differ (CMake handles the naming for you):

| | plugin file | finds `libhex9` at runtime via |
|---|---|---|
| **macOS**   | `libpdal_plugin_filter_hex9.dylib` | `INSTALL_RPATH` (set); `DYLD_LIBRARY_PATH` in-tree |
| **Linux**   | `libpdal_plugin_filter_hex9.so`    | `INSTALL_RPATH` (set); else `LD_LIBRARY_PATH` / `ldconfig` |
| **Windows** | `pdal_plugin_filter_hex9.dll` (no `lib` prefix) | **no rpath** — put `hex9.dll` next to the plugin or on `PATH` |

Notes:
- Get PDAL from your platform's packages (Linux: apt / conda-forge; Windows:
  conda-forge / OSGeo4W). `find_package(PDAL)` + `-DPDAL_DIR=<pdal>/lib/cmake/PDAL`
  works cross-platform; on Windows build the plugin with the **same MSVC toolchain**
  PDAL was built with.
- Windows has no rpath, so after `cmake --install` ensure `hex9.dll` is on `PATH`
  (or copy it beside the plugin) — the plugin won't load otherwise.
- Only macOS is CI-tested here; Linux/Windows are expected-to-work but unverified.

## Demo scripts

Input is a folder of Environment Agency `.laz` tiles (British National Grid,
EPSG:27700) — e.g. a 3×3 OS-grid cluster. The scripts reproject to WGS84 and
render Hex9 maps. They use the `hex9_ext` Python binding for the multi-layer
render; set `HEX9_BUILD` if the build dir isn't `../../build`.

```sh
# 1. merge + decimate + reproject every .laz in a dir  (drop a 3x3 cluster in)
python3 extract_cluster.py  test_data/2015-TQ28ne  cluster_pts.csv  200

# 2. annotated cover + elevation map at a fixed layer (bounds, projection, N arrow)
python3 render_map.py  cluster_pts.csv  cluster_map_L11.png  11

# 3. the DGGS story: fixed sweep L7..L11 + adaptive population-ceiling digest
python3 render_seq.py  cluster_pts.csv  cluster_seq.png
```

`extract_cluster.py` globs `*.laz`/`*.las`, so adding more tiles to the folder
just widens the map — no code change. Decimation `step` trades speed for
density (encode is the bottleneck; ~15 s / 470k pts on 8 threads).

`bin_to_las.json` is a pure-PDAL pipeline using `filters.hex9bin` directly
(no Python), emitting one LAS point per cell with the aggregates as dimensions.
