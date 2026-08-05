# hex9

**Deterministic hexagonal addressing for the whole planet.** hex9 is the
Python package for the Hex9 / H9 discrete global grid system (DGGS): an
aperture-9 hierarchical hexagonal grid on the AKW octahedral equal-area
projection, with 128-bit cell addresses that are **bit-identical on every
platform** — the same lon/lat mints the same uuid on your laptop, your
cluster, and your colleague's Windows machine, enforced by cross-platform
golden tests on every release.

```python
import numpy as np
import hex9

lon = np.array([-3.19, -0.1276])
lat = np.array([55.95, 51.5072])

uuids = hex9.encode(lon, lat)         # (n, 16) uint8 — full-depth addresses
hex9.decode(uuids)                    # cell centroids back to lon/lat
hex9.bin(uuids, 8)                    # layer-8 bins: join keys for analytics
```

All operations are numpy-vectorised and release the GIL.

## What this package is (and is not)

`hex9` is the compiled **libhex9** core — the *address authority* of the
Hex9 system, shared bit-for-bit with the PostGIS extension, the DuckDB
community extension (`INSTALL hex9 FROM community`), and the GeoPlegma Rust
adapter. It contains the full addressing surface: encode/decode, layer
binning, cell lineage and ancestry, k-rings, labels, WKB cell geometry,
grid enumeration, Hamiltonian-curve locality ordering, adaptive density
grids, the sphere-datum twins, and the face-coordinate
(bring-your-own-projection) research seam.

It contains **no code from `hhg9`**, the pure-Python research
implementation of Hex9. hhg9 is where the system is explored — rendering,
nets, composition, cartography, experiments; `hex9` is where addresses are
minted. If you store or share addresses, mint them here (or through any
libhex9-backed surface): pure-Python encoding is not bit-guaranteed across
platforms at full depth.

## The determinism contract

The encode chain's floating-point program is the *definition* of the
address (vendored deterministic math kernels; pinned FP contraction). Every
wheel, on every platform, must pass frozen full-depth golden tests before
it can be published. Practical consequences:

- Addresses are safe to store, share, join on, and deduplicate across
  heterogeneous fleets.
- Any future change to the chain is a versioned regime change, never a
  silent drift.

## Layers, addresses, bins

A full address identifies a ~50 nm cell at layer 30; `hex9.bin(uuid, L)`
coarsens it to the containing cell at layer L (a layer-8 cell is roughly
county-sized, layer-14 city-block-sized). Bins derived from a full address
are exact prefixes of the hierarchy — aggregation is a byte operation, not
a spatial query.

## History

One line per release; the full story is the
[CHANGELOG](https://github.com/MrBenGriffin/libhex9/blob/main/CHANGELOG.md).

- **2.3.0** — E4H aperture-4 tails (exact classifier, addresses to depth
  28) with hexagon binning (`e4h_hex` — one canonical key per fine
  hexagon, 4^d per host); grid verbs `hex9.verbs` (aim / walk_to /
  vision_cone) + true-bearing correction; universality made citable
  (`docs/universality.md`, `hex9.selftest()`); curve labels; marker
  guards (curve/E4H input rejected by h9 machinery).
- **2.2.1** — packaging/metadata only; no address moves.
- **2.2.0** — face-coordinate addressing (`encode_boct`/`decode_boct`):
  bring your own sphere→octahedron projection.
- **2.1.0** — the deterministic chain: bit-identical addresses on every
  platform (musl-derived kernels, FMA off); sphere-datum twins.
- **2.0.0** — via-sphere regime (REGIME CHANGE: addresses moved; re-derive
  from source geometry, never remap).

## Links

- Source, C ABI, and doctrine: <https://github.com/MrBenGriffin/libhex9>
- The projection (AKW) and the layer architecture:
  [docs/projection.md](https://github.com/MrBenGriffin/libhex9/blob/main/docs/projection.md)
- Changelog: [CHANGELOG.md](https://github.com/MrBenGriffin/libhex9/blob/main/CHANGELOG.md)

Apache-2.0. © Ben Griffin. Includes deterministic math kernels derived from
musl libc / fdlibm (see COPYRIGHT).
