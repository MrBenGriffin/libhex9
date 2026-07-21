# libhex9

The reference implementation of **Hex9 (H9)** — an equal-area hexagonal
Discrete Global Grid System with a 9-fold (3×3) hierarchy, from layer 0
(12 base cells) to layer 30 (~32 nm cell diameter).

libhex9 is a **geometry-free core** behind a stable C ABI (`hex9_c.h`). It
depends on no geometry library, no projection library and no database. Every
consumer — the PostGIS extension, the PDAL filters, the Python module, the
Rust crate — builds on that one surface, so the grid can never drift between
them.

```c
#include <hex9_c.h>

uint8_t uuid[16];
hex9_warp_init(NULL, 0);
hex9_encode(-3.19, 55.95, uuid);   /* lon, lat → cell address */
```

## What makes a Hex9 address

A cell identifier is a plain 16-byte UUID — it indexes, joins and `GROUP BY`s
with no custom type anywhere in the stack. It is *self-contained*: the address
carries its own layer and lineage, so a cell needs no external table to be
interpreted, and the layer-`n` bin of an address is a pure prefix operation.

The text form is `"<x_list>.<T>"`, e.g. `435878503.3`.

On the default L30 layout the old terminal nibble is reclaimed as a real body
digit, so a full UUID *is* its own deepest-layer bin. The legacy 29-layer
on-disk format remains available with `-DHEX9_USE_L29=ON`.

> **Addresses are not labels, and labels are not addresses.** The distinction
> is load-bearing and the failure modes are catalogued in
> [`docs/addressing-doctrine.md`](docs/addressing-doctrine.md). Read it before
> designing anything that persists a bin label.

## Capabilities

| Area | What it gives you |
|---|---|
| **Point addressing** | `lon/lat` → cell, single and batch (OpenMP where available) |
| **Hierarchy** | canonical parent / ancestor; the two exact relations — *lineage* and *ownership* |
| **Hamiltonian curve** | sortable curve-uuids, numeric indices, labels, exact prefix coarsening, curve-ordered descendants |
| **Geometry** | cell vertices and centroids; grid enumeration over a region |
| **Adaptive grids** | population-ceiling multi-layer digest |
| **Neighbours** | k-ring / k-disk |
| **`w_oct`** | warped octahedral cartesian CRS — the seamless 3D baseline for point-cloud storage |

The equal-area (authalic) warp holds cell areas uniform to ~0.014 % globally.

## Consumers

| | |
|---|---|
| `extension/postgis_hex9` | PostGIS extension — cells as native `uuid` columns |
| `pdal/` | `filters.hex9`, `filters.hex9bin` — LiDAR / point-cloud binning |
| `python/` | nanobind module |
| `geoplegma/` | `hex9-sys` Rust FFI crate, consumed by GeoPlegma |

## Build

```sh
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
cmake --install build
```

See [`INSTALL.md`](INSTALL.md) for the PostGIS extension, the PDAL plugins,
the Python module and the Rust crate — each builds *on top of* the installed
C library.

## Warp regimes — read before storing addresses

An address only means something in the context of the projection that produced
it, and **nothing in the 16 bytes records which projection that was.**

There is exactly one regime — 2.0.0 removed the second rather than let a point
have two indistinguishable addresses. But 2.0.0 *is* that change: addresses
written by 1.x and by 2.x differ from about layer 7 downward, and look
identical. Coarse bins often do not move at all, so a spot check proves
nothing.

If you hold addresses written by 1.x, **re-derive them from your source
geometry — never decode-and-re-encode.** Remapping is not idempotent, cannot
detect a second run, and silently invalidates anything already derived from
those addresses. See [`docs/warp-regimes.md`](docs/warp-regimes.md).

## Versioning

libhex9, `postgis_hex9`, and `hex9-sys` share one version number. `h9_version()`
reports the loaded build, stamped with its git revision.

## Licence

Apache-2.0 — see [`LICENSE`](LICENSE) and [`COPYRIGHT`](COPYRIGHT).
[`LLM.md`](LLM.md) records AI-tooling provenance.
