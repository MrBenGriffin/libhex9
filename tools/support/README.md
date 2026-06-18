# tools/support — warp/lattice artifact producers (staging for review)

These are the **build-time producer tools** that generate the artifacts libhex9
vendors in `core/`. They were migrated (copied) from the archived `hex9_cli`
project (`/Users/ben/Documents/Archive/hex9_cli/hex9_cli`) on 2026-06-18 so the
generators live alongside the artifacts they produce.

**Status: staged for review.** Not yet wired into the build. They require a
Python environment with the `hhg9` package on the path (and `scipy`/`numpy`) —
they import `hhg9.h9.*` and operate on `WGS84_l5_warp_data.npz`. They are *dev /
regeneration* tools, not part of the C/C++ library or its runtime.

## Why these were brought over (and the diagnostics were not)

libhex9 ships the **generated artifacts** but, until now, **none of the
generators that produce them**. Only the producers / dependencies were migrated;
the cli's `check_*`, `diag_*`, `diff_*`, `compare_grid*`, `trace_*`,
`verify_hexbin`, `warp_check`, `test_*` diagnostics were left in the archive.

## Artifact → producer map

| libhex9 vendored artifact            | Produced by                  |
|--------------------------------------|------------------------------|
| `core/h9_lattice_consts.h`           | `gen_h9_lattice_consts.py`   |
| `core/h9_uv_lattice.h`               | `export_uv_lattice.py`       |
| `core/WGS84_l5_warp.mirrored.f64.h9warp` (v2) | `export_warp_deltas.py` |
| `core/h9_warp_mesh.h`                | `export_ct_mesh.py`          |
| CT gradient arrays / `core/h9_ct.h`  | `export_ct_grads.py`         |
| warp validation header (test dep)    | `export_warp_test_vecs.py`   |
| `.npz` preprocessing (zero-equator)  | `warp_tidy.py`               |
| python reference reader for `core/h9_warp_io.h` | `h9warp_io.py`    |

`core/h9_lattice_gen.h` is a hand-written C++ port of `hhg9/h9/polygon.py`
(`tri_mesh`), not a generated file — no producer to migrate.

## Deliberately NOT migrated

- `gen_warp_embedded_header.py` — emits a C++ byte-array header from a `.h9warp`.
  **Obsolete here:** libhex9 embeds the raw `.h9warp` directly via CMake
  `.incbin` (`H9_WARP_BLOB`, see `core/h9_warp_embedded.cpp`), so no
  header-generation step is needed.
- `export_warp_header.py` — emits the 256×256 `h9_warp_data.h` grid for **PROJ**
  (`proj-9.8.0`), a different downstream consumer, not a libhex9 dependency.

## ⚠️ OPEN: the actually-shipped v3 blob has no producer here

The library actually embeds **`core/WGS84_l5_warp_f6.full.f64g.h9warp` — format
v3, with bundled per-vertex CT gradients** (6 × f64 per vertex; see the v3 path
in `core/h9_warp_io.h`). Verified 2026-06-18: its deltas are bit-identical
(max |diff| = 0.0 over 266,815 verts) to `hhg9/data/WGS84_l5_warp_data.npz`,
CRC valid.

But the migrated `export_warp_deltas.py` writes **v2 only (no gradients)**. A
search of the active python repo, the whole Archive tree, and libhex9's full git
history (`--all`) found **no v3 / f64g / gradient-bundling exporter** — the v3
blob was committed (`bd7c0be`) with no accompanying script.

**TODO (Ben to confirm):** locate the v3 (f64g) exporter — Ben believes the
script exists somewhere; the auto-memory note ("v3 blob ships hhg9 gradients")
may be stale on *where the producer lives*. If it's truly lost, reconstruct it
by merging `export_warp_deltas.py` (v2 deltas, tri_mesh order) with
`export_ct_grads.py` (gradients) into the 6-f64/vertex v3 layout that
`core/h9_warp_io.h` reads. Until then, libhex9 can regenerate the **v2** path but
not the blob it actually ships.
