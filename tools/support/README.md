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

## ✅ RESOLVED (2026-06-18): the v3 producer is `export_warp_v3.py`

The library embeds **`core/WGS84_l5_warp_f6.full.f64g.h9warp` — format v3, with
bundled per-vertex CT gradients** (6 × f64 per vertex; see the v3 path in
`core/h9_warp_io.h`). Its producer is **`export_warp_v3.py`** (in this dir).

That script was authored + run inline in a prior libhex9 Claude session and
produced the shipped blob, but was **never committed as a file** — which is why
the earlier search (active repo + whole Archive + libhex9 git `--all`) and
`grep "full.f64g"` found nothing: it lived only in the session transcript.
Recovered verbatim on 2026-06-18 and **confirmed as THE producer**:

- delta columns (blob `[:, 0:2]`) ≡ `target_pts − source_pts` — max diff `0.0`
- gradient columns (blob `[:, 2:6]`) ≡ `grad_dx` / `grad_dy` — max diff `0.0`

over all 266,815 verts (npz already in tri_mesh order → identity perm), CRC valid.

`export_warp_deltas.py` (v2, no gradients) is retained as the format ancestor /
mirrored-half reference, but `export_warp_v3.py` is what regenerates the shipped
blob.

## Full v3 pipeline (every step is now a tracked file)

```
(sinkhorn training, in hhg9)   -> hhg9/data/WGS84_l5_warp_data.npz     (deltas)
export_warp_grads.py           -> experimental/sinkhorn/output/WGS84_l5_warp_grads.npz
export_warp_v3.py              -> core/WGS84_l5_warp_f6.full.f64g.h9warp
CMake .incbin (H9_WARP_BLOB)   -> blob embedded in libhex9
```

- `export_warp_grads.py` — reads the FINAL per-vertex CT gradients off the hhg9
  `AuthalicWarp` loader (after ghost-orbit padding + x-mirror symmetrisation +
  edge-tangent projection). Recovered from the brief's snippet; it had never been
  a committed script.
- `warp-port-brief-f6-cside.md` — the F6 C-side port brief: full rationale,
  "Option A" (ship gradients) which both scripts implement, the construction
  steps (fallback Option B), spot-check values, and acceptance criteria.
  **Was git-ignored in `.claude/`** (would be lost on a clean checkout) — copied
  here on 2026-06-18 to co-locate it with the scripts that implement it and make
  it durable. Both scripts' docstrings reference this path.
