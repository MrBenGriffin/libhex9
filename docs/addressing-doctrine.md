# Hex9 addressing doctrine

Referenced from `core/h9_addressing.h`, `core/h9_kring.h`, `core/h9_grid.h`,
`hex9_c.h`, and `tools/diagnosis/fossil_probe.py`. This document records the
addressing contract, the three distinct roll-up operations, and the status of
the historical "fossil" caveats (F1–F4). Paper cross-references are to the
Hex9 paper (§10b, §12, §13) in the `hex9` repository.

Last full verification: 2026-07-08 (see *Verification artefacts* below).

## The contract: canonicalise-always

Every emitted address — full UUID or layer bin — presents its terminal cell
through the cell's **mode-0 parent** (the canonical fold). Consequences:

- A full UUID *is* the canonical bin at maximum depth; "address vs bin"
  dissolves into one decode path (the reclaimed layout).
- The key tail is a single nibble `(p_c2 << 1) | r_mo` with `p_mo = 0`
  implicit; backward passes over canonical input may seed `c_mo = 0`
  directly (the two-way c_mo recovery dance is a legacy accommodation).
- Canonical digit **bodies never collide**: no two cells at the same layer
  share a tail-stripped body — machine-verified exhaustively L0..L5
  (708,588 bodies at L5, all distinct; `hex9:experimental/body_census.py`).
- A body is nevertheless a **lineage path, not ancestry**: cutting digits
  walks lineage, which at split cells (terminal digit 6/7/8) differs from
  canonical ancestry. Worked example (paper §10b): cell `43585` cuts to
  `4358`, but its canonical L3 ancestor is the hexagon canonically named
  `4348` — `4358` canonically names a different hexagon elsewhere.

## The three roll-up operations — never conflate them

1. **Point roll-up** — *which layer-K cell contains this point.*
   `h9_bin_uuid` (address-arithmetic, from the FULL uuid) / `hex9_bin`.
   Exact; verified over 2M random points: a split cell's points resolve to
   exactly its two geometric parents, never more, no tie-break scatter.
2. **Canonical cell roll-up** — *which single layer-K cell is the canonical
   ancestor of this cell.* The canonical parent of a cell is the parent
   containing its **mode-0 d_cell**; multi-level ancestry is the *iterated*
   one-level parent (a single deep re-bin answers the point question
   instead, and differs on exactly 1/9 of cells).
   `h9kring::h9_cell_parent_uuid` / `h9_cell_ancestor_uuid`
   (`core/h9_kring.h`), C ABI `hex9_cell_parent{,_many}` /
   `hex9_cell_ancestor{,_many}`, ext `cell_parent` / `cell_ancestor`;
   Python reference `hhg9.h9.uuid_address.h9_cell_parent` /
   `h9_cell_ancestor`. Byte-identical across the two implementations for
   every cell L1..L4; exactly **9 canonical children per parent** for every
   cell at every layer pair through L5.
3. **Lineage cut** — raw digit truncation. A grouping by lineage, not by
   cell ancestry; the two coincide away from split lineages and diverge on
   the split-cell band, which thins as (1/6)·3^(1−k) with depth difference
   k (paper §12).

The recurring disease, three sightings, one cure: a split cell's centroid
lies exactly ON the boundary it straddles, so *descending or re-binning a
centroid* is classification-unstable (Python `h9_bin` on cell uuids
scattered 6–12 children per parent; `mesh.addr` had a sibling bug;
`uuid_from_cxcy` without the fold mis-parents 180/972 at L2). The cure is
always the same: operate from a point strictly interior to the **mode-0
half** (centroid nudged 0.10 toward the cell origin), or equivalently apply
the canonical fold. New code touching cell identity must use one of these.

## Fossil register

- **F1 — bare-label ambiguity** (`hex9_c.h` parse/common-ancestor caveats).
  *Status: UNDER REVIEW, likely retirable.* The caveat predates the body
  census. Canonical bodies are unique (above), and the six-tail parse with
  canonical re-encode verification empirically yields exactly one fixed
  point per body, split terminals included. Retirement needs a dedicated
  probe of `hex9_parse_label` / `hex9_common_ancestor` over all split
  bodies at a working layer; until then the caveat stands in the headers.
  2026-07-18 (Ben's ruling: the frozen L29 expectations are deprecated):
  the SQL regression's F1 pin now asserts the L30 reclaimed behaviour —
  `h9_parse_label('43')` returns Westminster's own L1 cell. The historical
  L29 behaviour survives only in `test/golden_l29.py` territory.
- **F2 — decode(bin) mislocation** (`hex9_c.h:hex9_decode`).
  *Status: STANDING as a caveat; does not reproduce on the L30 reclaimed
  layout* (the deepest address IS the canonical max-depth bin, so the
  identity machinery no longer guesses meta). The SQL regression's F2 pin
  asserts the L30 behaviour (Westminster's L1 bin decodes in place,
  2026-07-18 ruling). Decode remains guaranteed for full UUIDs only.
  Probe: `tools/diagnosis/fossil_probe.py`.
- **F3 — bin→coarser re-binning** (`h9_addressing.h`, `h9_kring.h`,
  `h9_grid.h`). *Status: RETIRED 2026-07-08.* The supported cell-level
  roll-up is `h9_cell_parent_uuid` / `h9_cell_ancestor_uuid` (operation 2
  above). Feeding a bin to `h9_bin_uuid` at a coarser layer remains
  unsupported — it was never the right question. `test/bin_prefix_guard.c`
  guards the old behaviour; the new functions carry the invariant tests.
  2026-07-18: the near-vertex garbage case no longer reproduces on L30
  (the SQL F3 pin asserts agreement now), but the path stays UNSUPPORTED —
  `hex9_bin(bin_L8, 1)` mislabelled GB's `43` ancestor as `65` on L30, and
  `hex9_common_ancestor` was fixed that day to label from the shared
  nibble prefix directly instead of re-coarsening through this path.
- **F4 — grid identity UUIDs** (`h9_grid.h:full_id_from_cell`).
  *Status: STANDING, by design.* Grid cells get full-depth identity UUIDs
  by descending a mode-0-half interior point (the same cure). decode(id)
  lands ~10% of a circumradius off the geometric centroid — use the
  geometric centroid for display, the UUID for the key.

## Marker registry (uuid kinds in one column)

Three uuid kinds coexist byte-distinguishably; the sentinels cannot occur
in each other's grammars, so every test below is positional/decisive:

| kind        | marker                        | operations                     |
|-------------|-------------------------------|--------------------------------|
| h9 cell     | body nibbles ≤ 11, pad `0xF`  | the h9 surface                 |
| curve       | nibble 0 = `0xC`              | `hex9_curve_*` (`hex9_is_curve`) |
| E4H tail    | any nibble = `0xE`            | `hex9_e4h_*` (`hex9_is_e4h`)   |

E4H tails are ADDRESSES (ruling 2026-08-04), not bins — the 2026-07-15
bin framing is superseded. Truncating an E4H tail IS binning (suffix-local,
exact) — the one place truncation is legitimate; the a9 body fossils
(F1–F3) still apply unchanged to the host part. Every h9/curve entry point
rejects `0xE`-marked input (the 2.3.0 marker-guard doctrine): tails have
their own operations, and the h9 machinery must never silently mis-read
one as a deep h9 body.

## Verification artefacts

- `hex9:experimental/body_census.py` — body uniqueness, L0..L5, all cells.
- `hex9:experimental/cell_ancestor_verify.py` — 9-per-parent invariant
  L1..L5 + composition identity (ancestor(L−2) = parent-of-parent).
- `hex9:tests/test_h9_uuid_address.py` — includes the London known-answer
  (`cell_parent(43585.1) = 4348.2`) and cell-ancestry property tests.
- `test/` (this repo) — ctest suite incl. `bin_prefix_guard.c`,
  `canonical_invariants.py`; byte-compare C↔Python is re-runnable from the
  cell_ancestor_verify harness.
- `tools/diagnosis/fossil_probe.py` — the historical failure-mode cases,
  including the near-vertex point (0.0091, −89.9863).
