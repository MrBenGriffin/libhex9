# tools/support/legacy — producers for the retired WGS84 warp regime

Everything in this directory produced, validated, or documented the
**WGS84-trained F6 v3 warp field** (`core/WGS84_l5_warp_f6.full.f64g.h9warp`
and its v2 mirrored ancestor). That field was libhex9's addressing regime up to
1.5.0 and was **removed in 2.0.0** along with the blobs it made — see
`docs/warp-regimes.md`.

Nothing here is on the build path. Nothing here produces an artifact the
library still ships.

## Why kept rather than deleted

The projection may change again, if something sufficiently compelling turns up.
These scripts are the worked record of how a trained field is taken from
sinkhorn output to an embedded blob — the gradient extraction, the v3 container
layout, the edge-tangent handling, and the acceptance checks that proved the
shipped blob matched its inputs bit-for-bit. That method survives the field it
was written for.

Git history would hold them either way; a directory is easier not to miss.

## Contents

| File | Was |
|---|---|
| `export_warp_v3.py` | THE producer of the shipped v3 blob |
| `export_warp_grads.py` | per-vertex CT gradients (v3 input) |
| `export_warp_deltas.py` | v2 format ancestor / mirrored-half reference |
| `export_ct_grads.py` | CT gradient arrays (`core/h9_ct.h`) |
| `export_ct_mesh.py` | `core/h9_warp_mesh.h` |
| `export_warp_test_vecs.py` | warp validation header |
| `warp_tidy.py` | `.npz` preprocessing (zero-equator) |
| `warp-port-brief-f6-cside.md` | the F6 C-side port brief — rationale, construction steps, acceptance criteria |

## The live equivalent

The current regime's field is produced by
[`../export_warp_fund_v4.py`](../export_warp_fund_v4.py) (Sphere-L6 v4 fund
wedge), with its reference generators `../gen_fund_warp_ref.py` and
`../gen_via_sphere_ref.py`. Those feed `test/fund_warp`, `test/authalic` and
`test/via_sphere` — the oracle the whole suite rests on.
