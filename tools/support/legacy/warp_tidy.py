#!/usr/bin/env python3
"""
warp_tidy.py — inspect and clean WGS84_l5_warp_data.npz.

Default action (no flags) is a structural report: vertex count, x/y range,
boundary-line populations, mirror-symmetry residual, |delta| histogram,
NaN/Inf check.

Actions:
    --zero-equator     zero the displacement at every point on y == Y_EQ
                       (target_pts <- source_pts). 730 points (L5); these
                       sit on the b_oct triangle's top edge and accumulate
                       trainer noise in dx that doesn't contribute usefully
                       to CT (dy is already 0 by symmetry). Preserves the
                       266,815-point shape so the tri_mesh premise stays.
    --verify           build AuthalicWarp before vs after the tidy, sample
                       across the b_oct domain, report worst-case drift.

Output:
    Writes to {stem}.tidied.npz by default; --out overrides; --dry-run
    skips writing.

Usage:
    python warp_tidy.py                                  # report only
    python warp_tidy.py --zero-equator --verify          # tidy + verify
    python warp_tidy.py WGS84_l5_warp_data_ot.npz        # report a variant
"""

import argparse
import os
import sys
import math
import numpy as np


Y_EQ = math.sqrt(6.0) / 6.0
R3   = math.sqrt(3.0)
W_EQ = math.sqrt(6.0) / 3.0   # |Ẇ| — lateral edges are y = R3·|x| - W_EQ

# 1 b_oct unit ≈ this many ground metres (mode-0 b_oct has area = 1 in its own
# units; one full octant is 1/8 of the sphere's surface).
GROUND_PER_BOCT = 2.5e7


def _hist(values, edges):
    """Counts in [edges[i], edges[i+1]); last bucket inclusive on right."""
    return np.histogram(values, bins=edges)[0]


def report(src, dst, label):
    n = src.shape[0]
    diff = dst - src
    print(f"\n=== {label} ===")
    print(f"  vertices       : {n:,}")
    print(f"  x range        : [{src[:,0].min():+.6f}, {src[:,0].max():+.6f}]")
    print(f"  y range        : [{src[:,1].min():+.6f}, {src[:,1].max():+.6f}]")
    print(f"  Y_EQ           : {Y_EQ:.6f}")

    nan_inf = ~(np.isfinite(src).all(axis=1) & np.isfinite(dst).all(axis=1))
    print(f"  NaN/Inf rows   : {int(nan_inf.sum())}")

    # Boundary populations
    on_eq = np.isclose(src[:, 1], Y_EQ, atol=1e-14)
    on_right = np.isclose(src[:, 1] - R3 * src[:, 0], -W_EQ, atol=1e-14)
    on_left  = np.isclose(src[:, 1] + R3 * src[:, 0], -W_EQ, atol=1e-14)
    on_apex  = np.isclose(src[:, 0], 0.0, atol=1e-14) & \
               np.isclose(src[:, 1], -W_EQ, atol=1e-14)
    print(f"  on equator y=Y_EQ : {int(on_eq.sum())}  "
          f"on right edge: {int(on_right.sum())}  "
          f"on left edge: {int(on_left.sum())}")

    # Delta magnitudes
    mag = np.linalg.norm(diff, axis=1)
    print(f"  |delta|        : max {mag.max():.3e}  p99 {np.percentile(mag,99):.3e}  "
          f"median {np.median(mag):.3e}")
    edges = np.array([0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
    h = _hist(mag, np.concatenate([edges, [np.inf]]))
    labels = ['<1e-6', '1e-6..1e-5', '1e-5..1e-4', '1e-4..1e-3', '1e-3..1e-2', '>=1e-2']
    pieces = [f"{l}: {c:,}" for l, c in zip(labels, h)]
    print(f"  |delta| hist   : {'  '.join(pieces)}")

    # Boundary-vertex delta magnitudes
    def grp(mask, name):
        if not np.any(mask):
            return
        m = mag[mask]
        d = diff[mask]
        print(f"  {name:15s}: n={int(mask.sum()):>5,}  "
              f"max|dx|={np.abs(d[:,0]).max():.3e}  "
              f"max|dy|={np.abs(d[:,1]).max():.3e}  "
              f"max|delta|={m.max():.3e}")

    grp(on_eq,   "  @equator")
    grp(on_right,"  @right-edge")
    grp(on_left, "  @left-edge")
    grp(on_apex, "  @apex")

    # Mirror-symmetry residual: for each point with x>0, find its mirror
    # at (-x, y) and compare deltas. Quantize via the same scheme as
    # _unique_rows_tol so the lookup is exact.
    tol = 1e-12
    k = int(np.ceil(-np.log2(tol)))
    qx = np.rint(np.ldexp(src[:, 0], k)).astype(np.int64)
    qy = np.rint(np.ldexp(src[:, 1], k)).astype(np.int64)
    key = {(int(qx[i]), int(qy[i])): i for i in range(n)}
    pos = np.where(qx > 0)[0]
    res_dx, res_dy = [], []
    missing = 0
    for i in pos:
        j = key.get((-int(qx[i]), int(qy[i])))
        if j is None:
            missing += 1
            continue
        # Expected: dx[i] == -dx[j], dy[i] == dy[j]
        res_dx.append(abs(diff[i, 0] + diff[j, 0]))
        res_dy.append(abs(diff[i, 1] - diff[j, 1]))
    if res_dx:
        rdx = np.asarray(res_dx)
        rdy = np.asarray(res_dy)
        print(f"  mirror residual: pairs={len(rdx):,}  missing={missing}  "
              f"max(|dx+dx_m|)={rdx.max():.3e}  max(|dy-dy_m|)={rdy.max():.3e}")
        print(f"                   p99 dx={np.percentile(rdx,99):.3e}  "
              f"p99 dy={np.percentile(rdy,99):.3e}")
    else:
        print("  mirror residual: no x>0 pairs found")


def zero_equator(src, dst):
    """Set dst[on_eq] = src[on_eq]. Returns (dst_new, n_changed, dx_was, dy_was)."""
    on_eq = np.isclose(src[:, 1], Y_EQ, atol=1e-14)
    n = int(on_eq.sum())
    if n == 0:
        return dst.copy(), 0, 0.0, 0.0
    dx_was = float(np.abs(dst[on_eq, 0] - src[on_eq, 0]).max())
    dy_was = float(np.abs(dst[on_eq, 1] - src[on_eq, 1]).max())
    dst_new = dst.copy()
    dst_new[on_eq] = src[on_eq]
    return dst_new, n, dx_was, dy_was


def verify_drift(npz_path, src, dst_old, dst_new, sample=4000, seed=42):
    """Build AuthalicWarp instances around two npz blobs and compare warp
    output on `sample` random b_oct interior points.

    AuthalicWarp loads from a file, so we materialise two temporary npz files.
    """
    import tempfile
    from hhg9.domains.octahedral_barycentric import AuthalicWarp

    rng = np.random.default_rng(seed)
    # Sample inside the mode-0 b_oct triangle:
    #   y in [-W_EQ, Y_EQ];  y >= R3·|x| - W_EQ
    xs = rng.uniform(-math.sqrt(2) / 2, math.sqrt(2) / 2, sample * 3)
    ys = rng.uniform(-W_EQ, Y_EQ, sample * 3)
    inside = (ys >= R3 * np.abs(xs) - W_EQ) & (ys <= Y_EQ)
    pts = np.column_stack([xs[inside][:sample], ys[inside][:sample]])
    print(f"\n[verify] sampling {len(pts)} interior b_oct points")

    with tempfile.TemporaryDirectory() as td:
        path_old = os.path.join(td, 'old.npz')
        path_new = os.path.join(td, 'new.npz')
        # Preserve `layer` etc. so AuthalicWarp can load without surprises.
        repo = np.load(npz_path, allow_pickle=True)
        extras = {k: repo[k] for k in repo.files
                  if k not in ('source_pts', 'target_pts')}
        repo.close()
        np.savez(path_old, source_pts=src, target_pts=dst_old, **extras)
        np.savez(path_new, source_pts=src, target_pts=dst_new, **extras)

        warp_old = AuthalicWarp(path_old)
        warp_new = AuthalicWarp(path_new)

        do_old = warp_old.do(pts, mo=0)
        do_new = warp_new.do(pts, mo=0)

        mask = np.isfinite(do_old).all(axis=1) & np.isfinite(do_new).all(axis=1)
        if not np.any(mask):
            print("  no finite warp outputs; cannot compare")
            return
        diff = np.linalg.norm(do_old[mask] - do_new[mask], axis=1)
        ground = diff * GROUND_PER_BOCT
        print(f"  do()  drift  : max {diff.max():.3e} b_oct  "
              f"({ground.max() * 1e3:.3f} mm on ground)")
        print(f"               : p99 {np.percentile(diff,99):.3e}  "
              f"median {np.median(diff):.3e}")

        # Verify the round trip stayed stable too: undo(do(pt)) ≈ pt
        rt_old = warp_old.undo(do_old, mo=0)
        rt_new = warp_new.undo(do_new, mo=0)
        rt_diff_old = np.linalg.norm(rt_old[mask] - pts[mask], axis=1)
        rt_diff_new = np.linalg.norm(rt_new[mask] - pts[mask], axis=1)
        print(f"  round-trip   : old max {rt_diff_old.max():.3e}  "
              f"new max {rt_diff_new.max():.3e} b_oct")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('npz', nargs='?', default='WGS84_l5_warp_data.npz',
                    help='Input npz (default: WGS84_l5_warp_data.npz)')
    ap.add_argument('--zero-equator', action='store_true',
                    help='Zero displacement at y==Y_EQ points')
    ap.add_argument('--verify', action='store_true',
                    help='Build AuthalicWarp before/after and report drift')
    ap.add_argument('--out', default=None,
                    help='Output path (default: {stem}.tidied.npz)')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report but do not write')
    args = ap.parse_args()

    in_path = args.npz
    if not os.path.exists(in_path):
        print(f"error: {in_path} not found", file=sys.stderr)
        return 2

    repo = np.load(in_path, allow_pickle=True)
    src = np.asarray(repo['source_pts'], dtype=np.float64)
    dst = np.asarray(repo['target_pts'], dtype=np.float64)
    extras = {k: repo[k] for k in repo.files
              if k not in ('source_pts', 'target_pts')}
    repo.close()

    report(src, dst, in_path)

    actions = []
    if args.zero_equator:
        actions.append('zero-equator')
    if not actions:
        print(f"\nNo --action flags; report-only mode. "
              f"Try --zero-equator [--verify].")
        return 0

    dst_new = dst.copy()
    if args.zero_equator:
        dst_new, n, dx_was, dy_was = zero_equator(src, dst_new)
        print(f"\n[zero-equator] {n} points at y=Y_EQ; "
              f"before: max |dx|={dx_was:.3e}, max |dy|={dy_was:.3e}; "
              f"after: 0, 0")

    if args.verify:
        verify_drift(in_path, src, dst, dst_new)

    if args.dry_run:
        print("\n[dry-run] not writing output.")
        return 0

    out_path = args.out
    if out_path is None:
        stem, ext = os.path.splitext(in_path)
        out_path = f"{stem}.tidied{ext}"

    np.savez(out_path, source_pts=src, target_pts=dst_new, **extras)
    sz = os.path.getsize(out_path)
    print(f"\nWrote {out_path}  ({sz:,} bytes)")
    report(src, dst_new, out_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
