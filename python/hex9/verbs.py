# Part of the Hex9 (H9) Project
# Copyright ©2026, Ben Griffin
# Licensed under the Apache License, Version 2.0

"""verbs — cell-first grid verbs on the symbolic hex9 grid.

Field primitives built on the core's validated neighbour algebra
(`k_disk` / `neighbors` — exact integer seam crossing, no geometry
sampled per cell). The API is CELL-FIRST: every verb takes and returns
hex KEYS (layer-scoped bin uuids as `bytes`); positions appear only
when a caller binds an agent to the grid. Internals are vectorised —
one encode/bin round per verb, not per cell.

    vision_cone(src_key, layer, bearing, half_angle, k, obstacles)
        source hex + bearing + cone angle + limit (k rings) -> the
        visible hex keys. Cells whose line of sight from the source
        centroid passes through an obstacle hex are occluded (the
        obstacle itself is visible — you can see the hedge, just not
        through it). Sight lines are great-circle, sampled at a third
        of the (self-calibrated) hex pitch and binned in one call.
    aim(src_key, layer, bearing)
        choose a bearing -> the neighbour hex it points at. Accepts a
        single key or a batch of keys with matching bearings.
    walk_to(src_key, dest_key, layer, obstacles)
        choose a destination hex -> the shortest hex path (A* over
        `neighbors`, obstacle hexes impassable), as a key list.

Moving-target pursuit is deliberately absent: that is the caller's
loop (re-aim / re-walk per tick), not a grid verb.

Bearings and distances here are advisory FP geometry over cell
centroids — they steer, they never mint. Every address these verbs
return comes from the core's canonical encode/bin/neighbour chain, so
none of this touches the universality contract (docs/universality.md).

Verbs cross octant seams like any other cells: the seam is a property
of the addressing, not of the field.
"""
import heapq

import numpy as np

try:
    from . import _core as _c
except ImportError:            # dev tree: the flat CMake module, same API
    import hex9_ext as _c

__all__ = ['gc_bearing', 'gc_angle', 'wrap180',
           'encode_keys', 'encode_keys_multi', 'centroids',
           'vision_cone', 'aim', 'walk_to']


# ── spherical helpers ────────────────────────────────────────────────

def gc_bearing(lon1, lat1, lon2, lat2):
    """Initial great-circle bearing, degrees clockwise from north."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(np.subtract(lon2, lon1))
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(y, x)) % 360.0


def gc_angle(lon1, lat1, lon2, lat2):
    """Central angle, degrees."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(np.subtract(lon2, lon1))
    h = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2)
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))))


def wrap180(a):
    """Fold angle differences into (-180, 180]."""
    return (np.asarray(a) + 180.0) % 360.0 - 180.0


def _xyz(lon, lat):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)], axis=-1)


def _lonlat(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return (np.degrees(np.arctan2(v[..., 1], v[..., 0])),
            np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))))


# ── key plumbing (full uuid in, bin KEY out — walk by re-encode) ─────

def encode_keys(lon, lat, layer):
    """Positions -> hex keys, one batch call. The caller's bind verb."""
    u = _c.encode(np.atleast_1d(np.asarray(lon, dtype=np.float64)),
                  np.atleast_1d(np.asarray(lat, dtype=np.float64)))
    b = _c.bin(np.asarray(u, dtype=np.uint8).reshape(-1, 16), layer)
    return [bytes(r) for r in np.asarray(b, dtype=np.uint8).reshape(-1, 16)]


def encode_keys_multi(lon, lat, layers):
    """Positions -> {layer: [keys]} for SEVERAL layers at once.

    The hierarchy leverage: binning is a re-bin of the SAME full uuid
    at each depth, so a coarse cell costs no extra encode — one call
    fixes a point on the sphere and every resolution follows. A caller
    can therefore run its processes at whatever grain each deserves
    (fine forage, coarse neighbourhoods) off a single bind."""
    u = np.asarray(
        _c.encode(np.atleast_1d(np.asarray(lon, dtype=np.float64)),
                  np.atleast_1d(np.asarray(lat, dtype=np.float64))),
        dtype=np.uint8).reshape(-1, 16)
    out = {}
    for layer in layers:
        b = np.asarray(_c.bin(u, layer), dtype=np.uint8).reshape(-1, 16)
        out[layer] = [bytes(r) for r in b]
    return out


def centroids(keys):
    """Keys -> (lon, lat) arrays, one batch call."""
    a = np.frombuffer(b''.join(keys), dtype=np.uint8).reshape(-1, 16)
    lon, lat = _c.decode(a)
    return np.asarray(lon), np.asarray(lat)


def _full(keys):
    """Bin keys -> (n, 16) full-depth uuids by re-encoding centroids —
    the guaranteed re-derivation path (bins are not addresses)."""
    lon, lat = centroids(keys)
    return np.asarray(_c.encode(lon, lat), dtype=np.uint8).reshape(-1, 16)


def _pitch_of(key, layer):
    """Median centroid step to the neighbours — the local hex pitch."""
    nbs, cnt = _c.neighbors(_full([key]), layer)
    nb = np.asarray(nbs, dtype=np.uint8)[0, :int(np.asarray(cnt).reshape(-1)[0])]
    lon0, lat0 = centroids([key])
    nlon, nlat = _c.decode(nb)
    return float(np.median(gc_angle(lon0, lat0, nlon, nlat)))


# ── the three verbs ──────────────────────────────────────────────────

def vision_cone(src_key, layer, bearing, half_angle, k,
                obstacles=frozenset()):
    """Visible hex keys from src_key: within k rings, within
    half_angle of bearing, and not occluded by an obstacle hex.
    The source hex is always visible. Returns a set of keys."""
    u = _full([src_key])
    slon, slat = (float(v[0]) for v in centroids([src_key]))
    disk = np.asarray(_c.k_disk(u[0], layer, k), dtype=np.uint8).reshape(-1, 16)
    keys = [bytes(r) for r in disk]
    dlon, dlat = _c.decode(disk)
    br = gc_bearing(slon, slat, dlon, dlat)
    dist = gc_angle(slon, slat, dlon, dlat)
    in_cone = np.abs(wrap180(br - bearing)) <= half_angle

    cand = [i for i in range(len(keys))
            if keys[i] != src_key and in_cone[i]]
    seen = {src_key}
    if not cand:
        return seen
    if not obstacles:
        seen.update(keys[i] for i in cand)
        return seen

    # occlusion: sample every sight line at pitch/3, bin all samples
    # in ONE encode/bin round, then test obstacle membership per line
    pitch = _pitch_of(src_key, layer)
    a = _xyz(slon, slat)
    ns = np.maximum(2, np.ceil(dist[cand] / (pitch / 3)).astype(int))
    pts, owner = [], []
    for j, i in enumerate(cand):
        t = (np.arange(ns[j]) + 0.5) / ns[j]
        seg = a[None, :] * (1 - t[:, None]) + _xyz(dlon[i], dlat[i]) * t[:, None]
        pts.append(seg)
        owner.extend([j] * ns[j])
    slon_s, slat_s = _lonlat(np.concatenate(pts))
    skeys = encode_keys(slon_s, slat_s, layer)
    blocked = np.zeros(len(cand), dtype=bool)
    for kk, j in zip(skeys, owner):
        if kk in obstacles and kk != src_key and kk != keys[cand[j]]:
            blocked[j] = True
    seen.update(keys[i] for j, i in enumerate(cand) if not blocked[j])
    return seen


def aim(src_key, layer, bearing):
    """The neighbour hex whose centroid bearing best matches.

    Vectorised: pass one key (bytes) and a scalar bearing for one
    neighbour back, or a sequence of keys with matching bearings (or
    one shared bearing) for a list — one neighbour/decode round for
    the whole batch."""
    single = isinstance(src_key, (bytes, bytearray))
    keys = [bytes(src_key)] if single else [bytes(r) for r in src_key]
    bear = np.atleast_1d(np.asarray(bearing, dtype=np.float64))
    if bear.size == 1 and len(keys) > 1:
        bear = np.repeat(bear, len(keys))
    nbs, cnt = _c.neighbors(_full(keys), layer)
    nbs = np.asarray(nbs, dtype=np.uint8)
    cnt = np.asarray(cnt, dtype=np.int64).reshape(-1)
    n, m = nbs.shape[0], nbs.shape[1]
    valid = np.arange(m)[None, :] < cnt[:, None]
    slon, slat = centroids(keys)
    nlon, nlat = _c.decode(nbs[valid])
    br = gc_bearing(np.repeat(slon, cnt), np.repeat(slat, cnt),
                    np.asarray(nlon), np.asarray(nlat))
    diff = np.full((n, m), np.inf)
    diff[valid] = np.abs(wrap180(br - np.repeat(bear, cnt)))
    j = np.argmin(diff, axis=1)
    out = [bytes(nbs[i, j[i]]) for i in range(n)]
    return out[0] if single else out


def walk_to(src_key, dest_key, layer, obstacles=frozenset(),
            max_expand=5000):
    """Shortest hex path src -> dest (A*; obstacle hexes impassable).
    Returns the key list (first = src, last = dest), or None if no
    path is found within max_expand node expansions."""
    dlon, dlat = (float(v[0]) for v in centroids([dest_key]))
    pitch = _pitch_of(src_key, layer)
    slon, slat = (float(v[0]) for v in centroids([src_key]))

    openq = [(float(gc_angle(slon, slat, dlon, dlat)) / pitch * 0.95,
              0, src_key)]
    came, g = {src_key: None}, {src_key: 0}
    for _ in range(max_expand):
        if not openq:
            break
        _, gc, cur = heapq.heappop(openq)
        if cur == dest_key:
            path = [cur]
            while came[path[-1]] is not None:
                path.append(came[path[-1]])
            return path[::-1]
        if gc > g[cur]:
            continue
        nbs, cnt = _c.neighbors(_full([cur]), layer)
        nb = np.asarray(nbs, dtype=np.uint8)[0, :int(np.asarray(cnt).reshape(-1)[0])]
        nlon, nlat = _c.decode(nb)
        hh = np.asarray(gc_angle(np.asarray(nlon), np.asarray(nlat),
                                 dlon, dlat)) / pitch * 0.95
        for i in range(nb.shape[0]):
            kk = bytes(nb[i])
            if kk in obstacles and kk != dest_key:
                continue
            if kk not in g or g[cur] + 1 < g[kk]:
                g[kk] = g[cur] + 1
                came[kk] = cur
                heapq.heappush(openq, (g[kk] + float(hh[i]), g[kk], kk))
    return None
