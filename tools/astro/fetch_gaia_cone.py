#!/usr/bin/env python3
"""Fetch a Gaia DR3 magnitude-laddered cone around a sky POI for the hex9 zoom.

The hex9 zoom gate is m_lim(t) = 6.5 + 1.9t and the viewport half-diagonal is
roughly 184°/3^t, so a magnitude band (g_lo, g_hi] first becomes visible at
t = (g_lo - 6.5)/1.9 and only ever needs a cone of that half-diagonal.  Tiers
are disjoint in magnitude (no cross-tier duplicates) and floored at G > 11.3,
below AT-HYG's completeness, so the all-sky catalogue keeps the bright end.

Queries the VizieR TAP mirror of gaiadr3 (I/355) — the ESA anonymous async
endpoint 500s intermittently; VizieR allows large synchronous pulls.

Usage:
    python3 fetch_gaia_cone.py                     # Baade's Window (default)
    python3 fetch_gaia_cone.py --ra 266.41683 --dec -29.00781 --tag sgr_a

Output: data/gaia_<tag>.csv with columns ra,dec,g,bp_rp (ra/dec degrees).
Credit line for renders: Gaia DR3 · ESA/Gaia/DPAC.
"""
import argparse
import os
import time

import pyvo

TAP_URL = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap'

# (radius_deg, g_lo, g_hi): half-diagonal at the tier's surfacing zoom, with
# margin.  Band width 1.9 mag = one zoom level; radius steps ~×3 to match.
TIERS = [
    (12.0, 11.3, 13.2),
    (4.2, 13.2, 15.1),
    (1.4, 15.1, 17.0),
    (0.48, 17.0, 18.9),
    (0.16, 18.9, 20.8),
]

ADQL = '''SELECT RA_ICRS, DE_ICRS, Gmag, "BP-RP"
FROM "I/355/gaiadr3"
WHERE 1=CONTAINS(POINT('ICRS', RA_ICRS, DE_ICRS),
                 CIRCLE('ICRS', {ra}, {dec}, {r}))
AND Gmag > {g_lo} AND Gmag <= {g_hi}'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ra', type=float, default=270.892,
                    help='POI RA, degrees (default: Baade\'s Window/NGC 6522)')
    ap.add_argument('--dec', type=float, default=-30.034,
                    help='POI dec, degrees')
    ap.add_argument('--tag', default='baade',
                    help='output name: data/gaia_<tag>.csv')
    args = ap.parse_args()

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'data', f'gaia_{args.tag}.csv')
    svc = pyvo.dal.TAPService(TAP_URL)

    def fetch_band(r, g_lo, g_hi, depth=0):
        """Fetch one (radius, mag-band) piece; bisect the band on maxrec.

        TAP truncation is in table (HEALPix) order — spatially biased — so a
        capped result must never be kept.
        """
        q = ADQL.format(ra=args.ra, dec=args.dec, r=r, g_lo=g_lo, g_hi=g_hi)
        for attempt in range(4):
            try:
                t = svc.search(q, maxrec=500_000).to_table()
                break
            except Exception as e:
                if attempt == 3:
                    raise
                print(f'  attempt {attempt + 1} failed ({e}); retrying …',
                      flush=True)
                time.sleep(15 * (attempt + 1))
        if len(t) >= 500_000:
            if g_hi - g_lo < 0.05:
                raise RuntimeError(f'band {g_lo}–{g_hi} at r={r}° still '
                                   f'caps at 0.05 mag — smaller radius needed')
            mid = round((g_lo + g_hi) / 2.0, 3)
            print(f'  {"  " * depth}band {g_lo}–{g_hi} capped; '
                  f'splitting at {mid}', flush=True)
            return (fetch_band(r, g_lo, mid, depth + 1)
                    + fetch_band(r, mid, g_hi, depth + 1))
        print(f'  {"  " * depth}{g_lo} < G ≤ {g_hi}: {len(t):,} rows',
              flush=True)
        return [t]

    tables = []
    for r, g_lo, g_hi in TIERS:
        print(f'tier r={r}° …', flush=True)
        tables.extend(fetch_band(r, g_lo, g_hi))

    os.makedirs(os.path.dirname(out), exist_ok=True)
    n = 0
    with open(out, 'w') as f:
        f.write('ra,dec,g,bp_rp\n')
        for t in tables:
            bp_col = 'BP-RP' if 'BP-RP' in t.colnames else 'BP_RP'
            for ra, dec, g, bp in zip(t['RA_ICRS'], t['DE_ICRS'],
                                      t['Gmag'], t[bp_col]):
                bps = '' if bp is None or bp != bp else f'{float(bp):.4f}'
                f.write(f'{float(ra):.7f},{float(dec):.7f},'
                        f'{float(g):.4f},{bps}\n')
                n += 1
    print(f'wrote {n:,} stars → {out}')


if __name__ == '__main__':
    main()
