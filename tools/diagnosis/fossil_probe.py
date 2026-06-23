# fossil_probe.py — ground-truth probe for docs/addressing-doctrine.md.
#
# Part of the Hex9 (H9) Project
# Copyright ©2026, Ben Griffin
# Licensed under the Apache License, Version 2.0
#
# Prints the Python reference (hhg9) answers for the failure-mode cases:
# full UUID, layer-L bin (from the point AND from the uuid — must agree),
# and the cell centre for the bin's hexagon. Compare against libhex9 /
# postgis_hex9 (h9_bin must be byte-identical; h9_decode(bin) is the
# fossil — see F2 in the doctrine doc).
#
# Run with cwd = the hhg9 repo (so `import hhg9` resolves):
#   cd /Users/ben/Documents/Projects/PyCharm/hex9 && python3 <this file>

import numpy as np
from hhg9 import Registrar, Points
from hhg9.h9.uuid_address import h9_enc, h9_bin_pts, h9_bin
import hhg9.h9.region as rg

reg = Registrar()
g_gcd = reg.domain('g_gcd')
b_oct = reg.domain('b_oct')

# (lat, lon, layer) — F2 cases from docs/addressing-doctrine.md.
cases = [
    (0.0, -90.0, 1),        # octahedron vertex
    (51.5074, -0.1276, 1),  # Westminster — fails L1,3,4,5,10; clean 2,6,7,8,9
    (51.5074, -0.1276, 4),
    (51.5074, -0.1276, 8),  # clean layer (but note the ~350 m F4 gap)
    (0.0, -30.0, 2),        # equator seam
    (-30.0, -150.0, 4),
    (0.0091, -89.9863, 1),  # F3 garbage-rebin point (1.5 km off the vertex)
    (0.0091, -89.9863, 2),
]

for lat, lon, L in cases:
    pos = Points(np.atleast_2d([lat, lon]), g_gcd)
    bry = reg.project(pos, [g_gcd, b_oct])
    oc, mo = bry.cm()
    full = h9_enc(bry)[0]
    bin_from_pt = h9_bin_pts(bry, L)[0]
    bin_from_uuid = h9_bin([full], L, reg)[0]
    cells = rg.xy_regions(bry.coords, mo, L)
    xym = rg.regions_xy(cells)
    ctr = Points(xym[:, :2], oid=oc, domain=b_oct)
    gcd = reg.project(ctr, [b_oct, g_gcd])
    clat, clon = float(gcd.coords[0, 0]), float(gcd.coords[0, 1])
    print(f"({lat},{lon}) L{L}")
    print(f"  full          : {full}")
    print(f"  bin(point)    : {bin_from_pt}")
    print(f"  bin(uuid)     : {bin_from_uuid}  same={bin_from_pt == bin_from_uuid}")
    print(f"  cell centre   : lat={clat:.6f} lon={clon:.6f}")
