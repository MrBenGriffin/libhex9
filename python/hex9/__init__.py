"""hex9 — the Hex9 / H9 discrete global grid system.

Deterministic aperture-9 hexagonal addressing on the AKW octahedral
projection: the same lon/lat yields the bit-identical cell uuid on every
platform (see libhex9's docs/projection.md and CHANGELOG 2.1.0
"universality"). This package is the compiled libhex9 core with numpy
array vectorisation; the pure-Python research implementation is `hhg9`.

    import numpy as np, hex9
    uuids = hex9.encode(np.array([-3.19]), np.array([55.95]))
    hex9.decode(uuids)

The bring-your-own-projection surface (encode_boct / decode_boct /
cell_ring_boct) admits a caller's own sphere→octahedron map to the full
cell machinery — such addresses are non-canonical: projection identity is
dataset metadata; never mix projections within one dataset.
"""

from ._core import *  # noqa: F401,F403 — the nanobind core IS the API
from .selftest import selftest  # noqa: F401 — runtime universality pins
from . import verbs  # noqa: F401 — cell-first grid verbs (aim / walk_to / vision_cone)

from ._core import version as _version

#: libhex9 version, e.g. "2.2.0" (the wheel version matches it).
__version__ = _version().split()[1]

del _version
