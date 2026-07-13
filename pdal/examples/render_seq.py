"""The DGGS story in one figure: fixed layer sweep L7..L11 + the adaptive
population-ceiling digest, coloured by classification (mode).

Usage: render_seq.py <pts.csv from extract_cluster.py> <out.png>
Set HEX9_BUILD if the libhex9 build dir isn't ../../build.
"""
import os, sys, numpy as np
sys.path.insert(0, os.environ.get('HEX9_BUILD',
                os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'build')))
import hex9_ext as h
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
from collections import Counter

CLS = {0: ('#dddddd', 'Created'), 1: ('#cfcfcf', 'Unclassified'), 2: ('#9c7a4d', 'Ground'),
       3: ('#c6e17a', 'Low veg'), 4: ('#82c341', 'Med veg'), 5: ('#2f7d24', 'High veg'),
       6: ('#e0473b', 'Building'), 7: ('#777777', 'Noise'), 9: ('#3a8ec9', 'Water')}

d = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1)
lon = np.ascontiguousarray(d[:, 0]); lat = np.ascontiguousarray(d[:, 1])
z = d[:, 2]; cls = d[:, 3].astype(int)
n = len(lon)
uu = h.encode(lon, lat)
classes = np.arange(0, int(cls.max()) + 1)


def modeclass(inv, ng):
    cc = np.zeros((ng, len(classes)))
    for c in classes:
        cc[:, c] = np.bincount(inv, weights=(cls == c), minlength=ng)
    return cc.argmax(1)


def keys(b):
    b = np.ascontiguousarray(b)
    return b.view(np.dtype((np.void, b.shape[1]))).ravel()


def colsof(mode):
    return [CLS.get(int(m), ('#000000', ''))[0] for m in mode]


panels = []; seen = set()
for L in [7, 8, 9, 10, 11]:
    b = h.bin(uu, L)
    _, idx, inv = np.unique(keys(b), return_index=True, return_inverse=True)
    rep = b[idx]; ng = len(idx)
    mode = modeclass(inv, ng)
    polys = [h.cell(rep[i], L, 1)[:, :2] for i in range(ng)]
    panels.append((f'Fixed L{L}  -  {ng} cells', polys, colsof(mode)))
    seen.update(int(m) for m in mode)

# adaptive digest; mixed layers overlap geometrically, so draw COARSEST first
# and finest last (fine detail on top). Tune `floor` for the LOD mix.
cu, lay, val, npts, assign = h.adaptive(uu, 7, 11, 100.0, 25.0)
dmode = modeclass(assign.astype(np.int64), len(cu))
order = np.argsort(lay, kind='stable')
dpolys = [h.cell(cu[i], int(lay[i]), 1)[:, :2] for i in order]
mix = ', '.join(f'L{k}:{v}' for k, v in sorted(Counter(lay.tolist()).items()))
panels.append((f'Adaptive L7-11  -  {len(cu)} cells  ({mix})', dpolys, colsof(dmode[order])))
seen.update(int(m) for m in dmode)

fig, axes = plt.subplots(2, 3, figsize=(22, 13), dpi=140)
for ax, (title, polys, cols) in zip(axes.ravel(), panels):
    ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors='#00000012', linewidths=0.08))
    ax.autoscale(); ax.set_aspect('equal'); ax.set_axis_off()
    ax.set_title(title, fontsize=12)
fig.legend(handles=[Patch(color=CLS[c][0], label=f'{c}: {CLS[c][1]}') for c in sorted(seen) if c in CLS],
           loc='lower center', ncol=len(seen), fontsize=10, frameon=False)
fig.suptitle('UK EA LiDAR → Hex9 classification bins - fixed L7..L11 + adaptive digest',
             fontsize=15, y=0.98)
fig.savefig(sys.argv[2], bbox_inches='tight', facecolor='white')
print('panels:', [p[0] for p in panels], '\nwrote', sys.argv[2])
