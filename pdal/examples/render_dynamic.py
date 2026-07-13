"""The adaptive (population-ceiling) digest as a map: one mixed-layer Hex9 cell
set, shown as mean elevation (DEM) + land cover (classification mode). Fine cells
where the cloud is dense, coarse where sparse; the sample is captured exactly.

Usage: render_dynamic.py <pts.csv> <out.png> [min_layer max_layer floor ceiling]
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

MINL = int(sys.argv[3]) if len(sys.argv) > 3 else 7
MAXL = int(sys.argv[4]) if len(sys.argv) > 4 else 11
FLOOR = float(sys.argv[5]) if len(sys.argv) > 5 else 25.0
CEIL = float(sys.argv[6]) if len(sys.argv) > 6 else 100.0
CLS = {0: ('#dddddd', 'Created'), 1: ('#cfcfcf', 'Unclassified'), 2: ('#9c7a4d', 'Ground'),
       3: ('#c6e17a', 'Low veg'), 4: ('#82c341', 'Med veg'), 5: ('#2f7d24', 'High veg'),
       6: ('#e0473b', 'Building'), 7: ('#777777', 'Noise'), 9: ('#3a8ec9', 'Water')}

d = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1)
lon = np.ascontiguousarray(d[:, 0]); lat = np.ascontiguousarray(d[:, 1])
z = d[:, 2]; cls = d[:, 3].astype(int)
classes = np.arange(0, int(cls.max()) + 1)
uu = h.encode(lon, lat)

cu, lay, val, npts, assign = h.adaptive(uu, MINL, MAXL, CEIL, FLOOR)
assign = assign.astype(np.int64); nc = len(cu)
cnt = np.bincount(assign, minlength=nc)
meanz = np.bincount(assign, weights=z, minlength=nc) / np.maximum(cnt, 1)
cc = np.zeros((nc, len(classes)))
for c in classes:
    cc[:, c] = np.bincount(assign, weights=(cls == c), minlength=nc)
mode = cc.argmax(1)

# mixed layers overlap geometrically -> draw COARSEST first, finest on top
order = np.argsort(lay, kind='stable')
polys = [h.cell(cu[i], int(lay[i]), 1)[:, :2] for i in order]
meanz_o = meanz[order]; mode_o = mode[order]
mix = ', '.join(f'L{k}:{v}' for k, v in sorted(Counter(lay.tolist()).items()))

lo0, lo1, la0, la1 = lon.min(), lon.max(), lat.min(), lat.max()
aspect = 1.0 / np.cos(np.radians((la0 + la1) / 2))


def deco(ax):
    ax.autoscale(); ax.set_aspect(aspect); ax.set_axis_off()
    ax.annotate('N', xy=(0.965, 0.97), xytext=(0.965, 0.86), xycoords='axes fraction',
                ha='center', va='center', fontsize=13, fontweight='bold',
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2))
    ns = 'N' if la0 >= 0 else 'S'; ew = 'E' if lo0 >= 0 else 'W'
    ax.text(0.015, 0.015, f'{abs(la0):.3f}-{abs(la1):.3f}°{ns}   '
            f'{abs(lo1):.3f}-{abs(lo0):.3f}°{ew}', transform=ax.transAxes,
            fontsize=8, va='bottom', bbox=dict(fc='white', alpha=0.75, ec='none', pad=2))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10.5), dpi=140)
pc = PolyCollection(polys, array=meanz_o, cmap='terrain', edgecolors='#00000010', linewidths=0.05)
ax1.add_collection(pc); deco(ax1)
fig.colorbar(pc, ax=ax1, shrink=0.6, label='mean elevation (m)')
ax1.set_title(f'Adaptive DEM — mean elevation  ·  {nc} cells', fontsize=12)

ax2.add_collection(PolyCollection(polys, facecolors=[CLS.get(int(m), ('#000', ''))[0] for m in mode_o],
                                  edgecolors='#00000010', linewidths=0.05))
deco(ax2)
present = sorted(set(int(m) for m in mode))
ax2.legend(handles=[Patch(color=CLS[c][0], label=f'{c}: {CLS[c][1]}') for c in present if c in CLS],
           loc='lower right', fontsize=8, framealpha=0.9, title='Classification')
ax2.set_title('Adaptive land cover — classification mode', fontsize=12)

fig.suptitle('UK EA LiDAR → Hex9 ADAPTIVE digest (L%d–%d, floor=%g)   ·   %s   ·   '
             'EPSG:27700 → EPSG:4326   ·   %d pts'
             % (MINL, MAXL, FLOOR, mix, len(lon)), fontsize=12, y=0.98)
fig.savefig(sys.argv[2], bbox_inches='tight', facecolor='white')
print(f'adaptive {nc} cells ({mix}), elev {meanz.min():.1f}-{meanz.max():.1f}m -> {sys.argv[2]}')
