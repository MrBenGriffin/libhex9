"""Annotated Hex9 map: land cover (classification mode) + mean elevation, with
north arrow, geographic bounds and projection.

Usage: render_map.py <pts.csv from extract_cluster.py> <out.png> [layer]
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

LAYER = int(sys.argv[3]) if len(sys.argv) > 3 else 11
CLS = {0: ('#dddddd', 'Created'), 1: ('#cfcfcf', 'Unclassified'), 2: ('#9c7a4d', 'Ground'),
       3: ('#c6e17a', 'Low veg'), 4: ('#82c341', 'Med veg'), 5: ('#2f7d24', 'High veg'),
       6: ('#e0473b', 'Building'), 7: ('#777777', 'Noise'), 9: ('#3a8ec9', 'Water')}

d = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1)
lon = np.ascontiguousarray(d[:, 0]); lat = np.ascontiguousarray(d[:, 1])
z = d[:, 2]; cls = d[:, 3].astype(int)
classes = np.arange(0, int(cls.max()) + 1)
uu = h.encode(lon, lat)
b = h.bin(uu, LAYER)
bv = np.ascontiguousarray(b).view(np.dtype((np.void, b.shape[1]))).ravel()
_, idx, inv = np.unique(bv, return_index=True, return_inverse=True)
rep = b[idx]; ng = len(idx)
cnt = np.bincount(inv, minlength=ng)
meanz = np.bincount(inv, weights=z, minlength=ng) / cnt
cc = np.zeros((ng, len(classes)))
for c in classes:
    cc[:, c] = np.bincount(inv, weights=(cls == c), minlength=ng)
mode = cc.argmax(1)
polys = [h.cell(rep[i], LAYER, 1)[:, :2] for i in range(ng)]

lo0, lo1, la0, la1 = lon.min(), lon.max(), lat.min(), lat.max()
aspect = 1.0 / np.cos(np.radians((la0 + la1) / 2))    # equal-distance lon/lat


def deco(ax):
    ax.autoscale(); ax.set_aspect(aspect); ax.set_axis_off()
    ax.annotate('N', xy=(0.965, 0.97), xytext=(0.965, 0.86), xycoords='axes fraction',
                ha='center', va='center', fontsize=13, fontweight='bold',
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2))   # true lon/lat: N is up
    ns = 'N' if la0 >= 0 else 'S'; ew = 'E' if lo0 >= 0 else 'W'
    ax.text(0.015, 0.015, f'{abs(la0):.3f}-{abs(la1):.3f}°{ns}   '
            f'{abs(lo1):.3f}-{abs(lo0):.3f}°{ew}', transform=ax.transAxes,
            fontsize=8, va='bottom', bbox=dict(fc='white', alpha=0.75, ec='none', pad=2))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10.5), dpi=140)
ax1.add_collection(PolyCollection(polys, facecolors=[CLS.get(int(m), ('#000', ''))[0] for m in mode],
                                  edgecolors='none'))
deco(ax1)
present = sorted(set(int(m) for m in mode))
ax1.legend(handles=[Patch(color=CLS[c][0], label=f'{c}: {CLS[c][1]}') for c in present if c in CLS],
           loc='lower right', fontsize=8, framealpha=0.9, title='Classification')
ax1.set_title(f'Land cover - classification mode  ·  Hex9 L{LAYER}, {ng} cells', fontsize=12)

pc = PolyCollection(polys, array=meanz, cmap='terrain', edgecolors='none')
ax2.add_collection(pc); deco(ax2)
fig.colorbar(pc, ax=ax2, shrink=0.6, label='mean elevation (m)')
ax2.set_title(f'Mean elevation  ·  Hex9 L{LAYER}', fontsize=12)

fig.suptitle('UK EA LiDAR (Environment Agency)  →  Hex9 L%d bins   ·   '
             'source EPSG:27700 (OSGB36 / British National Grid)  →  EPSG:4326 (WGS84)   ·   '
             '%d pts' % (LAYER, len(lon)), fontsize=13, y=0.98)
fig.savefig(sys.argv[2], bbox_inches='tight', facecolor='white')
print(f'L{LAYER}: {ng} cells, elev {meanz.min():.1f}-{meanz.max():.1f}m, '
      f'bounds lon[{lo0:.4f},{lo1:.4f}] lat[{la0:.4f},{la1:.4f}] -> {sys.argv[2]}')
