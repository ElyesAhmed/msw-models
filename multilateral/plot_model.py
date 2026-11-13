#!/usr/bin/env python3
"""Plot the multilateral stress case: map view of the wells over the oil column, and a J-K section."""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import make_multilateral as mm

deck = sys.argv[1] if len(sys.argv) > 1 else 'SECTOR01_MULTILAT.DATA'
geoms = np.load(deck.replace('.DATA', '_geom.npy'), allow_pickle=True)
act, depth = mm.load_props()

oil = (act > 0) & (depth > mm.GOC) & (depth < mm.OWC)
colors = plt.cm.tab10(np.arange(len(geoms)))

fig = plt.figure(figsize=(13, 5.4))
ax = [fig.add_subplot(1, 2, 1, projection='3d'), fig.add_subplot(1, 2, 2)]

# 3D view of all branches (I, J, K), K pointing down.
for g, c in zip(geoms, colors):
    m = np.array(g['mother'])
    ax[0].plot(m[:, 0], m[:, 1], m[:, 2], '-', color=c, lw=3, label='%s (K=%d)' % (g['name'], m[0, 2]))
    ax[0].scatter(m[0, 0], m[0, 1], m[0, 2], color=c, s=30)
    for lat in g['laterals']:
        l = np.array([g['mother'][[p[1] for p in g['mother']].index(lat[0][1])]] + lat)
        ax[0].plot(l[:, 0], l[:, 1], l[:, 2], '-', color=c, lw=2)
ax[0].plot([8] * 18, range(4, 22), [96] * 18, 'b--', lw=2, label='A3H injector (K=95-105)')
ax[0].set(xlabel='I', ylabel='J', zlabel='K', xlim=(1, 15), ylim=(1, 25), zlim=(105, 58),
          title='3D view of the wells')
ax[0].view_init(elev=22, azim=-60)
ax[0].legend(fontsize=7, loc='upper left')

# Section along J, through the injector column I=8: active cells and the oil column.
sec = np.where(oil[:, :, 7], 2, np.where(act[:, :, 7] > 0, 1, 0))
ax[1].imshow(sec, origin='upper', cmap=matplotlib.colors.ListedColormap(['white', '#cfe3f5', '#f2c27b']),
             extent=(0.5, mm.NY + 0.5, mm.NZ + 0.5, 0.5), aspect='auto')
for g, c in zip(geoms, colors):
    m = np.array(g['mother'])
    ax[1].plot(m[:, 1], m[:, 2], '-', color=c, lw=3)
ax[1].plot(range(4, 22), [96] * 18, 'b--', lw=2)
ax[1].set(xlabel='J', ylabel='K (layer)', ylim=(115, 50),
          title='Section I=8 (orange: oil column, blue: water zone); mother bores projected')

tot = lambda k: sum(g[k] for g in geoms)
fig.suptitle('SECTOR01 multilateral stress case: %d producers, %d branches, %d segments, %d connections'
             % (len(geoms), tot('nbranch'), tot('nseg'), tot('nconn')))
fig.tight_layout()
out = deck.replace('.DATA', '_wells.png')
fig.savefig(out, dpi=130)
print(out)
