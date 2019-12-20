#!/usr/bin/env python3

from __future__ import division, print_function
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation

usage = "usage: %s FILENAME" % sys.argv[0]

if len(sys.argv) != 2:
    print(usage)
    exit(1)

fname = sys.argv[1]

with h5py.File(fname, 'r') as f:
    list(f.keys())
    data = np.array(f['velocity'][()])

fig = plt.figure()
ax = plt.subplot(111)

vmin = data.min()
vmax = data.max()

im = plt.imshow(
        data[0].T,
        vmin=vmin,
        vmax=vmax,
        cmap="viridis",
        interpolation="none",
        animated=True
        )

def update(i):
    im.set_array(data[i].T)
    return [im]

anim = animation.FuncAnimation(
        fig,
        update,
        range(1, data.shape[0]),
        interval=42,
        blit=True
        )

plt.show()
