#!/usr/bin/env python3

from __future__ import division, print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

usage = "usage: %s LENGTH HEIGHT FILENAME" % sys.argv[0]

if len(sys.argv) != 4:
    print(usage)
    exit(1)

N = int(sys.argv[1])
M = int(sys.argv[2])
fname = sys.argv[3]

data = np.loadtxt(fname)

t = data[:,0]

fig = plt.figure()
ax = plt.subplot(111)

vmin = data[:,1:].min()
vmax = data[:,1:].max()

im = plt.imshow(
        np.reshape(data[0,1:], (N,M)).T,
        vmin=vmin,
        vmax=vmax,
        cmap="viridis",
        interpolation="none",
        animated=True
        )

def update(i):
    im.set_array(np.reshape(data[i,1:], (N,M)).T)
    return [im]

anim = animation.FuncAnimation(
        fig,
        update,
        range(1,len(t)),
        interval=42,
        blit=True
        )

plt.show()
