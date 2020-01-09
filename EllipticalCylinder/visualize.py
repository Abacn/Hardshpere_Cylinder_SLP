import sys

import numpy as np
from mayavi import mlab

from iomodule import readCoords


def visualize(filename):
    """Visualize the configuration"""
    coords, boxNow, diameter = readCoords(filename)
    xs = boxNow[0]*coords[:, 0]*np.cos(coords[:, 1])
    ys = boxNow[1]*coords[:, 0] * np.sin(coords[:, 1])
    zs = coords[:, 2]
    for x, y, z in zip(xs, ys, zs):
        mlab.points3d(x, y, z, scale_factor=diameter, resolution=100, color=(0.3, 0.3, 1))
    mlab.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'resultConfig_0.dat'
    visualize(filename)