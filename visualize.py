from iomodule import readCoords
import numpy as np
from mayavi import mlab
import sys

def visualize(filename):
    """Visualize the configuration"""
    coords, boxNow = readCoords(filename)
    xs = coords[:, 0] * np.cos(coords[:, 1])
    ys = coords[:, 0] * np.sin(coords[:, 1])
    zs = coords[:, 2]
    for x, y, z in zip(xs, ys, zs):
        mlab.points3d(x, y, z, resolution=100, color=(0.3, 0.3, 1))
    mlab.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'testInitialConfig.dat'
    visualize(filename)