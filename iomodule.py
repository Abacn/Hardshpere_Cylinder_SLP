import json
import re
import numpy as np


def readParameter(filename='config.json'):
    """Read in parameters from .json file
    {
  "NParticle": 50,        // Number of ball of a single box
  "diameter": 2.2,        // diameter of the cylinder
  "transMax": 0.002,   	  // translation max for sphere movements (given as fraction of largest sphere radii)
  "compMax":  0.002,      // compression max for change in strain matrix (given as fraction of largest sphere radii)
  "radMax":   0.002,      // maximal radial movement
  "rotMax":   0.002,      // max change of theta
  "maxIters": 20000,      // number of SLP steps to take before stopping
  "termTol":  1e-6        // volume-based termination criteria over 2 runs (i.e., if volume doesn't decrease at least this amount over 2 SLP steps, terminate)
}"""
    with open(filename) as data_file:
        data = json.load(data_file)
    return data


def writeCoords(coords, boxNow, filename=None):
    """Write coordinates to file or string (if filename is None)"""
    results = []
    results.append('Cylinder size: (%.4f, %.4f, %.4f)' % tuple(boxNow))
    results.append('Number of spheres: %d' % (coords.shape[0]))
    results.append('Corrdinates:')
    for i in range(coords.shape[0]):
        results.append(' %.4f %.4f %.4f' % tuple(coords[i, :]))
    resultstr = '\n'.join(results)
    if filename is not None:
        with open(filename, 'w') as fout:
            fout.write(resultstr)
    return resultstr


def readCoords(filename):
    """Read coordinates as coords"""
    with open(filename, 'r') as fin:
        line = fin.readline()  # Cylinder size
        boxNow = [float(i) for i in re.search(r'\((.+)\)', line).group(1).split(',')]
        N = int(fin.readline().split(':')[1])  # Number of particle
        _ = fin.readline()
        corrdstrlst = [np.fromstring(s.strip(), sep=' ') for s in fin.readlines()]
        coords = np.array(corrdstrlst)

    return coords, boxNow


