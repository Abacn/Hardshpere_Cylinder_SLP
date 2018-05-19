import numpy as np
from tools import *


def checkOverlap(coords, newcoord, idx):
    """Check overlap of particle idx (with coordinates newcoord) with former particles
    coords: coordinate matrix, idx: index of the particle need to be checked"""
    for i in range(idx):
        # distance
        dz = newcoord[2] - coords[i, 2]
        dz = pbc_dz(dz)
        if abs(np.abs(dz) < 1.0): # overlap possible, check the distance
            sqdis = sqdistance(newcoord, coords[i, :], dz)
            if sqdis < 1.0:
                return i
    return None


def genconfig(parameter):
    """Generate random configuration by given parameters"""
    N = parameter['NParticle']
    R = parameter['diameter']/2
    zstep = 5.0    # the box size increases zstep if overlap happens too much
    boxOrigin = np.array([R, 2*np.pi, zstep])
    boxNow = np.array(boxOrigin)
    update_zmax(boxNow)
    coords = np.zeros((N, 3))
    count = 0
    for i in range(N):
        ovlpflag = True
        while ovlpflag:
            if count < 1000:  # z component uniform distributes in the space
                newcoord = np.multiply(np.random.rand(3), boxNow)
            elif count < 2000: # z component is in the last grew space
                newcoord = np.multiply(np.random.rand(3), boxOrigin)
                newcoord[2] = boxNow[2] - newcoord[2]
            else:  # Too many trial is overlapped
                boxNow[2] += zstep
                count = 0
                update_zmax(boxNow)
                continue
            if parameter['radMax'] <= 0:
                # all particles are confined on the boundary
                newcoord[0] = R
            if checkOverlap(coords, newcoord, i) is None:
                coords[i, :] = newcoord
                ovlpflag = False
            else:
                count += 1
    return coords, boxNow


if __name__ == '__main__':
    from iomodule import *
    parameter = readParameter()
    coords, boxNow = genconfig(parameter)
    result = writeCoords(coords, boxNow, 'testInitialConfig.dat')
    print(result)