import numpy as np
from tools import pbc_dz, pbc_z, sqdistance

def checkOverlap(idx, attempt, coords):
    """Check if the new coordinate (attempt) will cause overlap"""
    for i in range(coords.shape[0]):
        if idx == i: continue
        dz = attempt[2] - coords[i, 2]
        dz = pbc_dz(dz)
        if np.abs(dz) < 1:
            sqrd = attempt[0]*attempt[0] + coords[i, 0]*coords[i, 0] - 2*attempt[0]*coords[i, 0]*np.cos(attempt[1]-coords[i, 1]) + dz*dz
            if sqrd < 1:
                return i
    return None

def mcmove(coords, boxNow, rOption=True):
    """Conduct Monte Carlo move"""
    N = coords.shape[0]
    mcmax = 0.1
    count = 0
    for rp in range(1000):
        i = np.random.randint(N)
        if rOption:
            icoord = np.random.randint(3)
        else:
            icoord = np.random.randint(1, 3)
        if rp%100==99 and 0==count:
            # maximum move step too high
            mcmax /= 10
        attempt = np.array(coords[i, :])
        attempt[icoord] += mcmax*(2*np.random.random()-1)
        if icoord == 0:
            # changing r
            if attempt[0] < 0:
                attempt[0] = -attempt[0]
            elif attempt[0] > boxNow[0]:
                attempt[0] = 2*boxNow[0] - attempt[0]
        elif icoord == 2:  # changing z need to check PBC
            attempt[2] = pbc_z(attempt[2])
        status = checkOverlap(i, attempt, coords)
        if status is None:
            # conduct movement
            coords[i, icoord] = attempt[icoord]
            count += 1
    return count
