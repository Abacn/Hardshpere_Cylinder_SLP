import numpy as np


def pbc_dz(dz):
    """Periodic boundary condition for pair distance"""
    global zmax, zmax_half

    if dz > zmax_half:
        dz = dz - zmax
    elif dz < -zmax_half:
        dz = dz + zmax
    return dz


def pbc_z(z):
    """Periodic boundary condition for z coordinate"""
    global zmax

    if z > zmax:
        z = z - zmax
    elif z < 0:
        z = z + zmax
    return z


def update_zmax(boxNow):
    global zmax, zmax_half
    zmax = boxNow[2]
    zmax_half = 0.5*zmax


def sqdistance(cor1, cor2, parameter, dz=None):
    """Calculate the square of distance by given coordinate pairs. Use dz as dz component"""
    if dz is None:
        dz = cor1[2]-cor2[2]
        dz = pbc_dz(dz)
    if parameter["radMax"] > 0:
        # general case
        dx = parameter["Ra"]*(cor1[0]*np.cos(cor1[1]) - cor2[0]*np.cos(cor2[1]))
        dy = parameter["Rb"]*(cor1[0]*np.sin(cor1[1]) - cor2[0]*np.sin(cor2[1]))
    else:
        # balls one tha boundary of the ellipse
        dx = parameter["Ra"]*(np.cos(cor1[1]) - np.cos(cor2[1]))
        dy = parameter["Rb"]*(np.sin(cor1[1]) - np.sin(cor2[1]))
    d = dx*dx + dy*dy + dz*dz
    return d
