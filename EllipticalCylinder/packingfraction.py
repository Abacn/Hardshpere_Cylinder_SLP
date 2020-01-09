import numpy as np
from scipy.special import ellipk as ellipk
from iomodule import readBoxSize

def getPackingFraction(fname):
    """Get packing fraction from a filename"""
    N, boxNow, diameter = readBoxSize(fname)
    rho = N/boxNow[2]
    rhoa = rho*boxNow[0]
    p = 2*boxNow[0]/diameter
    if boxNow[0] == boxNow[1]:
        # circular
        eta = rhoa/(p*(p+1.0)*(p+1.0))*4.0/3.0
    else:
        asr  = boxNow[0]/boxNow[1] # aspect ratio
        sqasr = asr*asr
        term1 = 4.0/(3.0*p)*np.pi*rhoa
        term2 = np.pi*p*p + 2.0*p*(ellipk(1.0 - 1.0/sqasr) + ellipk(1.0 - sqasr)*asr) + np.pi
        eta = term1/term2
    return eta

if __name__ == '__main__':
    eta = getPackingFraction('resultConfig_0.dat')
    print(eta)
