import neilist, randomconfig
import numpy as np
from scipy.optimize import linprog

class LPCore:
    """Linear programming core"""
    def __init__(self, parameter, coords=None, boxNow=None):
        if coords is None or boxNow is None:
            coords, boxNow = randomconfig.genconfig(parameter)
        self.neigh = neilist.NeiList(parameter)
        self.parameter = parameter
        self.coords = coords
        self.boxNow = boxNow
        self.N = coords.shape[0]
        # set variable constaints
        # \epsilon_z, d r, d \theta, d z
        self.bounds = [(0, parameter["compMax"])]
        for i in range(self.N):
            self.bounds.append((-parameter["rotMax"], parameter["rotMax"]))  # dtheta
        for i in range(self.N):
            self.bounds.append((-parameter["transMax"], parameter["transMax"]))  # dz
        for i in range(self.N):
            self.bounds.append((-parameter["radMax"], parameter["radMax"]))  # dr
        self.N_var = len(self.bounds)  # Number of variables

    def optimize(self):
        # check r
        thetaidx_start = 1
        dzidx_start = thetaidx_start + self.N
        dridx_start = dzidx_start + self.N
        for niter in range(self.parameter["maxIters"]):
            if niter%self.neigh.updInterval == 0:
                # update neighbor list
                self.neigh.genNeiList(self.coords)
            for i in range(self.N):
                r_ub = self.parameter["radMax"]
                r_lb = -r_ub
                if self.coords[i, 0] < r_ub:
                    r_lb = -self.coords[i, 0]
                elif self.coords[i, 0]+r_ub > self.boxNow[0]:
                    r_ub = self.boxNow[0] - self.coords[i, 0]
                self.bounds[dridx_start+i] = (r_lb, r_ub)
            # set |r_i - r_j| <= 1
            self.neigh.resetIter()
            for iidx, jidx in self.neigh:
                print(iidx, jidx)
            exit(1)

    def getCoords(self):
        return self.coords

    def getBox(self):
        return self.boxNow

    def getDensity(self):
        return self.N/self.boxNow[2]
