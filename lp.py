import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog
import neilist, randomconfig
from tools import pbc_dz, pbc_z, update_zmax
from mcmove import mcmove
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
        self.c_arr = np.zeros(self.N_var)
        self.c_arr[0] = -1  # minimize -epsilon=maximize epsilon

    def optimize(self):
        # check r
        thetaidx_start = 1
        dzidx_start = thetaidx_start + self.N
        dridx_start = dzidx_start + self.N
        returnstatus = 0
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
            ia = []
            ja = []
            ar = []
            b_ub = []
            pac = 0
            for i_idx in range(self.N):
                sini = np.sin(self.coords[i_idx, 1])
                cosi = np.cos(self.coords[i_idx, 1])
                sqri = self.coords[i_idx, 0]*self.coords[i_idx, 0]
                for j_idx in self.neigh.neighIdx[i_idx]:
                    sinj = np.sin(self.coords[j_idx, 1])
                    cosj = np.cos(self.coords[j_idx, 1])
                    sqrj = self.coords[j_idx, 0]*self.coords[j_idx, 0]
                    rirj = self.coords[i_idx, 0]*self.coords[j_idx, 0]
                    dz = self.coords[i_idx, 2] - self.coords[j_idx, 2]
                    dz = pbc_dz(dz)
                    cefc = sqri + sqrj - 2*rirj*np.cos(self.coords[i_idx, 1] - self.coords[j_idx, 1]) + dz*dz
                    if cefc < self.neigh.sqrNNLextraDist:
                        # Add to inequal constraint
                        # depsilon
                        deps = 2*dz*dz
                        # dtheta1, dtheta2
                        sinsub = sini*cosj-cosi*sinj  # sin(i-j)
                        dth2 = 2*rirj*sinsub
                        dth1 = -dth2
                        # dz1, dz2
                        dz2 = 2*dz
                        dz1 = -dz2
                        # dr1, dr2
                        cossub = cosi*cosj+sini*sinj  # cos(i-j)
                        dr1 = -2*(self.coords[i_idx, 0] + self.coords[j_idx, 0]*cossub)
                        dr2 = -2*(self.coords[j_idx, 0] + self.coords[i_idx, 0]*cossub)
                        # add to list
                        ia.append(pac); ja.append(0); ar.append(deps)
                        ia.append(pac); ja.append(thetaidx_start+i_idx); ar.append(dth1)
                        ia.append(pac); ja.append(thetaidx_start+j_idx); ar.append(dth2)
                        ia.append(pac); ja.append(dzidx_start+i_idx); ar.append(dz1)
                        ia.append(pac); ja.append(dzidx_start+j_idx); ar.append(dz2)
                        ia.append(pac); ja.append(dridx_start+i_idx); ar.append(dr1)
                        ia.append(pac); ja.append(dridx_start+j_idx); ar.append(dr2)
                        b_ub.append(cefc-1)
                        pac += 1
            if 0 == pac:
                epsz = self.parameter["compMax"]
                rescale = 1 - epsz
                for rp in range(self.N):
                    self.coords[rp, 2] *= rescale
                self.boxNow[2] *= 1 - epsz
            else:
                matdim = (pac, self.N_var)
                lpmat = coo_matrix((ar, (ia, ja)), shape=matdim)
                rst = linprog(c=self.c_arr, A_ub=lpmat.toarray(), b_ub=b_ub, bounds=self.bounds)

                if rst.status is 2:
                    fixstatus = self._tryfix()
                    if fixstatus is 0:
                        returnstatus = 2
                        break
                    else:
                        continue
                epsz = rst.x[0]

                if epsz < self.parameter["termTol"]:
                    break  # algorithm terminated
                else:
                    # update positions
                    rescale = 1 - epsz
                    self.boxNow[2] *= rescale
                    update_zmax(self.boxNow)
                    for i in range(self.N):
                        self.coords[i, 0] += rst.x[dridx_start+i]
                        self.coords[i, 1] += rst.x[thetaidx_start + i]
                        self.coords[i, 2] = pbc_z(rescale*(self.coords[i, 2] + rst.x[dzidx_start + i]))
            print(rescale, self.boxNow[2])
        return returnstatus

    def _tryfix(self):
        # trying to fix LP fail by randomly moving partcles
        status = mcmove(self.coords, self.boxNow)
        print("Trying to fix LP failure... Conducted %d movement" % status)
        self.neigh.genNeiList(self.coords)
        return status

    def getCoords(self):
        return self.coords

    def getBox(self):
        return self.boxNow

    def getDensity(self):
        return self.N/self.boxNow[2]
