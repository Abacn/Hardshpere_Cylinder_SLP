import neilist
import numpy as np
import randomconfig
import warnings
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from tools import pbc_dz, pbc_z, update_zmax, sqdistance

from mcmove import mcmove


class LPCore:
    """Linear programming core"""
    def __init__(self, parameter, coords=None, boxNow=None):
        self.N = int(parameter["NParticle"])
        self.rOption = parameter['radMax'] > 0
        self.mode = 0
        if parameter['mode'] == 'compression':
            # compressing axial length
            self.mode = 0
            self.diam = 1.0    # hard sphere diameter
        elif parameter['mode'] == 'expansion':
            # expanding sphere diameter
            self.mode = 1
            self.diam = float(parameter['z'])/self.N
        if coords is None or boxNow is None:
            coords, boxNow = randomconfig.genconfig(parameter)
        self.neigh = neilist.NeiList(parameter, diameter=self.diam)
        self.parameter = parameter
        self.coords = coords
        self.boxNow = boxNow
        self.sqRa = self.parameter["Ra"]*self.parameter["Ra"]
        self.sqRb = self.parameter["Rb"]*self.parameter["Rb"]
        # set variable constaints
        # \epsilon_z, d r, d \theta, d z
        self.bounds = [(0, parameter["compMax"])]
        for i in range(self.N):
            self.bounds.append((-parameter["rotMax"], parameter["rotMax"]))  # dtheta
        for i in range(self.N):
            self.bounds.append((-parameter["transMax"], parameter["transMax"]))  # dz
        if self.rOption:
            for i in range(self.N):
                self.bounds.append((-parameter["radMax"], parameter["radMax"]))  # dr
        self.N_var = len(self.bounds)  # Number of variables
        self.c_arr = np.zeros(self.N_var)
        self.c_arr[0] = -1  # minimize -epsilon=maximize epsilon

    def optimize(self):
        # check r
        thetaidx_start = 1
        dzidx_start = thetaidx_start + self.N
        dridx_start = dzidx_start + self.N  # Only works when self.parameter['radMax'] > 0
        returnstatus = 1
        lpfailcount = 0
        for niter in range(self.parameter["maxIters"]):
            if niter%self.neigh.updInterval == 0:
                # update neighbor list
                self.neigh.genNeiList(self.coords, self.diam)
            if self.rOption:
                for i in range(self.N):
                    r_ub = self.parameter["radMax"]
                    r_lb = -r_ub
                    if self.coords[i, 0] < r_ub:
                        # ensure the lower bound is no less than 0
                        r_lb = -self.coords[i, 0]
                    elif self.coords[i, 0]+r_ub > 1.0:
                        # ensure the upper bound is no more than 1
                        r_ub = 1.0 - self.coords[i, 0]
                    self.bounds[dridx_start+i] = (r_lb, r_ub)
            # set |r_i - r_j| <= 1
            ia = []  # index of the inequation
            ja = []  # index of the prefactor
            ar = []  # value of the prefactor
            b_ub = []
            pac = 0
            sqrdiam = self.diam*self.diam
            for i_idx in range(self.N):
                sini = np.sin(self.coords[i_idx, 1])
                cosi = np.cos(self.coords[i_idx, 1])
                # sqri = self.coords[i_idx, 0]*self.coords[i_idx, 0]
                for j_idx in self.neigh.neighIdx[i_idx]:
                    sinj = np.sin(self.coords[j_idx, 1])
                    cosj = np.cos(self.coords[j_idx, 1])
                    # sqrj = self.coords[j_idx, 0]*self.coords[j_idx, 0]
                    # rirj = self.coords[i_idx, 0]*self.coords[j_idx, 0]
                    dz = self.coords[i_idx, 2] - self.coords[j_idx, 2]
                    dz = pbc_dz(dz)
                    cefc = sqdistance(self.coords[i_idx, :], self.coords[j_idx, :], self.parameter, dz)
                    if cefc < self.neigh.sqrNNLextraDist*sqrdiam:
                        # Add to inequal constraint
                        # dtheta1, dtheta2
                        dx = self.coords[i_idx, 0]*cosi-self.coords[j_idx, 0]*cosj
                        dy = self.coords[i_idx, 0]*sini-self.coords[j_idx, 0]*sinj
                        # Attention to the signs. Note the inequation applies <=
                        dth2 = -2*self.coords[j_idx, 0]*(self.sqRa*sinj*dx - self.sqRb*cosj*dy)
                        dth1 = 2*self.coords[i_idx, 0]*(self.sqRa*sini*dx - self.sqRb*cosi*dy)
                        # dz1, dz2
                        dz2 = 2*dz
                        dz1 = -dz2
                        # add to list
                        if 0 == self.mode:
                            # depsilon_z
                            deps = 2*dz*dz
                            ia.append(pac); ja.append(0); ar.append(deps)
                        else:
                            # depsilon_diameter
                            deps = 2*sqrdiam
                            ia.append(pac); ja.append(0); ar.append(deps)
                        ia.append(pac); ja.append(thetaidx_start+i_idx); ar.append(dth1)
                        ia.append(pac); ja.append(thetaidx_start+j_idx); ar.append(dth2)
                        ia.append(pac); ja.append(dzidx_start+i_idx); ar.append(dz1)
                        ia.append(pac); ja.append(dzidx_start+j_idx); ar.append(dz2)
                        if self.rOption:
                            # dr1, dr2
                            dr1 = -2*(self.sqRa*cosi*dx + self.sqRb*sini*dy)
                            dr2 = 2*(self.sqRa*cosj*dx + self.sqRb*sinj*dy)
                            ia.append(pac); ja.append(dridx_start+i_idx); ar.append(dr1)
                            ia.append(pac); ja.append(dridx_start+j_idx); ar.append(dr2)
                        b_ub.append(cefc-sqrdiam)
                        pac += 1
            if 0 == pac:
                epsz = self.parameter["compMax"]
                if 0==self.mode:
                    rescale = 1 - epsz
                    for rp in range(self.N):
                        self.coords[rp, 2] *= rescale
                    self.boxNow[2] *= 1 - epsz
                    update_zmax(self.boxNow)
                else:
                    rescale = 1 + epsz
                    self.diam *= rescale
            else:
                matdim = (pac, self.N_var)
                # formalize prefactors into a sparse matrix
                lpmat = coo_matrix((ar, (ia, ja)), shape=matdim)
                try:
                    rst = linprog(c=self.c_arr, A_ub=lpmat.toarray(), b_ub=b_ub, bounds=self.bounds)
                except ValueError:
                    rst.status = 2
                if rst.status is 2:
                    fixstatus = self._tryfix()
                    if fixstatus is 0 or lpfailcount > 10:
                        returnstatus = 2
                        break
                    else:
                        lpfailcount += 1
                        continue
                else:
                    lpfailcount = 0
                epsz = rst.x[0]

                if epsz < self.parameter["termTol"]:
                    returnstatus = 0
                    break  # algorithm terminated
                else:
                    # update positions
                    if 0==self.mode:
                        rescale = 1 - epsz
                        rescalez = rescale
                        self.boxNow[2] *= rescale
                        update_zmax(self.boxNow)
                    else:
                        rescale = 1 + epsz
                        rescalez = 1
                        self.diam *= rescale
                    if self.rOption:
                        for i in range(self.N):
                            self.coords[i, 0] += rst.x[dridx_start + i]
                    for i in range(self.N):
                        self.coords[i, 1] += rst.x[thetaidx_start + i]
                        self.coords[i, 2] = pbc_z(rescalez*(self.coords[i, 2] + rst.x[dzidx_start + i]))
            if niter%2000==0:
                if 0==self.mode:
                    print(rescale, self.boxNow[2])
                else:
                    print(rescale, self.diam)
        return returnstatus

    def _tryfix(self):
        # trying to fix LP fail by randomly moving partcles
        status = mcmove(self.coords, self.boxNow, self.rOption, self.diam)
        print("Trying to fix LP failure... Conducted %d movement" % status)
        self.neigh.genNeiList(self.coords, self.diam)
        return status

    def getCoords(self):
        return self.coords

    def getBox(self):
        return self.boxNow

    def getDensity(self):
        return self.N/self.boxNow[2]
